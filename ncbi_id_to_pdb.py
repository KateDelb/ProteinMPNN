"""Map NCBI Gene / nucleotide / protein records to RCSB PDB structures.

Input: CSV (comma- or tab-delimited) with at least `sequence-id` and `sequence`.
Output: pdb_results/pdb_mapping.csv and unique structure files in
        pdb_results/structures/.

Alignment metrics are calculated locally with Biopython PairwiseAligner
(BLOSUM62, local alignment) against the canonical polymer-entity sequence
provided by the RCSB Data API. Percent identity is identical aligned columns
(including gap columns) / total alignment columns. Query/target coverage is
aligned residues / full query/target sequence length.
"""

import os
import re
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq


# =====================================================
# CONFIGURATION
# =====================================================

# Set NCBI_EMAIL in your environment, or replace the value below.
Entrez.email = os.getenv("NCBI_EMAIL", "katedelb@gmail.com")

# Optional: set NCBI_API_KEY in your environment.
if os.getenv("NCBI_API_KEY"):
    Entrez.api_key = os.environ["NCBI_API_KEY"]

INPUT_FILE = "sequence_db.csv"
OUTPUT_DIR = Path("pdb_results")
PDB_DIR = OUTPUT_DIR / "structures"
OUTPUT_CSV = OUTPUT_DIR / "pdb_mapping.csv"

# RCSB sequence search settings.
IDENTITY_CUTOFF = 0.30
MAX_PDB_RESULTS = 1  # Polymer entities retained per protein sequence.
E_VALUE_CUTOFF = 0.1

# Candidates pulled from RCSB and ranked before MAX_PDB_RESULTS are retained.
# RCSB orders hits by its own alignment score, which favours a long partial
# homologue over a short exact match; ranking locally applies this project's
# criteria instead. Ranking uses the per-hit statistics RCSB returns with the
# search, so widening the pool costs one request rather than one per candidate.
MAX_PDB_CANDIDATES = 200

# Candidates realigned locally before the final choice. RCSB's reported
# identity is measured over a tight BLAST HSP, while the identity this script
# reports comes from a local BLOSUM62 alignment that often extends further and
# scores lower, so the two disagree and ranking on one while reporting the
# other picks the wrong structure. The RCSB statistics are therefore used only
# to shortlist; the shortlist is realigned and reranked on the same metrics
# that reach the CSV.
RERANK_SHORTLIST = 10

# AlphaFold DB fallback. Many proteins have no experimental structure at all;
# a predicted model still gives ProteinMPNN a usable backbone. Models are
# accepted only above a mean-pLDDT floor, since a largely disordered model is
# not a sound template, and every such row is labelled predicted so it stays
# separable from experimental coordinates downstream.
ALPHAFOLD_FALLBACK = True
MIN_ALPHAFOLD_PLDDT = 70.0
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

SOURCE_EXPERIMENTAL = "RCSB PDB (experimental)"
SOURCE_PREDICTED = "AlphaFold DB (predicted)"

# RCSB search occasionally returns a transient 5xx. Without a retry a single
# blip silently costs that protein its structure, so retry with backoff.
SEARCH_ATTEMPTS = 4
SEARCH_RETRY_BACKOFF = 3.0

REQUEST_DELAY = 0.4
HTTP_TIMEOUT = 90

# Redundant-accession collapse. A row qualifies as the definitive mapping for a
# structure when it is a full-length exact match; competing accessions for the
# same sequence-id that resolve to the SAME pdb-id are then dropped.
PERFECT_IDENTITY_PCT = 100.0
PERFECT_QUERY_COVERAGE_PCT = 100.0
MIN_TARGET_COVERAGE_PCT = 95.0

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PDB_DIR.mkdir(parents=True, exist_ok=True)

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "NCBI-to-PDB-mapper/1.0"})

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
RCSB_FILE_URL = "https://files.rcsb.org/download/{pdb_id}.{extension}"

# Local alignment configuration. This is used to calculate identity and
# coverage against the PDB polymer entity sequence returned by RCSB.
ALIGNER = PairwiseAligner()
ALIGNER.mode = "local"
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.open_gap_score = -10.0
ALIGNER.extend_gap_score = -0.5

# Avoid repeatedly requesting the same PDB polymer-entity sequence.
ENTITY_SEQUENCE_CACHE = {}


# =====================================================
# GENERAL UTILITIES
# =====================================================

def clean_nucleotide_sequence(sequence):
    return re.sub(r"\s+", "", str(sequence)).upper()


def clean_protein_sequence(sequence):
    """Normalize protein sequence for alignment with BLOSUM62.

    Rare/ambiguous residues not represented in BLOSUM62 are converted to X.
    Terminal stop symbols are removed; internal stop symbols become X.
    """
    sequence = re.sub(r"\s+", "", str(sequence)).upper().rstrip("*")
    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZ")
    return "".join(aa if aa in allowed else "X" for aa in sequence)


def normalize_identifier(identifier):
    """Strip whitespace and a dangling accession separator from an input ID.

    NCBI rejects an accession with trailing punctuation: `X03179.` returns
    "Failed to understand id", which previously produced an empty GenBank
    record and silently dropped through to the unverified-frame fallback.
    The original string is still written to the output CSV, so joins back to
    the input file keep working.
    """
    identifier = str(identifier).strip().strip(",;|")
    return identifier.rstrip(".") or identifier


def is_gene_id(identifier):
    return identifier.isdigit()


def is_protein_accession(identifier):
    return bool(re.match(r"^(NP|XP|YP|WP|AP|ZP|YP)_", identifier, re.I))


def value_or_blank(row, column):
    value = row.get(column, "")
    return "" if pd.isna(value) else value


# =====================================================
# NCBI RETRIEVAL
# =====================================================

def get_ncbi_record(identifier, database):
    handle = Entrez.efetch(
        db=database,
        id=identifier,
        rettype="gb",
        retmode="text",
    )
    try:
        return list(SeqIO.parse(handle, "genbank"))
    finally:
        handle.close()


def get_protein_records(protein_ids):
    """Fetch FASTA records for NCBI protein Entrez UIDs."""
    if not protein_ids:
        return []

    proteins = []
    # Chunk requests to avoid excessively long EFetch URLs.
    for start in range(0, len(protein_ids), 100):
        chunk = protein_ids[start:start + 100]
        handle = Entrez.efetch(
            db="protein",
            id=",".join(chunk),
            rettype="fasta",
            retmode="text",
        )
        try:
            records = list(SeqIO.parse(handle, "fasta"))
        finally:
            handle.close()

        for record in records:
            proteins.append({
                "protein_accession": record.id,
                "protein_sequence": str(record.seq),
                "protein_name": record.description,
                "organism": "",
                "source": "NCBI Gene protein link",
                "gene_name": "",
                "cds_location": "",
            })

    return proteins


def get_proteins_from_gene(gene_id):
    """Retrieve RefSeq protein records linked to an NCBI Gene ID."""
    protein_ids = []

    for linkname in ("gene_protein_refseq", None):
        kwargs = {
            "dbfrom": "gene",
            "db": "protein",
            "id": gene_id,
        }
        if linkname:
            kwargs["linkname"] = linkname

        handle = Entrez.elink(**kwargs)
        try:
            linksets = Entrez.read(handle)
        finally:
            handle.close()

        for linkset in linksets:
            for link_db in linkset.get("LinkSetDb", []):
                protein_ids.extend(
                    item["Id"] for item in link_db.get("Link", [])
                )

        if protein_ids:
            break

    # Preserve order while removing duplicated UIDs.
    protein_ids = list(dict.fromkeys(protein_ids))
    return get_protein_records(protein_ids)


def extract_cds_from_accession(accession):
    """Extract all annotated CDS translations from a GenBank nucleotide record.

    A nucleotide accession may contain a single gene, a gene cluster, or a
    complete genome. This intentionally returns every translated CDS feature
    in the record, as in the earlier pipeline.
    """
    proteins = []
    records = get_ncbi_record(accession, "nuccore")

    for record in records:
        organism = record.annotations.get("organism", "Unknown")

        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers
            translation = qualifiers.get("translation", [None])[0]
            if not translation:
                continue

            protein_id = qualifiers.get("protein_id", [None])[0]
            gene_name = qualifiers.get(
                "gene", qualifiers.get("product", ["Unknown"])
            )[0]
            product = qualifiers.get("product", [gene_name])[0]

            proteins.append({
                "protein_accession": protein_id or "Not available",
                "protein_sequence": clean_protein_sequence(translation),
                "protein_name": product,
                "organism": organism,
                "source": "NCBI annotated CDS",
                "gene_name": gene_name,
                "cds_location": str(feature.location),
            })

    return proteins


def find_longest_orf(nucleotide_sequence):
    """Find the longest ORF across all six reading frames.

    Each frame is translated and split on stop codons. Within every stop-free
    stretch the first Met is located and the peptide from that Met onward is a
    candidate; the longest candidate across all frames wins. If no stretch
    contains a Met, the longest stop-free stretch is returned instead, flagged
    via has_start_codon.
    """
    sequence = Seq(nucleotide_sequence)
    best = {
        "protein_sequence": "",
        "frame": None,
        "has_start_codon": False,
        "location": "",
    }

    for strand, subject in (("+", sequence), ("-", sequence.reverse_complement())):
        for offset in range(3):
            usable = ((len(subject) - offset) // 3) * 3
            if usable <= 0:
                continue

            peptides = str(subject[offset:offset + usable].translate(table=1))

            # Track the codon offset of each stop-free stretch within the frame.
            cursor = 0
            for stretch in peptides.split("*"):
                start_aa = cursor
                cursor += len(stretch) + 1  # +1 for the stop codon consumed

                met = stretch.find("M")
                if met >= 0:
                    candidate, aa_offset, has_start = stretch[met:], start_aa + met, True
                else:
                    candidate, aa_offset, has_start = stretch, start_aa, False

                if len(candidate) <= len(best["protein_sequence"]):
                    continue

                # Report the location on the supplied (+) strand, 1-based.
                nt_start = offset + aa_offset * 3
                nt_end = nt_start + len(candidate) * 3
                if strand == "-":
                    nt_start, nt_end = len(sequence) - nt_end, len(sequence) - nt_start

                best = {
                    "protein_sequence": candidate,
                    "frame": f"{strand}{offset + 1}",
                    "has_start_codon": has_start,
                    "location": f"[{nt_start + 1}:{nt_end}]({strand})",
                }

    return best


def extract_proteins(row):
    """Resolve a Gene ID, protein accession, or nucleotide accession."""
    identifier = normalize_identifier(row["sequence-id"])
    nucleotide_sequence = clean_nucleotide_sequence(row["sequence"])

    if not re.fullmatch(r"[ACGTUNRYKMSWBDHV]+", nucleotide_sequence):
        raise ValueError("Input sequence contains invalid nucleotide characters")

    if is_gene_id(identifier):
        proteins = get_proteins_from_gene(identifier)
        if proteins:
            return proteins

        # A Gene ID with no linked RefSeq protein is normally a discontinued
        # record (esummary status 2) with no replacement ID, so there is
        # nothing to redirect to. The supplied nucleotide sequence is still
        # usable, so translate that rather than dropping the row.
        print(
            "No RefSeq protein linked to Gene ID (record is likely "
            "discontinued); translating the supplied sequence instead"
        )
        return [orf_fallback(row, identifier, nucleotide_sequence)]

    if is_protein_accession(identifier):
        handle = Entrez.efetch(
            db="protein",
            id=identifier,
            rettype="fasta",
            retmode="text",
        )
        try:
            records = list(SeqIO.parse(handle, "fasta"))
        finally:
            handle.close()

        if records:
            return [{
                "protein_accession": records[0].id,
                "protein_sequence": clean_protein_sequence(records[0].seq),
                "protein_name": records[0].description,
                "organism": "",
                "source": "NCBI protein record",
                "gene_name": "",
                "cds_location": "",
            }]

    # Nucleotide accession: retrieve GenBank CDS annotations.
    proteins = extract_cds_from_accession(identifier)
    if proteins:
        return proteins

    return [orf_fallback(row, identifier, nucleotide_sequence)]


def orf_fallback(row, identifier, nucleotide_sequence):
    """Translate the supplied sequence when NCBI yields no annotated protein.

    Frame 0 of the supplied sequence is rarely the coding frame, since genomic
    records carry UTRs, so all six frames are searched for the longest ORF.
    """
    if identifier.startswith("NC_"):
        raise ValueError("No annotated CDS found in the RefSeq genomic record")

    orf = find_longest_orf(nucleotide_sequence)
    if len(orf["protein_sequence"]) < 25:
        raise ValueError(
            "No annotated CDS found; longest ORF across all six reading "
            "frames is shorter than 25 amino acids"
        )

    return {
        "protein_accession": f"Predicted_{identifier}",
        "protein_sequence": clean_protein_sequence(orf["protein_sequence"]),
        "protein_name": value_or_blank(row, "sequence-name"),
        "organism": "Unknown",
        "source": (
            f"Fallback six-frame translation (frame {orf['frame']}, "
            f"{'start codon' if orf['has_start_codon'] else 'no start codon'})"
        ),
        "gene_name": "",
        "cds_location": orf["location"],
    }


# =====================================================
# BASIC PROTEIN VALIDATION
# =====================================================

def validate_protein(protein_sequence):
    sequence = clean_protein_sequence(protein_sequence)
    problems = []

    if len(sequence) < 25:
        problems.append("Protein shorter than 25 aa; RCSB sequence search requires >=25 aa")

    if not sequence:
        problems.append("Empty protein sequence")

    return problems


# =====================================================
# RCSB SEARCH + PDB ENTITY SEQUENCE RETRIEVAL
# =====================================================

def get_pdb_entity_sequence(pdb_id, entity_id):
    """Get canonical amino-acid sequence for one RCSB polymer entity."""
    cache_key = (pdb_id.upper(), str(entity_id))
    if cache_key in ENTITY_SEQUENCE_CACHE:
        return ENTITY_SEQUENCE_CACHE[cache_key]

    url = RCSB_ENTITY_URL.format(
        pdb_id=pdb_id,
        entity_id=entity_id,
    )
    response = SESSION.get(url, timeout=HTTP_TIMEOUT)
    response.raise_for_status()
    data = response.json()

    entity_poly = data.get("entity_poly", {})
    sequence = entity_poly.get("pdbx_seq_one_letter_code_can") or ""
    if not sequence:
        sequence = entity_poly.get("pdbx_seq_one_letter_code") or ""

    sequence = clean_protein_sequence(sequence)
    ENTITY_SEQUENCE_CACHE[cache_key] = sequence
    return sequence


def calculate_alignment_metrics(query_sequence, target_sequence):
    """Local-align query protein to PDB entity sequence and return metrics.

    query = protein sequence extracted from NCBI.
    target = canonical protein sequence for the PDB polymer entity.

    Percent identity is identical aligned columns divided by all alignment
    columns (gap columns included). Coverage is the count of aligned residues
    divided by the full length of the corresponding sequence.
    Coordinates are 1-based and inclusive.
    """
    query = clean_protein_sequence(query_sequence)
    target = clean_protein_sequence(target_sequence)

    if not query or not target:
        raise ValueError("Cannot align an empty protein sequence")

    alignment = ALIGNER.align(target, query)[0]
    indices = alignment.indices
    target_indices = [int(i) for i in indices[0]]
    query_indices = [int(i) for i in indices[1]]

    alignment_columns = len(target_indices)
    if alignment_columns == 0:
        raise ValueError("Alignment contains no columns")

    paired_columns = [
        (t_idx, q_idx)
        for t_idx, q_idx in zip(target_indices, query_indices)
        if t_idx >= 0 and q_idx >= 0
    ]

    identical = sum(
        target[t_idx] == query[q_idx]
        for t_idx, q_idx in paired_columns
    )

    query_positions = [i for i in query_indices if i >= 0]
    target_positions = [i for i in target_indices if i >= 0]

    query_covered = len(set(query_positions))
    target_covered = len(set(target_positions))

    return {
        "sequence_identity_pct": round(100.0 * identical / alignment_columns, 2),
        "query_coverage_pct": round(100.0 * query_covered / len(query), 2),
        "target_coverage_pct": round(100.0 * target_covered / len(target), 2),
        "alignment_length": alignment_columns,
        "aligned_residue_pairs": len(paired_columns),
        "query_start": min(query_positions) + 1 if query_positions else None,
        "query_end": max(query_positions) + 1 if query_positions else None,
        "target_start": min(target_positions) + 1 if target_positions else None,
        "target_end": max(target_positions) + 1 if target_positions else None,
        "pdb_entity_sequence_length": len(target),
        "alignment_score": round(float(alignment.score), 2),
        "alignment_method": "local BLOSUM62; PairwiseAligner",
    }


def summarize_match_context(item):
    """Provisional identity/coverage for one RCSB hit, from the search reply.

    The sequence service reports BLAST statistics per hit when the search is
    run with verbose verbosity, so candidates can be ranked without fetching
    and aligning every polymer-entity sequence. These figures are used only to
    choose which candidates to keep; the retained ones are then realigned
    locally so the reported columns stay consistent with the rest of the file.
    """
    best = (0.0, 0.0, 0.0)

    for service in item.get("services", []):
        for node in service.get("nodes", []):
            for context in node.get("match_context", []):
                identity = (context.get("sequence_identity") or 0.0) * 100.0

                query_length = context.get("query_length") or 0
                subject_length = context.get("subject_length") or 0

                query_span = (
                    context.get("query_end", 0) - context.get("query_beg", 0) + 1
                )
                subject_span = (
                    context.get("subject_end", 0) - context.get("subject_beg", 0) + 1
                )

                query_coverage = (
                    query_span / query_length * 100.0 if query_length else 0.0
                )
                target_coverage = (
                    subject_span / subject_length * 100.0 if subject_length else 0.0
                )

                best = max(best, (identity, query_coverage, target_coverage))

    return {
        "provisional_identity_pct": best[0],
        "provisional_query_coverage_pct": best[1],
        "provisional_target_coverage_pct": best[2],
    }


def match_rank(match):
    """Order realigned matches by the metrics actually reported in the CSV.

    Mirrors candidate_rank, but reads the locally computed alignment metrics.
    A match whose alignment failed sorts last rather than winning on nulls.
    """
    if match.get("sequence_identity_pct") is None:
        return (1, 0.0, 0.0, 0.0, match["pdb_entity"])

    return (
        0,
        -match["sequence_identity_pct"],
        -match["query_coverage_pct"],
        -match["target_coverage_pct"],
        match["pdb_entity"],
    )


def candidate_rank(candidate):
    """Shortlist candidate structures by this project's match-quality criteria.

    Sequence identity ranks first, then query coverage, then target coverage.
    This favours an exact match over a longer but divergent homologue, which
    is what RCSB's own alignment score does not do.

    Target coverage is deliberately the last numeric key rather than a
    threshold applied first. Gating on MIN_TARGET_COVERAGE_PCT up front lets a
    divergent hit that just clears the threshold beat a perfect match that
    just misses it -- measured on this dataset, that cost eight rows their
    100% identity match, one dropping to 63%. Ranking target coverage
    descending already places any hit above MIN_TARGET_COVERAGE_PCT ahead of
    any hit below it among otherwise equal candidates, so the threshold is
    honoured as a tie-break without overriding identity.

    The entity ID is the final key so repeated runs agree.
    """
    return (
        -candidate["provisional_identity_pct"],
        -candidate["provisional_query_coverage_pct"],
        -candidate["provisional_target_coverage_pct"],
        candidate["pdb_entity"],
    )


def post_rcsb_search(payload):
    """POST a search to RCSB, retrying transient server-side failures.

    Returns the decoded response, or None when RCSB reports no content.
    Raises the last error if every attempt fails, so the row is recorded as a
    search failure rather than being mistaken for a genuine absence of hits.
    """
    last_error = None

    for attempt in range(1, SEARCH_ATTEMPTS + 1):
        try:
            response = SESSION.post(
                RCSB_SEARCH_URL,
                json=payload,
                timeout=HTTP_TIMEOUT,
            )

            if response.status_code == 204:
                return None

            # 4xx means the query itself is wrong; retrying cannot help.
            if 400 <= response.status_code < 500:
                response.raise_for_status()

            if response.status_code >= 500:
                raise requests.HTTPError(
                    f"{response.status_code} Server Error from RCSB search"
                )

            return response.json()

        except (requests.RequestException, ValueError) as error:
            last_error = error
            if attempt == SEARCH_ATTEMPTS:
                break
            delay = SEARCH_RETRY_BACKOFF * attempt
            print(
                f"RCSB search attempt {attempt}/{SEARCH_ATTEMPTS} failed "
                f"({error}); retrying in {delay:.0f}s"
            )
            time.sleep(delay)

    raise last_error


def search_pdb(protein_sequence):
    """Search RCSB for protein polymer entities and calculate alignment stats."""
    query_sequence = clean_protein_sequence(protein_sequence)
    if len(query_sequence) < 25:
        return []

    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "target": "pdb_protein_sequence",
                "value": query_sequence,
                "identity_cutoff": IDENTITY_CUTOFF,
                "evalue_cutoff": E_VALUE_CUTOFF,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max(MAX_PDB_CANDIDATES, MAX_PDB_RESULTS),
            },
            "results_verbosity": "verbose",
        },
    }

    data = post_rcsb_search(payload)
    if data is None:
        return []

    candidates = []
    seen_entities = set()

    for item in data.get("result_set", []):
        entity_identifier = item.get("identifier", "")
        if "_" not in entity_identifier or entity_identifier in seen_entities:
            continue

        seen_entities.add(entity_identifier)
        candidate = {
            "pdb_entity": entity_identifier,
            "rcsb_search_score": item.get("score"),
            "rcsb_search_rank": len(seen_entities),
        }
        candidate.update(summarize_match_context(item))
        candidates.append(candidate)

    candidates.sort(key=candidate_rank)

    matches = []
    for candidate in candidates[:max(RERANK_SHORTLIST, MAX_PDB_RESULTS)]:
        entity_identifier = candidate["pdb_entity"]
        pdb_id, entity_id = entity_identifier.split("_", 1)

        match = {
            "pdb_id": pdb_id.upper(),
            "entity_id": entity_id,
            "structure_source": SOURCE_EXPERIMENTAL,
            "pdb_entity": entity_identifier,
            "rcsb_search_score": candidate["rcsb_search_score"],
            "rcsb_search_rank": candidate["rcsb_search_rank"],
            "candidates_considered": len(candidates),
            "sequence_identity_pct": None,
            "query_coverage_pct": None,
            "target_coverage_pct": None,
            "alignment_length": None,
            "aligned_residue_pairs": None,
            "query_start": None,
            "query_end": None,
            "target_start": None,
            "target_end": None,
            "pdb_entity_sequence_length": None,
            "alignment_score": None,
            "alignment_method": "",
            "alignment_error": "",
        }

        try:
            target_sequence = get_pdb_entity_sequence(pdb_id, entity_id)
            if not target_sequence:
                raise ValueError("RCSB returned an empty polymer-entity sequence")
            match.update(
                calculate_alignment_metrics(query_sequence, target_sequence)
            )
        except Exception as error:
            match["alignment_error"] = str(error)

        matches.append(match)

    # Rerank on the locally computed metrics, then keep the best.
    matches.sort(key=match_rank)
    return matches[:MAX_PDB_RESULTS]


# =====================================================
# ALPHAFOLD DB FALLBACK
# =====================================================

def find_uniprot_accession(protein_accession, protein_sequence):
    """Resolve a protein to a UniProt accession with an identical sequence.

    The accession is looked up rather than the sequence searched, then the
    candidate's sequence is compared with the query. Requiring an exact match
    avoids attaching a homologue's predicted model to this protein.
    """
    base = str(protein_accession or "").split(".")[0]
    if not base or base.startswith("Predicted_"):
        return None

    try:
        response = SESSION.get(
            UNIPROT_SEARCH_URL,
            params={
                "query": base,
                "fields": "accession,sequence",
                "format": "json",
                "size": 5,
            },
            timeout=HTTP_TIMEOUT,
        )
        response.raise_for_status()
        results = response.json().get("results", [])
    except (requests.RequestException, ValueError):
        return None

    query = clean_protein_sequence(protein_sequence)
    for entry in results:
        candidate = entry.get("sequence", {}).get("value", "")
        if clean_protein_sequence(candidate) == query:
            return entry.get("primaryAccession")

    return None


def mean_ca_plddt(pdb_text):
    """Mean pLDDT over CA atoms, read from the B-factor column."""
    values = [
        float(line[60:66])
        for line in pdb_text.splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    ]
    return sum(values) / len(values) if values else 0.0


def fetch_alphafold_model(uniprot_accession, folder):
    """Download the AlphaFold model for a UniProt accession.

    Returns (path, mean_plddt), or (None, mean_plddt) when the model falls
    below MIN_ALPHAFOLD_PLDDT, or (None, 0.0) when none is available.
    """
    try:
        response = SESSION.get(
            ALPHAFOLD_API_URL.format(accession=uniprot_accession),
            timeout=HTTP_TIMEOUT,
        )
        if response.status_code != 200:
            return None, 0.0
        entries = response.json()
        if not entries:
            return None, 0.0

        url = entries[0].get("pdbUrl")
        if not url:
            return None, 0.0

        model = SESSION.get(url, timeout=HTTP_TIMEOUT)
        model.raise_for_status()
        text = model.text
    except (requests.RequestException, ValueError, KeyError, IndexError):
        return None, 0.0

    confidence = mean_ca_plddt(text)
    if confidence < MIN_ALPHAFOLD_PLDDT:
        return None, confidence

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"AF-{uniprot_accession}.pdb"
    path.write_text(text, encoding="utf-8")
    return str(path), confidence


def search_alphafold(protein):
    """Build a predicted-structure match for a protein with no PDB entry."""
    accession = find_uniprot_accession(
        protein.get("protein_accession"),
        protein.get("protein_sequence", ""),
    )
    if not accession:
        return None

    path, confidence = fetch_alphafold_model(accession, PDB_DIR)
    if not path:
        if confidence:
            print(
                f"AlphaFold model for {accession} rejected: mean pLDDT "
                f"{confidence:.1f} is below {MIN_ALPHAFOLD_PLDDT:.0f}"
            )
        return None

    sequence_length = len(clean_protein_sequence(protein.get("protein_sequence", "")))

    # The model is built from this exact sequence, so it covers it completely.
    return {
        "pdb_id": f"AF-{accession}",
        "entity_id": "1",
        "pdb_entity": f"AF-{accession}_1",
        "structure_source": SOURCE_PREDICTED,
        "uniprot_accession": accession,
        "model_confidence_plddt": round(confidence, 2),
        "rcsb_search_score": None,
        "rcsb_search_rank": None,
        "candidates_considered": None,
        "sequence_identity_pct": 100.0,
        "query_coverage_pct": 100.0,
        "target_coverage_pct": 100.0,
        "alignment_length": sequence_length,
        "aligned_residue_pairs": sequence_length,
        "query_start": 1,
        "query_end": sequence_length,
        "target_start": 1,
        "target_end": sequence_length,
        "pdb_entity_sequence_length": sequence_length,
        "alignment_score": None,
        "alignment_method": "Not aligned; model predicted from this sequence",
        "alignment_error": "",
        "downloaded_file": path,
    }


# =====================================================
# DOWNLOAD STRUCTURE ONCE PER PDB ID
# =====================================================

def download_structure(pdb_id, folder):
    """Download PDB; if unavailable, download mmCIF.

    Existing valid files are reused. Since all files go into the shared
    structures folder, the same PDB ID is not downloaded repeatedly for
    different proteins/isoforms.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.upper()

    pdb_path = folder / f"{pdb_id}.pdb"
    cif_path = folder / f"{pdb_id}.cif"

    if pdb_path.is_file() and pdb_path.stat().st_size > 0:
        return str(pdb_path)
    if cif_path.is_file() and cif_path.stat().st_size > 0:
        return str(cif_path)

    # Try legacy PDB format first.
    try:
        response = SESSION.get(
            RCSB_FILE_URL.format(pdb_id=pdb_id, extension="pdb"),
            timeout=HTTP_TIMEOUT,
        )
        if response.status_code == 200:
            text = response.text
            if re.search(r"^(ATOM  |HETATM|HEADER)", text, re.MULTILINE):
                pdb_path.write_text(text, encoding="utf-8")
                return str(pdb_path)
    except requests.RequestException:
        pass

    # Try modern mmCIF format when PDB format is not available.
    try:
        response = SESSION.get(
            RCSB_FILE_URL.format(pdb_id=pdb_id, extension="cif"),
            timeout=HTTP_TIMEOUT,
        )
        if response.status_code == 200:
            text = response.text
            if text.lstrip().startswith("data_") and "_atom_site." in text:
                cif_path.write_text(text, encoding="utf-8")
                return str(cif_path)
    except requests.RequestException:
        pass

    return None


# =====================================================
# OUTPUT HELPERS
# =====================================================

BASE_OUTPUT_FIELDS = [
    "sequence-id",
    "sequence-name",
    "sequence-risk",
    "sequence-type",
    "sequence-length",
    "organism",
    "gene-name",
    "protein-accession",
    "protein-name",
    "protein-source",
    "cds-location",
    "protein-length",
    "pdb-id",
    "pdb-entity",
    "structure-source",
    "uniprot-accession",
    "model-confidence-plddt",
    "rcsb-search-score",
    "rcsb-search-rank",
    "candidates-considered",
    "sequence-identity-pct",
    "query-coverage-pct",
    "target-coverage-pct",
    "alignment-length",
    "aligned-residue-pairs",
    "query-start",
    "query-end",
    "target-start",
    "target-end",
    "pdb-entity-sequence-length",
    "alignment-score",
    "alignment-method",
    "downloaded-file",
    "error",
    "status",
]


def make_result(row, identifier, protein=None, match=None,
                downloaded_file=None, error="", status=""):
    protein = protein or {}
    match = match or {}

    return {
        "sequence-id": identifier,
        "sequence-name": value_or_blank(row, "sequence-name"),
        "sequence-risk": value_or_blank(row, "sequence-risk"),
        "sequence-type": value_or_blank(row, "sequence-type"),
        "sequence-length": value_or_blank(row, "sequence-length"),
        "organism": protein.get("organism", ""),
        "gene-name": protein.get("gene_name", ""),
        "protein-accession": protein.get("protein_accession", ""),
        "protein-name": protein.get("protein_name", ""),
        "protein-source": protein.get("source", ""),
        "cds-location": protein.get("cds_location", ""),
        "protein-length": len(protein.get("protein_sequence", "")) or None,
        "pdb-id": match.get("pdb_id", ""),
        "pdb-entity": match.get("pdb_entity", ""),
        "structure-source": match.get("structure_source", ""),
        "uniprot-accession": match.get("uniprot_accession", ""),
        "model-confidence-plddt": match.get("model_confidence_plddt"),
        "rcsb-search-score": match.get("rcsb_search_score"),
        "rcsb-search-rank": match.get("rcsb_search_rank"),
        "candidates-considered": match.get("candidates_considered"),
        "sequence-identity-pct": match.get("sequence_identity_pct"),
        "query-coverage-pct": match.get("query_coverage_pct"),
        "target-coverage-pct": match.get("target_coverage_pct"),
        "alignment-length": match.get("alignment_length"),
        "aligned-residue-pairs": match.get("aligned_residue_pairs"),
        "query-start": match.get("query_start"),
        "query-end": match.get("query_end"),
        "target-start": match.get("target_start"),
        "target-end": match.get("target_end"),
        "pdb-entity-sequence-length": match.get("pdb_entity_sequence_length"),
        "alignment-score": match.get("alignment_score"),
        "alignment-method": match.get("alignment_method", ""),
        "downloaded-file": downloaded_file or "",
        "error": error,
        "status": status,
    }


# =====================================================
# REDUNDANT ACCESSION COLLAPSE
# =====================================================

def redundancy_rank(record):
    """Sort key for competing accessions that hit the same structure.

    Ranks on target coverage descending, then prefers a curated RefSeq
    accession over a predicted (XP_/ZP_) model, then the lowest accession so
    the choice is deterministic across runs.
    """
    accession = str(record.get("protein-accession") or "")
    predicted = 1 if re.match(r"^(XP|ZP)_", accession, re.I) else 0

    target_coverage = pd.to_numeric(record.get("target-coverage-pct"), errors="coerce")
    if pd.isna(target_coverage):
        target_coverage = -1.0

    return (-float(target_coverage), predicted, accession)


def collapse_redundant_entities(output):
    """Drop redundant accessions that resolve to an already-covered structure.

    A row qualifies when identity and query coverage are both 100% and target
    coverage exceeds MIN_TARGET_COVERAGE_PCT. Within each
    (sequence-id, pdb-id) group that contains at least one qualifying row, the
    best qualifying row is kept and the remainder of that group is dropped.

    Grouping on pdb-id as well as sequence-id is deliberate. For polyprotein
    inputs the additional accessions are distinct mature proteins that map to
    different structures rather than redundant isoforms, so they survive.
    Groups with no qualifying row are left untouched.
    """
    if output.empty:
        return output

    metrics = {
        column: pd.to_numeric(output[column], errors="coerce")
        for column in (
            "sequence-identity-pct",
            "query-coverage-pct",
            "target-coverage-pct",
        )
    }

    qualifies = (
        (metrics["sequence-identity-pct"] >= PERFECT_IDENTITY_PCT)
        & (metrics["query-coverage-pct"] >= PERFECT_QUERY_COVERAGE_PCT)
        & (metrics["target-coverage-pct"] > MIN_TARGET_COVERAGE_PCT)
        & (output["status"] == "Downloaded")
    )

    mapped = output[
        (output["status"] == "Downloaded")
        & output["pdb-id"].notna()
        & (output["pdb-id"].astype(str) != "")
    ]

    dropped = []
    for _, group in mapped.groupby(["sequence-id", "pdb-id"], dropna=False):
        if len(group) < 2:
            continue

        winners = list(group.index[qualifies.loc[group.index]])
        if not winners:
            continue

        best = min(winners, key=lambda index: redundancy_rank(output.loc[index]))
        dropped.extend(index for index in group.index if index != best)

    if dropped:
        print(f"\nCollapsed {len(dropped)} redundant accession row(s):")
        for index in sorted(dropped):
            record = output.loc[index]
            print(
                f"  dropped {record['sequence-id']} {record['protein-accession']} "
                f"-> {record['pdb-id']} "
                f"(identity {record['sequence-identity-pct']}, "
                f"query cov {record['query-coverage-pct']}, "
                f"target cov {record['target-coverage-pct']})"
            )

    return output.drop(index=dropped)


def drop_unmapped_rows(output):
    """Drop structure-less rows for inputs that mapped elsewhere.

    A protein that produced no usable structure leaves a row with an empty
    pdb-id. Where another protein from the same sequence-id did map, that row
    carries no information and is dropped.

    A structure-less row is kept when it is the only row for its sequence-id,
    since dropping it would remove the input from the mapping altogether and
    silently shrink the record of what was attempted.
    """
    if output.empty:
        return output

    has_structure = (
        output["pdb-id"].notna()
        & (output["pdb-id"].astype(str).str.strip() != "")
    )
    mapped_ids = set(output.loc[has_structure, "sequence-id"])

    removable = ~has_structure & output["sequence-id"].isin(mapped_ids)

    if removable.any():
        print(f"\nDropped {int(removable.sum())} structure-less row(s) "
              f"superseded by a mapped protein from the same input:")
        for _, record in output[removable].iterrows():
            print(f"  dropped {record['sequence-id']} "
                  f"{record['protein-accession']} ({record['status']})")

    return output[~removable]


# =====================================================
# MAIN PIPELINE
# =====================================================

def process_row(row, identifier=None):
    """Resolve a single input row into output records.

    Shared by main() and by targeted reruns of previously failed rows. The
    identifier is reported verbatim from the input so the output CSV still
    joins back to the input file, while NCBI lookups use the normalized form.
    """
    if identifier is None:
        identifier = str(row["sequence-id"]).strip()

    results = []

    try:
        proteins = extract_proteins(row)
        print(f"Extracted {len(proteins)} protein(s)")
    except Exception as error:
        print(f"PROTEIN EXTRACTION FAILED: {error}")
        results.append(make_result(
            row, identifier,
            error=str(error),
            status="Extraction failed",
        ))
        time.sleep(REQUEST_DELAY)
        return results

    for protein in proteins:
        protein_sequence = protein.get("protein_sequence", "")
        validation_errors = validate_protein(protein_sequence)

        if validation_errors:
            error_text = "; ".join(validation_errors)
            results.append(make_result(
                row, identifier,
                protein=protein,
                error=error_text,
                status="Protein validation failed",
            ))
            continue

        try:
            matches = search_pdb(protein_sequence)
        except Exception as error:
            print(f"PDB SEARCH FAILED for {protein.get('protein_accession')}: {error}")
            results.append(make_result(
                row, identifier,
                protein=protein,
                error=str(error),
                status="PDB search failed",
            ))
            continue

        if not matches:
            # No experimental structure. A predicted model still gives
            # ProteinMPNN a backbone, provided it is confident enough.
            predicted = search_alphafold(protein) if ALPHAFOLD_FALLBACK else None

            if predicted:
                print(f"Using AlphaFold model {predicted['pdb_id']} "
                      f"(mean pLDDT {predicted['model_confidence_plddt']})")
                results.append(make_result(
                    row, identifier,
                    protein=protein,
                    match=predicted,
                    downloaded_file=predicted["downloaded_file"],
                    error="",
                    status="Downloaded",
                ))
            else:
                results.append(make_result(
                    row, identifier,
                    protein=protein,
                    error="No matching PDB polymer entity found; "
                          "no usable AlphaFold model either",
                    status="No PDB match",
                ))
            continue

        print(f"Found {len(matches)} PDB match(es) for {protein.get('protein_accession')}")

        for match in matches:
            pdb_id = match["pdb_id"]

            # IMPORTANT: shared PDB_DIR, not a protein-specific subfolder.
            # Existing PDB/CIF files are reused, so each PDB ID is saved once.
            filepath = download_structure(
                pdb_id,
                PDB_DIR,
            )

            errors = []
            if match.get("alignment_error"):
                errors.append(
                    "Alignment metrics unavailable: " + match["alignment_error"]
                )

            if filepath:
                status = "Downloaded"
                print(f"Downloaded/reused {pdb_id}: {filepath}")
            else:
                status = "Download failed"
                errors.append("PDB and mmCIF downloads failed")
                print(f"FAILED to download {pdb_id}")

            results.append(make_result(
                row,
                identifier,
                protein=protein,
                match=match,
                downloaded_file=filepath,
                error="; ".join(errors),
                status=status,
            ))

            time.sleep(REQUEST_DELAY)

    time.sleep(REQUEST_DELAY)
    return results


def main():
    if Entrez.email == "YOUR_EMAIL@example.com":
        raise ValueError(
            "Please set Entrez.email in the script or define NCBI_EMAIL. "
            "NCBI requests require a valid contact email."
        )

    # Detect comma- or tab-separated CSV/TSV.
    df = pd.read_csv(INPUT_FILE, sep=None, engine="python", dtype=str)
    df.columns = [str(column).strip().lower().rstrip(":") for column in df.columns]

    for required in ("sequence-id", "sequence"):
        if required not in df.columns:
            raise ValueError(f"Missing required input column: {required}")

    results = []
    input_has_download = {index: False for index in df.index}
    input_download_failures = {index: 0 for index in df.index}

    total = len(df)
    print(f"\nProcessing {total} input rows.\n")

    for index, row in df.iterrows():
        identifier = str(row["sequence-id"]).strip()
        print(f"\n[{index + 1}/{total}] Processing {identifier}")

        row_results = process_row(row, identifier)
        results.extend(row_results)

        statuses = [record["status"] for record in row_results]
        input_has_download[index] = "Downloaded" in statuses
        input_download_failures[index] = statuses.count("Download failed")

    # Save consistent columns even if there are no results.
    output = pd.DataFrame(results, columns=BASE_OUTPUT_FIELDS)
    output = collapse_redundant_entities(output)
    output = drop_unmapped_rows(output)
    output.to_csv(OUTPUT_CSV, index=False)

    successful_inputs = sum(input_has_download.values())
    failed_inputs = total - successful_inputs
    total_download_failures = sum(input_download_failures.values())

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED")
    print("=" * 60)
    print(f"Input rows processed: {total}")
    print(f"Input rows with at least one downloaded structure: {successful_inputs}")
    print(f"Input rows without any downloaded structure: {failed_inputs}")
    print(f"PDB candidate download failures: {total_download_failures}")
    if total:
        print(f"Input-row failure percentage: {failed_inputs / total * 100:.2f}%")
    print(f"\nMapping CSV: {OUTPUT_CSV}")
    print(f"Unique structure files: {PDB_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
