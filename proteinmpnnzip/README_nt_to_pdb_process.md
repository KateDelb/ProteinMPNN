# NCBI sequence to PDB mapping pipeline

This project takes a CSV of NCBI identifiers and nucleotide sequences, retrieves the associated protein sequence(s), searches the RCSB Protein Data Bank (PDB) for similar protein structures, downloads the best matching structure file, and records the results in a CSV.

The pipeline is intended for small mixed-organism datasets (for example, bacterial, viral, and eukaryotic records). It uses NCBI annotations where available rather than assuming that every supplied nucleotide sequence is a single, uninterrupted coding sequence.

> **Important scientific limitation:** the current script does not yet align the nucleotide sequence in your CSV to individual CDS features in the NCBI record. For a nucleotide accession, it extracts **all translated CDS features in that GenBank record**. This is appropriate for a gene-cluster accession if you want every encoded protein, but it can return many unrelated proteins when the accession represents a whole genome. Review the `gene-name`, `protein-name`, and `cds-location` columns before interpreting the hits. See [What the script extracts](#what-the-script-extracts-and-important-limitations) below.

## Contents

- [Requirements](#requirements)
- [Input CSV](#input-csv)
- [Installation](#installation)
- [Configuration](#configuration)
- [Run the pipeline](#run-the-pipeline)
- [How the pipeline works](#how-the-pipeline-works)
- [Output files](#output-files)
- [Understanding the alignment metrics](#understanding-the-alignment-metrics)
- [Statuses and error reporting](#statuses-and-error-reporting)
- [Important scientific considerations](#important-scientific-considerations)
- [Troubleshooting](#troubleshooting)
- [Useful documentation](#useful-documentation)

## Requirements

- Python 3.9 or later is recommended.
- An internet connection is required to contact NCBI and RCSB PDB.
- A valid email address is required by NCBI Entrez. The script reads it from the `NCBI_EMAIL` environment variable, or it can be set directly in the script.

Python packages:

- `biopython`
- `pandas`
- `requests`

## Input CSV

By default, the script looks for a file named:

```text
input_sequences.csv
```

The script automatically detects comma-separated or tab-separated input. At minimum, the file must contain these two columns:

| Column | Required? | Description |
|---|---|---|
| `sequence-id` | Yes | NCBI Gene ID, nucleotide accession, or supported protein accession. |
| `sequence` | Yes | Nucleotide sequence using IUPAC DNA/RNA characters. |
| `sequence-name` | No | Short description of the gene or sequence. |
| `sequence-risk` | No | User-provided classification, such as `benign` or `hazardous`. Preserved in the output; not used by the script to make biological judgments. |
| `sequence-type` | No | User-provided classification, such as `known` or `unknown`. Preserved in the output. |
| `sequence-length` | No | User-provided nucleotide length. Preserved in the output; the script does not use it to choose a CDS. |

Example:

```csv
sequence-id,sequence-name,sequence-risk,sequence-type,sequence-length,sequence
51738,Homo sapiens Ghrelin (GHRL peptide),benign,known,354,ATGGCT...
AF183379.1,Microcystin gene cluster,benign,known,12345,ATG...
NC_014395.1,Example genome,benign,known,2500000,ATG...
```

The nucleotide strings above are placeholders; replace them with the complete sequences in your file. Column names are normalized to lowercase and a trailing colon is removed, so headers such as `sequence-id:` are accepted. Column names should otherwise match the names above.

### Identifier types handled

The script identifies records using the `sequence-id` value:

| Identifier example | How it is handled |
|---|---|
| `51738` | Treated as an NCBI Gene ID. Retrieves proteins linked to that Gene record. |
| `NP_415477.1` | Recognized as a common RefSeq protein accession and retrieves that protein record. The script recognizes accessions beginning with `NP_`, `XP_`, `YP_`, `WP_`, `AP_`, or `ZP_`. |
| `AF183379.1`, `M16593.1`, `NC_014395.1` | Treated as a nucleotide accession. Retrieves the GenBank record and extracts its annotated translated CDS features. |

For Gene IDs, the script uses the NCBI Gene-to-protein links; the nucleotide sequence in the CSV is checked for valid nucleotide characters but is not used to choose among the linked protein isoforms. For nucleotide accessions, the script fetches the accession's GenBank record from NCBI; it does not use a local alignment to select just the CDS matching the supplied CSV sequence.

## Installation

Open a terminal or command prompt in the project directory and install the required packages:

```bash
python -m pip install biopython pandas requests
```

If your system uses `python3` instead of `python`, use:

```bash
python3 -m pip install biopython pandas requests
```

## Configuration

Place these files in the same project directory:

```text
project/
├── input_sequences.csv
└── ncbi_to_pdb_updated.py
```

### Set your NCBI email address

The recommended method is to set an environment variable.

**macOS / Linux:**

```bash
export NCBI_EMAIL="you@example.com"
```

**Windows PowerShell:**

```powershell
$env:NCBI_EMAIL = "you@example.com"
```

Alternatively, edit the script and replace `YOUR_EMAIL@example.com` with your email address. The script stops with an error if the placeholder has not been replaced and `NCBI_EMAIL` is not set.

An NCBI API key is optional. If you have one, set `NCBI_API_KEY` in your environment:

**macOS / Linux:**

```bash
export NCBI_API_KEY="your_api_key"
```

**Windows PowerShell:**

```powershell
$env:NCBI_API_KEY = "your_api_key"
```

### Check the settings in the script

Near the top of `ncbi_to_pdb_updated.py`:

```python
INPUT_FILE = "input_sequences.csv"

IDENTITY_CUTOFF = 0.30
MAX_PDB_RESULTS = 1
E_VALUE_CUTOFF = 0.1
```

- `INPUT_FILE`: input CSV/TSV path, relative to the directory from which the script is run, or an absolute path.
- `IDENTITY_CUTOFF`: minimum sequence identity threshold used by the RCSB sequence search (`0.30` means 30%).
- `MAX_PDB_RESULTS`: maximum number of matching PDB polymer entities returned for each protein sequence. `1` means up to one PDB entity per extracted protein, **not** one protein per input Gene ID.
- `E_VALUE_CUTOFF`: E-value threshold passed to the RCSB sequence search.

For a small first run, keeping `MAX_PDB_RESULTS = 1` is reasonable. Increasing it can return more candidate structures per protein.

## Run the pipeline

From the project directory, run:

```bash
python ncbi_to_pdb_updated.py
```

On systems that use `python3`:

```bash
python3 ncbi_to_pdb_updated.py
```

The script prints progress as it processes each input row. It contacts NCBI to resolve the identifier and obtain protein sequences, searches RCSB for PDB polymer entities, retrieves the sequence of each matched PDB entity, calculates local alignment metrics, and downloads the corresponding structure file.

## How the pipeline works

### 1. Validate and classify the input identifier

The script checks that the supplied nucleotide sequence contains recognized nucleotide characters. It then classifies the `sequence-id`:

- Numeric ID: NCBI Gene ID.
- Common RefSeq protein prefix: protein accession.
- Other accession: treated as a nucleotide accession.

### 2. Retrieve protein sequence(s) from NCBI

The route depends on the identifier:

- **Gene ID:** NCBI Gene-to-protein links are used to retrieve linked protein records. A gene may have multiple RefSeq protein products, for example because it has multiple transcript isoforms.
- **Protein accession:** the protein FASTA record is retrieved directly.
- **Nucleotide accession:** the GenBank record is retrieved and every CDS feature with an annotated `/translation` qualifier is extracted. The record may therefore produce one protein, several proteins, or many proteins.
- **Fallback translation:** if no annotated CDS is found for a non-`NC_` nucleotide accession, the script attempts a simple translation of the CSV sequence in reading frame 1 using the standard genetic code. This is explicitly marked as an unverified fallback. It does not predict intron/exon structure or search all six reading frames.

### 3. Basic protein validation

The script checks that the extracted protein sequence is at least 25 amino acids long, because the RCSB sequence search requires a sequence of that minimum length. Ambiguous or uncommon amino-acid symbols are normalized to `X` for alignment. This is a basic screening step, not a full validation that a nucleotide sequence is a biologically correct CDS.

### 4. Search RCSB PDB by protein sequence

The protein sequence is submitted to the RCSB Search API using its protein sequence-similarity search. The search returns PDB polymer entities, for example `7F9Y_3`:

- `7F9Y` is the PDB entry ID.
- `3` is the polymer entity number within that entry.

A polymer entity is not necessarily the same thing as a chain ID. One polymer entity may be represented by one or more chains in the deposited structure.

### 5. Retrieve the PDB entity sequence and calculate alignment metrics

For each returned PDB entity, the script retrieves its canonical amino-acid sequence from the RCSB Data API. It then locally aligns the NCBI protein sequence to that PDB entity sequence with Biopython `PairwiseAligner`, using local alignment and the BLOSUM62 substitution matrix.

The script records sequence identity, query coverage, target coverage, alignment coordinates, alignment length, aligned residue-pair count, and alignment score. See [Understanding the alignment metrics](#understanding-the-alignment-metrics).

### 6. Download the structure once

The script tries to download the legacy `.pdb` file first. If that format is unavailable, it tries the modern `.cif` (mmCIF) format.

All structures are saved in one shared directory:

```text
pdb_results/
└── structures/
    ├── 7F9Y.pdb
    ├── 2JSH.pdb
    └── ...
```

If the same PDB ID is matched by several proteins or isoforms, the existing file is reused instead of downloading another copy. The complete deposited structure is downloaded; the script does not extract only the matching chain or polymer entity.

## Output files

After the run, the project directory will look like this:

```text
project/
├── input_sequences.csv
├── ncbi_to_pdb_updated.py
└── pdb_results/
    ├── pdb_mapping.csv
    └── structures/
        ├── 7F9Y.pdb
        ├── 2JSH.pdb
        └── ...
```

### `pdb_results/pdb_mapping.csv`

The output has one row per extracted protein/PDB entity match. A single input row can therefore produce multiple output rows if the Gene ID has multiple linked protein isoforms, or if a nucleotide record contains multiple CDS features.

The output columns are:

| Output column | Meaning |
|---|---|
| `sequence-id` | Original input identifier. |
| `sequence-name` | Original descriptive name, if present. |
| `sequence-risk` | Original user-provided risk label, if present. It is carried through as metadata only. |
| `sequence-type` | Original user-provided type label, if present. |
| `sequence-length` | Original sequence length value from the input file, if present. |
| `organism` | Organism annotation returned with the protein, when available. This may be blank for some retrieval routes. |
| `gene-name` | Gene qualifier from the GenBank CDS feature, when available. |
| `protein-accession` | NCBI protein accession, or a generated `Predicted_...` label for fallback translations. |
| `protein-name` | NCBI protein description or GenBank CDS product. |
| `protein-source` | Where the protein sequence came from (Gene link, protein record, annotated CDS, or fallback translation). |
| `cds-location` | CDS feature location in the GenBank record, when applicable. Joined locations can indicate a CDS assembled from multiple nucleotide segments. |
| `protein-length` | Length of the extracted amino-acid sequence. |
| `pdb-id` | PDB entry identifier, for example `7F9Y`. |
| `pdb-entity` | RCSB polymer entity identifier, for example `7F9Y_3`. This is not a chain ID. |
| `rcsb-search-score` | Score returned by the RCSB search service. It is not a percentage identity. |
| `sequence-identity-pct` | Percentage of identical residues in the local alignment, with alignment gap columns included in the denominator. |
| `query-coverage-pct` | Percentage of the full NCBI protein sequence covered by the local alignment. The NCBI protein is the query. |
| `target-coverage-pct` | Percentage of the full PDB polymer-entity sequence covered by the local alignment. The PDB entity is the target. |
| `alignment-length` | Number of columns in the local alignment, including columns containing gaps. |
| `aligned-residue-pairs` | Number of alignment columns in which both sequences have a residue (neither position is a gap). |
| `query-start`, `query-end` | 1-based inclusive start and end positions in the NCBI protein sequence covered by the alignment. |
| `target-start`, `target-end` | 1-based inclusive start and end positions in the PDB entity sequence covered by the alignment. |
| `pdb-entity-sequence-length` | Full amino-acid length of the PDB polymer entity sequence. |
| `alignment-score` | Local BLOSUM62 alignment score. This is a score, not a percentage. |
| `alignment-method` | Method used for the locally calculated metrics. |
| `downloaded-file` | Relative path to the downloaded `.pdb` or `.cif` file. Blank if download did not succeed. |
| `error` | Explanation of extraction, validation, alignment, or download problems. Empty when no error was recorded. |
| `status` | Processing outcome for that output row. |

## Understanding the alignment metrics

The two sequences being compared are:

- **Query:** the protein sequence obtained from NCBI.
- **Target:** the canonical protein sequence for the matched PDB polymer entity.

The three most useful percentages answer different questions:

| Metric | Question | Example interpretation |
|---|---|---|
| `sequence-identity-pct` | Among the aligned columns, how many contain the same amino acid? | 95% means the aligned region is highly similar. |
| `query-coverage-pct` | How much of my complete NCBI protein is represented in the alignment? | 25% means only about one quarter of the NCBI protein aligns to the PDB entity. |
| `target-coverage-pct` | How much of the complete PDB entity sequence is represented in the alignment? | 100% means the entire PDB entity sequence is covered by the alignment. |

For example, a 117-amino-acid precursor that aligns over 28 amino acids to a 28-amino-acid PDB peptide could have approximately:

- 100% sequence identity;
- 23.9% query coverage;
- 100% target coverage.

That is an exact match to a short region, not a structure of the complete 117-amino-acid precursor. This can happen when a mature peptide is processed from a longer precursor.

### How to use the metrics when reviewing hits

As an initial screening guide (not universal biological cutoffs):

- **More convincing full-protein representation:** high sequence identity and high query coverage, with substantial target coverage.
- **Potential domain or mature-peptide match:** high identity but low query coverage.
- **Match to a region of a larger protein:** high identity and high query coverage, but low target coverage.
- **More distant candidate:** lower identity, even if coverage is high.

A useful exploratory starting point is at least 70% identity and 80% query coverage for a relatively close, broad match. Lower thresholds can still be useful for homolog searches, but need more careful biological review. Coverage and identity do not by themselves prove that the PDB structure has the same function, activity, processing state, or biological context as the query protein.

## Statuses and error reporting

Common `status` values include:

| Status | Meaning |
|---|---|
| `Downloaded` | The PDB or mmCIF file was downloaded successfully, or an existing file was reused. |
| `Extraction failed` | The identifier could not be resolved to a protein sequence, or NCBI retrieval failed. |
| `Protein validation failed` | The extracted protein was empty or shorter than the minimum sequence-search length. |
| `No PDB match` | RCSB returned no matching PDB polymer entity at the configured search thresholds. |
| `PDB search failed` | The RCSB request or response processing failed. |
| `Download failed` | A PDB and mmCIF file could not be downloaded for a returned PDB ID. |

The script writes errors to the `error` column and prints a summary at the end:

- Number of input rows processed.
- Number of input rows with at least one successfully downloaded/reused structure.
- Number of input rows without any successfully downloaded structure.
- Number of candidate structure downloads that failed.
- Percentage of input rows without a downloaded structure.

An input row is counted as successful if at least one structure was downloaded or reused for any protein extracted from that row. An input row can still have some failed candidate downloads and be counted as successful if another structure for that row was retrieved.

**No PDB match and download failure are different outcomes.** `No PDB match` means the search did not return a candidate; `Download failed` means a candidate was returned but neither file format could be retrieved.

## What the script extracts and important limitations

Read this section carefully when interpreting results from mixed bacterial, viral, and eukaryotic inputs.

### A Gene ID can yield multiple protein records

One Gene ID can link to multiple RefSeq protein accessions. These may be alternative transcript products, isoforms, or other annotated products. The number appended to a folder or row in earlier versions of the script was only a processing index, not a biological isoform number.

### A nucleotide accession can yield multiple CDS proteins

For nucleotide accessions, the script extracts every GenBank CDS feature that has a `/translation` annotation. This is expected for a gene-cluster record containing several genes. It is also expected for a complete genome, which may contain many genes. The current script does not determine which CDS corresponds to the descriptive label in `sequence-name`, nor does it compare each CDS nucleotide segment against the `sequence` column in the CSV.

### Introns and exons

A GenBank CDS location containing joined segments (for example, `join(...)`) indicates that the annotated CDS is assembled from multiple segments of the record. In eukaryotic genomic records, these often correspond to coding exons separated by introns. The script records the feature location in `cds-location`, but it does not independently classify the biological reason for a joined location.

### Fallback translations are provisional

If NCBI provides no annotated CDS for a non-`NC_` nucleotide accession, the script may translate the supplied sequence directly in frame 1 using the standard genetic code. That fallback assumes the input starts at the first base of a CDS and contains no introns. It is not a gene predictor and should be verified before using the resulting protein sequence for biological conclusions.

### The downloaded structure may contain more than the matched protein

The script downloads the full PDB/mmCIF entry. That file may contain multiple protein chains, nucleic acids, ligands, cofactors, or binding partners. `pdb-entity` identifies the matched polymer entity in the entry; it does not mean that the file has been reduced to that entity or chain. The script does not currently extract an individual chain.

### PDB does not guarantee an experimental structure of the full query protein

A high-identity hit may cover only a short peptide or one domain. Always interpret sequence identity together with query coverage and target coverage. A PDB hit from another species may be a useful homologous structure, but it is not automatically an experimentally determined structure of the exact query protein.

### PDB versus mmCIF

The script tries `.pdb` first and falls back to `.cif` (mmCIF) when legacy PDB format is unavailable. Some large structures are only available in mmCIF. The `downloaded-file` column records which file was saved.

## Troubleshooting

**Error: `Please set Entrez.email...`**

Set the `NCBI_EMAIL` environment variable or edit `Entrez.email` in the script.

**Error: `Missing required input column`**

Check that the input header includes `sequence-id` and `sequence`. Other metadata columns are optional.

**A single Gene ID produces many proteins.**

This can happen when the gene has multiple protein products. For a nucleotide accession, it can also happen because the record contains multiple CDS features. Inspect the protein accessions, protein descriptions, and CDS locations in the output.

**A genome accession produces hundreds or thousands of proteins.**

That is expected with the current all-CDS extraction behavior. If you only want the gene represented by the CSV sequence, the next improvement is to match the supplied nucleotide sequence to the relevant CDS feature(s) before translation, rather than extracting every CDS from the record.

**A sequence has no PDB match.**

This means no polymer entity met the configured RCSB search criteria. It does not prove that no related structure exists under any search condition. Consider reviewing the sequence, lowering the identity cutoff cautiously, or running a separate BLAST/sequence search.

**A structure downloads as `.cif`, not `.pdb`.**

This is expected for entries not available in legacy PDB format. Use an mmCIF-compatible molecular visualization or analysis tool, or convert it if your downstream software requires PDB format (be aware that PDB format has size and format limitations).

**The same PDB ID appears on multiple output rows.**

Different proteins or isoforms can match the same PDB entry. The structure file is stored only once in `pdb_results/structures/`; several CSV rows may point to that same file.

## Useful documentation

- [NCBI Entrez Programming Utilities (E-utilities)](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [NCBI GenBank format and feature table documentation](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
- [RCSB PDB sequence similarity search](https://www.rcsb.org/docs/search-and-browse/advanced-search/sequence-similarity-search)
- [RCSB Data API](https://data.rcsb.org/)
- [RCSB structure file download services](https://www.rcsb.org/docs/programmatic-access/file-download-services)
- [Biopython documentation](https://biopython.org/wiki/Documentation)
