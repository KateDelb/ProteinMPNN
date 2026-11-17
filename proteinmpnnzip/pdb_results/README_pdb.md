Here is how each one is calculated in ncbi_id_to_pdb.py:350-390. The script aligns the NCBI protein (the query) against the PDB entity's sequence (the target) with a local BLOSUM62 alignment.

sequence-identity-pct

Identical aligned residues ÷ total alignment columns × 100.
Gap columns count in the denominator, so a gap lowers identity.
Example: the bovine insulin row (6ZI8) has identity 97.14% and alignment length 105. That's 102 identical columns out of 105, so 3 columns differ.
query-coverage-pct

Number of distinct query (NCBI protein) residues that appear in the alignment ÷ the full query length × 100.
Measures how much of your protein the structure covers.
Example: 9A2M covers 100% of the 346-residue ompA protein.
target-coverage-pct

Number of distinct PDB-entity residues in the alignment ÷ the full PDB entity sequence length × 100.
Measures how much of the structure's chain the match covers.
It's usually lower than query coverage when the structure is a fragment or has extra residues, such as tags or disordered regions that aren't in the alignment.