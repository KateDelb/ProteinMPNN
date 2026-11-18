# ProteinMPNN
![ProteinMPNN](https://docs.google.com/drawings/d/e/2PACX-1vTtnMBDOq8TpHIctUfGN8Vl32x5ISNcPKlxjcQJF2q70PlaH2uFlj2Ac4s3khnZqG1YxppdMr0iTyk-/pub?w=889&h=358)
Read [ProteinMPNN paper](https://www.biorxiv.org/content/10.1101/2022.06.03.494563v1).

To run ProteinMPNN clone this github repo and install Python>=3.0, PyTorch, Numpy. 

Full protein backbone models: `vanilla_model_weights/v_48_002.pt, v_48_010.pt, v_48_020.pt, v_48_030.pt`, `soluble_model_weights/v_48_010.pt, v_48_020.pt`.

CA only models: `ca_model_weights/v_48_002.pt, v_48_010.pt, v_48_020.pt`. Enable flag `--ca_only` to use these models.

Helper scripts: `helper_scripts` - helper functions to parse PDBs, assign which chains to design, which residues to fix, adding AA bias, tying residues etc.

Code organization:
* `protein_mpnn_run.py` - the main script to initialialize and run the model.
* `protein_mpnn_utils.py` - utility functions for the main script.
* `examples/` - simple code examples.
* `inputs/` - input PDB files for examples
* `outputs/` - outputs from examples
* `colab_notebooks/` - Google Colab examples
* `training/` - code and data to retrain the model
-----------------------------------------------------------------------------------------------------

### Notes on plotting addition:

To run the plots:
``` 
python plot_folding_scores.py 
```

## Creation of unknown protein sequence creation with ProteinMPNN

A pdb equivalent for all sequences was found and stored under pdb_results/ (see README_nt_to_pdb_process) (note: some multi-chain or deprecated (often hazardous) sequences did not have an exact multi-chain equivalent, or had multiple equivalents, we filtered on high sequence identity and high coverage of the original sequence), some sequences have multiple backbones, resulting in a total of 88 backbones.

We designed 5 variants to all these sequences (and all 440 designs can be found in data/proteinmpnn_designs_all_full.csv). For each structure we chose the top 3 based on MPNN-score to limit database size. The resulting sequences can be found in proteinmpnn_designs_all (264 sequence). 

These resultign sequences were folded with ESMFold and then filtered for pLDDT > 80, dpLDDT > -10 and TM-score > 0.5. 
The included sequences are flagged in the included column as True (138/264 included).