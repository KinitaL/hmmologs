# HMMologs

## Installation
1. Clone the repository:
```shell
git clone https://github.com/KinitaL/hmmologs.git
```
2. Create the conda env:
```shell
conda env create --file=env.yaml
```
3. Activate the env:
```shell
conda activate hmmologs
```

## Configuration
Fill the config.yaml file:
- domain -> interpro_id: InterPro ID of the domain homologs of which you want to find
- output_dir: directory for the results
- pdb_filter -> resolution_threshold: PDB structures with the resolution lower that this threshold will be discarded 
- pdb_filter -> max_residues_count: maximum allowed number of residues in the chain
- max_amount_of_proteins: limit in the query to Uniprot

## Execution
```shell
python hmmologs.py 
```