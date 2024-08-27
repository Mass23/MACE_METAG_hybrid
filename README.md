-> Python wrapper that uses several tools (trim_galore, megahit, etc.) to assemble metagenomes and obtain MAGs

-> Author: Massimo Bourquin, 2024

# 1. Metagenome assembly

  a. Clone the Github repo `git clone https://github.com/Mass23/MACE_METAG_hybrid`, place the data folder (reads) somewhere and write down the folder path.
  
  b. Create the conda environment for metagenome assembly `conda env create -f MACE_METAG_assemble.yml \ conda activate MACE_METAG_assemble`
  
  c. Run example (in the example file) `python3 assemble_metagenome.py -f path_to_reads/ -n run_name -m metadata.tsv -t 24` from within the repository, -f is the data folder (reads), -n is the name for the output folder, metadata.tsv is the metadata file (look at the one in the repo for guidance, the Subset column gives the information on what samples should be assembled together), -t is the numbers of threads to use.

# 2. Bins creation and analysis

  a. ...

  b. lll

# 3. Results/ folder
to do...
