-> Python wrapper that uses several tools (trim_galore, megahit, etc.) to assemble metagenomes and obtain MAGs

-> Author: Massimo Bourquin, 2024

# 1. Installation & Usage

  a. Clone the Github repo `git clone https://github.com/Mass23/MACE_METAG_hybrid`, place the data folder (reads) somewhere and write down the folder.
  
  b. Create the conda environment for metagenome assembly `conda env create -f MACE_METAG_assemble.yml \ conda activate MACE_METAG_assemble`
  
  c. Run example (in the example file) `python3 process_16S_nanopore.py -f alpine_soil/ -n alpine_soil -m metadata.tsv -t 24` from within the repository, -f is the data folder (reads), -n is the name for the output, metadata.tsv is the metadata file (look at the one in the repo for guidance, needs the Barcode and Well columns), -t is the numbers of threads to use.

# 2. Description

  a. Porechop: removes adapters
  
  b. Chopper: quality control (q > 9) and size (only keeps reads from 1300bp to 1800bp)
  
  c. Qiime2: dereplicates sequences, creates abundance table, and taxonokmy based on greengenes2

# 3. Results/ folder
to do...
