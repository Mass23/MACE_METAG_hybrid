import os
import argparse
import subprocess
import pandas as pd # type: ignore
import datetime
import glob
from Bio import SeqIO

################################################################################
#################           FUNCTIONS           ################################
################################################################################

def create_binning_folder(binning_folder_name):
    if not os.path.exists(binning_folder_name):
        os.makedirs(binning_folder_name)
        os.makedirs(f'{binning_folder_name}/coverm')
        os.makedirs(f'{binning_folder_name}/concoct')
        os.makedirs(f'{binning_folder_name}/metabat2')
    
    with open(f'{binning_folder_name}/log.txt', 'w') as log:
        log.write(f"LStarting binning, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

def rename_filter_contigs(results_folder_name, assembly_path):
    assembly = f'{results_folder_name}/{assembly_path}'
    assembly_filtered = assembly.replace('.fa','.filtered.fa')
    assembly_renamed = assembly.replace('.fa','.filtered.renamed.fa')

    args1 = f"seqkit seq -j 10 --remove-gaps -o {assembly_filtered} -m 1500 {assembly}"
    args2 = f"perl -pe s/ .*//g {assembly_filtered} > {assembly_renamed}"

    subprocess.call(args1, shell = True)
    subprocess.call(args2, shell = True)

# Coverage computation
def run_coverm(results_folder_name, assembly_path, threads):
    """
    Run coverm on the assembly and trimmed reads to get coverage information for 
    1) metabat and 2) mean coverage for binning and contig abundance information.
    """

    assembly = assembly_path.replace('.fa', '.filtered.renamed.fa')
    reads_1 = f'{results_folder_name}/trimmed_reads/*val_1.fq.gz'
    reads_2 = f'{results_folder_name}/trimmed_reads/*val_2.fq.gz'
    coverm_output_path = f'{results_folder_name}/binning/coverm'

    args1 = f'coverm contig -m metabat -r {assembly} -t {threads}  -1 {reads_1} -2 {reads_2} -o {coverm_output_path}/coverm_metabat.tsv'
    args2 = f'coverm contig -m mean -r {assembly} -t {threads}  -1 {reads_1} -2 {reads_2} -o {coverm_output_path}/coverm_mean.tsv'

    subprocess.call(args1, shell = True)
    subprocess.call(args2, shell = True)

# Binning softwares
def run_concoct(results_folder_name, assembly_path, threads):
    """ Function to run concoct. """
    assembly = assembly_path.replace('.fa', '.filtered.renamed.fa')
    coverm_output_path = f'{results_folder_name}/binning/coverm'

    args = f'concoct --coverage_file {coverm_output_path}/coverm_mean.tsv --composition_file {assembly} -t {threads} -c 800 -b {results_folder_name}/binning/concoct'
    subprocess.call(args, shell = True)

def run_rosella(results_folder_name, assembly_path, threads):
    """ Function to run rosella. """
    assembly = assembly_path.replace('.fa', '.filtered.renamed.fa')
    coverm_output_path = f'{results_folder_name}/binning/coverm'

    args = f'rosella recover --coverage-file {coverm_output_path}/coverm_metabat.tsv --assembly {assembly} --output-directory {results_folder_name}/binning/rosella_out/ --threads {threads}'
    subprocess.call(args, shell = True)

def run_metabat2(results_folder_name, assembly_path, threads):
    """ Function to run metabat2. """
    assembly = assembly_path.replace('.fa', '.filtered.renamed.fa')
    coverm_output_path = f'{results_folder_name}/binning/coverm'

    args1 = f'gzip -c {assembly} > {assembly}.gz' # metabat2 requires gzipped fasta file for the assembly
    args2 = f'metabat2 -i {assembly}.gz -l -o {results_folder_name}/binning/metabat2/bins_metabat --saveCls -a {coverm_output_path}/coverm_metabat.tsv -t {threads} --minContig 1500 --seed 23'
    
    subprocess.call(args1, shell = True)
    subprocess.call(args2, shell = True)

def summarise_rosella(results_folder_name):
    fasta_files = glob.glob(f'{results_folder_name}/binning/rosella_out/*.fna')
    with open(f'{results_folder_name}/binning/rosella_bins.tsv', 'w') as out_bins:
        for file in fasta_files:
            fasta_bins = SeqIO.parse(open(file),'fasta')
            for contig in fasta_bins:
                out_bins.write(f'{contig.id}\t{file.split('/')[-1]}\n')

def summarise_concoct(results_folder_name):
    clustering_file = f'{results_folder_name}/binning/concoct/clustering_gt1000.csv' 
    with open(f'{results_folder_name}/binning/concoct_bins.tsv', 'w') as out_bins:
        opened_clustering_file = open(clustering_file,'r').readlines()
        for i, line in enumerate(opened_clustering_file):
            if i > 0:
                out_bins.write(f"{line.replace(',', '\tconcoct_bin.')}\n")

def summarise_metabat2(results_folder_name):
    contig_files = glob.glob(f'{results_folder_name}/binning/metabat2/*.*')
    with open(f'{results_folder_name}/binning/metabat2_bins.tsv', 'w') as out_bins:
        for file in contig_files:
            contigs = open(file, 'r').readlines()
            for contig in contigs:
                out_bins.write(f'{contig}\t{file}\n')


def run_dastool(results_folder_name):
    os.makedirs(f'{results_folder_name}/binning/dastool')
    args = f'DAS_Tool'


################################################################################
#################             MAIN             #################################
################################################################################

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="List files in a folder")

    # Add the folder path argument
    parser.add_argument("-a", "--assembly", type=str,
                        help="Path to the file as a string", required=True)
    parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end)", required=True)
    #parser.add_argument("-m", "--metadata_file", type=str,
    #                    help="Path to the metadata tsv file", required=True)
    parser.add_argument("-t", "--threads", type=str,
                        help="Number of threads to use for multiprocessing-compatible tasks", required=True)

    # Parse arguments
    args = parser.parse_args()
    out_folder = f'{args.name}_results'

    create_binning_folder(f'{out_folder}/binning')
    rename_filter_contigs(out_folder, args.assembly)
    
    run_coverm(out_folder, args.assembly, args.threads)

    run_concoct(out_folder, args.assembly, args.threads)
    run_rosella(out_folder, args.assembly, args.threads)
    run_metabat2(out_folder, args.assembly, args.threads)

    summarise_rosella(out_folder)
    summarise_metabat2(out_folder)
    summarise_concoct(out_folder)

    run_dastool(out_folder)
    

if __name__ == "__main__":
    main()
