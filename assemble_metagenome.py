import os
import argparse
import subprocess
import pandas as pd # type: ignore
import datetime
import glob
import multiprocessing

################################################################################
#################           FUNCTIONS           ################################
################################################################################

# Preprocessing part

def create_result_folder(results_folder_name):
    if not os.path.exists(results_folder_name):
        os.makedirs(results_folder_name)
    
    with open(f'{results_folder_name}/log.txt', 'w') as log:
        log.write(f"Log file for the run {results_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

def print_env_summary(results_folder_name):
    subprocess.call(f'conda list > {results_folder_name}/list_conda.txt', shell = True)
    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write('Creating the list_conda.txt file that summarises the conda env.' + '\n\n')

def load_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep='\t', header=0)
    return(metadata)

def list_samples(folder_path):
    """
    Lists all folders in the given folder (= samples fastq files).
    """
    try:
        # Check if the path is a directory
        if not os.path.isdir(folder_path):
            print(f"The path '{folder_path}' is not a directory.")
            return

        # List all entries in the directory
        files = glob.glob(f'{folder_path}*/*.fastq.gz')
        files = [file.replace('_R1_001.fastq.gz', '').replace('_R2_001.fastq.gz', '') for file in files]
        files = [file for file in files if 'control' not in file]
        samples = list(set(files))
        print(samples)
        return(samples)

    except Exception as e:
        print(f"An error occurred: {e}")

def check_metadata_samples(metadata, samples, results_folder_name):
    """
    Check that the metadata agrees with the sample names listed.
    """
    metadata_samples = set(metadata['Barcode']) # samples in the metadata
    files_samples = set(samples) # samples in the reads' files
    print(files_samples)

    print('There are %d samples in the metadata file'%(len(metadata_samples)))
    print("There are %d samples in the reads' files"%(len(files_samples)))
    print("Barcode in the metadata but not in the reads' files:")
    print(metadata_samples.difference(files_samples))

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write('There are %d samples in the metadata file'%(len(metadata_samples)) + '\n\n')
        log.write("There are %d samples in the reads' files"%(len(files_samples)) + '\n\n')
        log.write("Barcode in the metadata but not in the reads' files:" + '\n\n')
        log.write(str(metadata_samples.difference(files_samples)) + '\n\n')

    out_list = [sample for sample in list(metadata_samples) if sample in files_samples]
    return metadata.loc[metadata['Barcode'].isin(out_list)], out_list



# Reads preprocessing part
def sample_trim_galore(sample, results_folder_name):
    reads1_in = f'{sample}_R1_001.fastq.gz'
    reads2_in = f'{sample}_R2_001.fastq.gz'
        
    args = ['trim_galore --fastqc --paired --length 50 -j 4',
            f'-o {results_folder_name}/trimmed_reads', reads1_in, reads2_in]
    subprocess.call(' '.join(args), shell = True)

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(' '.join(args) + '\n\n')

def run_trimming(samples, results_folder_name):
    os.makedirs(f'{results_folder_name}/trimmed_reads')

    pool = multiprocessing.Pool(10)
    pool.starmap(sample_trim_galore, zip(samples, [results_folder_name for sample in samples])) 

# Metagenome assembly
def run_megahit(file1, file2, results_folder_name, threads):
    args = [f'megahit --presets meta-large --k-min 27 --k-max 87 --k-step 10 --min-contig-len 1000 -t {threads}', 
            f'-1 {file1} -2 {file2} -o {results_folder_name}/full_coassembly_megahit']
    subprocess.call(' '.join(args), shell = True)

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(' '.join(args) + '\n\n')

def run_metaspades(file1, file2, results_folder_name, threads):
    args = [f'metaspades.py -t {threads} -m 400 -k 27,37,47,57,67,77,87', 
            f'-1 {file1} -2 {file2} -o {results_folder_name}/full_coassembly_metaspades']
    subprocess.call(' '.join(args), shell = True)

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(' '.join(args) + '\n\n')

def full_coassembly(results_folder_name, threads, software = 'both'):
    all_r1 = f'{results_folder_name}/trimmed_reads/all_trimmed_R1.fq.gz'
    all_r2 = f'{results_folder_name}/trimmed_reads/all_trimmed_R2.fq.gz'

    subprocess.call(f'cat {results_folder_name}/trimmed_reads/*_val_1.fq.gz > {all_r1}', shell = True)
    subprocess.call(f'cat {results_folder_name}/trimmed_reads/*_val_2.fq.gz > {all_r2}', shell = True)
    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(f'cat {results_folder_name}/trimmed_reads/*_val_1.fq.gz > {all_r1}' + '\n\n')
        log.write(f'cat {results_folder_name}/trimmed_reads/*_val_1.fq.gz > {all_r2}' + '\n\n')

    if software in ['both', 'megahit']:
        run_megahit(all_r1, all_r2, results_folder_name, threads)
    
    if software in ['both', 'metaspades']:
        run_metaspades(all_r1, all_r2, results_folder_name, threads)

def assemblies_stats(results_folder_name, threads, software = 'both'):
    if software == 'both':
        assembly_files = ''
    elif software == 'metaspades':
        assembly_files = ''
    elif software == 'megahit':
        assembly_files = '{results_folder_name}/full_coassembly_megahit/final.contigs.fa'
    else:
        print('The method set for software is not recognised, skipping the metaquast run.')
        return()


    args = f'metaquast {assembly_files} -t {threads} -o {results_folder_name}/metaquast_results'
    subprocess.call(args, shell = True)

################################################################################
#################             MAIN             #################################
################################################################################

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="List files in a folder")

    # Add the folder path argument
    parser.add_argument("-illu", "--illuminafolder", type=str,
                        help="Path to the folder as a string", required=True)
    parser.add_argument("-nano", "--nanoporefolder", type=str,
                        help="Path to the folder as a string")
    parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end)", required=True)
    #parser.add_argument("-m", "--metadata_file", type=str,
    #                    help="Path to the metadata tsv file", required=True)
    parser.add_argument("-t", "--threads", type=str,
                        help="Number of threads to use for multiprocessing-compatible tasks", required=True)
    parser.add_argument("--skippreprocessing", action='store_true',
                        help="To add if you want to skip preprocessing")


    # Parse arguments
    args = parser.parse_args()
    out_folder = f'{args.name}_results'
    
    #create_result_folder(out_folder)
    #print_env_summary(out_folder)
    
    #samples = list_samples(args.illuminafolder)
    #run_trimming(samples, out_folder)

    full_coassembly(out_folder, args.threads, software = 'both')
    assemblies_stats(out_folder, args.threads, software = 'both')
    

if __name__ == "__main__":
    main()
