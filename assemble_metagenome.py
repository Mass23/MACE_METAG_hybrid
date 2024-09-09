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
        samples = list(set(files))
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
            '-o {results_folder_name}/trimmed_reads', reads1_in, reads2_in]
    subprocess.call(' '.join(args), shell = True)

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(' '.join(args) + '\n\n')

def run_trimming(samples, results_folder_name):
    os.makedirs(f'{results_folder_name}/trimmed_reads')

    pool = multiprocessing.Pool(10)
    pool.starmap(sample_trim_galore, zip(samples, [results_folder_name for sample in samples])) 


def run_megahit(file1, file2, out_folder):
    args = ['megahit --presets meta-large', '--k-min 27', '--k-step 10',
            '-1', file1, '-2', file2, '-o', out_folder]
    subprocess.call(' '.join(args), shell = True)

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(' '.join(args) + '\n\n')


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
    
    create_result_folder(out_folder)
    print_env_summary(out_folder)
    
    samples = list_samples(args.illuminafolder)
    run_trimming(samples, out_folder)

    #full_coassembly(out_folder, metadata, samples_names, args.threads, software = 'megahit')
    #full_coassembly(out_folder, metadata, samples_names, args.threads, software = 'metaspades')

    #sub_coassemblies(out_folder, metadata, samples_names, args.threads, software = 'megahit')
    #sub_coassemblies(out_folder, metadata, samples_names, args.threads, software = 'metaspades')
    

if __name__ == "__main__":
    main()
