import os
import argparse
import subprocess
import pandas as pd # type: ignore
import datetime
import glob

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

def list_subfolders(folder_path):
    """
    Lists all folders in the given folder (= samples fastq files).
    """
    try:
        # Check if the path is a directory
        if not os.path.isdir(folder_path):
            print(f"The path '{folder_path}' is not a directory.")
            return

        # List all entries in the directory
        entries = os.listdir(folder_path)

        # Filter out files, keep only directories
        subfolders = [entry for entry in entries if os.path.isdir(os.path.join(folder_path, entry))]
        filtered_subfolders = [subfolder for subfolder in subfolders if subfolder != 'unclassified']
        return(filtered_subfolders)

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

def concatenate_files(folder_path, metadata, samples, results_folder_name):
    """
    Takes folder paths, metadata files and samples list and concatenate the data
    for each sample in the list, + renames it to the sample name using metadata.
    """
    os.makedirs('%s/raw_data' % (results_folder_name))
    new_samples = []
    for sample in samples:
        new_sample = metadata.loc[metadata['Barcode'] == sample, '#SampleID'].values[0]
        new_path = f'{results_folder_name}/raw_data/{new_sample}.fastq.gz'
        args = ['cat', folder_path + sample + '/*.fastq.gz', '>', new_path]
        subprocess.call(' '.join(args), shell = True)
        new_samples.append(new_sample)
    return(new_samples)


# Reads preprocessing part

def run_trim_galore(samples, threads, results_folder_name):
    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write(f'Running porechop...' + '\n\n')

    for sample in samples:
        reads_in = f'{results_folder_name}/raw_data/{sample}.fastq.gz'
        reads_out = f'{results_folder_name}/raw_data/{sample}_porechopped.fastq.gz'
        args = ['trim_galore --paired --illumina ',
                *.clock_UMI.R1.fq.gz *.clock_UMI.R2.fq.gz']
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
    parser.add_argument("-illu", "--illumina-folder", type=str,
                        help="Path to the folder as a string", required=True)
    parser.add_argument("-nano", "--nanopore-folder", type=str,
                        help="Path to the folder as a string")
    parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end)", required=True)
    parser.add_argument("-m", "--metadata_file", type=str,
                        help="Path to the metadata tsv file", required=True)
    parser.add_argument("-t", "--threads", type=str,
                        help="Number of threads to use for multiprocessing-compatible tasks", required=True)
    parser.add_argument("--skippreprocessing", action='store_true',
                        help="To add if you want to skip preprocessing")
    parser.add_argument("--skipqiime2", action='store_true',
                        help="To add if you want to skip the qiime2 part (only preprocessing)")


    # Parse arguments
    args = parser.parse_args()
    out_folder = f'{args.name}_results'
    
    create_result_folder(out_folder)
    print_env_summary(out_folder)
    metadata = load_metadata(args.metadata_file)
    
    samples = list_subfolders(args.folder)
    metadata, samples_to_process = check_metadata_samples(metadata, samples, out_folder)

        # Concatenate files belonging to the same sample in the new directory,
        # run porechop and chopper
    run_trim_galore(samples, threads, out_folder)
    samples_names = concatenate_files(args.folder, metadata, samples_to_process, out_folder)

    full_coassembly(out_folder, metadata, samples_names, args.threads, software = 'megahit')
    full_coassembly(out_folder, metadata, samples_names, args.threads, software = 'metaspades')

    sub_coassemblies(out_folder, metadata, samples_names, args.threads, software = 'megahit')
    sub_coassemblies(out_folder, metadata, samples_names, args.threads, software = 'metaspades')
    

if __name__ == "__main__":
    main()
