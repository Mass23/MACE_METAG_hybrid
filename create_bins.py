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

def create_binning_folder(binning_folder_name):
    if not os.path.exists(binning_folder_name):
        os.makedirs(binning_folder_name)
    
    with open(f'{binning_folder_name}/log.txt', 'w') as log:
        log.write(f"Log file for the run {binning_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

# Coverage computation
def run_coverm(results_folder_name, threads, assembly_path):
    args = ['coverm contig', '-t', str(threads), '-r', assembly_path, '-m', 'metabat mean variance'
            '-1', f'{results_folder_name}/trimmed_reads/*val_1.fq.gz',
            '-2', f'{results_folder_name}/trimmed_reads/*val_2.fq.gz',
            '-o', f'{results_folder_name}/binning/coverm_depth.tsv']
    subprocess.call(' '.join(args), shell = True)
    

    

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
    
    run_coverm(out_folder, args.threads, args.assembly)
    run_metabat2()
    run_maxbin2()

    run_dastool()
    

if __name__ == "__main__":
    main()
