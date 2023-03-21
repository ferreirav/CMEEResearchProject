#!/usr/bin/env python3

"""
This pipeline calls `emp_data_wrang` script which loads Earth Microbiome Project (EMP) data, RefSeq Release 201 and 
combine these into a data frame to be loaded here with the relevant information to sample microbial communities per environmental material.
This is achieved by generating a set of unique sample identifiers and fed into multiprocessing processes, sample by sample.
Once generated the data for unique microbial communities, it is piped into bash scripts that runs CarveMe, to generate metabolic models, and 
SMETANA that calculates the interaction potential of the communities. Throughout the process, temporary files and folders are generated per 
sample which are removed once each processed is completed.

The final output files can be found in output folder:
- "{specific_sample}_global.tsv" for whole community interaction metrics.
- "{specific_sample}_detailed.tsv" with the individual interaction values between organisms and metabolites produced.

"""

__appname__ = 'get_interactions.py'
__author__ = 'Vitor Ferreira (f.ferreira22@imperial.ac.uk)'
__version__ = '0.0.1'


### IMPORTS ###
import os
import time
import sys
import subprocess
import pandas as pd
from multiprocessing import Pool
from emp_data_wrang import get_data, get_genome


### Functions ###
def processing_interactions():
    """
    This function gets the EMP data for the respective sample, downloads the genomes and stores the information in a csv file for the community organisms.
    Within this function BASH is called to create a copy of a generic template that will handle individual samples and runs CarveMe and SMETANA functions within.
    """

    start = time.time()

    ##### Getting input data #####
    sample_list, samples_metadata = get_data()
    
    #####################################################
    ##### UNCOMMENT BELOW FOR MULTIPROCESSING RUNS ######
    #####          AND COMMENT LINE ABOVE          ######
    #####################################################
    #samples_metadata = get_data()[1]
    
    temp_id = sample_list

    print(f"Started processing the sample {temp_id}...")
    
    temp_id_sample = samples_metadata[samples_metadata['sample'] == temp_id]
    org_ids = set(temp_id_sample['org_id'])

    ### Loading REF-SEQ REFERENCE TABLE ###
    # Currently, using the same reference table as in the Polarisation paper by Machado et al. (2021) - https://www.nature.com/articles/s41559-020-01353-4
    # However, a new RefSeq Release 216 is available for download!
    refseq = pd.read_csv("../data/refseq_release_201.tsv", sep='\t', index_col=[0])

    ### Downloading organisms genomes for the sampled data
    get_genome(org_ids, temp_id_sample, temp_id, refseq)

    ### Saving CSV file to be accessed in BASH
    temp_id_sample.to_csv(f"../results/output/{temp_id}.csv", index=False,
                          header=False, columns=['assembly_accession', 'org_id'])
    
    ### Copying BASH template to sample_specific template and running it accordingly
    subprocess.Popen(f"cp get_acessions.sh {temp_id}.sh", shell=True).wait()
    subprocess.Popen(f"bash ./{temp_id}.sh", shell=True).wait()
    # Delete temporary sample_specific template
    os.system(f"rm {temp_id}.sh")
    
    print(f"the sample {temp_id} has taken {time.time() - start} s to process.")


if __name__ == '__main__':
    processing_interactions()

    #####################################################
    ##### UNCOMMENT BELOW FOR MULTIPROCESSING RUNS ######
    #####            AND COMMENT LINE ABOVE        ######
    #####################################################

    # sample_list = get_data()[0]

    # with Pool() as p:
    #     p.map(processing_interactions, sample_list)