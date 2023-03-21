#!/usr/bin/env python3

"""T"""

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
from emp_data_wrang_profiling import get_data, get_genome


### Functions ###
def processing_interactions(comm_size):

    print("Community Size is: " + comm_size)
    start = time.time()
    
    
    sample_list, samples_metadata = get_data(size = comm_size)
    
    temp_id = sample_list

    print(f"Started processing the sample {temp_id}...")
    
    temp_id_sample = samples_metadata[samples_metadata['sample'] == temp_id]
    org_ids = set(temp_id_sample['org_id'])

    ### Loading REF-SEQ REFERENCE TABLE ###

    # Currently, using the same reference table as in the Polarisation paper by Machado et al. (2021) - https://www.nature.com/articles/s41559-020-01353-4
    # However, a new RefSeq Release 216 is available for download!
    refseq = pd.read_csv("../../data/refseq_release_201.tsv", sep='\t', index_col=[0])

    get_genome(org_ids, temp_id_sample, temp_id, refseq)

    print("done subset")

    temp_id_sample.to_csv(f"../../results/output/{temp_id}.csv", index=False,
                          header=False, columns=['assembly_accession', 'org_id'])
    
    #subprocess.call(['bash', './get_acessions.sh'])
    subprocess.Popen(f"cp get_acessions_profiling.sh {temp_id}.sh", shell=True).wait()
    subprocess.Popen(f"bash ./{temp_id}.sh", shell=True).wait()

    #os.system(f"cp get_acessions.sh {temp_id}.sh")
    #os.system(f"bash ./{temp_id}.sh")

    os.system(f"rm {temp_id}.sh")
    
    print(f"the sample {temp_id} has taken {time.time() - start} s to process.")



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Provides the community size for scaling computation")
    parser.add_argument('--size', metavar='size', required=True, help="The size of community to be assessed")
    args = parser.parse_args()
    processing_interactions(comm_size=args.size)

    # sample_list = get_data()[0]

    # with Pool() as p:
    #     p.map(processing_interactions, sample_list)