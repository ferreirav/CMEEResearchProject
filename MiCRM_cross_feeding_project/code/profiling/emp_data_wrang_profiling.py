#!/usr/bin/env python3

"""This script loads and subsets data from the Earth Microbiome Project by sample to be 
fed into a BASH script that will run CarveME to produce individual metabolic models and SMETANA to 
calculate community interactions"""

__appname__ = 'emp_data_wrang_profiling.py'
__author__ = 'Vitor Ferreira (f.ferreira22@imperial.ac.uk)'
__version__ = '0.0.1'

###IMPORTS###
import os
import pandas as pd
import urllib.request
# import sys

# Adding to system path
# sys.path.insert(0, '/Users/vitorferreira/Desktop/Imperial_College/ResearchProject/Tutorials/Metabolic_Pathway/cooccurrence-master/pipeline/code')


### Function to generate data ###
def get_data(size):

    species = int(size)

    ### Loading Metadata ###
    # Data can be found at the EMP GitHub Repository - https://github.com/biocore/emp/tree/master/data/mapping-files
    metadata = pd.read_csv("../../data/emp_qiime_mapping_qc_filtered.tsv", sep="\t")
    metadata.rename(columns={'#SampleID': 'sample'}, inplace=True)


    ### Loading Sample Data ###
    # Data that has been filtered and used in the Polarisation paper from Machado et al. (2021) and can be downloaded directly from https://oc.embl.de/index.php/s/SbhoQa9YJoa748V/download
    samples = pd.read_csv("../../data/emp_150bp_filtered.tsv", sep="\t")


    ### Genome-Scale Metabolic Models ###
    # A collection of Genome-Scale Models for bacterial species can be found at the repository - https://github.com/cdanielmachado/embl_gems
    metabolic_models = pd.read_csv("../../data/model_list.tsv", sep="\t", usecols=[0,4])
    # Generate organisms ID from the file path column
    metabolic_models['org_id'] = metabolic_models['file_path'].apply(lambda x: os.path.basename(x)[:-7])
    metabolic_models = metabolic_models.drop('file_path', axis=1)


    ### Merging data sets ###
    # Step 1: 'samples' and 'metabolic_models' by 'org_id'
    samples = pd.merge(samples, metabolic_models, on='org_id')

    # Step 2: Selecting features of interest to merge
    features_interest = ['sample', 'Description', 'title', 'principal_investigator', 'ebi_accession',
                        'target_gene', 'sample_taxid', 'sample_scientific_name', 'host_taxid', 
                        'host_common_name_provided', 'host_common_name', 'host_scientific_name', 'host_kingdom',
                        'host_phylum', 'host_class', 'host_order', 'host_family', 'host_genus', 'host_species', 
                        'env_biome', 'env_feature', 'env_material', 'envo_biome_0', 'envo_biome_1', 'envo_biome_2', 
                        'envo_biome_3', 'envo_biome_4', 'envo_biome_5', 'empo_0', 'empo_1', 'empo_2', 'empo_3']

    # Step3: Merging 'samples' and 'metadata'
    samples_envo_metadata = pd.merge(samples, metadata[features_interest], on='sample')


    ### FINALLY: ###
    # I will store the samples environmental entology (ENVO) metadata and a sample list
    #sample_list = sorted(set(samples_envo_metadata["sample"]))



    # *********************************************** #
    #                Testing for Speed                #
    # Communities of 2 species
    if species == 2:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # Communities of 4 species
    elif species == 4:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # Communities of 8 species
    elif species == 8:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # Communities of 16 species
    elif species == 16:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # Communities of 32 species
    elif species == 32:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # Communities of 64 species

    elif species == 64:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(species,species)]
        sample_list = comms_2_size.iloc[0][0]

    #-------------------------------------------------#
    # ELSE...get the smallest
    else:
        samples_comms_sizes = samples_envo_metadata.groupby(['sample']).size().reset_index(name='community_size')
        comms_2_size = samples_comms_sizes[samples_comms_sizes['community_size'].between(2,2)]
        sample_list = comms_2_size.iloc[0][0]

    # *********************************************** #



    # to_csv
    #sample_list.to_csv("../results/sample_list.csv")
    #samples_envo_metadata.to_csv("../results/samples_envo_metadata.csv", index = False, header = False)

    return sample_list, samples_envo_metadata



def download_ncbi_genome_temp(accession, refseq_table, sample, prefer_protein=True, overwrite=False):

    # Consider to add Subdirectory xxx/sample
    genomes_folder = '../genomes_folder/'

    if not os.path.exists(genomes_folder + sample):
        os.makedirs(genomes_folder + sample)

    sample_path = genomes_folder + sample + '/'

    if accession not in refseq_table.index:
        print('Invalid accession code')
        return

    entry = refseq_table.loc[accession, :]

    
    downloaded = False

    if prefer_protein:
        url = 'https://{}/{}_protein.faa.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = sample_path + '{}.faa.gz'.format(accession)

        if os.path.exists(outputfile) and not overwrite:
            print('File exists, skipping.')
            return outputfile

        _, result = urllib.request.urlretrieve(url, outputfile)

        if result.get_content_type() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if not downloaded:
        url = 'https://{}/{}_genomic.fna.gz'.format(
            entry['ftp_path'][6:], entry['ftp_path'].split('/')[-1])

        outputfile = sample_path + '{}.fna.gz'.format(accession)

        if os.path.exists(outputfile) and not overwrite:
            print('File exists, skipping.')
            return outputfile

        _, result = urllib.request.urlretrieve(url, outputfile)

        if result.get_content_type() != 'application/x-gzip':
            os.remove(outputfile)
        else:
            downloaded = True

    if downloaded:
        return outputfile
    


def get_genome(org_list, temp_id_sample, temp_id, refseq):
    for org in range(len(org_list)):
        accession = temp_id_sample.iloc[org]['assembly_accession']
        download_ncbi_genome_temp(accession, refseq, temp_id)

