#############################################################
# Exploring the EMP project and Polarization paper data     #
#############################################################

require(tidyverse)

# LOAD METADATA & SAMPLES

metadata <- read_tsv("../data/emp_qiime_mapping_qc_filtered.tsv")
samples <- read_tsv("../data/emp_150bp_filtered.tsv")

# Filter relevant metadata for community assembly

meta_filtered <- metadata %>% 
  select("#SampleID", "Description", "host_subject_id", "study_id", "title",
         "principal_investigator", "sample_taxid", "sample_scientific_name",
         "host_taxid", "host_common_name", "host_scientific_name", "host_family", 
         "host_genus", "host_species", "env_biome", "env_feature", "env_material",
         "envo_biome_0", "envo_biome_1", "envo_biome_2", "envo_biome_3",
         "envo_biome_4", "envo_biome_5", "empo_0", "empo_1", "empo_2", "empo_3")


# Community Sizes per sample
community_sizes <- samples %>% 
  group_by(sample) %>% 
  summarise(Community_size = n())

hist(community_sizes$Community_size, breaks = 50)


# Exploring the data

# Per Principal Investigator
pi_metadata <- metadata %>% 
  filter(principal_investigator == "Catherine Pfister")


empo_by_pi <- metadata %>% group_by(principal_investigator) %>% 
  summarise(empo_3 = length(unique(empo_3)))



# Unique values for ontologies (envo, empo, "extra")

# envo
unique(metadata$envo_biome_0)
unique(metadata$envo_biome_1)
unique(metadata$envo_biome_2)
unique(metadata$envo_biome_3)
unique(metadata$envo_biome_4)


# empo
unique(meta_df$empo_0)
unique(meta_df$empo_1)
unique(meta_df$empo_2)
unique(meta_df$empo_3)


# "extra"
unique(meta_df$env_biome)
unique(meta_df$env_feature)
unique(meta_df$env_material)




