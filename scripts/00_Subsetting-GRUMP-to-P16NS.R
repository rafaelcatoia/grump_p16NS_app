#####################################################################
## Script that filters only P16NS Cruise & Saves 3 different datasets
## grump_P16NS.csv stores all ASVs that were observed
## grump_P16NS_phyto.csv stores all Phytoplankton ASVs
## grump_P16NS_zoop.csv stores all Zoop ASVs
#####################################################################

## Packages ########################################################
library(dplyr) ; library(tidyr)

## Directories & functions #####################################################
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")
path_grump = root$find_file('data/grump_asv_long_May21_24.csv')
path_complete_finnest <- root$find_file('functions/complete_finnest_taxa.R')
source(path_complete_finnest)


df_grump_all <- data.table::fread(path_grump) %>% 
  filter(Cruise %in% c('P16N','P16S')) %>% 
  mutate(Raw.Sequence.Counts = Corrected_sequence_counts) %>% 
  filter(Depth<600)

## Saving All ASVs but only for P16SN
data.table::fwrite(x = df_grump_all,file = paste0(datadir,'/','grump_P16NS.csv'))
## Summary #####################################################

## Part 01 -- Filter Only Phyto

## Part 02 -- Filter Only Zoop

##############################################################
## Begin #####################################################
## Part 01 -- Filter Only Phyto ##############################
## This filtering was provided by Yubin ######################
##############################################################

df_grump_all[, c(2:11,54)][is.na(df_grump_all[, c(2:11,54)])] <- "Not Applicable" 
cyano <- rbind(df_grump_all[df_grump_all$Eco_relevant_plank_groups=="Prochlorococcus",],
               df_grump_all[df_grump_all$Eco_relevant_plank_groups=="Synechococcus",],
               df_grump_all[grepl("Richelia", df_grump_all$Genus),],
               df_grump_all[grepl("Trichodesmium", df_grump_all$Genus),],
               df_grump_all[grepl("UCYN-A", df_grump_all$Genus),],
               df_grump_all[grepl("Crocosphaera", df_grump_all$Genus),],
               df_grump_all[grepl("UCYN-C", df_grump_all$Genus),])
## Subset out the chloroplast 16S ASVs from the main dataframe
chl.16s <- df_grump_all[df_grump_all$Sequence_Type=="Chloroplast_16S",]
## Remove the following ASVs from that subset
chl.16s <- chl.16s[!chl.16s$Supergroup=="Rhizaria",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Excavata",]
chl.16s <- chl.16s[!chl.16s$Supergroup=="Alveolata",]
chl.16s <- chl.16s[!chl.16s$Division=="Rhodophyta",]
## Final "phytoplankton" community we will use for the time being

phyto.all <- rbind(cyano, chl.16s)
## From here, subset out the P16 S/N cruises to focus on that <- this is more or less the 
## transect that I have been using for a lot of other analysis

## Saving file:: 
data.table::fwrite(x = phyto.all,file = paste0(datadir,'/','grump_P16NS_phyto.csv'))

##############################################################
## Begin #####################################################
## Part 02 -- Filter Only Zoop ###############################
## The OTUs were provided by Yubin ###########################
##############################################################

## Re-loading all grump
df_grump_all <- data.table::fread(path_grump) %>% 
  filter(Cruise %in% c('P16N','P16S')) %>% 
  mutate(Raw.Sequence.Counts = Corrected_sequence_counts) %>% 
  filter(Depth<600)

dat_tax = data.table::fread(
  'https://raw.githubusercontent.com/rafaelcatoia/zoop_16N/main/treated_taxonomy_dat.csv') %>%
    as_tibble()

grump_zoop <- df_grump_all %>% filter(ASV_hash %in% dat_tax$ASV_ID)

## Saving file::
data.table::fwrite(x = grump_zoop,file = paste0(datadir,'/','grump_P16NS_zoop.csv'))