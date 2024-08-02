########################################################
## Script that generates all the files of the shiny app
## We load grump_P16NS.csv load grump_P16NS_phyto.csv and load grump_P16NS_zoop.csv
## that were generated in 00_Subsetting-GRUMP...
########################################################

## Packages ########################################################
library(dplyr) ; library(tidyr)

## Directories and dataframes #####################################################
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")
path_grump_all <- root$find_file("data/grump_P16NS.csv")
path_grump_phyto <- root$find_file("data/grump_P16NS_phyto.csv")
path_grump_zoop <- root$find_file("data/grump_P16NS_zoop.csv")

df_grump_all <- data.table::fread(path_grump_all) %>% as.data.frame()
df_grump_phyto <- data.table::fread(path_grump_phyto) %>% as.data.frame()
df_grump_zoop <- data.table::fread(path_grump_zoop) %>% as.data.frame()


## Functions #####################################################

normalizeMatrix <- function(X){
  normMat = norm(X,type='2')
  return(X/normMat)
}

path_complete_finnest <- root$find_file('functions/complete_finnest_taxa.R')
source(path_complete_finnest)


path_create_biotic_dist_matrices <- root$find_file('functions/create_biotic_dist_matrices')
source(path_create_biotic_dist_matrices)

##############################################################################################
## 01 - Creating df_abiotic of each subset ###################################################
##############################################################################################
vet_abiotic = c(
  "Temperature","Salinity",
  "Oxygen","Silicate","NO2","NO3","PO4")

##############################################################################################
## 01 - Creating df_abiotic of each subset ###################################################
##############################################################################################

## All ----------------------------------------------------
df_geo_abiotics_all <- df_grump_all %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID)

## Fixing missing
df_geo_abiotics_all$Oxygen[df_geo_abiotics_all$SampleID=='P16N-S40-N10'] <- 212.9
df_geo_abiotics_all$NO3[df_geo_abiotics_all$SampleID=='P16S-S05-N15'] <- 1.375
df_geo_abiotics_all$Oxygen[df_geo_abiotics_all$SampleID=='P16S-S96-N28'] <- 346.9

## Phyto ----------------------------------------------------
df_geo_abiotics_phyto <- df_grump_phyto %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID)

## Fixing missing
df_geo_abiotics_phyto$Oxygen[df_geo_abiotics_phyto$SampleID=='P16N-S40-N10'] <- 212.9
df_geo_abiotics_phyto$NO3[df_geo_abiotics_phyto$SampleID=='P16S-S05-N15'] <- 1.375
df_geo_abiotics_phyto$Oxygen[df_geo_abiotics_phyto$SampleID=='P16S-S96-N28'] <- 346.9

## Zoop ----------------------------------------------------
df_geo_abiotics_zoop <- df_grump_zoop %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID)

## Fixing missing
df_geo_abiotics_zoop$Oxygen[df_geo_abiotics_zoop$SampleID=='P16N-S40-N10'] <- 212.9
df_geo_abiotics_zoop$NO3[df_geo_abiotics_zoop$SampleID=='P16S-S05-N15'] <- 1.375
df_geo_abiotics_zoop$Oxygen[df_geo_abiotics_zoop$SampleID=='P16S-S96-N28'] <- 346.9

## Storing ----------------------------------------------------
list_geo_abiotics<-list()
list_geo_abiotics$all <- df_geo_abiotics_all %>% arrange(SampleID)
list_geo_abiotics$phyto <- df_geo_abiotics_phyto %>% arrange(SampleID)
list_geo_abiotics$zoop <- df_geo_abiotics_zoop %>% arrange(SampleID)
saveRDS(list_geo_abiotics,paste0(savingdir,'/list_geo_abiotcs'))


##############################################################################################
## 02 - Creating biotic_distance matrices of each subset #####################################
##############################################################################################

## 02.1 Identifing finest taxonomic level -----------------------------------------------------------------------------------
## Fixing Taxonomy 
vet_tax_names <- c('Domain','Phylum','Class','Order','Family','Genus','Species')

df_grump_tax_wide = 
  df_grump_all %>% select(ASV_hash,Supergroup) %>% distinct() %>% 
  left_join(df_grump_all %>% select(ASV_hash,Division) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Domain) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Phylum) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Class) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Order) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Family) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Genus) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Species) %>% distinct()) %>% 
  left_join(df_grump_all %>% select(ASV_hash,Eco_relevant_plank_groups) %>% distinct()) 
  
# df_grump_tax_wide %>% select(-ASV_hash) %>% apply(.,2,function(x){sum(x=='')})

df_grump_tax_wide = df_grump_tax_wide %>% select(ASV_hash,Supergroup,Division,Eco_relevant_plank_groups) %>% 
  bind_cols(complete_finnest_taxa(df_grump_tax_wide %>% select(-ASV_hash,-Supergroup,-Division,-Eco_relevant_plank_groups)))


vet_tax_names <- c('Supergroup','Division','Domain','Phylum','Class','Order','Family','Genus','Species','Eco_relevant_plank_groups')

df_grump_all = df_grump_all %>% select(-any_of(vet_tax_names)) %>% 
  left_join(df_grump_tax_wide)

df_grump_phyto = df_grump_phyto %>% select(-any_of(vet_tax_names)) %>% 
  left_join(df_grump_tax_wide)

df_grump_zoop = df_grump_phyto %>% select(-any_of(vet_tax_names)) %>% 
  left_join(df_grump_tax_wide)

#################################################################################################################
## Creating distance matrices -----------------------------------------------------------------------------------
#################################################################################################################

list_dist_matrices_normalized <- list()
list_dist_matrices <- list()

df_geo_abiotics_all = df_geo_abiotics_all %>% arrange(SampleID)
df_geo_abiotics_phyto = df_geo_abiotics_phyto %>% arrange(SampleID)
df_geo_abiotics_zoop = df_geo_abiotics_zoop %>% arrange(SampleID)

#################################################################################################################
#################################################################################################################
## Creating distance matrices ----------------------------- All Plankton ----------------------------------------
#################################################################################################################
#################################################################################################################
abiotic_list_dist_matrices <- list()

abiotic_list_dist_matrices$all_grump<-list()

abiotic_list_dist_matrices$all_grump$geo = df_geo_abiotics_all %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

abiotic_list_dist_matrices$all_grump$geo_normalized = df_geo_abiotics_all %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

abiotic_list_dist_matrices$all_grump$abiotic <- df_geo_abiotics_all %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()

abiotic_list_dist_matrices$all_grump$abiotic_normalized <- df_geo_abiotics_all %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2), 
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() %>%  normalizeMatrix()

#################################################################################################################
#################################################################################################################
## Creating distance matrices ----------------------------- Phyto ----------------------------------------
#################################################################################################################
#################################################################################################################

abiotic_list_dist_matrices$phyto<-list()

abiotic_list_dist_matrices$phyto$geo = df_geo_abiotics_phyto %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

abiotic_list_dist_matrices$phyto$geo_normalized = df_geo_abiotics_phyto %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

abiotic_list_dist_matrices$phyto$abiotic <- df_geo_abiotics_phyto %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()

abiotic_list_dist_matrices$phyto$abiotic_normalized <- df_geo_abiotics_phyto %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2), 
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() %>%  normalizeMatrix()

#################################################################################################################
#################################################################################################################
## 03 - Creating distance matrices ----------------------------- Zoop ----------------------------------------
#################################################################################################################
#################################################################################################################

abiotic_list_dist_matrices$zoop<-list()

abiotic_list_dist_matrices$zoop$geo = df_geo_abiotics_zoop %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

abiotic_list_dist_matrices$zoop$geo_normalized = df_geo_abiotics_zoop %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

abiotic_list_dist_matrices$zoop$abiotic <- df_geo_abiotics_zoop %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()

abiotic_list_dist_matrices$zoop$abiotic_normalized <- df_geo_abiotics_zoop %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2), 
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() %>%  normalizeMatrix()


saveRDS(abiotic_list_dist_matrices,paste0(savingdir,'/geo-abiotic-list-dist-matrices'))
## In case I prefer separated files ::
saveRDS(abiotic_list_dist_matrices$all_grump,paste0(savingdir,'/geo-abiotic-list-dist-matrices-all-grump'))
saveRDS(abiotic_list_dist_matrices$phyto,paste0(savingdir,'/geo-abiotic-list-dist-matrices-all-phyto'))
saveRDS(abiotic_list_dist_matrices$zoop,paste0(savingdir,'/geo-abiotic-list-dist-matrices-all-zoop'))

#################################################################################################################
#################################################################################################################
## 04 - Creating Biotic Distance Matrices
#################################################################################################################
#################################################################################################################

## 04.1 Identifying rare ASVs, i.e. which ASVs are present in only one samples ###############
df_asvs_appearence = df_grump_all %>% 
  select(SampleID,ASV_hash) %>% distinct() %>% 
  group_by(ASV_hash) %>% 
  summarise(n_samples = n()) %>% arrange(n_samples)

not_rare_asvs <- df_asvs_appearence %>% filter(n_samples>1)

min_raw_count<-list()
min_raw_count$all_grump <- df_grump_all %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count$phyto <- df_grump_phyto %>% select(Raw.Sequence.Counts) %>% min()
min_raw_count$zoop <- df_grump_zoop %>% select(Raw.Sequence.Counts) %>% min()

## Since all of them have the same minimun small count:
min_raw_count = min_raw_count$zoop/100

#################################################################################################################
## 04.2 - Creating Biotic Distance Matrices - ALL Grump
#################################################################################################################
biotic_list_dist_matrices_all_grump <- create_biotic_dist_matrices(
  df_grump_all,min_raw_count = 0.01,
  vet_asv_at_least_two_samples = not_rare_asvs)

biotic_list_dist_matrices_phyto <- create_biotic_dist_matrices(
  df_grump_phyto,min_raw_count = 0.01,
  vet_asv_at_least_two_samples = not_rare_asvs)

biotic_list_dist_matrices_zoop <- create_biotic_dist_matrices(
  df_grump_zoop,min_raw_count = 0.01,
  vet_asv_at_least_two_samples = not_rare_asvs)

## So here we have two different list of distance matrices. 
saveRDS(biotic_list_dist_matrices_all_grump$biotic_list_dist_all_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_all_grump'))
saveRDS(biotic_list_dist_matrices_phyto$biotic_list_dist_all_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_phyto'))
saveRDS(biotic_list_dist_matrices_zoop$biotic_list_dist_all_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_zoop'))

saveRDS(biotic_list_dist_matrices_all_grump$biotic_list_dist_not_rare_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_all_grump_not_sparse'))
saveRDS(biotic_list_dist_matrices_phyto$biotic_list_dist_not_rare_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_phyto_not_sparse'))
saveRDS(biotic_list_dist_matrices_zoop$biotic_list_dist_not_rare_ASVs,file = paste0(savingdir,'/','biotic_list_dist_matrices_zoop_sparse'))

##################################################################################################
# MDS_Coordinates ################################################################################
##################################################################################################
biotic_list_dist_matrices_all_grump = readRDS(file = paste0(savingdir,'/','biotic_list_dist_matrices_all_grump'))
biotic_list_dist_matrices_phyto = readRDS(file = paste0(savingdir,'/','biotic_list_dist_matrices_phyto'))
biotic_list_dist_matrices_zoop = readRDS(file = paste0(savingdir,'/','biotic_list_dist_matrices_zoop'))



### Future work
# cmdscale(d = as.dist(justin_list$dist_matrices$biotic),k = 3,eig = T,add=T)
# pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)

list_mds_coord_all_grump <- list(
  ASV = cmdscale(biotic_list_dist_matrices_all_grump$ASV,k = 3),
  Species = cmdscale(biotic_list_dist_matrices_all_grump$Species,k = 3),
  Genus = cmdscale(biotic_list_dist_matrices_all_grump$Genus,k = 3),
  Family = cmdscale(biotic_list_dist_matrices_all_grump$Family,k = 3),
  Order = cmdscale(biotic_list_dist_matrices_all_grump$Order,k = 3),
  Class = cmdscale(biotic_list_dist_matrices_all_grump$Class,k = 3),
  Phylum = cmdscale(biotic_list_dist_matrices_all_grump$Phylum,k = 3),
  Domain = cmdscale(biotic_list_dist_matrices_all_grump$Domain,k = 3),
  Division = cmdscale(biotic_list_dist_matrices_all_grump$Division,k = 3),
  Supergroup = cmdscale(biotic_list_dist_matrices_all_grump$Supergroup,k = 3)
)

saveRDS(list_mds_coord_all_grump,file = paste0(savingdir,'/','list_mds_coord_all_grump'))

list_mds_coord_phyto <- list(
  ASV = cmdscale(biotic_list_dist_matrices_phyto$ASV,k = 3),
  Species = cmdscale(biotic_list_dist_matrices_phyto$Species,k = 3),
  Genus = cmdscale(biotic_list_dist_matrices_phyto$Genus,k = 3),
  Family = cmdscale(biotic_list_dist_matrices_phyto$Family,k = 3),
  Order = cmdscale(biotic_list_dist_matrices_phyto$Order,k = 3),
  Class = cmdscale(biotic_list_dist_matrices_phyto$Class,k = 3),
  Phylum = cmdscale(biotic_list_dist_matrices_phyto$Phylum,k = 3),
  Domain = cmdscale(biotic_list_dist_matrices_phyto$Domain,k = 3),
  Division = cmdscale(biotic_list_dist_matrices_phyto$Division,k = 3),
  Supergroup = cmdscale(biotic_list_dist_matrices_phyto$Supergroup,k = 3)
)

saveRDS(list_mds_coord_phyto,file = paste0(savingdir,'/','list_mds_coord_phyto'))

list_mds_coord_zoop <- list(
  ASV = cmdscale(biotic_list_dist_matrices_zoop$ASV,k = 3),
  Species = cmdscale(biotic_list_dist_matrices_zoop$Species,k = 3),
  Genus = cmdscale(biotic_list_dist_matrices_zoop$Genus,k = 3),
  Family = cmdscale(biotic_list_dist_matrices_zoop$Family,k = 3),
  Order = cmdscale(biotic_list_dist_matrices_zoop$Order,k = 3),
  Class = cmdscale(biotic_list_dist_matrices_zoop$Class,k = 3),
  Phylum = cmdscale(biotic_list_dist_matrices_zoop$Phylum,k = 3),
  Domain = cmdscale(biotic_list_dist_matrices_zoop$Domain,k = 3),
  Division = cmdscale(biotic_list_dist_matrices_zoop$Division,k = 3),
  Supergroup = cmdscale(biotic_list_dist_matrices_zoop$Supergroup,k = 3)
)

saveRDS(list_mds_coord_zoop,file = paste0(savingdir,'/','list_mds_coord_zoop'))


##################################################################################################
# Ecologic Groups ################################################################################
##################################################################################################
## All
eco_groups_all_grump_tax_lat_depth = df_grump_all %>% 
  mutate(Eco_relevant_plank_groups=ifelse(Eco_relevant_plank_groups=='','NotAssigned',Eco_relevant_plank_groups)) %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Raw.Sequence.Counts=sum(Raw.Sequence.Counts)) %>% 
  select(SampleID,Eco_relevant_plank_groups,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Eco_relevant_plank_groups ,
                     values_from = Sum_RawCounts,
                     values_fill = 0.01) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  left_join(list_geo_abiotics$all %>% select(SampleID,Latitude,Depth)) 

## Phyto
eco_groups_phyto_tax_lat_depth = df_grump_phyto %>% 
  mutate(Eco_relevant_plank_groups=ifelse(Eco_relevant_plank_groups=='','NotAssigned',Eco_relevant_plank_groups)) %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Raw.Sequence.Counts=sum(Raw.Sequence.Counts)) %>% 
  select(SampleID,Eco_relevant_plank_groups,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Eco_relevant_plank_groups ,
                     values_from = Sum_RawCounts,
                     values_fill = 0.01) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  left_join(list_geo_abiotics$phyto %>% select(SampleID,Latitude,Depth)) 

## zoop
eco_groups_zoop_tax_lat_depth = df_grump_zoop %>% 
  mutate(Eco_relevant_plank_groups=ifelse(Eco_relevant_plank_groups=='','NotAssigned',Eco_relevant_plank_groups)) %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Raw.Sequence.Counts=sum(Raw.Sequence.Counts)) %>% 
  select(SampleID,Eco_relevant_plank_groups,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Eco_relevant_plank_groups) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Eco_relevant_plank_groups ,
                     values_from = Sum_RawCounts,
                     values_fill = 0.01) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  left_join(list_geo_abiotics$zoop %>% select(SampleID,Latitude,Depth)) 

saveRDS(eco_groups_zoop_tax_lat_depth,paste0(savingdir,'/eco_groups_zoop_tax_lat_depth'))
saveRDS(eco_groups_phyto_tax_lat_depth,paste0(savingdir,'/eco_groups_phyto_tax_lat_depth'))
saveRDS(eco_groups_all_grump_tax_lat_depth,paste0(savingdir,'/eco_groups_all_grump_tax_lat_depth'))

