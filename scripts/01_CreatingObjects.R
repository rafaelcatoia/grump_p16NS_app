


normalizeMatrix <- function(X){
  normMat = norm(X,type='2')
  return(X/normMat)
}

#files_vec <- list.files(funsdir)
#
#for( i in 1:length(files_vec)){
#  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
#}

# dat_tax = data.table::fread('https://raw.githubusercontent.com/rafaelcatoia/zoop_16N/main/treated_taxonomy_dat.csv') %>%
#   as_tibble()

### Loading the new data 

## Use this if you are using the new GRUMP Data Set
datapath = root$find_file(paste0(datadir,'/','grump_phyto.csv'))

dframe = data.table::fread(input = datapath) %>%
  filter(Cruise %in% c('P16N','P16S')) %>% 
  #filter(Raw.Sequence.Counts>0) %>% 
  #filter(Domain!='Unassigned') %>% 
  mutate(Raw.Sequence.Counts = Corrected_sequence_counts) %>% 
  #### New filter by Depth!
  filter(Depth<600)

### -----------------------------------------------------------
## Creating abiotics dataframe  -------------------------------
### -----------------------------------------------------------
vet_abiotic = c(
  "Temperature",
  "Salinity",
  "Oxygen",
  "Silicate",
  "NO2",
  "NO3",#this causes duplicates
  #"NH3",#this is empty
  "PO4"
)

## Filtering the abiotics per sample
df_geo_abiotics <- dframe %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  mutate(Oxygen=ifelse(Oxygen<0,NA,Oxygen),NO3=ifelse(NO3<0,NA,NO3)) %>% 
  distinct() %>% arrange(SampleID)


## Solving the missing problem -------------------------------------------------

df_geo_abiotics <- dframe %>%
  select(SampleID,one_of(vet_abiotic),Latitude,Longitude,Depth,Longhurst_Short) %>%
  distinct() %>% arrange(SampleID)


## Here we have the sample with missing values on Oxygen and NO3
df_geo_abiotics %>% 
  filter(NO3<0 |Oxygen<0)


## First inputation = 212.9
df_geo_abiotics %>% filter(Latitude ==18) %>% arrange(Depth)
df_geo_abiotics$Oxygen[df_geo_abiotics$SampleID=='P16N-S40-N10'] <- 212.9

## Second inputation = 0.13 + 2.49 / 2 = 1.375
df_geo_abiotics %>% filter(Latitude ==-18.0) %>% arrange(Depth)
df_geo_abiotics$NO3[df_geo_abiotics$SampleID=='P16S-S05-N15'] <- 1.375

## Third inputation = 346.9
df_geo_abiotics %>% filter(Latitude ==-63.5) %>% arrange(Depth)
df_geo_abiotics$Oxygen[df_geo_abiotics$SampleID=='P16S-S96-N28'] <- 346.9

## Saving the abiotics df
saveRDS(df_geo_abiotics,file = paste0(savingdir,'/','df_geo_abiotics'))
dim(df_geo_abiotics)

## Saving the filtered Grump stacked
saveRDS(dframe,file = paste0(savingdir,'/','dfGrump_longer_filtered'))

## Creating distance matrices ----------------------------------------------------------------
list_dist_matrices_normalized <- list()
list_dist_matrices <- list()

df_geo_abiotics = df_geo_abiotics %>% arrange(SampleID)

## First the normalized ones 

list_dist_matrices$geo = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix()

list_dist_matrices_normalized$geo = df_geo_abiotics %>% 
  mutate(Latitude = abs(Latitude)) %>% 
  transmute(lat_scaled = (Latitude -mean(Latitude)) /sd(Latitude),
            depht_scaled = (Depth -mean(Depth)) /sd(Depth)) %>%
  as.matrix() %>% dist() %>% as.matrix() %>% normalizeMatrix()

list_dist_matrices$abiotic <- df_geo_abiotics %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2),
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() #%>%  normalizeMatrix()

list_dist_matrices_normalized$abiotic <- df_geo_abiotics %>% 
  transmute(Temperature = scale(Temperature),
            Salinity = scale(Salinity),
            Oxygen = scale(Oxygen),
            Silicate = scale(Silicate),
            NO2 = scale(NO2), 
            NO2 = scale(NO3),
            PO4 = scale(PO4)) %>% 
  as.matrix() %>% dist() %>% as.matrix() %>%  normalizeMatrix()

#######################################################################################
# Now lets create the biotic part.
###########################################################################################################s

## First defining the ASVs that are " rare "
df_asvs_appearence = dframe %>% 
  select(SampleID,ASV_hash) %>% distinct() %>% 
  group_by(ASV_hash) %>% 
  summarise(n_samples = n()) %>% arrange(n_samples)

vet_asv_at_least_two_samples <- df_asvs_appearence %>% filter(n_samples>1) %>% 
  select(ASV_hash) %>% pull()

# vet_asv_at_least_three_samples <- df_asvs_appearence %>% filter(n_samples>2) %>% 
#   select(ASV_hash) %>% pull()


## adding a small value to make sure we can calculate aitchison distances

min_raw_count=dframe %>%
  select(Raw.Sequence.Counts) %>% min()

min_raw_count = min_raw_count/100

biotic_list_dist <- list()
biotic_list_dist_at_least_two <- list()
##########################################################################################
# ASV ####################################################################################
##########################################################################################

biotic_list_dist$ASV <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ASV_hash,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,ASV_hash) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ASV_hash ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$ASV_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ASV_hash,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,ASV_hash) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ASV_hash ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$ASV <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ASV_hash,Raw.Sequence.Counts) %>% distinct() %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  group_by(SampleID,ASV_hash) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ASV_hash ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$ASV_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,ASV_hash,Raw.Sequence.Counts) %>% distinct() %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  group_by(SampleID,ASV_hash) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = ASV_hash ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


##########################################################################################
# Species ################################################################################
##########################################################################################


biotic_list_dist$Species <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Species,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Species) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Species ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Species_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Species,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Species) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Species ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Species <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Species,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Species) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Species ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Species_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Species,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Species) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Species ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()



##########################################################################################
# Genus ################################################################################
##########################################################################################

biotic_list_dist$Genus <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Genus,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Genus) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Genus ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Genus_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Genus,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Genus) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Genus ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Genus <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Genus,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Genus) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Genus ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Genus_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Genus,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Genus) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Genus ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()

##########################################################################################
# Family ################################################################################
##########################################################################################

biotic_list_dist$Family <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Family,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Family) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Family ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Family_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Family,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Family) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Family ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Family <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Family,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Family) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Family ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Family_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Family,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Family) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Family ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()

##########################################################################################
# Order ################################################################################
##########################################################################################

biotic_list_dist$Order <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Order,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Order) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Order ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Order_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Order,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Order) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Order ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Order <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Order,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Order) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Order ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Order_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Order,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Order) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Order ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()

##########################################################################################
# Class ################################################################################
##########################################################################################

biotic_list_dist$Class <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Class,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Class) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Class ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Class_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Class,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Class) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Class ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Class <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Class,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Class) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Class ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Class_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Class,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Class) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Class ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()

##########################################################################################
# Division ################################################################################
##########################################################################################

biotic_list_dist$Division <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Division,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Division) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Division ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Division_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Division,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Division) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Division ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Division <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Division,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Division) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Division ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Division_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Division,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Division) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Division ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()

##########################################################################################
# Supergroup ################################################################################
##########################################################################################

biotic_list_dist$Supergroup <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Supergroup,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Supergroup) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Supergroup ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist$Supergroup_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  select(SampleID,Supergroup,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Supergroup) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Supergroup ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


biotic_list_dist_at_least_two$Supergroup <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Supergroup,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Supergroup) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Supergroup ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()

biotic_list_dist_at_least_two$Supergroup_normalized <- dframe %>%
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
  filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
  select(SampleID,Supergroup,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(SampleID,Supergroup) %>% 
  summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
  tidyr::pivot_wider(id_cols = SampleID,names_from = Supergroup ,
                     values_from = Sum_RawCounts,
                     values_fill = min_raw_count) %>%
  data.frame() %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
  relocate(SampleID,sort(names(.))) %>% 
  arrange(SampleID) %>% 
  select(-SampleID) %>% 
  vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()


## So here we have two different list of distance matrices. 
saveRDS(list_dist_matrices,file = paste0(savingdir,'/','list_dist_matrices'))
saveRDS(list_dist_matrices_normalized,file = paste0(savingdir,'/','list_dist_matrices_normalized'))
saveRDS(biotic_list_dist,file = paste0(savingdir,'/','biotic_list_dist'))
saveRDS(biotic_list_dist_at_least_two,file = paste0(savingdir,'/','biotic_list_dist_at_least_two'))

##################################################################################################
# MDS_Coordinates ################################################################################
##################################################################################################

### Future work
# cmdscale(d = as.dist(justin_list$dist_matrices$biotic),k = 3,eig = T,add=T)
# pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)

list_mds_coord <- list(
  ASV = cmdscale(biotic_list_dist$ASV,k = 3),
  Species = cmdscale(biotic_list_dist$Species,k = 3),
  Genus = cmdscale(biotic_list_dist$Genus,k = 3),
  Family = cmdscale(biotic_list_dist$Family,k = 3),
  Order = cmdscale(biotic_list_dist$Order,k = 3),
  Class = cmdscale(biotic_list_dist$Class,k = 3),
  Division = cmdscale(biotic_list_dist$Division,k = 3),
  Supergroup = cmdscale(biotic_list_dist$Supergroup,k = 3)
)

saveRDS(list_mds_coord,file = paste0(savingdir,'/','list_mds_coord'))

##################################################################################################
# Taxonomic ################################################################################
##################################################################################################

####-TODO

