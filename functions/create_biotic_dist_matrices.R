create_biotic_dist_matrices <- function(
    grump_df_aux,
    min_raw_count=0.01,
    vet_asv_at_least_two_samples){
  
  biotic_list_dist <- list()
  biotic_list_dist_at_least_two <- list()
  grump_df_aux = grump_df_aux %>% mutate(Division=ifelse(Division=='','NotAssigned',Division))
  grump_df_aux = grump_df_aux %>% mutate(Supergroup=ifelse(Supergroup=='','NotAssigned',Supergroup))
  
  ##########################################################################################
  # ASV ####################################################################################
  ##########################################################################################
  cat('ASV \n')
  
  biotic_list_dist$ASV <- grump_df_aux %>%
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
  
  biotic_list_dist$ASV_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$ASV <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$ASV_normalized <- grump_df_aux %>%
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
  cat('Species \n')
  
  biotic_list_dist$Species <- grump_df_aux %>%
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
  
  biotic_list_dist$Species_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Species <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Species_normalized <- grump_df_aux %>%
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
  cat('Genus \n')
  
  
  biotic_list_dist$Genus <- grump_df_aux %>%
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
  
  biotic_list_dist$Genus_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Genus <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Genus_normalized <- grump_df_aux %>%
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
  cat('Family \n')
  
  biotic_list_dist$Family <- grump_df_aux %>%
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
  
  biotic_list_dist$Family_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Family <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Family_normalized <- grump_df_aux %>%
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
  cat('Order \n')
  
  biotic_list_dist$Order <- grump_df_aux %>%
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
  
  biotic_list_dist$Order_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Order <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Order_normalized <- grump_df_aux %>%
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
  cat('Class \n')
  
  biotic_list_dist$Class <- grump_df_aux %>%
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
  
  biotic_list_dist$Class_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Class <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Class_normalized <- grump_df_aux %>%
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
  # Phylum ################################################################################
  ##########################################################################################
  cat('Phylum \n')
  
  biotic_list_dist$Phylum <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,Phylum,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Phylum) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Phylum ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()
  
  biotic_list_dist$Phylum_normalized <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,Phylum,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Phylum) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Phylum ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()
  
  
  biotic_list_dist_at_least_two$Phylum <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
    select(SampleID,Phylum,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Phylum) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Phylum ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()
  
  biotic_list_dist_at_least_two$Phylum_normalized <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
    select(SampleID,Phylum,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Phylum) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Phylum ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()
  
  ##########################################################################################
  # Domain ################################################################################
  ##########################################################################################
  cat('Domain \n')
  
  biotic_list_dist$Domain <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,Domain,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Domain) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Domain ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()
  
  biotic_list_dist$Domain_normalized <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,Domain,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Domain) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Domain ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() %>% normalizeMatrix()
  
  
  biotic_list_dist_at_least_two$Domain <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
    select(SampleID,Domain,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Domain) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Domain ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()
  
  biotic_list_dist_at_least_two$Domain_normalized <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    filter(ASV_hash%in%vet_asv_at_least_two_samples) %>% 
    select(SampleID,Domain,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Domain) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,names_from = Domain ,
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
  cat('Division \n')
  
  biotic_list_dist$Division <- grump_df_aux %>%
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts + min_raw_count) %>% 
    select(SampleID,Division,Raw.Sequence.Counts) %>% distinct() %>% 
    group_by(SampleID,Division) %>% 
    summarise(Sum_RawCounts = sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    tidyr::pivot_wider(id_cols = SampleID,
                       names_from = Division ,
                       values_from = Sum_RawCounts,
                       values_fill = min_raw_count) %>%
    data.frame() %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    relocate(SampleID,sort(names(.))) %>% 
    arrange(SampleID) %>% 
    select(-SampleID) %>% 
    vegan::vegdist(method = 'aitchison') %>% as.matrix() #%>% normalizeMatrix()
  
  biotic_list_dist$Division_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Division <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Division_normalized <- grump_df_aux %>%
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
  
  biotic_list_dist$Supergroup <- grump_df_aux %>%
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
  
  biotic_list_dist$Supergroup_normalized <- grump_df_aux %>%
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
  
  
  biotic_list_dist_at_least_two$Supergroup <- grump_df_aux %>%
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
  
  biotic_list_dist_at_least_two$Supergroup_normalized <- grump_df_aux %>%
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
  
  
  
  
  
  return(list(
    biotic_list_dist_all_ASVs=biotic_list_dist,
    biotic_list_dist_not_rare_ASVs=biotic_list_dist_at_least_two
  ))
  
}
