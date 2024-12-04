############################################################################
#########Filter samples and Merge the corresponding filtered 16S and ITS datasets
###########################################################################
library(readxl)
library(here)
library(tidyverse)
library(dplyr)



#######Read the data (didn't consider the prune OTU dataset here)------------------
bact.n.otu = read_xlsx(here("PMI/data/TOX2_16S_TSS_meta_OTUtable.xlsx"))
bact.n.class = read_xlsx(here("PMI/data/TOX2_16S_TSS_meta_Class_table.xlsx"))
bact.n.order = read_xlsx(here("PMI/data/TOX2_16S_TSS_meta_Order_table.xlsx"))
bact.n.phylum = read_xlsx(here("PMI/data/TOX2_16S_TSS_meta_Phylum_table.xlsx"))



ITS.n.otu = read_xlsx(here("PMI/data/TOX2_ITS_TSS_meta_OTU_table.xlsx"))
ITS.n.class = read_xlsx(here("PMI/data/TOX2_ITS_TSS_meta_Class_table.xlsx"))
ITS.n.order = read_xlsx(here("PMI/data/TOX2_ITS_TSS_meta_Order_table.xlsx"))
ITS.n.phylum = read_xlsx(here("PMI/data/TOX2_ITS_TSS_meta_Phylum_table.xlsx"))




list_of_df <- list(bact.n.otu, bact.n.class, bact.n.order, bact.n.phylum, ITS.n.otu, ITS.n.class, ITS.n.order, ITS.n.phylum)
names(list_of_df) <- c("bact.n.otu", "bact.n.class", "bact.n.order", "bact.n.phylum", "ITS.n.otu", "ITS.n.class", "ITS.n.order", "ITS.n.phylum")



# calculat the number of features(Otu,ITS) in each dataframe
get_n_features <- function(df) {
  
  n_features <- df %>%
    select(matches(c("Otu", "ITS"))) %>%
    summarise(features = ncol(.)) %>%
    pull()
  
  return(n_features)
}

n_features <- list_of_df %>% map_df(get_n_features)




#######Rename feature names based on taxtable-------------------------
bact.otu.tax = read_xlsx(here("PMI/data/TOX2_16S_TSS_OTU_taxtable.xlsx"))
bact.class.tax = read_xlsx(here("PMI/data/TOX2_16S_TSS_Class_taxtable.xlsx"))
bact.order.tax = read_xlsx(here("PMI/data/TOX2_16S_TSS_Order_taxtable.xlsx"))
bact.phylum.tax = read_xlsx(here("PMI/data/TOX2_16S_TSS_Phylum_taxtable.xlsx"))


TSS.otu.tax = read_xlsx(here("PMI/data/TOX2_ITS_TSS_OTU_taxtable.xlsx"))
TSS.class.tax = read_xlsx(here("PMI/data/TOX2_ITS_TSS_Class_taxtable.xlsx"))
TSS.order.tax = read_xlsx(here("PMI/data/TOX2_ITS_TSS_Order_taxtable.xlsx"))
TSS.phylum.tax = read_xlsx(here("PMI/data/TOX2_ITS_TSS_Phylum_taxtable.xlsx"))

# rename colnames of bact.n.class, bact.n.order, bact.n.phylum,
rename_tax_fun <- function(df, mapping_df, tax) {
  # Ensure the target_column exists in mapping_df
  if (!tax %in% colnames(mapping_df)) {
    stop("The specified target column does not exist in the mapping data frame.")
  }
  
  # Create the mapping vector using the specified column
  mapping_vector <- setNames(mapping_df$OTU,mapping_df[[tax]])
  # Make the new names unique by appending suffixes to duplicates
  unique_mapping_vector <- make.unique(names(mapping_vector))
  # Update the mapping vector with unique names
  names(mapping_vector) <- unique_mapping_vector
  
  # Rename the columns in df using dplyr
  df <- df %>% rename(any_of(mapping_vector))
  
  return(df)
}

colnames(bact.n.phylum)
colnames(ITS.n.phylum)
bact.n.class <- rename_tax_fun(bact.n.class,bact.class.tax,'Class')
bact.n.order <- rename_tax_fun(bact.n.order,bact.order.tax,'Order')
bact.n.phylum  <-  rename_tax_fun(bact.n.phylum,bact.phylum.tax,'Phylum')


ITS.n.class = rename_tax_fun(ITS.n.class,TSS.class.tax,'Class')
ITS.n.order = rename_tax_fun(ITS.n.order,TSS.order.tax,'Order')
ITS.n.phylum = rename_tax_fun(ITS.n.phylum,TSS.phylum.tax,'Phylum')
colnames(bact.n.phylum)
colnames(ITS.n.phylum)

colnames(bact.n.class)
colnames(ITS.n.class)


list_of_df <- list(bact.n.otu, bact.n.class, bact.n.order, bact.n.phylum, ITS.n.otu, ITS.n.class, ITS.n.order, ITS.n.phylum)
names(list_of_df) <- c("bact.n.otu", "bact.n.class", "bact.n.order", "bact.n.phylum", "ITS.n.otu", "ITS.n.class", "ITS.n.order", "ITS.n.phylum")









#######Filtering-------------------------
filter_AHD10 <- function(df, adh_cutoff) {
  df_filter <- df %>% filter(ADH_10_actual < adh_cutoff)
  
  return(df_filter)
}

df_list_5000 <- list_of_df %>% map(filter_AHD10, 5000)




#######merge-----------------------
merge_16S_ITS_df <- function(bact_its_pair, df_list) {
  bact_df_name <- as.character(bact_its_pair[1])
  its_df_name <- as.character(bact_its_pair[2])
  
  
  ITS_df_sub <- df_list[[its_df_name]] %>% select(-Sample , -Type,-Donor,-Trt,-Sample_Type,-ADH,-Actual_ADH,-Timepoint,-Per_ADH,-Total_ADH,-Total_Days,-Season,-Samp_Season,-Sector,-Age,-Sex,-Stature,-Weight,-BMI,-BMI_level,-Temperature,-ADH_10_actual,-Diabetes,-Cancer,-Cardio,-Respiratory,-Neuro,-Pneumonia,-moisture,-pH,-pH_LRR,-EC,-EC_LRR,-BG_avg,-BG_LRR,-NAG_avg,-NAG_LRR,-PHOS_avg,-PHOS_LRR,-LAP_avg,-LAP_LRR,-Evolved_CO2,-LRR_resp)
  
  merged_df <- df_list[[bact_df_name]] %>% left_join(ITS_df_sub, by = "group")
  return(merged_df)
}



merge_list <- list("bact.ITS.n.otu" = c("bact.n.otu", "ITS.n.otu"),
                   "bact.ITS.n.phylum" = c("bact.n.phylum", "ITS.n.phylum"),
                   "bact.ITS.n.class" = c("bact.n.class", "ITS.n.class"),
                   "bact.ITS.n.order" = c("bact.n.order", "ITS.n.order"))

merged_df_5000 <- merge_list %>% map(merge_16S_ITS_df, df_list_5000)

merged_df_5000 %>% map_df(get_n_features)



#######combie bact,ITS and merged dataframe-----------------------
df_list_5000 <- append(df_list_5000, merged_df_5000, after = length(df_list_5000))
n_features <- df_list_5000 %>% map_df(get_n_features)
df_list_5000
colnames(df_list_5000$bact.n.otu)


####### No env: keep only observations from donor assigned to the training set and the response variable-------
filter_noenv <- function(df) {
  df_filtered <- df %>% select(-group,-Sample,-Type,-Donor,-Trt,-Sample_Type,-ADH,-Actual_ADH,
                                  -Timepoint,-Per_ADH,-Total_ADH,-Total_Days,-Season,-Samp_Season,
                                  -Sector,-Age,-Sex,-Stature,-Weight,-BMI,-BMI_level,-Temperature,
                                  -Diabetes,-Cancer,-Cardio,-Respiratory,-Neuro,-Pneumonia,-moisture,
                                  -pH,-pH_LRR,-EC,-EC_LRR,-BG_avg,-BG_LRR,-NAG_avg,-NAG_LRR,-PHOS_avg,
                                  -PHOS_LRR,-LAP_avg,-LAP_LRR,-Evolved_CO2,-LRR_resp)
  return(df_filtered)
}



####### env:Keep response variable (ADH_actual) and numeric environmental data-------------------
filter_yesenv <- function(df) {
  df_filtered <- df %>% select(-group,-Sample,-Type,-Donor,-Trt,-Sample_Type,-ADH,-Actual_ADH,
                               -Timepoint,-Per_ADH,-Total_ADH,-Total_Days,-Season,-Samp_Season,
                               -Sector,-Age,-Sex,-Stature,-Weight, -BMI, -BMI_level, 
                               -Diabetes,-Cancer,-Cardio,-Respiratory,-Neuro,-Pneumonia,
                               -pH,-EC,-BG_avg,-NAG_avg,-PHOS_avg,-LAP_avg,-Evolved_CO2,-LRR_resp)
  return(df_filtered)
}


df_list_5000_noenv <- lapply(df_list_5000, filter_noenv)
df_list_5000_env <- lapply(df_list_5000, filter_yesenv)




for (i in 1:length(df_list_5000_noenv)) {
  print(paste0("data_final/", names(df_list_5000_noenv)[[i]], ".noenv.csv")
  )
}


#######save the data: without env variables and with env variables, save them in the SelectMicro folder
for (i in 1:length(df_list_5000_noenv)) {
  write.csv(df_list_5000_noenv[[i]], paste0(here("SelectMicro_24new/Analysis/PMI/data/count_table"), names(df_list_5000_noenv)[[i]], ".noenv.csv"), row.names = FALSE)
}
for (i in 1:length(df_list_5000_env)) {
  write.csv(df_list_5000_env[[i]], paste0(here("SelectMicro_24new/Analysis/PMI/data/count_table"), names(df_list_5000_env)[[i]], ".env.csv"), row.names = FALSE)
}
