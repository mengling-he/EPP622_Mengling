set.seed(6)

library(readxl)
library(here)
library(phyloseq)
library(tidyverse)
library(writexl)

######################################
######## load metadata ###############
######################################
treatments <- read_excel(here("PMI/rawdata/treatments_pub.xlsx"))# metadata
treatments$ADH=as.numeric(treatments$ADH)
treatments$Type=as.factor(treatments$Type)
treatments$Trt=as.factor(treatments$Trt)
treatments$Sample_Type=as.factor(treatments$Sample_Type)
treatments$Timepoint=as.factor(treatments$Timepoint)
treatments$Season=as.factor(treatments$Season)
treatments$Sector=as.factor(treatments$Sector)
treatments$Sex=as.factor(treatments$Sex)
treatments$BMI_level = factor(treatments$BMI_level, ordered = TRUE, levels = c("Underweight","Normal","Overweight","Obese"))

TOXdata = read_excel(here("PMI/rawdata/TOX_data_all_pub.xlsx"), sheet = "Data")#other meta and environmental infor

TOXdata$Diabetes=as.factor(TOXdata$Diabetes)
TOXdata$Cancer=as.factor(TOXdata$Cancer)
TOXdata$Cardio=as.factor(TOXdata$Cardio)
TOXdata$Respiratory=as.factor(TOXdata$Respiratory)
TOXdata$Neuro=as.factor(TOXdata$Neuro)
TOXdata$Pneumonia=as.factor(TOXdata$Pneumonia)
TOXdata = TOXdata %>% select("Sample","Temperature","ADH_10_actual", "Diabetes", "Cancer","Cardio","Respiratory","Neuro","Pneumonia","Gravimetric Moisture","pH","pH_LRR","EC","EC_LRR","BG_avg","BG_LRR","NAG_avg","NAG_LRR","PHOS_avg","PHOS_LRR","LAP_avg", "LAP_LRR","Evolved_CO2","LRR_resp")
colnames(TOXdata)[10] = "moisture"

metadata = treatments %>%
  left_join(TOXdata, by = "Sample")

#prep metadata for phyloseq
map=metadata
head(map)
map=sample_data(map) #convert metadata into phyloseq format
rownames(map)=map$group #assign rownames to be group

levels=c("TOX001CON", "TOX001", "TOX001_1500", "TOX001_3000", "TOX001_4500", "TOX001_8000", "TOX001_11500", "TOX001_15500","TOX001DF", "TOX002CON", "TOX002", "TOX002_1500A", "TOX002_1500B", "TOX002_2750", "TOX002_3750","TOX002DF", "TOX003CON", "TOX003", "TOX003_1500", "TOX003_3000", "TOX003_4500","TOX003DF", "TOX004CON", "TOX004", "TOX004_1500", "TOX004_3000", "TOX004_3500","TOX004DF", "TOX005CON", "TOX005", "TOX005_1500", "TOX005_3000", "TOX005_4500", "TOX005_6500", "TOX005_8500","TOX005DF", "TOX006CON", "TOX006", "TOX006_1500", "TOX006_3000", "TOX006_4500","TOX006DF", "TOX007CON", "TOX007", "TOX007_1500", "TOX007_3000", "TOX007_4500", "TOX007_6000","TOX007DF", "TOX008CON", "TOX008", "TOX008_1500", "TOX008_3000", "TOX008_4500", "TOX008_6500","TOX008DF", "TOX009CON", "TOX009", "TOX009_1500", "TOX009_3000", "TOX009_4500","TOX009DF", "TOX010CON", "TOX010", "TOX010_1500", "TOX010_2500", "TOX010_4500", "TOX010_9000", "TOX010_13000", "TOX010_17500", "TOX010DF", "TOX011CON", "TOX011", "TOX011_1500", "TOX011_3000", "TOX011_3500", "TOX011_5500", "TOX011_7000", "TOX012CON", "TOX012", "TOX012_1500", "TOX012_3000", "TOX012_4000", "TOX012_6000","TOX012DF", "TOX013CON", "TOX013", "TOX013_1500", "TOX013_3000", "TOX013_4000", "TOX013_7500", "TOX013DF","TOX015CON","TOX015","TOX015_250","TOX015_500","TOX015_1000","TOX015_1500", "TOX016CON","TOX016","TOX016_500","TOX016_1000","TOX016_1500","TOX016DF","TOX017CON","TOX017","TOX017_1500","TOX017_3000","TOX017_4500","TOX017_8500","TOX017_13000","TOX017_18500","TOX017DF","TOX018CON","TOX018","TOX018_1500","TOX018_2000","TOX018_3000","TOX018_4000", "TOX019CON","TOX019","TOX019_1500","TOX019_3000","TOX019_3500","TOX019_4500","TOX019_5500","TOX019DF","TOX020CON","TOX020","TOX020_1500","TOX020_3500","TOX020_4500","TOX020_6000","TOX020DF")
map$group=ordered(map$group, levels=levels)

######################################
###   Where do we cut off ADH?     ###
######################################

dist <- metadata %>%
  filter(Trt == "donor") %>%
  ggplot() +
  aes(y = Donor, x = ADH_10_actual) +
  geom_point() +
  geom_vline(xintercept = 5000, linetype="dotted", 
                color = "blue", size=1.5) +
  labs(x =  bquote(paste("ADH (base 10", degree, " C)"))) +
  theme_bw()
  
tiff(here::here("PMI/figures/samples_by_donor_5000ADH.tiff"), units = "in", width = 7.5, height = 5.5, res = 300)
dist
dev.off()

#tiff(here::here("Mason_SoilMicrobe-PMI_PLoS_2024/figures/S1Fig.tiff"), units = "in", width = 7.5, height = 5.5, res = 300)
dist
dev.off()


#####################################################
#### obtain 16S OTU tables from phyloseq object #####
######################################################
mothur_data= import_mothur(mothur_shared_file = here("PMI/rawdata/TOX2_16S.shared"), mothur_constaxonomy_file = here("PMI/rawdata/TOX2_16S.taxonomy"))
moth_merge = merge_phyloseq(mothur_data, map) #merge mothurdata with the metadata
## Mothur is an open-source platform for analyzing 16S rRNA gene sequences, commonly used in microbiome research to assess microbial community structure and diversity.
colnames(tax_table(moth_merge)) #print column names of tax file
colnames(tax_table(moth_merge)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") #rename tax column names
moth_sub=subset_samples(moth_merge, Type == "sample") #remove check samples
moth_soil=subset_samples(moth_merge, Trt == "donor") #remove fluid samples and control soil samples
moth_soil=subset_samples(moth_soil, group != "TOX002_1500A") #remove TOX002-1500
moth_soil=subset_samples(moth_soil, group != "TOX002_1500B") #remove TOX002-1500

min(sample_sums(moth_soil)) #29510
sum(sample_sums(moth_soil)) # a total of 7858040 reads across all samples

bact_phylo = moth_soil

####### 16S - OTU level
#################################################
bact_phylo_n=transform_sample_counts(bact_phylo, function(x) {x/sum(x)*10000})#need to do total sum scaling: aka take relative abundance and then multiply by fixed library size of 10,000
bact_phylo_n=prune_taxa(taxa_sums(bact_phylo_n) > 10, bact_phylo_n)#remove normalized OTUs with less than 10 reads across all samples 
#note went from 50,720 taxa to 5,195 taxa




bact_phylo_n_5000 = subset_samples(bact_phylo_n, ADH <= 5000)

sum(sample_sums(bact_phylo_n_5000)) # a total of 745467.8 reads across all samples <= 5000 ADH
sample_data(bact_phylo_n_5000)$Sample_Type

#mean number of timepoints per donor # there are 19 donor with 78 samples with 745467.8 reads
sample_data(bact_phylo_n_5000) %>%
  data.frame()%>%
  count(Donor) %>%
  summarise(mean_timepoints = mean(n)) %>%
  pull(mean_timepoints)
  





bact_OTU_n=as(otu_table(bact_phylo_n), "matrix")
if(taxa_are_rows(bact_phylo_n)){bact_OTU_n=t(bact_OTU_n)}
bact_OTU_n_df=as.data.frame(bact_OTU_n)
bact_OTU_n_df$group = row.names(bact_OTU_n_df)

bact.n_meta_otu = metadata %>% left_join(bact_OTU_n_df, by = "group")%>%filter(Otu00001 != "NA")

bact.n_OTU_tax_table = as.data.frame(tax_table(bact_phylo_n))
bact.n_OTU_tax_table$OTU = row.names(bact.n_OTU_tax_table) #make a new row names "OTU" with otu names
#write_xlsx(bact.n_meta_otu, here("Mason_SoilMicrobe-PMI_PLoS_2024/data/TOX2_16S_TSS_meta_OTUtable.xlsx"))
#write_xlsx(bact.n_OTU_tax_table, here("Mason_SoilMicrobe-PMI_PLoS_2024/data/TOX2_16S_TSS_OTU_taxtable.xlsx")) 


##### OTU level - remove singletons and doubletons
bact_phylo_prune = prune_taxa(taxa_sums(bact_phylo) > 2, bact_phylo) 
bact_phylo_prune # this removed 26585 taxa

bact_phylo_prune_n = transform_sample_counts(bact_phylo_prune, function(x) {x/sum(x)*10000})#need to do total sum scaling: aka take relative abundance and then multiply by fixed library size of 10,000
bact_phylo_prune_n = prune_taxa(taxa_sums(bact_phylo_prune_n) > 10, bact_phylo_prune_n)#remove normalized OTUs with less than 10 reads across all samples 
#note went from 24135 taxa to 5,205 taxa

bact_OTU_prune_n = as(otu_table(bact_phylo_prune_n), "matrix")
if(taxa_are_rows(bact_phylo_prune_n)){bact_OTU_prune_n = t(bact_OTU_prune_n)}
bact_OTU_prune_n_df = as.data.frame(bact_OTU_prune_n)
bact_OTU_prune_n_df$group = row.names(bact_OTU_prune_n_df)

bact.prune.n_meta_otu = metadata %>% left_join(bact_OTU_prune_n_df, by = "group") %>% filter(Otu00001 != "NA")

bact.prune.n_OTU_tax_table = as.data.frame(tax_table(bact_phylo_prune_n))
bact.prune.n_OTU_tax_table$OTU = row.names(bact.prune.n_OTU_tax_table) #make a new row names "OTU" with otu names
write_xlsx(bact.prune.n_meta_otu, here("PMI/data/TOX2_16S_TSS_meta_prune_OTUtable.xlsx"))
write_xlsx(bact.prune.n_OTU_tax_table, here("PMI/data/TOX2_16S_TSS_prune_OTU_taxtable.xlsx")) 


####### 16S - Phylum level
################################
bact_phylum = bact_phylo %>%
  tax_glom(taxrank="Phylum")
bact_phylum_n=transform_sample_counts(bact_phylum, function(x) {x/sum(x)*10000})
bact_phylum_n=prune_taxa(taxa_sums(bact_phylum_n) > 10, bact_phylum_n)
#note: went from 46 to 35 taxa
bact_phylum_n_5000 = prune_samples(ADH <= 5000, bact_phylum_n)


bact_phylum_n_table=as(otu_table(bact_phylum_n), "matrix")
if(taxa_are_rows(bact_phylum_n)){bact_phylum_n_table=t(bact_phylum_n_table)}
bact_phylum_n_df=as.data.frame(bact_phylum_n_table)
bact_phylum_n_df$group = row.names(bact_phylum_n_df)

bact.n_meta_p = metadata %>% left_join(bact_phylum_n_df, by = "group")%>%filter(Otu00002 != "NA")

bact.n_p_tax_table = as.data.frame(tax_table(bact_phylum_n))
bact.n_p_tax_table$OTU = row.names(bact.n_p_tax_table) #make a new row names "OTU" with otu names
write_xlsx(bact.n_meta_p, here("PMI/data/TOX2_16S_TSS_meta_Phylum_table.xlsx") )
write_xlsx(bact.n_p_tax_table, here("PMI/data/TOX2_16S_TSS_Phylum_taxtable.xlsx") )

####### 16S - Class level
######################################
bact_class = bact_phylo %>%
  tax_glom(taxrank="Class")
bact_class_n=transform_sample_counts(bact_class, function(x) {x/sum(x)*10000})
bact_class_n=prune_taxa(taxa_sums(bact_class_n) > 10, bact_class_n)
#note: went from 144 to 111 taxa

bact_class_n_table=as(otu_table(bact_class_n), "matrix")
if(taxa_are_rows(bact_class_n)){bact_class_n_table=t(bact_class_n_table)}
bact_class_n_df=as.data.frame(bact_class_n_table)
bact_class_n_df$group = row.names(bact_class_n_df)

bact.n_meta_c = metadata %>% left_join(bact_class_n_df, by = "group")%>%filter(Otu00002 != "NA")

bact.n_c_tax_table = as.data.frame(tax_table(bact_class_n))
bact.n_c_tax_table$OTU = row.names(bact.n_c_tax_table) #make a new row names "OTU" with otu names
write_xlsx(bact.n_meta_c, here("PMI/data/TOX2_16S_TSS_meta_Class_table.xlsx") )
write_xlsx(bact.n_c_tax_table, here("PMI/data/TOX2_16S_TSS_Class_taxtable.xlsx"))

####### 16S - Order level
######################################
bact_order = bact_phylo %>%
  tax_glom(taxrank="Order")
bact_order_n=transform_sample_counts(bact_order, function(x) {x/sum(x)*10000})
bact_order_n=prune_taxa(taxa_sums(bact_order_n) > 10, bact_order_n)
#note: went from 377 to 264 taxa

bact_order_n_table=as(otu_table(bact_order_n), "matrix")
if(taxa_are_rows(bact_order_n)){bact_order_n_table=t(bact_order_n_table)}
bact_order_n_df=as.data.frame(bact_order_n_table)
bact_order_n_df$group = row.names(bact_order_n_df)

bact.n_meta_o = metadata %>% left_join(bact_order_n_df, by = "group")%>%filter(Otu00002 != "NA")

bact.n_o_tax_table = as.data.frame(tax_table(bact_order_n))
bact.n_o_tax_table$OTU = row.names(bact.n_o_tax_table) #make a new row names "OTU" with otu names
write_xlsx(bact.n_meta_o, here("PMI/data/TOX2_16S_TSS_meta_Order_table.xlsx") )
write_xlsx(bact.n_o_tax_table, here("PMI/data/TOX2_16S_TSS_Order_taxtable.xlsx"))











