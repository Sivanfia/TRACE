# ------------------------------------------------------------
# Purpose: Reproduce the GTEx subset used in the paper from
#          GTEx bulk tissue expression (TPM) files.
# Inputs:
#   - module_137.csv (columns: ENTREZID, Name, P.value, Description)
#   - GTEx TPM by-tissue files: gene_tpm_2017-06-05_v8_<tissue>.gct.gz
#       Download page: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
#       Put in the folder: /Linear/realdata/
# Outputs:
#   - GTEx.RData  (contains: data_all, tissue_n, tissue_cons_id)
# ------------------------------------------------------------

library(caret)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(here)

module_137<-read.csv(here("Linear","realdata","module_137.csv"),header = TRUE)
colnames(module_137)<-c("ENTREZID","Name","P.value","Description")

symbol<-bitr(module_137[,1], fromType = "ENTREZID",
             toType = "SYMBOL", OrgDb = "org.Hs.eg.db",drop = FALSE)

module<-merge(module_137,symbol,by="ENTREZID")


tissue.names <- c("adipose_subcutaneous", "adipose_visceral_omentum", "adrenal_gland", 
                 "artery_aorta", "artery_coronary", "artery_tibial", 
                 "brain_amygdala", "brain_anterior_cingulate_cortex_ba24",
                 "brain_caudate_basal_ganglia", "brain_cerebellar_hemisphere",
                 "brain_cerebellum", "brain_cortex", "brain_frontal_cortex_ba9",
                 "brain_hippocampus", "brain_hypothalamus", 
                 "brain_nucleus_accumbens_basal_ganglia", "brain_putamen_basal_ganglia",
                 "brain_spinal_cord_cervical_c-1", "brain_substantia_nigra",
                 "breast_mammary_tissue", "cells_cultured_fibroblasts",
                 "cells_ebv-transformed_lymphocytes", "colon_sigmoid", "colon_transverse",
                 "esophagus_gastroesophageal_junction", "esophagus_mucosa",
                 "esophagus_muscularis", "heart_atrial_appendage", 
                 "heart_left_ventricle", "kidney_cortex", "liver", 
                 "lung", "minor_salivary_gland", "muscle_skeletal", "nerve_tibial", 
                 "ovary", "pancreas",  "pituitary", "prostate", 
                 "skin_not_sun_exposed_suprapubic", "skin_sun_exposed_lower_leg", 
                 "small_intestine_terminal_ileum", "spleen", "stomach", "testis", 
                 "thyroid", "uterus", "vagina", "whole_blood")

data_all <- data.frame(matrix(ncol = 1620, nrow = 0))
tissue.n <- c()
tissue.cons.id <- list()
k <- 1
for (tissue in tissue.names){
  loc.tissue.temp <- c(here("Linear","realdata"), "/GTEx_TPM_by_tissue/gene_tpm_2017-06-05_v8_",tissue,".gct.gz")
  loc.tissue <- paste(loc.tissue.temp, collapse = "")
  x<-read.table(gzfile(loc.tissue),
               skip = 2, header = TRUE, sep = "\t")
  data.gene<-x[which(x$Description %in% module$SYMBOL),]
  data <- as.data.frame(t(data.gene[,4:ncol(x)]))
  colnames(data)<- data.gene$Description
  #### standardization ####
  data <- predict(preProcess(data,method = c("center", "scale")),data)
  data$tissue<-tissue
  tissue.n<-c(tissue.n,ncol(x)-3)
  
  gene.cons.id <- which(apply(data[,1:(ncol(data)-1)],2,sd)==0)
  tissue.cons.id[[k]] <- gene.cons.id
  k<-k+1
  data_all<-rbind(data_all,data)
}
## two genes expressing in constant levels are excluded:
## DEFB104A (637) and DAZ1 (1619)
cons.id <- c(637,1619)
data_all<-data_all[,-cons.id]


keep_vars <- c("data_all", "tissue.n", "tissue.names")
rm(list = setdiff(ls(), keep_vars))

save(list = ls(), file = here("Linear", "realdata", "GTEx", "GTEx.RData"))

