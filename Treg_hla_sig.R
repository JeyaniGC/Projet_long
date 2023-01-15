rm(list = ls())
setwd("~/M2_BI/projet_long_Jeyani")

#Load packages
suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(readr)
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# Load dataset
hla_genes <- read_csv("data/HLA_genes.csv")
hla_genotypes <- read_csv("data/HLA_genotypes.csv")
tcr_hv <- read_csv("data/TCR_data_HV.csv")
tcr_ra <- read_csv("data/TCR_data_RA.csv")
tcr_t1d <- read_csv("data/TCR_data_T1D.csv")
tcr_meta <- read_csv("data/TCR_metadata.csv")
tcr_sig <- read_csv("data/TCR_signature.csv")

# For Treg analysis
# Analysis RAvsHV
# In tcr_sig subset with Treg, RA and enriched data
ra_sig <- tcr_sig %>% 
  dplyr::filter(cell_subset == "CD4_Treg", disease == "RA",
                signature_type == "enriched")

# match tcr_ra patients with these signatures
tcr_ra_sig <- tcr_ra %>% 
  dplyr::filter(cdr3aa %in% ra_sig$cdr3aa)

# match tcr_hv patients with these signatures
tcr_hv_sig <- tcr_hv %>% 
  dplyr::filter(cdr3aa %in% ra_sig$cdr3aa)

# Check that they do not have the same id
length(intersect(tcr_hv_sig$sample_id, tcr_ra_sig$sample_id))

# Merge 
ra_hv_sig <- bind_rows(tcr_ra_sig, tcr_hv_sig)

# Subset sample_id and subject_id
meta_sig <- tcr_meta %>% 
  dplyr::select(subject_id, sample_id)

ra_hv_sig <- left_join(ra_hv_sig, meta_sig, by= "sample_id")

# Keep only patients in hla_genotypes
ra_hv_sig <- ra_hv_sig %>% 
  dplyr::filter(subject_id %in% hla_genotypes$subject_id) 

hm_RAvsHV <- data.frame(table(ra_hv_sig$cdr3aa, ra_hv_sig$subject_id))
hm_1_RAvsHV <- reshape(hm_RAvsHV,direction="wide",timevar="Var2",idvar="Var1") %>% 
  column_to_rownames("Var1")

names(hm_1_RAvsHV) <- gsub(x = names(hm_1_RAvsHV), pattern = "Freq\\.", replacement = "")

# Verification 
test_1 <- ra_hv_sig %>% 
  filter(subject_id == "FD47") %>% 
  dplyr::select(cdr3aa) %>% 
  table() %>% 
  as.data.frame() %>% 
  column_to_rownames(var="cdr3aa")

test_1_bis <- hm_1_RAvsHV %>% 
  dplyr::select(FD47) %>% 
  dplyr::filter(FD47 != 0)

identical(rownames(test_1), rownames(test_1_bis))
all(test_1$Freq == test_1_bis$FD47)

test_2 <- ra_hv_sig %>% 
  dplyr::filter(subject_id == "32D0") %>% 
  dplyr::select(cdr3aa) %>% 
  table() %>% 
  as.data.frame() %>% 
  column_to_rownames(var="cdr3aa")

# if suject_id start with a number
subject <- rownames(hm_1_RAvsHV)
filter <- subject[hm_1_RAvsHV$"32D0" != 0]
test_2_bis <- hm_1_RAvsHV[filter,]["32D0"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$`32D0`)

rm(test_1, test_1_bis, test_2, test_2_bis,
   subject,filter)

# Change values greater than 1 to 1 
for(i in 1:ncol(hm_1_RAvsHV)) {       
  hm_1_RAvsHV[,i] <- ifelse(hm_1_RAvsHV[,i] > 1, 1, hm_1_RAvsHV[, i])
}

hm_2_RAvsHV <- as.matrix(hm_1_RAvsHV)

hla_RAvsHV_label <- hla_genotypes %>% 
  dplyr::filter(subject_id %in% as.character(hm_RAvsHV$Var2)) %>% 
  dplyr::select(subject_id, DRB11, DPA11) %>% 
  arrange(subject_id)

disease_RAvsHV_label <- ra_hv_sig %>% 
  dplyr::filter(subject_id %in% as.character(hm_RAvsHV$Var2)) %>% 
  dplyr::select(subject_id, disease) %>% 
  unique() %>% 
  arrange(subject_id)

length(intersect(hla_RAvsHV_label$subject_id, disease_RAvsHV_label$subject_id))

# Heatmap de RAvsHV

ha <-  HeatmapAnnotation(disease= disease_RAvsHV_label$disease, 
                         DPA11 = hla_RAvsHV_label$DPA11, 
                         DRB11 = hla_RAvsHV_label$DRB11, 
                         col = list(disease = c("RA" = "red", "HV" = "blue"), 
                                    DPA11 = c("DPA1*01:03:01" = "sea green", "DPA1*02:01:01" = "coral")))


pdf("Heatmap_RAvsHV_Treg.pdf",width =50, height =50)
hm_plot_RAvsHV <-Heatmap(hm_2_RAvsHV,
                         name="Signature",
                         show_column_names=FALSE,
                         top_annotation=ha,
                         col = colorRamp2(c(0, 1), c("azure3", "darkcyan")),
                         show_row_dend=TRUE,
                         show_column_dend=TRUE,
                         cluster_rows=TRUE,
                         cluster_columns=TRUE, 
                         column_km = 2,
                         row_km = 1,
                         row_names_gp = gpar(fontsize = 10))
hm_plot_RAvsHV
dev.off()

disease_RA <- ra_hv_sig %>% 
  dplyr::filter(subject_id %in% as.character(hm_RAvsHV$Var2)) %>% 
  dplyr::filter(disease == "RA") %>% 
  unique() %>% 
  arrange(subject_id)

hla_RA <- hla_genotypes %>% 
  dplyr::filter(subject_id %in% disease_RA$subject_id) %>% 
  dplyr::select(DPA11, DRB11)

# Analysis T1DvsHV
t1d_sig <- tcr_sig %>% 
  dplyr::filter(cell_subset == "CD4_Treg", disease == "T1D",
                signature_type == "enriched")

tcr_t1d_sig <- tcr_t1d %>% 
  dplyr::filter(cdr3aa %in% t1d_sig$cdr3aa)

hv_t1d_sig <- tcr_hv %>% 
  dplyr::filter(cdr3aa %in% t1d_sig$cdr3aa)

length(intersect(hv_t1d_sig$sample_id, tcr_t1d_sig$sample_id))

t1d_hv_sig <- bind_rows(tcr_t1d_sig, hv_t1d_sig)

meta_sig <- tcr_meta %>% 
  dplyr::select(subject_id, sample_id)

t1d_hv_sig <- left_join(t1d_hv_sig, meta_sig, by= "sample_id")

t1d_hv_sig <- t1d_hv_sig %>% 
  dplyr::filter(subject_id %in% hla_genotypes$subject_id) 

hm_T1DvsHV <- data.frame(table(t1d_hv_sig$cdr3aa, t1d_hv_sig$subject_id))
hm_1_T1DvsHV <- reshape(hm_T1DvsHV,direction="wide",timevar="Var2",idvar="Var1") %>% 
  column_to_rownames("Var1")

names(hm_1_T1DvsHV) <- gsub(x = names(hm_1_T1DvsHV), pattern = "Freq\\.", replacement = "")

test_1 <- t1d_hv_sig %>% 
  filter(subject_id == "FFC0") %>% 
  dplyr::select(cdr3aa) %>% 
  table() %>% 
  as.data.frame() %>% 
  column_to_rownames(var="cdr3aa")

test_1_bis <- hm_1_T1DvsHV %>% 
  dplyr::select(FFC0) %>% 
  dplyr::filter(FFC0 != 0)

identical(rownames(test_1), rownames(test_1_bis))
all(test_1$Freq == test_1_bis$FD47)

test_2 <- t1d_hv_sig %>% 
  dplyr::filter(subject_id == "71B1") %>% 
  dplyr::select(cdr3aa) %>% 
  table() %>% 
  as.data.frame() %>% 
  column_to_rownames(var="cdr3aa")

subject <- rownames(hm_1_T1DvsHV)
filter <- subject[hm_1_T1DvsHV$"71B1" != 0]
test_2_bis <- hm_1_T1DvsHV[filter,]["71B1"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$`71B1`)

for(i in 1:ncol(hm_1_T1DvsHV)) {       
  hm_1_T1DvsHV[,i] <- ifelse(hm_1_T1DvsHV[,i] > 1, 1, hm_1_T1DvsHV[, i])
}

hm_2_T1DvsHV <- as.matrix(hm_1_T1DvsHV)

hla_T1D_label <- hla_genotypes %>% 
  dplyr::filter(subject_id %in% as.character(hm_T1DvsHV$Var2)) %>% 
  dplyr::select(subject_id, DRB11, DPA11) %>% 
  arrange(subject_id)

disease_T1D_label <- t1d_hv_sig %>% 
  dplyr::filter(subject_id %in% as.character(hm_T1DvsHV$Var2)) %>% 
  dplyr::select(subject_id, disease) %>% 
  unique() %>% 
  arrange(subject_id)

length(intersect(hla_T1D_label$subject_id, disease_T1D_label$subject_id))

# Heatmap de T1DvsHV

ha_T1DvsHV <-  HeatmapAnnotation(disease= disease_T1D_label$disease, 
                                 DPA11 = hla_T1D_label$DPA11, 
                                 DRB11 = hla_T1D_label$DRB11, 
                                 col = list(disease = c("T1D" = "red", "HV" = "blue"), 
                                            control = c("DPA1*01:03:01" = "sea green", "DPA1*02:01:01" = "coral",
                                                        "DPA1*02:01:02" = "purple3")))


pdf("Heatmap_T1DvsHV_Treg.pdf",width =50, height =50)
hm_plot_T1DvsHV <-Heatmap(hm_2_T1DvsHV,
                          name="Signature",
                          show_column_names=FALSE,
                          top_annotation=ha_T1DvsHV,
                          col = colorRamp2(c(0, 1), c("azure3", "darkcyan")),
                          show_row_dend=TRUE,
                          show_column_dend=TRUE,
                          cluster_rows=TRUE,
                          cluster_columns=TRUE, 
                          column_km = 2,
                          row_km = 1,
                          row_names_gp = gpar(fontsize = 10))
hm_plot_T1DvsHV
dev.off()

disease_T1D <- t1d_hv_sig %>% 
  dplyr::filter(subject_id %in% as.character(hm_T1DvsHV$Var2)) %>% 
  dplyr::filter(disease == "T1D") %>% 
  unique() %>% 
  arrange(subject_id)

hla_T1D <- hla_genotypes %>% 
  dplyr::filter(subject_id %in% disease_T1D$subject_id) %>% 
  dplyr::select(DPA11, DRB11)





