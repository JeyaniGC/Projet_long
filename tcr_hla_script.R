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

# Pour les Teff
# select alleles HLA class II DRB11, control : DA11 

# Analyse RAvsHV
# tcr_sig prendre les Teff, RA, enriched
ra_sig <- tcr_sig %>% 
  dplyr::filter(cell_subset == "CD4_Teff", disease == "RA",
                signature_type == "enriched")

# faire correspondre les patients tcr_ra ayant ces sig
tcr_ra_sig <- tcr_ra %>% 
  dplyr::filter(cdr3aa %in% ra_sig$cdr3aa)

# faire la même chose avec patient HV
tcr_hv_sig <- tcr_hv %>% 
  dplyr::filter(cdr3aa %in% ra_sig$cdr3aa)

# voir si tcr_hv_sig et tcr_ra_sig ont les mêmes patients 
length(intersect(tcr_hv_sig$sample_id, tcr_ra_sig$sample_id))
# on voit que les patients sont bien uniques en fonction de la maladie

# merge les deux dtf pour la hm 
ra_hv_sig <- bind_rows(tcr_ra_sig, tcr_hv_sig)

# faire le lien sample_id et subject_id
meta_sig <- tcr_meta %>% 
  dplyr::select(subject_id, sample_id)

ra_hv_sig <- left_join(ra_hv_sig, meta_sig, by= "sample_id")

# récupérer les patients uniquemen dans hla_genes
ra_hv_sig <- ra_hv_sig %>% 
  dplyr::filter(subject_id %in% hla_genotypes$subject_id) 

hm_RAvsHV <- data.frame(table(ra_hv_sig$cdr3aa, ra_hv_sig$subject_id))
hm_1_RAvsHV <- reshape(hm_RAvsHV,direction="wide",timevar="Var2",idvar="Var1") %>% 
  column_to_rownames("Var1")

names(hm_1_RAvsHV) <- gsub(x = names(hm_1_RAvsHV), pattern = "Freq\\.", replacement = "")

# Verification des count sur 2 colonnes 

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

# cas si nom des gènes commence par des chiffres (ou si on aime pas dplyr) 
subject <- rownames(hm_1_RAvsHV)
filter <- subject[hm_1_RAvsHV$"32D0" != 0]
test_2_bis <- hm_1_RAvsHV[filter,]["32D0"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$`32D0`)


# modifier hm_1 
# changer toute les valeurs > 1 par des 1
for(i in 1:ncol(hm_1_RAvsHV)) {       
  hm_1_RAvsHV[,i] <- ifelse(hm_1_RAvsHV[,i] > 1, 1, hm_1_RAvsHV[, i])
}

# récupérer le nom des sample_id dans tcr_ra et tcr_hv 
#  -> les mettre dans une dtf et regarder les colnames et voir si 
# discriminer en fonction des patient ou des sig ?  

# faire la matrice pour la heatmap
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
# on a bien les mêmes patients dans les deux labels
# pour l'annotation 19 genotype pour DRB11 et 
# 2 genotypes pour DPA11

# Heatmap de RAvsHV

ha <-  HeatmapAnnotation(disease= disease_RAvsHV_label$disease, 
                         control = hla_RAvsHV_label$DPA11, 
                         target = hla_RAvsHV_label$DRB11, 
                         col = list(disease = c("RA" = "red", "HV" = "blue"), 
                                    control = c("DPA1*01:03:01" = "sea green", "DPA1*02:01:01" = "coral")))


pdf("Heatmap_RAvsHV.pdf",width =50, height =50)
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


# faire la même chose pour T1D
t1d_sig <- tcr_sig %>% 
  dplyr::filter(cell_subset == "CD4_Teff", disease == "T1D",
                signature_type == "enriched")

# faire correspondre les patients tcr_t1d ayant ces sig
tcr_t1d_sig <- tcr_t1d %>% 
  dplyr::filter(cdr3aa %in% t1d_sig$cdr3aa)

# faire la même chose avec patient HV
hv_t1d_sig <- tcr_hv %>% 
  dplyr::filter(cdr3aa %in% t1d_sig$cdr3aa)

# voir si tcr_hv_sig et tcr_ra_sig ont les mêmes patients 
length(intersect(hv_t1d_sig$sample_id, tcr_t1d_sig$sample_id))
# on voit que les patients sont bien uniques en fonction de la maladie

# merge les deux dtf pour la hm 
t1d_hv_sig <- bind_rows(tcr_t1d_sig, hv_t1d_sig)

# faire le lien sample_id et subject_id
meta_sig <- tcr_meta %>% 
  dplyr::select(subject_id, sample_id)

t1d_hv_sig <- left_join(t1d_hv_sig, meta_sig, by= "sample_id")

# récupérer les patients uniquement dans hla_genes
t1d_hv_sig <- t1d_hv_sig %>% 
  dplyr::filter(subject_id %in% hla_genotypes$subject_id) 

hm_T1DvsHV <- data.frame(table(t1d_hv_sig$cdr3aa, t1d_hv_sig$subject_id))
hm_1_T1DvsHV <- reshape(hm_T1DvsHV,direction="wide",timevar="Var2",idvar="Var1") %>% 
  column_to_rownames("Var1")

names(hm_1_T1DvsHV) <- gsub(x = names(hm_1_T1DvsHV), pattern = "Freq\\.", replacement = "")

# Verification des count sur 2 colonnes 

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

# cas si nom des gènes commence par des chiffres (ou si on aime pas dplyr) 
subject <- rownames(hm_1_T1DvsHV)
filter <- subject[hm_1_T1DvsHV$"71B1" != 0]
test_2_bis <- hm_1_T1DvsHV[filter,]["71B1"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$`71B1`)

# modifier hm_1 
# changer toute les valeurs > 1 par des 1
for(i in 1:ncol(hm_1_T1DvsHV)) {       
  hm_1_T1DvsHV[,i] <- ifelse(hm_1_T1DvsHV[,i] > 1, 1, hm_1_T1DvsHV[, i])
}

# faire la matrice pour la heatmap
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
# on a bien les mêmes patients dans les deux labels
# pour l'annotation 19 genotype pour DRB11 et 
# 2 genotypes pour DPA11

# Heatmap de T1DvsHV

ha_T1DvsHV <-  HeatmapAnnotation(disease= disease_T1D_label$disease, 
                         control = hla_T1D_label$DPA11, 
                         target = hla_T1D_label$DRB11, 
                         col = list(disease = c("T1D" = "red", "HV" = "blue"), 
                                    control = c("DPA1*01:03:01" = "sea green", "DPA1*02:01:01" = "coral",
                                                "DPA1*02:01:02" = "purple3")))


pdf("Heatmap_T1DvsHV.pdf",width =50, height =50)
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


# Pour les Treg

# Visualiser distribution de HLA
# PCA 

train_pca <- pca(hla_genes[, -1])
plotIndiv(train_pca, var.names = F)
plotVar(train_pca, var.names = F)

# lier data de tcr et de hla
## le faire par maladie ou HV

# T1D
tcr_t1d_link <- tcr_t1d %>% 
  select(sample_id, cdr3aa)

tcr_meta_link <- tcr_meta %>% 
  select(sample_id, subject_id)

tcr_meta_t1d_link <- tcr_t1d_link %>% 
  left_join(tcr_meta_link, by = "sample_id")

tcr_hla_t1d_link <- tcr_meta_t1d_link %>% 
  left_join(hla_genes, by = "subject_id")


tcr_hla_r1d_link_clean <- tcr_hla_t1d_link %>% 
  na.omit()

# RA
tcr_ra_link <- tcr_ra %>% 
  select(sample_id, cdr3aa)

tcr_meta_ra_link <- tcr_ra_link %>% 
  left_join(tcr_meta_link, by = "sample_id")

tcr_hla_ra_link <- tcr_meta_ra_link %>% 
  left_join(hla_genes, by = "subject_id")


tcr_hla_ra_link_clean <- tcr_hla_ra_link %>% 
  na.omit()

# HV
tcr_hv_link <- tcr_hv %>% 
  select(sample_id, cdr3aa)

tcr_meta_hv_link <- tcr_hv_link %>% 
  left_join(tcr_meta_link, by = "sample_id")

tcr_hla_hv_link <- tcr_meta_hv_link %>% 
  left_join(hla_genes, by = "subject_id")


tcr_hla_link_clean <- tcr_hla_link %>% 
  na.omit()
