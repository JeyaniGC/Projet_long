rm(list = ls())
setwd("~/M2_BI/projet_long_Jeyani")

#Load packages
suppressPackageStartupMessages(library(tidyverse))
library(dplyr)
library(readr)
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(ComplexHeatmap))

# Load dataset
hla_genes <- read_csv("data/HLA_genes.csv")
hla_genotypes <- read_csv("data/HLA_genotypes.csv")
tcr_hv <- read_csv("data/TCR_data_HV.csv")
tcr_ra <- read_csv("data/TCR_data_RA.csv")
tcr_t1d <- read_csv("data/TCR_data_T1D.csv")
tcr_meta <- read_csv("data/TCR_metadata.csv")
tcr_sig <- read_csv("data/TCR_signature.csv")

# select alleles HLA class II DRB11, control : DA11 

# Analyse RAvsHV
# dans tcr_sig prendre les Teff, RA, enriched
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

intersect(tcr_hv_sig$sample_id, tcr_ra_sig$sample_id)
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

hm <- data.frame(table(ra_hv_sig$cdr3aa, ra_hv_sig$subject_id))
hm_1 <- reshape(hm,direction="wide",timevar="Var2",idvar="Var1") %>% 
  column_to_rownames("Var1")

names(hm_1) <- gsub(x = names(hm_1), pattern = "Freq\\.", replacement = "")

# Verification des count sur 2 colonnes 

test_1 <- ra_hv_sig %>% 
  filter(subject_id == "FD47") %>% 
  dplyr::select(cdr3aa) %>% 
  table() %>% 
  as.data.frame() %>% 
  column_to_rownames(var="cdr3aa")

test_1_bis <- hm_1 %>% 
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

subject <- rownames(hm_1)
filter <- subject[hm_1$"32D0" != 0]
test_2_bis <- hm_1[filter,]["32D0"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$FD47)


# cas si nom des gènes commence par des chiffres (ou si on aime pas dplyr) 
subject <- rownames(hm_1)
filter <- subject[hm_1$"32D0" != 0]
test_2_bis <- hm_1[filter, ]["32D0"]

identical(rownames(test_2), rownames(test_2_bis))
all(test_2$Freq == test_2_bis$FE48)

# modifier hm_1 
# changer toute les valeurs > 1 par des 1
for(i in 1:ncol(hm_1)) {       
  hm_1[,i] <- ifelse(hm_1[,i] > 1, 1, hm_1[, i])
}

# récupérer le nom des sample_id dans tcr_ra et tcr_hv 
#  -> les mettre dans une dtf et regarder les colnames et voir si 
# discriminer en fonction des patient ou des sig ?  

# faire la matrice pour la heatmap
hm_2 <- as.matrix(hm_1)

hla_label <- hla_genotypes %>% 
  dplyr::filter(subject_id %in% as.character(hm$Var2)) %>% 
  dplyr::select(subject_id, DRB11, DPA11)

disease_label <- ra_hv_sig %>% 
  dplyr::filter(subject_id %in% as.character(hm$Var2)) %>% 
  dplyr::select(subject_id, disease)

# pour l'annotation 19 genotype pour DRB11 et 
# 2 genotypes pour DPA11
dplyr::select(subject_id, DRB11, DPA11)
# commencer par créer une dtf qui sera ensuite transformé en matrix


test <- sapply(unique(hm$cdr3aa), function(x) str_count(hm$cdr3aa,x))
test<-as.data.frame(test)

hm_count <- hm %>% 
  count(cdr3aa, subject_id, name = "count")

# preparer la matrice de donnes pour lheatmap()
hm_mat <- matrix(0, length(unique(ra_hv_sig$cdr3aa)),
                 length(unique(ra_hv_sig$subject_id)))

rownames(hm_mat) <- unique(ra_hv_sig$cdr3aa)
colnames(hm_mat) <- unique(ra_hv_sig$subject_id)
count_mat <- ra_hv_sig %>% 
  group_by(cdr3aa, subject_id) %>% 
  mutate(count = n())


county <- ra_hv_sig %>% 
  mutate(count_sig= count(cdr3aa))




t1d_sig <- tcr_sig %>% 
  dplyr::filter(cell_subset == "CD4_Teff", disease == "T1D",
                signature_type == "enriched")

ra_sig_meta <- tcr_meta %>% 
  left_join(ra_sig)

t1d_sig_meta <- tcr_meta %>% 
  left_join(t1d_sig)

ra_sig_hla <- ra_sig_meta %>% 
  left_join(hla_genes)

ra_sig_hla$label_hla <- ifelse(ra_sig_hla$subject_id %in% hla_genes$subject_id, 1, 0)

ra_sig_hm <- ra_sig_hla %>% 
  dplyr::select(subject_id, DRA1, DRA2, label_hla) %>% 
  as.matrix()

# pour tout les genes
ra_sig_hm <- ra_sig_hla %>% 
  dplyr::select(subject_id, 20:80) %>% 
  as.matrix()

rownames(ra_sig_hm) <- ra_sig_hla$subject_id
ra_sig_hm <-ra_sig_hm[, -1]

t1d_sig_hla <- t1d_sig_meta %>% 
  left_join(hla_genes_II)

# ra_sig_real <-  subset(ra_sig_hla, ra_sig_hla$subject_id %in% hla_genes_II$subject_id) 
# t1d_sig_real <-  subset(t1d_sig_hla, t1d_sig_hla$subject_id %in% hla_genes_II$subject_id) 
# il faut garder tout les sujet ou tout les sample_id 

# Heatmap
ha = HeatmapAnnotation(Category_HLA= ra_sig_hla$label_hla,col = list(Category_HLA = c("1" =  "red", "0" = "blue")))
jpeg(filename="Heatmap.jpeg",width =1000, height =1000,quality=100)
hm <-Heatmap(hm_2,
             name="expression",
             show_column_names=FALSE,
             # top_annotation=ha,
             show_row_dend=TRUE,
             show_column_dend=TRUE,
             cluster_rows=TRUE,
             cluster_columns=TRUE,
             column_title="Signature")
hm
dev.off()


# faire le lien HLA/TCR
ra_meta <- tcr_meta %>% 
  left_join(tcr_ra, by = "sample_id")

hv_meta <- tcr_meta %>% 
  left_join(tcr_hv, by = "sample_id")

t1d_meta <- tcr_meta %>% 
  left_join(tcr_t1d, by = "sample_id")

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
