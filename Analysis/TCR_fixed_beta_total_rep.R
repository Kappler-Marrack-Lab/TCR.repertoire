library(tidyr)
library(immunarch)
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(stringr)
library(vroom)
library(plyr)
library(ggpubr)
library(rstatix)
options(dplyr.summarise.inform = FALSE) # Suppress summarise info

# ---------------------------------------------------------------------------------
# Save Workspace file.
# save.image(file = "TCR_fixed_beta_total_rep.RData")
# # To restore your workspace, type this:
# load("TCR_fixed_beta_total_rep.RData")

# ---------------------------------------------------------------------------------
# Load in repertoire dataset via immune arch. Subset CD4 and CD8 data and export this data into a separate file to be loaded into all subsequent files.
# Set working directory to source file location then save target file path:
file_path = paste0(getwd(), "/total_repertoire_data")
TCR_DoVb_TR <- repLoad(file_path)
print(names(TCR_DoVb_TR))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# filter out the sequences where V.name

# Subset the TCR Vα DOβ samples in ImmuneArch format based on CD4 or CD8 data.
# Define filter condition names:
CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'
ff = 'ff'

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define functions:
filter_nonfunctional <- function(x){
  x <- x %>%
    # mutate(V.name=recode(V.name,  `TRAV05-4DN`="TRAV05-4A_NEW", `TRAV05-4A`="TRAV05-4DN_NEW")) %>%
    # mutate(V.name=recode(V.name,  `TRAV05-4A_NEW` = "TRAV05-4A", `TRAV05-4DN_NEW` = "TRAV05-4A")) %>%
    filter((V.name %in% c("TRAV14-1D", "TRAV15-1D")) == FALSE)
    # filter(Clones > 1)
  return(x)
} 
Remove_CDR3NT <- function(x){
  x2 <- x %>% 
    dplyr::select(Clones, Proportion, CDR3.aa, V.name, J.name) %>%
    mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    group_by(CDR3.aa, V.name, J.name) %>%
    dplyr::summarise(Clones = sum(Clones), Proportion = sum(Proportion)) %>% # add up the clonal count and proportional value the different CDR3NT sequences had
    as.data.frame()
  x3  <- x2 %>%
    # dplyr::mutate(Proportion = Clones/sum(Clones)) %>%
    # filter(Proportion >= 5e-06) %>%
    add_column(CDR3.nt = NA) %>% add_column(D.name = NA) %>% add_column(V.end = NA) %>%
    add_column(D.start = NA) %>% add_column(D.end = NA) %>% add_column(J.start = NA) %>%
    add_column(VJ.ins = NA) %>% add_column(VD.ins = NA) %>% add_column(DJ.ins = NA) %>% add_column(Sequence = NA) %>%
    select(Clones, Proportion, CDR3.nt, CDR3.aa, V.name, D.name, J.name, V.end,  D.start, D.end,  J.start, VJ.ins, VD.ins, DJ.ins, Sequence)
  return(x3)
} # function to combine the CDR3NT infromation to CDR3AA only
filter_low_frequency <- function(x){
  x2 <- x %>%
    select(Clones, Proportion, CDR3.aa, V.name, J.name) %>%
    mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    group_by(CDR3.aa, V.name, J.name) %>%
    dplyr::summarise(Clones = sum(Clones), Proportion = sum(Proportion)) %>% # add up the clonal count and proportional value the different CDR3NT sequences had
    as.data.frame()
  x3  <- x2 %>%
    # dplyr::mutate(Proportion = Clones/sum(Clones)) %>%
    filter(Proportion >= 5e-06) %>%
    add_column(CDR3.nt = NA) %>% add_column(D.name = NA) %>% add_column(V.end = NA) %>%
    add_column(D.start = NA) %>% add_column(D.end = NA) %>% add_column(J.start = NA) %>%
    add_column(VJ.ins = NA) %>% add_column(VD.ins = NA) %>% add_column(DJ.ins = NA) %>% add_column(Sequence = NA) %>%
    select(Clones, Proportion, CDR3.nt, CDR3.aa, V.name, D.name, J.name, V.end,  D.start, D.end,  J.start, VJ.ins, VD.ins, DJ.ins, Sequence)
  return(x3)
} # this also combines the CDR3NT informaiton to CDR3AA only but it also removes sequences below a frequency cutoff of: 5e-06.
subset_ImmuArchList <- function(TCR_DoVb_TR, condition){
  C_meta <- filter(TCR_DoVb_TR$meta, `Co-receptor` == condition | MHC == condition)
  C_samples <- C_meta$Sample
  subset.data <- TCR_DoVb_TR$data %>%
    magrittr::extract(C_samples)
  NEW_DoVb_TR <- list(
    data = subset.data,
    meta = C_meta
  )
  # NEW_DoVb_TR$data <- lapply(NEW_DoVb_TR$data, filter_nonfunctional)
  return(NEW_DoVb_TR)
} # Function takes in the ImmuneArch dataframe list and then subsets based on the co-receptor value.

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# First apply a function to remove nonfunctional TRAV genes.
TCR_DoVb_TR$data <- lapply(TCR_DoVb_TR$data, filter_nonfunctional) # ones are retained here. Just removes non functional TRAV genes

# create a column in the metadata which holds the total number of productive sequences for each run
TCR_DoVb_TR$meta$Total_Productive_Reads <- sapply(TCR_DoVb_TR$data, function(i) sum(i$Clones))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Next copy data to remove sequences which only occur once and copy these datasets
TCR_DoVb_TR_NT.ones.removed <- repFilter(TCR_DoVb_TR, .method = "by.clonotype", .query = list(Clones = morethan(1))) # ones removed from the TRAV-CDR3NT-TRAJ unique sequences.


TCR_DoVb_TR_AA.ones.removed <- TCR_DoVb_TR # copy the original data
TCR_DoVb_TR_AA.ones.removed$data <- lapply(TCR_DoVb_TR$data, Remove_CDR3NT) # Removes CDR3NT information so we can use the CDR3AA as the idenfifyer
TCR_DoVb_TR_AA.ones.removed <- repFilter(TCR_DoVb_TR_AA.ones.removed, .method = "by.clonotype", .query = list(Clones = morethan(1))) # ones removed from the TRAV-CDR3AA-TRAJ unique sequences.

# Summary:
# TCR_DoVb_TR_NT.ones.removed : Contains original data minus sequences which occur once in TRAV-CDR3NT-TRAJ
# TCR_DoVb_TR_AA.ones.removed : Contains original data minus sequences which occur once in TRAV-CDR3AA-TRAJ

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Next apply a frequency cutoff to the data only containing the TRAV-CDR3AA-TRAJ unique sequences.
# THEN also remove the sequences which only show up once. Applying the frequency cut off may not always remove these sequecnces which appear once.
# But we cannot trust that these sequences which appear once are not there due to some PCR/sequencing error so we need to remove them regardless.
# Note: the order of operations of these two steps does matter.
TCR_DoVb_TR_AA.ones.removed.LFSR <- TCR_DoVb_TR # copy the original data
TCR_DoVb_TR_AA.ones.removed.LFSR$data <- lapply(TCR_DoVb_TR$data, filter_low_frequency) # Removes CDR3NT information so we can use the CDR3AA as the idenfifyer, and low frequency sequences.
TCR_DoVb_TR_AA.ones.removed.LFSR <- repFilter(TCR_DoVb_TR_AA.ones.removed.LFSR, .method = "by.clonotype", .query = list(Clones = morethan(1)))  # ones removed from the TRAV-CDR3AA-TRAJ unique sequences.

# Summary:
# TCR_DoVb_TR_AA.ones.removed.LFSR : CDR3NT sequences are aggrigated, frequency is recalcualted after agrigation. a frequency cutoff is applied. ones are then removed.

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# separate datasets into CD4, CD8 and DP
TCR_DoVb_CD4.ones.retained <- repFilter(TCR_DoVb_TR, .method = "by.meta", .query = list(`Co-receptor` = include("CD4"))) # CD4
TCR_DoVb_CD8.ones.retained <- repFilter(TCR_DoVb_TR, .method = "by.meta", .query = list(`Co-receptor` = include("CD8"))) # CD8
TCR_DoVb_DP.ones.retained <- repFilter(TCR_DoVb_TR, .method = "by.meta", .query = list(`Co-receptor` = include("DP"))) # DP

TCR_DoVb_CD4_NT.ones.removed <- repFilter(TCR_DoVb_TR_NT.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("CD4"))) # CD4
TCR_DoVb_CD8_NT.ones.removed <- repFilter(TCR_DoVb_TR_NT.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("CD8"))) # CD8
TCR_DoVb_DP_NT.ones.removed  <- repFilter(TCR_DoVb_TR_NT.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("DP"))) # DP

TCR_DoVb_CD4_AA.ones.removed <- repFilter(TCR_DoVb_TR_AA.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("CD4"))) # CD4
TCR_DoVb_CD8_AA.ones.removed <- repFilter(TCR_DoVb_TR_AA.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("CD8"))) # CD8
TCR_DoVb_DP_AA.ones.removed  <- repFilter(TCR_DoVb_TR_AA.ones.removed, .method = "by.meta", .query = list(`Co-receptor` = include("DP"))) # DP

TCR_DoVb_CD4_AA.ones.removed.LFSR <- repFilter(TCR_DoVb_TR_AA.ones.removed.LFSR, .method = "by.meta", .query = list(`Co-receptor` = include("CD4"))) # CD4
TCR_DoVb_CD8_AA.ones.removed.LFSR <- repFilter(TCR_DoVb_TR_AA.ones.removed.LFSR, .method = "by.meta", .query = list(`Co-receptor` = include("CD8"))) # CD8
TCR_DoVb_DP_AA.ones.removed.LFSR  <- repFilter(TCR_DoVb_TR_AA.ones.removed.LFSR, .method = "by.meta", .query = list(`Co-receptor` = include("DP"))) # DP

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SAVE DATA
save(TCR_DoVb_TR, TCR_DoVb_CD4.ones.retained, TCR_DoVb_CD8.ones.retained, TCR_DoVb_DP.ones.retained, file = "TCR_fixed_beta_total_rep.ones.retained.RData") # Save objects with ones retained
save(TCR_DoVb_TR_NT.ones.removed, TCR_DoVb_CD4_NT.ones.removed, TCR_DoVb_CD8_NT.ones.removed, TCR_DoVb_DP_NT.ones.removed, file = "TCR_fixed_beta_total_rep.CDR3NT.ones.removed.RData") # Save objects with CDR3NT ones removed
save(TCR_DoVb_TR_AA.ones.removed, TCR_DoVb_CD4_AA.ones.removed, TCR_DoVb_CD8_AA.ones.removed, TCR_DoVb_DP_AA.ones.removed,  file = "TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData") # Save objects with CDR3AA ones removed
save(TCR_DoVb_TR_AA.ones.removed.LFSR, TCR_DoVb_CD4_AA.ones.removed.LFSR, TCR_DoVb_CD8_AA.ones.removed.LFSR, TCR_DoVb_DP_AA.ones.removed.LFSR,  file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # Save objects with CDR3NT collapsed into CDR3AA, low frequency sequences removed and ones removed.

print("task finished.")
print("Data stored in:")
print("TCR_fixed_beta_total_rep.ones.retained.RData, TCR_fixed_beta_total_rep.CDR3NT.ones.removed.RData, TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData, TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData")


# # Most analysis will be done with this data.
# TCR.CD4_DoVb_TR <- subset_ImmuArchList(TCR_DoVb_TR_tetra_separated_one.mouse, CD4)
# TCR_DoVb_TR_tetra_grouped_CD4 <- subset_ImmuArchList(TCR_DoVb_TR_tetra_grouped, CD4)
# 
# TCR.CD8_DoVb_TR <- subset_ImmuArchList(TCR_DoVb_TR_tetra_separated_one.mouse, CD8)
# TCR.DP_DoVb_TR <- subset_ImmuArchList(TCR_DoVb_TR_tetra_separated_one.mouse, DP)
# 
# # P2M = plus two mice, DMR = Divided mice removed
# TCR.CD4_DoVb_TR_P2M.DMR <- subset_ImmuArchList(TCR_DoVb_TR_tetra_separated_divided_mice.removed, CD4)
# TCR.CD8_DoVb_TR_P2M.DMR <- subset_ImmuArchList(TCR_DoVb_TR_tetra_separated_divided_mice.removed, CD8)
# 
# # Save objects
# save(TCR_DoVb_TR, TCR_DoVb_TR_tetra_separated, TCR.CD4_DoVb_TR, TCR.CD8_DoVb_TR, TCR.DP_DoVb_TR, TCR_DoVb_TR_tetra_grouped_CD4,
#      file = "TCR_fixed_beta_total_rep.RData")
# 
# #----------------------------------------------------------------------------------
# # Save objects
# save(TCR.CD4_DoVb_TR_P2M.DMR, TCR.CD8_DoVb_TR_P2M.DMR,
#      file = "TCR_fixed_beta_total_rep_plus two mice_and Divided mice removed.RData")
# 
# #----------------------------------------------------------------------------------
