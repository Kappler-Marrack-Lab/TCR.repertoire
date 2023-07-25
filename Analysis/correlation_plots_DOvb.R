# Author: Alex Brown
# License: MIT
#-------------------STARTING NOTES-------------------------------------------------
# 2020-10-27.
# Make sure the source data is correct and the path is set before running
# 2020-1-11
# Simplified and added functions to make more generalizeable for large repertoire datasets.
# Old code is commented out at the end which could still be useful.
#----------------------------------------------------------------------------------
# Load Libraries

library(devtools)
library(tidyr)
library(immunarch)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggforce)
library(scales)
library(dplyr)
library(stringr)

library(vroom)
library(plyr)
library(ggpubr)
library(rstatix)
library(data.table)

library(Rmisc)
library(ggsignif)
library(stringr)
library(openxlsx)
#----------------------------------------------------------------------------------
# Save Workspace file.
# save.image(file = "RepOverlapCorrelation.RData")

# To restore your workspace, type this:
# load("RepOverlapCorrelation.RData")
# Save multiple objects:
# save(
#   melted_cormat_j.CD4, melted_cormat_j.CD8,
#   melted_cormat_cj.CD4, melted_cormat_cj.CD8,
#   melted_cormat_j.CD4.simX, melted_cormat_cj.CD4.simX,
#   file = "RepOverlapCorrelation_matrix.RData"
#   )
# To load the data again
# load("RepOverlapCorrelation_matrix.RData")
#----------------------------------------------------------------------------------

# load the filtered data objects
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed
load(file = "counts tables filtered.LFSR.ones.removed.RData") # counts table

CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"
IAbHet = "b+/-"
tetra_bb.ff = "bb&ff"
tetra_bb.g7g7 = "bb&g7g7"

# Define MHC group colors:
ff_color ="#c5944e"
ss_color = "#fbff00"
bb_color = "#54a0fb" 
g7g7_color = "#9f104c"
bxf_color = "#808000"
fxs_color = "#fa5f00"
bxs_color = "#1afe02"
bxg7_color = "#f74ed6"
IAbHet_color = "#bddbff"
tetra_bb.ff_color = "#EEC0C4"
tetra_bb.g7g7_color = "#335566"


tetra_bb.ff_color.A = "#54a0fb"
tetra_bb.ff_color.B = "#c5944e"
tetra_bb.g7g7_color.A = "#54a0fb"
tetra_bb.g7g7_color.B = "#9f104c"

# Define the MHC groups and input elements for the current plot
F.MHC_groups <- c(ff, ss, bb, g7g7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color)
Jac = "Jaccard"
Chao = "Chao Jaccard"
#------------- Functions to calculate value matrix: ---------------------------------------------------------------------
# This re-assigns the databases lazily to add the tetraprental sequences I am interested in
#Data in use:
# TCR_DoVb_CD4_AA.ones.removed.LFSR
# TCR_DoVb_CD8_AA.ones.removed.LFSR
TCR_DoVb_CD4_AA.ones.removed.LFSR$meta$Tetraparental_Cell_Source <- NA
TCR_DoVb_CD8_AA.ones.removed.LFSR$meta$Tetraparental_Cell_Source <- NA
# comment the above lines to return back to the sourced TTCR database without the added tetraparental sequences.

# Jaccard Index Function    (Correct function based on Pippas Data)
Jaccard_Matrix <- function(x, y) {
  x_merged <- within(x, Merge <- paste(V.name, J.name, D.name, CDR3.aa, sep = "_"))
  y_merged <- within(y, Merge <- paste(V.name, J.name, D.name, CDR3.aa, sep = "_"))
  x_merged$Origin <- 1
  y_merged$Origin <- 1
  x_merged <- x_merged %>% uncount(Clones)
  y_merged <- y_merged %>% uncount(Clones)  
  
  x_merged <- x_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  y_merged <- y_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  
  J.numerator <- length(dplyr::intersect(x_merged$Merge, y_merged$Merge))
  x.row <- nrow(x_merged)
  y.row <- nrow(y_merged)
  Jaccard_index <- J.numerator/(x.row+y.row-J.numerator)
  # This is equivalent to:
  # (|A∩B|)/(|A|+|B|- |A∩B|)
}
Chao_Jaccard_Matrix <- function(x, y) {
  x_merged <- within(x, Merge <- paste(V.name, J.name, D.name, CDR3.aa, sep = "_"))
  y_merged <- within(y, Merge <- paste(V.name, J.name, D.name, CDR3.aa, sep = "_"))
  x_merged$Origin <- 1
  y_merged$Origin <- 1
  x_merged <- x_merged %>% uncount(Clones)
  y_merged <- y_merged %>% uncount(Clones)  
  
  x_merged <- x_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  y_merged <- y_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  
  x_merged$CountX <- x_merged$Origin
  y_merged$CountY <- y_merged$Origin
  x_merged = subset(x_merged, select = -c(Origin))
  y_merged = subset(y_merged, select = -c(Origin))
  test <- inner_join(x_merged, y_merged, by = "Merge") #this is identical to the intersection and keeps the count columns
  total.X.count <- sum(x_merged$CountX)
  total.Y.count <- sum(y_merged$CountY)
  test$CountX.proportion.of.Xtotal <- test$CountX / total.X.count
  test$CountY.proportion.of.Ytotal <- test$CountY / total.Y.count
  U <- sum(test$CountX.proportion.of.Xtotal)
  V <- sum(test$CountY.proportion.of.Ytotal)
  ChaoJaccard_index <- (U*V)/(U+V-(U*V))
}

jaccard.CD4 <- apply_symm(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix
jaccard.CD8 <- apply_symm(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix

# Chao Jaccard Index Function  # Correct function based on Pippas Data
chao_jaccard.CD4 <- apply_symm(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, Chao_Jaccard_Matrix)
chao_jaccard.CD8 <- apply_symm(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, Chao_Jaccard_Matrix)
# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  diag(cormat) <- 1
  return(cormat)
}
# Overlap Statistics Function
Overlap.Stats <- function(melted_cormat, TCR_DoVb_DB){
  melted_cormat_tib <- as_tibble(melted_cormat) # this is what I want to keep and compare
  melted_cormat_tib
  
  # Grab the MHC and +/- Tetraparental MHC information from the metadata
  TCR_DoVb_DB$meta[TCR_DoVb_DB$meta == "N/A"] <- NA
  
  merge_meta <- within(TCR_DoVb_DB$meta,  Merge <- paste(MHC, Tetraparental_Cell_Source, sep = "_"))
  merge_meta$Merge <- gsub("_NA.*", "\\1", merge_meta$Merge)
  
  # add in the meta data MHC to the Var1 as var1_MHC and to Var2 as var2_MHC then sort by vector name [var1_MHC] then [var2_MHC]
  melted_cormat_tib$var1_MHC <- merge_meta$Merge[match(melted_cormat_tib$Var1,merge_meta$Sample)]
  melted_cormat_tib$var2_MHC <- merge_meta$Merge[match(melted_cormat_tib$Var2,merge_meta$Sample)]
  
  
  melted_cormat_tib <- melted_cormat_tib[with(melted_cormat_tib, order(var1_MHC, var2_MHC)),]
  melted_cormat_tib
  # Next Create a Comparator and Comparatee Column in a new dataset.
  MHC_sets  <- unique(melted_cormat_tib$var1_MHC) # list of the unique MHC sets
  
  # Make new dataframe and give it the info for final column names:
  Overlap_Statistics.DF <- data.frame(matrix(ncol = 8, nrow = 0))
  x <- c("Var1", "Var2", "value", "rnd_value", "var1_MHC", "var2_MHC", "Comparator", "Comparatee")
  colnames(Overlap_Statistics.DF) <- x
  
  for(i in 1:length(MHC_sets)){
    print(paste0("Currently running ", MHC_sets[i], "  to subset MHC sets"))
    df <- melted_cormat_tib[melted_cormat_tib$var1_MHC == MHC_sets[i] | melted_cormat_tib$var2_MHC == MHC_sets[i],] #subset data
    df$Comparator <- MHC_sets[i]
    df$Comparatee<-ifelse(df$var1_MHC != MHC_sets[i], df$var1_MHC, df$var2_MHC)
    df <- df[with(df, order(Comparator, Comparatee)),]
    #merge all to one tibble
    Overlap_Statistics.DF <- rbind(Overlap_Statistics.DF, df)
  }
  # For plotting:
  Overlap_Statistics.DF = subset(Overlap_Statistics.DF, select = -c(var1_MHC,var2_MHC))
  
  # ADD P-VALUES TO GGPLOT FACETS
  Overlap_Statistics.DF$Comparator <- as.factor(Overlap_Statistics.DF$Comparator)
  Overlap_Statistics.DF$Comparatee <- as.factor(Overlap_Statistics.DF$Comparatee)
  # Facet by the Comparator variable and compare the levels of the Comparatee variable on the x-axis.
  
  return(Overlap_Statistics.DF)
}
# 
#-----------------------Calculate Jaccard correlation matrix-----------------------
Value_Matrix <- function(correlation_MX){
  lower_tri_j <- get_lower_tri(correlation_MX)# grab lower triangle of matrix
  melted_cormat <- reshape2::melt(lower_tri_j, na.rm = TRUE) # reshape to dataframe
  # This will restore the whole square and make the scale bars look the same relative to other plots.
  # melted_cormat <- reshape2::melt(correlation_MX, na.rm = FALSE) # uncomment this and comment the top to to ge the square matrix
  melted_cormat[is.na(melted_cormat)] <- 1
  melted_cormat$rnd_value <- scientific(melted_cormat$value, digits = 2)
  return(melted_cormat)
}

melted_cormat_j.CD4 <- Value_Matrix(jaccard.CD4)
melted_cormat_j.CD8 <- Value_Matrix(jaccard.CD8)
melted_cormat_cj.CD4 <- Value_Matrix(chao_jaccard.CD4)
melted_cormat_cj.CD8 <- Value_Matrix(chao_jaccard.CD8)

center_Jaccard.CD4 <- (max(melted_cormat_j.CD4$value)+min(melted_cormat_j.CD4$value))/2
center_Jaccard.CD8 <- (max(melted_cormat_j.CD8$value)+min(melted_cormat_j.CD8$value))/2
center_cj.CD4 <- (max(melted_cormat_cj.CD4$value)+min(melted_cormat_cj.CD4$value))/2
center_cj.CD8 <- (max(melted_cormat_cj.CD8$value)+min(melted_cormat_cj.CD8$value))/2
# 
# Value_Matrix_Plot_Jaccard <- function(melted_cormat, center){
#   cormat_HeatMap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#     geom_tile(colour="white", size=1.50, stat="identity") +
#     scale_x_discrete(position = "top")+
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          midpoint = center, space = "Lab", 
#                          name="Jaccard\nIndex") +
#     labs(title="Repertoire overlap", 
#          subtitle="Jaccard index", x = "", y = ""
#          ) +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(angle = 90),
#           panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()
#           ) +
#     coord_fixed()
#   return(cormat_HeatMap)
# }
# Value_Matrix_Plot_ChaoJaccard <- function(melted_cormat, center){
#   cormat_HeatMap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#     geom_tile(colour="white", size=1.50, stat="identity") +
#     scale_x_discrete(position = "top")+
#     scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                          midpoint = center, space = "Lab", 
#                          name="ChaoJaccard\nIndex") +
#     labs(title="Repertoire overlap", 
#          subtitle="ChaoJaccard index", x = "", y = "") +
#     theme_minimal()+ 
#     theme(axis.text.x = element_text(angle = 90),
#           panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()
#     ) +
#     coord_fixed()
#   return(cormat_HeatMap)
# }
# 
# 
# 
# # # Export PDF captures of the Jaccard and Chao Jaccard Correlation Matrix
# # pdf("Figure Export/Repertoire Overlap/Va_CD4_jaccard_minimal.pdf", width = 12, height = 12)
# # Value_Matrix_Plot_Jaccard(melted_cormat_j.CD4, center_Jaccard.CD4)
# # dev.off()
# # pdf("Figure Export/Repertoire Overlap/Va_CD8_jaccard_minimal.pdf", width = 12, height = 12)
# # Value_Matrix_Plot_Jaccard(melted_cormat_j.CD8, center_Jaccard.CD8)
# # dev.off()
# # pdf("Figure Export/Repertoire Overlap/Va_CD4_ChaoJaccard_minimal.pdf", width = 12, height = 12)
# # Value_Matrix_Plot_ChaoJaccard(melted_cormat_cj.CD4, center_cj.CD4)
# # dev.off()
# # pdf("Figure Export/Repertoire Overlap/Va_CD8_ChaoJaccard_minimal.pdf", width = 12, height = 12)
# # Value_Matrix_Plot_ChaoJaccard(melted_cormat_cj.CD8, center_cj.CD8)
# # dev.off()
# 
# # Pairwise correlation matrix of parental sequences only separated by CD4 and CD8
# #----------------------------------------------------------------------------------
# # Pairwise correlation matrix of parental sequences only separated by CD4 and CD8
# # This will rename the samples based on MHC and assiggns them an ordered index for the pairwise coorilation matrixs
# subset_ImmuArchList <- function(TCR_DoVb_DB, condition){
#   C_meta <- filter(TCR_DoVb_DB$meta, MHC %in% condition)
#   C_meta <- C_meta %>%
#     arrange(MHC) %>%
#     mutate(Index = row_number()) %>%
#     mutate(Sample_New = paste(Index, MHC))
#   C_samples <- C_meta$Sample
#   
#   subset.data <- TCR_DoVb_DB$data %>%
#     magrittr::extract(C_samples)
#   TCR.NEW_DoVb_DB <- list(
#     data = subset.data,
#     meta = C_meta
#   )
#   names(TCR.NEW_DoVb_DB$data) <- C_meta$Sample_New
#   return(TCR.NEW_DoVb_DB)
# }
# Parental_MHC <- c(bb, ff, g7g7, ss)
# TCR.CD4_Parental_MHC <- subset_ImmuArchList(TCR_DoVb_CD4_AA.ones.removed.LFSR, Parental_MHC)
# TCR.CD8_Parental_MHC <- subset_ImmuArchList(TCR_DoVb_CD8_AA.ones.removed.LFSR, Parental_MHC)
# # Generate Value Matrix
# TCR.CD4_Parental_MHC_jaccard <- apply_symm(TCR.CD4_Parental_MHC$data, Jaccard_Matrix)
# TCR.CD8_Parental_MHC_jaccard <- apply_symm(TCR.CD8_Parental_MHC$data, Jaccard_Matrix)
# TCR.CD4_Parental_MHC_chao.jaccard <- apply_symm(TCR.CD4_Parental_MHC$data, Chao_Jaccard_Matrix)
# TCR.CD8_Parental_MHC_chao.jaccard <- apply_symm(TCR.CD8_Parental_MHC$data, Chao_Jaccard_Matrix)
# # Melt value Matrix to extract data
# TCR.CD4_Parental_MHC_jaccard_melted_cormat <- Value_Matrix(TCR.CD4_Parental_MHC_jaccard)
# TCR.CD8_Parental_MHC_jaccard_melted_cormat <- Value_Matrix(TCR.CD8_Parental_MHC_jaccard)
# TCR.CD4_Parental_MHC_chao.jaccard_melted_cormat <- Value_Matrix(TCR.CD4_Parental_MHC_chao.jaccard)
# TCR.CD8_Parental_MHC_chao.jaccard_melted_cormat <- Value_Matrix(TCR.CD8_Parental_MHC_chao.jaccard)
# # Find the center of the values to make the plot scaling make sense
# TCR.CD4_Parental_center_Jaccard <- (max(TCR.CD4_Parental_MHC_jaccard_melted_cormat$value)+min(TCR.CD4_Parental_MHC_jaccard_melted_cormat$value))/2
# TCR.CD8_Parental_center_Jaccard <- (max(TCR.CD8_Parental_MHC_jaccard_melted_cormat$value)+min(TCR.CD8_Parental_MHC_jaccard_melted_cormat$value))/2
# TCR.CD4_Parental_center_Chao.Jaccard <- (max(TCR.CD4_Parental_MHC_chao.jaccard_melted_cormat$value)+min(TCR.CD4_Parental_MHC_chao.jaccard_melted_cormat$value))/2
# TCR.CD8_Parental_center_Chao.Jaccard <- (max(TCR.CD8_Parental_MHC_chao.jaccard_melted_cormat$value)+min(TCR.CD8_Parental_MHC_chao.jaccard_melted_cormat$value))/2
# # Generate ggPlots
# A <- Value_Matrix_Plot_Jaccard(TCR.CD4_Parental_MHC_jaccard_melted_cormat, TCR.CD4_Parental_center_Jaccard)
# B <- Value_Matrix_Plot_Jaccard(TCR.CD8_Parental_MHC_jaccard_melted_cormat, TCR.CD8_Parental_center_Jaccard)
# C <- Value_Matrix_Plot_ChaoJaccard(TCR.CD4_Parental_MHC_chao.jaccard_melted_cormat, TCR.CD4_Parental_center_Chao.Jaccard)
# D <- Value_Matrix_Plot_ChaoJaccard(TCR.CD8_Parental_MHC_chao.jaccard_melted_cormat, TCR.CD8_Parental_center_Chao.Jaccard)
# # Export PDF captures of the Jaccard and Chao Jaccard Correlation Matrix
# pdf("Figure Export/Repertoire Overlap/Matrix Plots/Parrents_CD4_CD8_Jaccard_minimal.pdf", width = 12, height = 6)
# A + B
# dev.off()
# pdf("Figure Export/Repertoire Overlap/Matrix Plots/Parrents_CD4_CD8_ChaoJaccard_minimal.pdf", width = 12, height = 6)
# C + D
# dev.off()
# #----------------------------------------------------------------------------------
# # 
# subset_ImmuArchList_2 <- function(TCR_DoVb_DB, condition){
#   C_meta <- filter(TCR_DoVb_DB$meta, MHC %in% condition)
#   C_meta <- C_meta %>%
#     arrange(`Co-receptor`) %>%
#     mutate(Index = row_number()) %>%
#     mutate(Sample_New = paste(`Co-receptor`, MHC, Sample))
#   C_samples <- C_meta$Sample
#   
#   subset.data <- TCR_DoVb_DB$data %>%
#     magrittr::extract(C_samples)
#   TCR.NEW_DoVb_DB <- list(
#     data = subset.data,
#     meta = C_meta
#   )
#   names(TCR.NEW_DoVb_DB$data) <- C_meta$Sample_New
#   return(TCR.NEW_DoVb_DB)
# }
# Parental_MHC <- c(ff) # see what is going on with the ff CD4 and CD8 samples with high overlap.
# TCR_DoVb_DB_CD4_CD8_overlap <- subset_ImmuArchList_2(TCR_DoVb_TR, Parental_MHC)
# TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard <- apply_symm(TCR_DoVb_DB_CD4_CD8_overlap$data, Chao_Jaccard_Matrix)
# TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_melted <- Value_Matrix(TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard)
# TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_center <- (max(TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_melted$value)+min(TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_melted$value))/2
# E <- Value_Matrix_Plot_ChaoJaccard(TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_melted, TCR_DoVb_DB_CD4_CD8_overlap_chao.jaccard_center)
# pdf("Figure Export/Repertoire Overlap/Matrix Plots/ff_CD4vsCD8_ChaoJaccard_minimal.pdf", width = 10, height = 10)
# E
# dev.off()
#------------- Calculate Overlap Statistics Separately for CD4 and CD8 ---------------------------------------------------------------------
Jaccard_Overlap_Statistics.CD4 <- Overlap.Stats(melted_cormat_j.CD4, TCR_DoVb_CD4_AA.ones.removed.LFSR)
Jaccard_Overlap_Statistics.CD8 <- Overlap.Stats(melted_cormat_j.CD8, TCR_DoVb_CD8_AA.ones.removed.LFSR)
ChaoJaccard_Overlap_Statistics.CD4 <- Overlap.Stats(melted_cormat_cj.CD4, TCR_DoVb_CD4_AA.ones.removed.LFSR)
ChaoJaccard_Overlap_Statistics.CD8 <- Overlap.Stats(melted_cormat_cj.CD8, TCR_DoVb_CD8_AA.ones.removed.LFSR)
# Rename tetrapartental groups. Name is too long
Jaccard_Overlap_Statistics.CD4 <- Jaccard_Overlap_Statistics.CD4 %>% mutate(across(Comparator:Comparatee,~ gsub("tetraparental", "tetra.", .))) %>% mutate(across(Comparator:Comparatee, as.factor))
ChaoJaccard_Overlap_Statistics.CD4 <- ChaoJaccard_Overlap_Statistics.CD4 %>% mutate(across(Comparator:Comparatee,~ gsub("tetraparental", "tetra.", .))) %>% mutate(across(Comparator:Comparatee, as.factor))


#  ---------------------------------- Define statistical tests to use ----------------------------------
# These make just the comparisons we want to make and adds t-test satistical bars
sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

# Grouped One-way ANOVA test
#:::::::::::::::::::::::::::::::::::::::::
# one-way ANOVA is generally followed up by Tukey post-hoc tests
ANOVA.Tu <- function(df.in){
  tukey.stats <- df.in %>% 
    select(value, Comparator, Comparatee) %>%
    group_by(Comparator) %>%
    tukey_hsd(value ~ Comparatee)
  return(tukey.stats)
}



Simple.T.test <- function(df.in){
    t.stats <- df.in %>% select(value, Comparator, Comparatee) %>%
      group_by(Comparator) %>%
      t_test(data =., value ~ Comparatee) %>%
      add_significance("p")
    return(t.stats)
  }
# -------------- Define mininal overlap plots ----------------------------------------------------------------------

minimally_Defined.Overlap.plots <- function(Overlap_Statistics, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, type){
  MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct statistical  comparisons
  filtered.df <- Overlap_Statistics %>%
    filter(Comparator %in% Comparator_group) %>%
    filter(Comparatee %in% F.MHC_groups) %>%
    filter(value != 1)
  
  filtered.df$Comparatee <- factor(filtered.df$Comparatee, levels = rev(F.MHC_groups))
  
  OL_plot <- ggplot(data=filtered.df, aes(x=Comparatee, y=value, fill=Comparatee))+
    geom_bar(stat="summary", fun=mean)+
    stat_summary(fun.data=mean_se, geom = "linerange")+
    # geom_point(alpha = 0.1 )+
    geom_jitter(shape = 1, width = 0.1, size = 0.5)+
    # geom_signif(comparisons = MHC_list,
    #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
    # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    coord_flip()+
    guides(fill = guide_legend(reverse = TRUE))+
    theme_classic()+
    theme(legend.position = "none")+
    facet_wrap(~Comparator, ncol=1, strip.position = "left") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face="bold", size = 12, color = "black"),
          axis.text.x = element_text(color = "black"),
          strip.text.y.left = element_text(angle=0, size = 20, face="bold"),
          strip.placement = "outside",
          panel.spacing.y=unit(-1.0, "pt"),
          panel.border=element_rect(color="black", fill = NA, size=1)
          ) +
    # theme(aspect.ratio = 0.3) + # This can make the plots all constantly narrower without playing around with the export size.
    # ylim(0, ylimit)+
    scale_y_continuous(breaks=seq(0,1,by=.25), expand = expansion(mult = c(0, .1)), limits = c(0, ylimit))+
    labs(caption=paste0(type), x="Comparatee", y=paste0(type, " similarity coefficient"))
  return(OL_plot)
}
# ------------------------------------------------------------------------------------
# Plot all CD4 together then plot CD4 and CD8 for the MHC hom comparisons
F.MHC_groups <- c(bb, ff, g7g7, ss)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color)
Comparator_group <- c(bb, ff, g7g7, ss)
ylimit = 0.5
MHChom_OL_CD4_Jac.plot <- minimally_Defined.Overlap.plots(Jaccard_Overlap_Statistics.CD4, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
MHChom_OL_CD4_CJ.plot <- minimally_Defined.Overlap.plots(ChaoJaccard_Overlap_Statistics.CD4, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)

ANOVA.Tu.Jaccard.CD4 <- Jaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD4 <- ChaoJaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()
# ------------------------
F.MHC_groups <- c(bb, ff, g7g7)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color)
Comparator_group <- c(bb, ff, g7g7)
MHChom_OL_CD8_Jac.plot <- minimally_Defined.Overlap.plots(Jaccard_Overlap_Statistics.CD4, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
MHChom_OL_CD8_CJ.plot <- minimally_Defined.Overlap.plots(ChaoJaccard_Overlap_Statistics.CD8, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
ANOVA.Tu.Jaccard.CD8 <- Jaccard_Overlap_Statistics.CD8 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD8 <- ChaoJaccard_Overlap_Statistics.CD8 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()
# -------------------------
Title_CD4 <- ggtitle("TCR Alpha Repertoire Overlap:CD4")
Title_CD8 <- ggtitle("TCR Alpha Repertoire Overlap:CD8")
compiled_MHC_hom.CD4 <- MHChom_OL_CD4_Jac.plot + Title_CD4 + MHChom_OL_CD4_CJ.plot
compiled_MHC_hom.CD8 <- MHChom_OL_CD8_Jac.plot + Title_CD8 + MHChom_OL_CD8_CJ.plot
compiled_MHC_hom <- (MHChom_OL_CD4_Jac.plot + Title_CD4 + MHChom_OL_CD4_CJ.plot) /
  (MHChom_OL_CD8_Jac.plot + Title_CD8 + MHChom_OL_CD8_CJ.plot)

# make plots
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.MHChom.pdf", compiled_MHC_hom.CD4, width=190,height=80, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD8.MHChom.pdf", compiled_MHC_hom.CD8, width=190,height=80, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.MHChom_nostats.pdf", compiled_MHC_hom.CD4, width=245,height=100, units = "mm")

# Export Statistical test tables
dataset_names <- list('ANOVA.Tu.Jaccard.CD4' = ANOVA.Tu.Jaccard.CD4, 'ANOVA.Tu.ChaoJaccard.CD4' = ANOVA.Tu.ChaoJaccard.CD4)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.ChaoJaccardJaccard Tu.xlsx')

dataset_names <- list('ANOVA.Tu.Jaccard.CD8' = ANOVA.Tu.Jaccard.CD8, 'ANOVA.Tu.ChaoJaccard.CD8' = ANOVA.Tu.ChaoJaccard.CD8)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD8.ChaoJaccardJaccard Tu.xlsx')
rm(dataset_names)
# ---------------- Define plotting function for making clean plots --------------------
# Difficult to make the plots look perfect with stats. Add the statistical bars in Illustrator later I think.
minimally_Defined.Overlap.plots_het <- function(Overlap_Statistics, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, type){
  # MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct statistical  comparisons
  filtered.df <- Overlap_Statistics %>%
    filter(Comparator %in% Comparator_group) %>%
    filter(Comparatee %in% F.MHC_groups) %>%
    filter(value != 1)
  
  filtered.df$Comparatee <- factor(filtered.df$Comparatee, levels = rev(F.MHC_groups))
  
  OL_plot <- ggplot(data=filtered.df, aes(x=Comparatee, y=value, fill=Comparatee))+
    geom_bar(stat="summary", fun=mean)+
    stat_summary(fun.data=mean_se, geom = "linerange")+
    # geom_point(alpha = 0.1 )+
    geom_jitter(shape = 1, width = 0.1, size = 0.5)+
    # geom_signif(comparisons = MHC_list,
    #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
    # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    coord_flip()+
    guides(fill = guide_legend(reverse = TRUE))+
    theme_classic()+
    theme(legend.position = "none")+
    facet_grid(rows = vars(Comparator), scales="free_y", switch = "y", space='free') +
    # facet_grid(. ~ SP, scales = "free", space='free') +
    # facet_wrap(~Comparator, ncol=1, strip.position = "left", scales="free_y") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face="bold", size = 12, color = "black"),
          axis.text.x = element_text(color = "black"),
          strip.text.y.left = element_text(angle=0, size = 20, face="bold"),
          strip.placement = "outside",
          panel.spacing.y=unit(-1.0, "pt"),
          panel.border=element_rect(color="black", fill = NA, size=1)
    ) +
    # theme(aspect.ratio = 0.3) + # This can make the plots all constantly narrower without playing around with the export size.
    # ylim(0, ylimit)+
    scale_y_continuous(breaks=seq(0,1,by=.25), expand = expansion(mult = c(0, .1)), limits = c(0, ylimit))+
    labs(caption=paste0(type), x="Comparatee", y=paste0(type, " similarity coefficient"))
  return(OL_plot)
}
# --------------------------- pretty  het groups (CD4) --------------------
F.MHC_groups <- c(bb, bxf, ff)
Comparator_group <- c(bxf)
bb.ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(bb, bxg7, g7g7)
Comparator_group <- c(bxg7)
bb.g7g7_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(bb, bxs, ss)
Comparator_group <- c(bxs)
ss.bb_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss.bb_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(ff, fxs, ss)
Comparator_group <- c(fxs)
ss.ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss.ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <-c(bb, bxf, ff, bxg7, g7g7, bxs, fxs, ss)
C.MHC_groups <- c(bb_color, bxf_color, ff_color, bxg7_color, g7g7_color, bxs_color, fxs_color, ss_color)
Comparator_group <- c(bxf, bxg7, bxs, fxs)
Jac_hets <- rbind(bb.ff_Jac_CD4_Filt, bb.g7g7_Jac_CD4_Filt, ss.bb_Jac_CD4_Filt, ss.ff_Jac_CD4_Filt)
CJ_hets <- rbind(bb.ff_CJ_CD4_Filt, bb.g7g7_CJ_CD4_Filt, ss.bb_CJ_CD4_Filt, ss.ff_CJ_CD4_Filt)
ylimit = 0.5
Jac_hets.plot <- minimally_Defined.Overlap.plots_het(Jac_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_hets.plot <- minimally_Defined.Overlap.plots_het(CJ_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
compiled_hets.plots <- Jac_hets.plot + CJ_hets.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.MHChet_nostats.pdf", compiled_hets.plots, width=244,height=77, units = "mm")

# Create statistical tables for MHC heterozygous plots
F.MHC_groups <-c(bb, bxf, ff, bxg7, g7g7, bxs, fxs, ss)
Comparator_group <- c(bxf, bxg7, bxs, fxs)
ANOVA.Tu.Jaccard.CD4.het <- Jaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD4.het <- ChaoJaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

dataset_names <- list('ANOVA.Tu.Jaccard.CD4.het' = ANOVA.Tu.Jaccard.CD4.het, 'ANOVA.Tu.ChaoJaccard.CD4.het' = ANOVA.Tu.ChaoJaccard.CD4.het)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.het.ChaoJaccard.Jaccard Tu.xlsx')
rm(dataset_names)
# --------------------------- pretty  het groups (CD4) with added comparator groups  --------------------
F.MHC_groups <- c(bb, bxf, ff)
Comparator_group <- c(bb, bxf, ff)
bb.ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(bb, bxg7, g7g7)
Comparator_group <- c(bb, bxg7, g7g7)
bb.g7g7_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(bb, bxs, ss)
Comparator_group <- c(bb, bxs, ss)
ss.bb_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss.bb_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c(ff, fxs, ss)
Comparator_group <- c(ff, fxs, ss)
ss.ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss.ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

# F1 Plots: bb, bxf, ff
ylimit = 0.5
F.MHC_groups <-c(bb, bxf, ff)
C.MHC_groups <- c(bb_color, bxf_color, ff_color)
Comparator_group <- c(bb, bxf, ff)
Jac_bxf.plot <- minimally_Defined.Overlap.plots_het(bb.ff_Jac_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_bxf.plot <- minimally_Defined.Overlap.plots_het(bb.ff_CJ_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
het.plot.bxf <- Jac_bxf.plot+CJ_bxf.plot

# F1 Plots: bb, bxg7, g7g7
ylimit = 0.5
F.MHC_groups <-c(bb, bxg7, g7g7)
C.MHC_groups <- c(bb_color, bxg7_color, g7g7_color)
Comparator_group <- c(bb, bxg7, g7g7)
Jac_bxg7.plot <- minimally_Defined.Overlap.plots_het(bb.g7g7_Jac_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_bxg7.plot <- minimally_Defined.Overlap.plots_het(bb.g7g7_CJ_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
het.plot.bxg7 <- Jac_bxg7.plot+CJ_bxg7.plot

# F1 Plots: bb, bxs, ss
ylimit = 0.5
F.MHC_groups <-c(bb, bxs, ss)
C.MHC_groups <- c(bb_color, bxs_color, ss_color)
Comparator_group <- c(bb, bxs, ss)
Jac_bxs.plot <- minimally_Defined.Overlap.plots_het(ss.bb_Jac_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_bxs.plot <- minimally_Defined.Overlap.plots_het(ss.bb_CJ_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
het.plot.bxs <- Jac_bxs.plot+CJ_bxs.plot


# F1 Plots: ff, fxs, ss
ylimit = 0.5
F.MHC_groups <-c(ff, fxs, ss)
C.MHC_groups <- c(ff_color, fxs_color, ss_color)
Comparator_group <- c(ff, fxs, ss)
Jac_fxs.plot <- minimally_Defined.Overlap.plots_het(ss.ff_Jac_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_fxs.plot <- minimally_Defined.Overlap.plots_het(ss.ff_CJ_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
het.plot.fxs <- Jac_fxs.plot+CJ_fxs.plot

# Export new F1 plots
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bxf_nostats.pdf", het.plot.bxf, width=226,height=63, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bxg7_nostats.pdf", het.plot.bxg7, width=244,height=63, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bxs_nostats.pdf", het.plot.bxs, width=230,height=63, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.fxs_nostats.pdf", het.plot.fxs, width=224,height=63, units = "mm")

# Create new statistical tables for detailed MHC heterozygous plots
F.MHC_groups <-c(bb, bxf, ff, bxg7, g7g7, bxs, fxs, ss)
Comparator_group <- c(bxf, bxg7, bxs, fxs)





dataset_names <- list(
  'ANOVA.Tu.Jaccard.CD4.bxf' = bb.ff_Jac_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  'ANOVA.Tu.ChaoJaccard.CD4.bxf' = bb.ff_CJ_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  
  'ANOVA.Tu.Jaccard.CD4.bxg7' = bb.g7g7_Jac_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  'ANOVA.Tu.ChaoJaccard.CD4.bxg7' = bb.g7g7_CJ_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  
  'ANOVA.Tu.Jaccard.CD4.bxs' = ss.bb_Jac_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  'ANOVA.Tu.ChaoJaccard.CD4.bxs' = ss.bb_CJ_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),

  'ANOVA.Tu.Jaccard.CD4.fxs' = ss.ff_Jac_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu(),
  'ANOVA.Tu.ChaoJaccard.CD4.fxs' = ss.ff_CJ_CD4_Filt %>% filter(value != 1) %>% ANOVA.Tu()
  )


write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.het_full.ChaoJaccard.Jaccard Tu.xlsx')
rm(dataset_names)


# ------------------ make pretty CD4 tetra included plots -----------------
# bb&ff tetra plots
F.MHC_groups <- c(bb, bxf, ff, tetra_bb.ff)
Comparator_group <- c(bb)
bb.ff_bb_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_bb_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(bxf)
bb.ff_bxf_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_bxf_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(ff)
bb.ff_ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(tetra_bb.ff)
bb.ff_tetra_bb.ff_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.ff_tetra_bb.ff_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
#
F.MHC_groups <-c(bb, bxf, ff, tetra_bb.ff)
C.MHC_groups <- c(bb_color, bxf_color, ff_color, tetra_bb.ff_color)
Comparator_group <- c(bb, bxf, ff, tetra_bb.ff)
Jac_hets <- rbind(bb.ff_bb_Jac_CD4_Filt, bb.ff_bxf_Jac_CD4_Filt, bb.ff_ff_Jac_CD4_Filt, bb.ff_tetra_bb.ff_Jac_CD4_Filt)
CJ_hets <- rbind(bb.ff_bb_CJ_CD4_Filt, bb.ff_bxf_CJ_CD4_Filt, bb.ff_ff_CJ_CD4_Filt, bb.ff_tetra_bb.ff_CJ_CD4_Filt)
ylimit = 0.5
Jac_hets.plot <- minimally_Defined.Overlap.plots_het(Jac_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = 1.0
CJ_hets.plot <- minimally_Defined.Overlap.plots_het(CJ_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
compiled_bb.ff_tetra.plots <- Jac_hets.plot + CJ_hets.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bb.ff_tetra_nostats.pdf", compiled_bb.ff_tetra.plots, width=250,height=93, units = "mm")

# bb&g7g7 tetra plots
F.MHC_groups <- c(bb, bxg7, g7g7, tetra_bb.g7g7)
Comparator_group <- c(bb)
bb.g7g7_bb_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_bb_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(bxg7)
bb.g7g7_bxg7_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_bxg7_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(g7g7)
bb.g7g7_g7g7_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_g7g7_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
Comparator_group <- c(tetra_bb.g7g7)
bb.g7g7_tetra_bb.g7g7_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb.g7g7_tetra_bb.g7g7_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
#
F.MHC_groups <-c(bb, bxg7, g7g7, tetra_bb.g7g7)
C.MHC_groups <- c(bb_color, bxg7_color, g7g7_color, tetra_bb.g7g7_color)
Comparator_group <- c(bb, bxg7, g7g7, tetra_bb.g7g7)
Jac_hets <- rbind(bb.g7g7_bb_Jac_CD4_Filt, bb.g7g7_bxg7_Jac_CD4_Filt, bb.g7g7_g7g7_Jac_CD4_Filt, bb.g7g7_tetra_bb.g7g7_Jac_CD4_Filt)
CJ_hets <- rbind(bb.g7g7_bb_CJ_CD4_Filt, bb.g7g7_bxg7_CJ_CD4_Filt, bb.g7g7_g7g7_CJ_CD4_Filt, bb.g7g7_tetra_bb.g7g7_CJ_CD4_Filt)
ylimit = 0.5
Jac_hets.plot <- minimally_Defined.Overlap.plots_het(Jac_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = 1.0
CJ_hets.plot <- minimally_Defined.Overlap.plots_het(CJ_hets, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
compiled_bb.g7g7_tetra.plots <- Jac_hets.plot + CJ_hets.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bb.g7g7_tetra_nostats.pdf", compiled_bb.g7g7_tetra.plots, width=288,height=93, units = "mm")

# Create statistical tables for tetraparental containing plots
F.MHC_groups <-c(bb, bxf, ff, tetra_bb.ff)
Comparator_group <- c(bb, bxf, ff, tetra_bb.ff)
ANOVA.Tu.Jaccard.CD4.tetra_bb.ff <- Jaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD4.tetra_bb.ff <- ChaoJaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

F.MHC_groups <-c(bb, bxg7, g7g7, tetra_bb.g7g7)
Comparator_group <- c(bb, bxg7, g7g7, tetra_bb.g7g7)
ANOVA.Tu.Jaccard.CD4.tetra_bb.g7g7 <- Jaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD4.tetra_bb.g7g7 <- ChaoJaccard_Overlap_Statistics.CD4 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

dataset_names <- list('Tu.tJac.CD4.tetra_bb.ff' = ANOVA.Tu.Jaccard.CD4.tetra_bb.ff, 'Tu.tCJ.CD4.tetra_bb.ff' = ANOVA.Tu.ChaoJaccard.CD4.tetra_bb.ff,
                      'Tu.tJac.CD4.tetra_bb.g7g7' = ANOVA.Tu.Jaccard.CD4.tetra_bb.g7g7, 'Tu.tCJ.CD4.tetra_bb.g7g7' = ANOVA.Tu.ChaoJaccard.CD4.tetra_bb.g7g7
                      )
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.tetra.ChaoJaccard.Jaccard Tu.xlsx')
rm(dataset_names)
# ------------------- make pretty  CD4 bb vs b+/- plots -----------------
F.MHC_groups <- c(bb, IAbHet)
C.MHC_groups <- c(bb_color, IAbHet_color)
Comparator_group <- c(bb, IAbHet)
bb_IAbHet_Jac_CD4_Filt <- Jaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb_IAbHet_CJ_CD4_Filt <- ChaoJaccard_Overlap_Statistics.CD4 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ylimit = 0.5
Jac_hets.plot <- minimally_Defined.Overlap.plots(bb_IAbHet_Jac_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_hets.plot <- minimally_Defined.Overlap.plots(bb_IAbHet_CJ_CD4_Filt, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
compiled_bb.b.plots <- Jac_hets.plot + CJ_hets.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.bb.vs.b_nostats.pdf", compiled_bb.b.plots, width=244,height=60, units = "mm")

ANOVA.Tu.Jaccard.CD4.bb.vs.IAbhet <- bb_IAbHet_Jac_CD4_Filt %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.ChaoJaccard.CD4.bb.vs.IAbhet <- bb_IAbHet_CJ_CD4_Filt %>%
  filter(value != 1) %>% ANOVA.Tu()

dataset_names <- list('Tu.JacCD4.bb.IAbhet' = ANOVA.Tu.Jaccard.CD4.bb.vs.IAbhet, 'Tu.CJCD4.bb.IAbhet' = ANOVA.Tu.ChaoJaccard.CD4.bb.vs.IAbhet)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.bb.vs.IAbhet.ChaoJaccard.Jaccard Tu.xlsx')
rm(dataset_names)
# --------------- Comparisons of CD4 vs CD8 -------------------------------------------------------------------

# comparison of the CD4 and CD8 samples to show they are not similar to identical haplotypes
# select CD4 and CD8 data of the desired MHC type and remove samples with "one_mouse" from the total rep data
# Combine CD4 and CD8 data
TCR.CD4.CD8_DoVb_TR_P2M.DMR.LFSR <- list(
  data = c(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, TCR_DoVb_CD8_AA.ones.removed.LFSR$data),
  meta = rbind(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta, TCR_DoVb_CD8_AA.ones.removed.LFSR$meta)
)

CD4_CD8_Selected <- repFilter(TCR.CD4.CD8_DoVb_TR_P2M.DMR.LFSR, .method = "by.meta", .query = list(MHC = include("bb", "ff", "g7g7"), Mouse_Count = exclude("1/2 of 2"), `Co-receptor` = exclude("DP")))

CD4_CD8_Selected$meta <- CD4_CD8_Selected$meta %>% mutate(MHC = paste(MHC, `Co-receptor`))

CD4_CD8_Selected_MHC_jaccard <- apply_symm(CD4_CD8_Selected$data, Jaccard_Matrix)
CD4_CD8_Selected_MHC_chao.jaccard <- apply_symm(CD4_CD8_Selected$data, Chao_Jaccard_Matrix)

CD4_CD8_Selected_jaccard_melted_cormat <- Value_Matrix(CD4_CD8_Selected_MHC_jaccard)
CD4_CD8_Selected_chao.jaccard_melted_cormat <- Value_Matrix(CD4_CD8_Selected_MHC_chao.jaccard)

CD4_CD8_Selected_Jaccard_Overlap_Stat <- Overlap.Stats(CD4_CD8_Selected_jaccard_melted_cormat, CD4_CD8_Selected)
CD4_CD8_Selected_ChaoJaccard_Overlap_Stat <- Overlap.Stats(CD4_CD8_Selected_chao.jaccard_melted_cormat, CD4_CD8_Selected)

# Set up data for plots
F.MHC_groups <- c("bb CD4", "bb CD8")
C.MHC_groups <- c(bb_color, "#2A507E")
Comparator_group <- c("bb CD4", "bb CD8")
# ylimit = 0.6
# bb_CD4_CD8_Selected_jac.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_Jaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
# ylimit = 1.2
# bb_CD4_CD8_Selected_CJ.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_ChaoJaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
# # Set up pretty plot
bb_CD4_CD8_Jac_CD4_Filt <- CD4_CD8_Selected_Jaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb_CD4_CD8_CJ_CD4_Filt <- CD4_CD8_Selected_ChaoJaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

F.MHC_groups <- c("ff CD4", "ff CD8")
C.MHC_groups <- c(ff_color, "#835D06")
Comparator_group <- c("ff CD4", "ff CD8")
# ylimit = 0.6
# ff_CD4_CD8_Selected_jac.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_Jaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
# ylimit = 1.2
# ff_CD4_CD8_Selected_CJ.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_ChaoJaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
# # Set up pretty plot
ff_CD4_CD8_Jac_CD4_Filt <- CD4_CD8_Selected_Jaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ff_CD4_CD8_CJ_CD4_Filt <- CD4_CD8_Selected_ChaoJaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)

# make plots
F.MHC_groups <- c("g7g7 CD4", "g7g7 CD8")
C.MHC_groups <- c(g7g7_color, "#5F0A2E")
Comparator_group <- c("g7g7 CD4", "g7g7 CD8")
# ylimit = 0.6
# g7g7_CD4_CD8_Selected_jac.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_Jaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
# ylimit = 1.2
# g7g7_CD4_CD8_Selected_CJ.plot <- minimally_Defined.Overlap.plots(CD4_CD8_Selected_ChaoJaccard_Overlap_Stat, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
# # Set up pretty plot
g7g7_CD4_CD8_Jac_CD4_Filt <- CD4_CD8_Selected_Jaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
g7g7_CD4_CD8_CJ_CD4_Filt <- CD4_CD8_Selected_ChaoJaccard_Overlap_Stat %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)


# # make stats  plots
# pdf("Figure Export/Repertoire Overlap/Bar Graphs/CD4_CD8_Stats.pdf", width = 10, height = 6)
# (bb_CD4_CD8_Selected_jac.plot + bb_CD4_CD8_Selected_CJ.plot) /
# (ff_CD4_CD8_Selected_jac.plot + ff_CD4_CD8_Selected_CJ.plot) /
# (g7g7_CD4_CD8_Selected_jac.plot + g7g7_CD4_CD8_Selected_CJ.plot)
# dev.off()

# make pretty plots
F.MHC_groups <- c("bb CD4", "bb CD8", "ff CD4", "ff CD8", "g7g7 CD4", "g7g7 CD8")
C.MHC_groups <- c(bb_color, "#2A507E", ff_color, "#835D06", g7g7_color, "#5F0A2E")
Comparator_group <- c("bb CD4", "bb CD8", "ff CD4", "ff CD8", "g7g7 CD4", "g7g7 CD8")
Jac_CD4.CD8 <- rbind(bb_CD4_CD8_Jac_CD4_Filt, ff_CD4_CD8_Jac_CD4_Filt, g7g7_CD4_CD8_Jac_CD4_Filt)
CJ_CD4.CD8 <- rbind(bb_CD4_CD8_CJ_CD4_Filt, ff_CD4_CD8_CJ_CD4_Filt, g7g7_CD4_CD8_CJ_CD4_Filt)

ylimit = 0.5
Jac_CD4.CD8.plot <- minimally_Defined.Overlap.plots_het(Jac_CD4.CD8, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = NA
CJ_CD4.CD8.plot <- minimally_Defined.Overlap.plots_het(CJ_CD4.CD8, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Chao)
compiled_CD4.CD8.plots <- Jac_CD4.CD8.plot + CJ_CD4.CD8.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4_CD8_Stats_nostats.pdf", compiled_CD4.CD8.plots, width=296,height=75, units = "mm")

Ttest.Jaccard.hom.CD4.vs.CD8 <- Jac_CD4.CD8 %>%
  filter(value != 1) %>% Simple.T.test()

Ttest.ChaoJaccard.hom.CD4.vs.CD8 <- CJ_CD4.CD8 %>%
  filter(value != 1) %>% Simple.T.test()

dataset_names <- list('Ttest Jac CD4 vs CD8' = Ttest.Jaccard.hom.CD4.vs.CD8, 'Ttest CJ CD4 vs CD8' = Ttest.ChaoJaccard.hom.CD4.vs.CD8)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4 vs CD8.ChaoJaccard.Jaccard Ttest.xlsx')
rm(dataset_names)

# ----------------------------- Calculate Jaccard for CD4 raw CDR3 and VJ only -----------------------------
# Data input:
# TCR_DoVb_CD4_AA.ones.removed.LFSR


# NEW Jaccard Index Functions for CDR3 and VJ only   (Correct function based on Pippas Data)
Jaccard_Matrix_CDR3 <- function(x, y) {
  x_merged <- within(x, Merge <- paste(CDR3.aa))
  y_merged <- within(y, Merge <- paste(CDR3.aa))
  x_merged$Origin <- 1
  y_merged$Origin <- 1
  x_merged <- x_merged %>% uncount(Clones)
  y_merged <- y_merged %>% uncount(Clones)  
  
  x_merged <- x_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  y_merged <- y_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  
  J.numerator <- length(dplyr::intersect(x_merged$Merge, y_merged$Merge))
  x.row <- nrow(x_merged)
  y.row <- nrow(y_merged)
  Jaccard_index <- J.numerator/(x.row+y.row-J.numerator)
  # This is equivalent to:
  # (|A∩B|)/(|A|+|B|- |A∩B|)
}
Jaccard_Matrix_VJ <- function(x, y) {
  x_merged <- within(x, Merge <- paste(V.name, J.name, sep = "_"))
  y_merged <- within(y, Merge <- paste(V.name, J.name, sep = "_"))
  x_merged$Origin <- 1
  y_merged$Origin <- 1
  x_merged <- x_merged %>% uncount(Clones)
  y_merged <- y_merged %>% uncount(Clones)  
  
  x_merged <- x_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  y_merged <- y_merged  %>%
    as.data.table() %>%
    .[,lapply(.SD,sum),by=Merge,.SDcols=c("Origin")] %>% # Calculate sum of Origin in .SD grouped by Merge. 20x faster than dplyr
    as_tibble() %>% # convert back to tibble and drop values = 1
    filter(Origin > 1)
  
  J.numerator <- length(dplyr::intersect(x_merged$Merge, y_merged$Merge))
  x.row <- nrow(x_merged)
  y.row <- nrow(y_merged)
  Jaccard_index <- J.numerator/(x.row+y.row-J.numerator)
  # This is equivalent to:
  # (|A∩B|)/(|A|+|B|- |A∩B|)
}

jaccard.CD4.CDR3 <- apply_symm(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, Jaccard_Matrix_CDR3)
jaccard.CD4.VJ <- apply_symm(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, Jaccard_Matrix_VJ)
melted_cormat_j.CD4.CDR3 <- Value_Matrix(jaccard.CD4.CDR3)
melted_cormat_j.CD4.VJ <- Value_Matrix(jaccard.CD4.VJ)
Jaccard_Overlap_Statistics.CD4.CDR3 <- Overlap.Stats(melted_cormat_j.CD4.CDR3, TCR_DoVb_CD4_AA.ones.removed.LFSR)
Jaccard_Overlap_Statistics.CD4.VJ <- Overlap.Stats(melted_cormat_j.CD4.VJ, TCR_DoVb_CD4_AA.ones.removed.LFSR)

Jaccard_Overlap_Statistics.CD4.CDR3 <- Jaccard_Overlap_Statistics.CD4.CDR3 %>% mutate(across(Comparator:Comparatee,~ gsub("tetraparental", "tetra.", .))) %>% mutate(across(Comparator:Comparatee, as.factor))
Jaccard_Overlap_Statistics.CD4.VJ <- Jaccard_Overlap_Statistics.CD4.VJ %>% mutate(across(Comparator:Comparatee,~ gsub("tetraparental", "tetra.", .))) %>% mutate(across(Comparator:Comparatee, as.factor))

# Plot CD4 VJ for the MHC hom comparisons
F.MHC_groups <- c(bb, ff, g7g7, ss)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color)
Comparator_group <- c(bb, ff, g7g7, ss)
ylimit = 0.50
MHChom_OL_CD4_Jac.V.CDR3.J.plot <- minimally_Defined.Overlap.plots(Jaccard_Overlap_Statistics.CD4, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac) + labs(caption = "Jaccard: V-CDR3-J", title = "TCR Alpha Repertoire Overlap:CD4")
ylimit = 0.50
MHChom_OL_CD4_Jac.CDR3.plot <- minimally_Defined.Overlap.plots(Jaccard_Overlap_Statistics.CD4.CDR3, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Jac) + labs(caption = "Jaccard: CDR3")
ylimit = 1.0
MHChom_OL_CD4_Jac.VJ.plot <- minimally_Defined.Overlap.plots(Jaccard_Overlap_Statistics.CD4.VJ, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Jac) + labs(caption = "Jaccard: V-J")

compiled_MHC_hom.CD4.versus <- MHChom_OL_CD4_Jac.CDR3.plot + MHChom_OL_CD4_Jac.VJ.plot
ggsave("Figure Export/Repertoire Overlap/Jaccard V.CDR3.J vs CDR3 vs. VJ/CD4.Jaccard.MHChom_versus_no stats.pdf", compiled_MHC_hom.CD4.versus, width=245,height=95, units = "mm")

Title_CD4 <- ggtitle("TCR Alpha Repertoire Overlap:CD4")
compiled_MHC_hom.CD4 <- MHChom_OL_CD4_Jac.plot + Title_CD4 + MHChom_OL_CD4_CJ.plot
compiled_MHC_hom.CD8 <- MHChom_OL_CD8_Jac.plot + Title_CD8 + MHChom_OL_CD8_CJ.plot
compiled_MHC_hom <- (MHChom_OL_CD4_Jac.plot + Title_CD4 + MHChom_OL_CD4_CJ.plot) /
  (MHChom_OL_CD8_Jac.plot + Title_CD8 + MHChom_OL_CD8_CJ.plot)

# make plots
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.MHChom.pdf", compiled_MHC_hom.CD4, width=190,height=80, units = "mm")

# Create statistical tables for MHC heterozygous plots for CD4 Jaccard VJ and Jaccard CDR3
F.MHC_groups <- c(bb, ff, g7g7, ss)
Comparator_group <- c(bb, ff, g7g7, ss)
ANOVA.Tu.Jac_hom.CDR3.CD4 <- Jaccard_Overlap_Statistics.CD4.CDR3 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.Jac_hom.VJ.CD4 <- Jaccard_Overlap_Statistics.CD4.VJ %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

dataset_names <- list('ANOVA.Tu.Jac_hom.CDR3.CD4' = ANOVA.Tu.Jac_hom.CDR3.CD4, 'ANOVA.Tu.Jac_hom.VJ.CD4' = ANOVA.Tu.Jac_hom.VJ.CD4)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.hom.CDR3_Jaccard.and.VJ_Jaccard Tu.xlsx')
rm(dataset_names)


# --------------- Het data VJ, CDR3----------
F.MHC_groups <- c(bb, bxf, ff)
Comparator_group <- c(bxf)
bb_ff_Jac_CD4.CDR3_Filt <- Jaccard_Overlap_Statistics.CD4.CDR3 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb_ff_Jac_CD4.VJ_Filt <- Jaccard_Overlap_Statistics.CD4.VJ %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)


F.MHC_groups <- c(bb, bxg7, g7g7)
Comparator_group <- c(bxg7)
bb_g7g7_Jac_CD4.CDR3_Filt <- Jaccard_Overlap_Statistics.CD4.CDR3 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
bb_g7g7_Jac_CD4.VJ_Filt <- Jaccard_Overlap_Statistics.CD4.VJ %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)


F.MHC_groups <- c(bb, bxs, ss)
Comparator_group <- c(bxs)
ss_bb_Jac_CD4.CDR3_Filt <- Jaccard_Overlap_Statistics.CD4.CDR3 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss_bb_Jac_CD4.VJ_Filt <- Jaccard_Overlap_Statistics.CD4.VJ %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)


F.MHC_groups <- c(ff, fxs, ss)
Comparator_group <- c(fxs)
ss_ff_Jac_CD4.CDR3_Filt <- Jaccard_Overlap_Statistics.CD4.CDR3 %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)
ss_ff_Jac_CD4.VJ_Filt <- Jaccard_Overlap_Statistics.CD4.VJ %>% filter(Comparator %in% Comparator_group) %>% filter(Comparatee %in% F.MHC_groups)


F.MHC_groups <-c(bb, bxf, ff, bxg7, g7g7, bxs, fxs, ss)
C.MHC_groups <- c(bb_color, bxf_color, ff_color, bxg7_color, g7g7_color, bxs_color, fxs_color, ss_color)
Comparator_group <- c(bxf, bxg7, bxs, fxs)
Jac_hets.CDR3 <- rbind(bb_ff_Jac_CD4.CDR3_Filt, bb_g7g7_Jac_CD4.CDR3_Filt, ss_bb_Jac_CD4.CDR3_Filt, ss_ff_Jac_CD4.CDR3_Filt)
Jac_hets.VJ <- rbind(bb_ff_Jac_CD4.VJ_Filt, bb_g7g7_Jac_CD4.VJ_Filt, ss_bb_Jac_CD4.VJ_Filt, ss_ff_Jac_CD4.VJ_Filt)
ylimit = 0.5
Jac_hets.CDR3.plot <- minimally_Defined.Overlap.plots_het(Jac_hets.CDR3, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac) + labs(caption = "Jaccard: CDR3")
ylimit = 0.95
Jac_hets.VJ.plot <- minimally_Defined.Overlap.plots_het(Jac_hets.VJ, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac) +  labs(caption = "Jaccard: V-J")
compiled_hets.plots_VJ_CDR3 <- Jac_hets.CDR3.plot + Jac_hets.VJ.plot
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_hets.plots_VJ_CDR3.pdf", compiled_hets.plots_VJ_CDR3, width=244,height=77, units = "mm")

# Create statistical tables for MHC heterozygous plots for CD4 Jaccard VJ and Jaccard CDR3
F.MHC_groups <-c(bb, bxf, ff, bxg7, g7g7, bxs, fxs, ss)
Comparator_group <- c(bxf, bxg7, bxs, fxs)
ANOVA.Tu.Jac_hets.CDR3.CD4 <- Jac_hets.CDR3 %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

ANOVA.Tu.Jac_hets.VJ.CD4 <- Jac_hets.VJ %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

dataset_names <- list('Tu.Jac_hets.CDR3.CD4' = ANOVA.Tu.Jac_hets.CDR3.CD4, 'Tu.Jac_hets.VJ.CD4' = ANOVA.Tu.Jac_hets.VJ.CD4)
write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4.het.CDR3_Jaccard.and.VJ_Jaccard Tu.xlsx')
rm(dataset_names)



### WORKING SECTION OF CODE/PLOTS BUT NOT IN USE FOR MANUSCRIPT
### # -------------------------- MHChom1 ∩ MHChet1 vs. MHChom1 ∩ MHChet2  (Intersection)----------------------------------------
 # first identify the sequences in MHChom that are in the MHC het and compare these intersected MHC hets vs. jaccard and chao jaccard.
 # If these sequences are largely not the same between the two hets then it suggests that negative seleciton is the dominat force in the presence
 # of a competing secondary MHC molecule.
 # Function takes in the ImmuneArch dataframe list and then subsets based on the co-receptor value.
 filter_intersect <- function(x, TCR.query){
   x <- x %>%
     dplyr::mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
     filter(V.CDR3.J %in% TCR.query) %>%
     select(-V.CDR3.J)
   return(x)
 }
 intersected.TCR <- function(matrix.pr, filtered.rep, MHC_hom.Samples, MHC_het.Samples){
   # identify the sequences shared in the target MHC_hom and MHC_het haplotypes
   TCR.query <- matrix.pr %>%
     select(V.name, CDR3.aa, J.name, {{MHC_hom.Samples}}, {{MHC_het.Samples}}) %>%
     filter({{MHC_hom.Samples}} > 0) %>%
     filter({{MHC_het.Samples}} > 0) %>%
     mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
     pull(V.CDR3.J)
   
   # Filter data for the specific clones identified as shared in TCR.query
   filtered.rep$data <- lapply(filtered.rep$data, filter_intersect,TCR.query=TCR.query)
   
   return(filtered.rep)
 }
 
 # pull out the haplotypes of interest using repFilter()
 bxf.in.bb_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 bxf.in.bb_rep <- repFilter(bxf.in.bb_rep, .method = "by.meta", .query = list(MHC = include("bb", "bxf")))
 bxg7.in.bb_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 bxg7.in.bb_rep <- repFilter(bxg7.in.bb_rep, .method = "by.meta", .query = list(MHC = include("bb", "bxg7")))
 bxs.in.bb_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 bxs.in.bb_rep <- repFilter(bxs.in.bb_rep, .method = "by.meta", .query = list(MHC = include("bb", "bxs")))
 # 
 bxf.in.ff_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 bxf.in.ff_rep <- repFilter(bxf.in.ff_rep, .method = "by.meta", .query = list(MHC = include("ff", "bxf")))
 fxs.in.ff_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 fxs.in.ff_rep <- repFilter(fxs.in.ff_rep, .method = "by.meta", .query = list(MHC = include("ff", "fxs")))
 
 # Apply intersection filter functions to target haplotypes
 bxf.in.bb_intersected <- intersected.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxf.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxf.Samples)
 bxg7.in.bb_intersected <- intersected.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxg7.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxg7.Samples)
 bxs.in.bb_intersected <- intersected.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxs.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxs.Samples)
 bxf.in.ff_intersected <- intersected.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxf.in.ff_rep, MHC_hom.Samples = ff.Samples, MHC_het.Samples = bxf.Samples)
 fxs.in.ff_intersected <- intersected.TCR(pr_CD4s_filtered.LFSR.ones.removed, fxs.in.ff_rep, MHC_hom.Samples = ff.Samples, MHC_het.Samples = fxs.Samples)
 
 # Drop MHC homozygous samples
 bxf.in.bb_intersected <- repFilter(bxf.in.bb_intersected, .method = "by.meta", .query = list(MHC = include("bxf")))
 bxg7.in.bb_intersected <- repFilter(bxg7.in.bb_intersected, .method = "by.meta", .query = list(MHC = include("bxg7")))
 bxs.in.bb_intersected <- repFilter(bxs.in.bb_intersected, .method = "by.meta", .query = list(MHC = include("bxs")))
 
 bxf.in.ff_intersected <- repFilter(bxf.in.ff_intersected, .method = "by.meta", .query = list(MHC = include("bxf")))
 fxs.in.ff_intersected <- repFilter(fxs.in.ff_intersected, .method = "by.meta", .query = list(MHC = include("fxs")))
 
 # Also extract the un modified MHC homozygous references
 bb_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 bb_rep <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb")))
 
 ff_rep <- TCR_DoVb_CD4_AA.ones.removed.LFSR
 ff_rep <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("ff")))
 
 
 # Bind up data into two new data lists
 # Data list 1
 `TCR.CD4_bb.∩.bxf&bxg7&bxs` <- list(
   data = c(
     bb_rep$data,
     bxf.in.bb_intersected$data,
     bxg7.in.bb_intersected$data,
     bxs.in.bb_intersected$data),
   meta = rbind(
     bb_rep$meta,
     bxf.in.bb_intersected$meta,
     bxg7.in.bb_intersected$meta,
     bxs.in.bb_intersected$meta)
 )
 # Data list 2
 `TCR.CD4_ff.∩.bxf&fxs` <- list(
   data = c(
     ff_rep$data,
     bxf.in.ff_intersected$data,
     fxs.in.ff_intersected$data),
   meta = rbind(
     ff_rep$meta,
     bxf.in.ff_intersected$meta,
     fxs.in.ff_intersected$meta)
 )
 TCR.CD4_bb.F1.intersections <- `TCR.CD4_bb.∩.bxf&bxg7&bxs`
 TCR.CD4_ff.F1.intersections <- `TCR.CD4_ff.∩.bxf&fxs`
 
 save(TCR.CD4_bb.F1.intersections, TCR.CD4_ff.F1.intersections, file = "Intersection.CD4.rep.RData")
 
 
### 
### # Run Jaccard and Chao Jaccard transformations on these lists
### # Generate Jaccard anc Chao Jaccard Value Matrix:
### `jaccard.CD4.bb.∩.bxf&bxg7&bxs` <- apply_symm(`TCR.CD4_bb.∩.bxf&bxg7&bxs`$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix
### `chao_jaccard.CD4.bb.∩.bxf&bxg7&bxs` <- apply_symm(`TCR.CD4_bb.∩.bxf&bxg7&bxs`$data, Chao_Jaccard_Matrix)
### `jaccard.CD4.ff.∩.bxf&fxs` <- apply_symm(`TCR.CD4_ff.∩.bxf&fxs`$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix
### `chao_jaccard.CD4.ff.∩.bxf&fxs` <- apply_symm(`TCR.CD4_ff.∩.bxf&fxs`$data, Chao_Jaccard_Matrix)
### 
### `melted_jaccard.CD4.bb.∩.bxf&bxg7&bxs` <- Value_Matrix(`jaccard.CD4.bb.∩.bxf&bxg7&bxs`)
### `melted_chao_jaccard.CD4.bb.∩.bxf&bxg7&bxs` <- Value_Matrix(`chao_jaccard.CD4.bb.∩.bxf&bxg7&bxs`)
### `melted_jaccard.CD4.ff.∩.bxf&fxs` <- Value_Matrix(`jaccard.CD4.ff.∩.bxf&fxs`)
### `melted_chao_jaccard.CD4.ff.∩.bxf&fxs` <- Value_Matrix(`chao_jaccard.CD4.ff.∩.bxf&fxs`)
### 
### `Jaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs` <- Overlap.Stats(`melted_jaccard.CD4.bb.∩.bxf&bxg7&bxs`, `TCR.CD4_bb.∩.bxf&bxg7&bxs`)
### `ChaoJaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs` <- Overlap.Stats(`melted_chao_jaccard.CD4.bb.∩.bxf&bxg7&bxs`, `TCR.CD4_bb.∩.bxf&bxg7&bxs`)
### 
### `Jaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` <- Overlap.Stats(`melted_jaccard.CD4.ff.∩.bxf&fxs`, `TCR.CD4_ff.∩.bxf&fxs`)
### `ChaoJaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` <- Overlap.Stats(`melted_chao_jaccard.CD4.ff.∩.bxf&fxs`, `TCR.CD4_ff.∩.bxf&fxs`)
### 
### # generate plots:
### bxf_color_intersect_bb = "#FF847F"
### bxs_color_intersect = "#A9FF13"
### bxg7_color_intersect = "#BC62FC"
### 
### bxf_color_intersect_ff = "#DBBA19"
### fxs_color_intersect = "#CB2D15"
### 
### 
### # F.MHC_groups <- c(bb, bxf, bxg7, bxs)
### # C.MHC_groups <- c(bb_color, bxf_color, bxg7_color, bxs_color)
### # Comparator_group <- c(bb, bxf, bxg7, bxs)
### F.MHC_groups <- c(bxf, bxg7, bxs)
### C.MHC_groups <- c(bxf_color_intersect_bb, bxg7_color_intersect, bxs_color_intersect)
### Comparator_group <- c(bxf, bxg7, bxs)
### ylimit = 0.5
### `intersect.bb.∩.bxf&bxg7&bxs_CD4_Jac.plot` <- minimally_Defined.Overlap.plots(`Jaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs`, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
### ylimit = 1.0
### `intersect.bb.∩.bxf&bxg7&bxs_CD4_CJ.plot` <- minimally_Defined.Overlap.plots(`ChaoJaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs`, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
### # 
### # F.MHC_groups <- c(ff, bxf, fxs)
### # C.MHC_groups <- c(ff_color, bxf_color, fxs_color)
### # Comparator_group <- c(ff, bxf, fxs)
### F.MHC_groups <- c(bxf, fxs)
### C.MHC_groups <- c(bxf_color_intersect_ff, fxs_color_intersect)
### Comparator_group <- c(bxf, fxs)
### ylimit = 0.5
### `intersect.ff.∩.bxf&fxs_CD4_Jac.plot` <- minimally_Defined.Overlap.plots(`Jaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs`, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
### ylimit = 1.0
### `intersect.ff.∩.bxf&fxs_CD4_CJ.plot` <- minimally_Defined.Overlap.plots(`ChaoJaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs`, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
### 
### # Display Plots:
### compiled_intersection_bb_hets <- `intersect.bb.∩.bxf&bxg7&bxs_CD4_Jac.plot`+`intersect.bb.∩.bxf&bxg7&bxs_CD4_CJ.plot`
### compiled_intersection_ff_hets <- `intersect.ff.∩.bxf&fxs_CD4_Jac.plot`+`intersect.ff.∩.bxf&fxs_CD4_CJ.plot`
### 
### # Plot for including the MHC homozygous raw data
### # ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_intersection_bb_hets.pdf", compiled_intersection_bb_hets, width=245,height=93, units = "mm")
### # ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_intersection_ff_hets.pdf", compiled_intersection_ff_hets, width=225,height=62, units = "mm")
### 
### # MHC homozygous excluded
### ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_intersection_bb_hets_v2.pdf", compiled_intersection_bb_hets, width=245,height=67, units = "mm")
### ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_intersection_ff_hets_v2.pdf", compiled_intersection_ff_hets, width=225,height=40, units = "mm")
### 
### 
### # Create statistical tables for MHChom1 ∩ MHChet1 vs. MHChom1 ∩ MHChet2  Intersections
### 
### # bb.∩.bxf&bxg7&bxs
### # MHC homozygous data is commented out here in statistical comparisons
### F.MHC_groups <- c(bxf, bxg7, bxs)
### Comparator_group <- c(bxf, bxg7, bxs)
### `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs` %>%
###   filter(Comparator %in% Comparator_group) %>%
###   filter(Comparatee %in% F.MHC_groups) %>%
###   filter(value != 1) %>% ANOVA.Tu()
### 
### `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.bb.∩.bxf&bxg7&bxs` %>%
###   filter(Comparator %in% Comparator_group) %>%
###   filter(Comparatee %in% F.MHC_groups) %>%
###   filter(value != 1) %>% ANOVA.Tu()
### 
### # # ff.∩.bxf&fxs
### # F.MHC_groups <- c(ff, bxf, fxs)
### # Comparator_group <- c(ff, bxf, fxs)
### # `ANOVA.Tu.intersect.ff.∩.bxf&fxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` %>%
### #   filter(Comparator %in% Comparator_group) %>%
### #   filter(Comparatee %in% F.MHC_groups) %>%
### #   filter(value != 1) %>% ANOVA.Tu()
### # 
### # `ANOVA.Tu.intersect.ff.∩.bxf&fxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` %>%
### #   filter(Comparator %in% Comparator_group) %>%
### #   filter(Comparatee %in% F.MHC_groups) %>%
### #   filter(value != 1) %>% ANOVA.Tu()
### 
### F.MHC_groups <- c(bxf, fxs)
### Comparator_group <- c(bxf, fxs)
### `Ttest.intersect.ff.∩.bxf&fxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` %>%
###   filter(Comparator %in% Comparator_group) %>%
###   filter(Comparatee %in% F.MHC_groups) %>%
###   filter(value != 1) %>% Simple.T.test()
### 
### `Ttest.intersect.ff.∩.bxf&fxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.ff.∩.bxf&fxs` %>%
###   filter(Comparator %in% Comparator_group) %>%
###   filter(Comparatee %in% F.MHC_groups) %>%
###   filter(value != 1) %>% Simple.T.test()
### 
### dataset_names <- list('Ttest Jac CD4 vs CD8' = Ttest.Jaccard.hom.CD4.vs.CD8, 'Ttest CJ CD4 vs CD8' = Ttest.ChaoJaccard.hom.CD4.vs.CD8)
### write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/CD4 vs CD8.ChaoJaccard.Jaccard Ttest.xlsx')
### rm(dataset_names)
### 
### 
### # dataset_names <- list('Tu.Jac.bb.∩.bxf&bxg7&bxs' = `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_Jac`, 'Tu.Jac.ff.∩.bxf&fxs' = `ANOVA.Tu.intersect.ff.∩.bxf&fxs_CD4_Jac`,
### #                       'Tu.CJ.bb.∩.bxf&bxg7&bxs' = `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_CJ`, 'Tu.CJ.ff.∩.bxf&fxs' = `ANOVA.Tu.intersect.ff.∩.bxf&fxs_CD4_CJ`)
### 
### dataset_names_v2 <- list('Tu.Jac.bb.∩.bxf&bxg7&bxs' = `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_Jac`, 'Ttest.Jac.ff.∩.bxf&fxs' = `Ttest.intersect.ff.∩.bxf&fxs_CD4_Jac`,
###                       'Tu.CJ.bb.∩.bxf&bxg7&bxs' = `ANOVA.Tu.intersect.bb.∩.bxf&bxg7&bxs_CD4_CJ`, 'Ttest.CJ.ff.∩.bxf&fxs' = `Ttest.intersect.ff.∩.bxf&fxs_CD4_CJ`)
### 
### 
### # write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/Intersections_Jaccard.and.VJ_Jaccard Tu.xlsx') # MHC homozygous data included
### write.xlsx(dataset_names_v2, file = 'Figure Export/Repertoire Overlap/Statistical Tables/Intersections_Jaccard.and.VJ_Jaccard Tu_v2.xlsx') # MHC homozygous data excluded
### rm(dataset_names)
### rm(dataset_names_v2)
### 
### # Conclusions
### # the ff ∩ b×f and ff ∩ f×s intersected sequences are quite different. ff ⋒ are different
### # the bb ∩ b×f and bb ∩ b×g7 and bb ∩ b×s are simmilar. bb ⋒ are similar


# -------------------------- MHChom1-MHChet1 vs. MHChom1-MHChet2 (Relative Complement A-B) ----------------------------------------
# Compare sequences in b absent in bxf VS. sequences in b absent in bxg7
# This will tell us how different negative selection is because the f and g7 molecules should be deleting different b selected sequences.
# create a new function which excludeds intersected sequences from the bb samples to find the relative complement
# Function takes in the ImmuneArch dataframe list and then subsets based on the co-receptor value.
exclude_intersect <- function(x, TCR.query){
  x <- x %>%
    dplyr::mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    filter(V.CDR3.J %in% TCR.query == FALSE) %>%
    select(-V.CDR3.J)
  return(x)
}
Relative_Complement.TCR <- function(matrix.pr, filtered.rep, MHC_hom.Samples, MHC_het.Samples){
  # identify the sequences shared in the target MHC_hom and MHC_het haplotypes
  TCR.query <- matrix.pr %>%
    select(V.name, CDR3.aa, J.name, {{MHC_hom.Samples}}, {{MHC_het.Samples}}) %>%
    filter({{MHC_hom.Samples}} > 0) %>%
    filter({{MHC_het.Samples}} > 0) %>%
    mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    pull(V.CDR3.J)
  
  # exclude the clones which are in the TCR.query column as these are shared between the F1 and parent and we don't want those this time.
  filtered.rep$data <- lapply(filtered.rep$data, exclude_intersect,TCR.query=TCR.query)
  
  return(filtered.rep)
}

# Apply intersection filter functions to target haplotypes
bxf.in.bb_Relative_Complement <- Relative_Complement.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxf.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxf.Samples)
bxg7.in.bb_Relative_Complement <- Relative_Complement.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxg7.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxg7.Samples)
bxs.in.bb_Relative_Complement <- Relative_Complement.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxs.in.bb_rep, MHC_hom.Samples = bb.Samples, MHC_het.Samples = bxs.Samples)
bxf.in.ff_Relative_Complement <- Relative_Complement.TCR(pr_CD4s_filtered.LFSR.ones.removed, bxf.in.ff_rep, MHC_hom.Samples = ff.Samples, MHC_het.Samples = bxf.Samples)
fxs.in.ff_Relative_Complement <- Relative_Complement.TCR(pr_CD4s_filtered.LFSR.ones.removed, fxs.in.ff_rep, MHC_hom.Samples = ff.Samples, MHC_het.Samples = fxs.Samples)

# # Check if this did anything.
# bxf.in.bb_intersected$data$`388 postJ` %>% head() #1586 x 15
# bxf.in.bb_Relative_Complement$data$`388 postJ` %>% head() #1320 x 15
# TCR_DoVb_CD4_AA.ones.removed.LFSR$data$`388 postJ`%>% head() #2906 x 15

# Drop MHC heterozygous samples
bxf.in.bb_Relative_Complement <- repFilter(bxf.in.bb_Relative_Complement, .method = "by.meta", .query = list(MHC = include("bb")))
bxg7.in.bb_Relative_Complement <- repFilter(bxg7.in.bb_Relative_Complement, .method = "by.meta", .query = list(MHC = include("bb")))
bxs.in.bb_Relative_Complement <- repFilter(bxs.in.bb_Relative_Complement, .method = "by.meta", .query = list(MHC = include("bb")))
bxf.in.ff_Relative_Complement <- repFilter(bxf.in.ff_Relative_Complement, .method = "by.meta", .query = list(MHC = include("ff")))
fxs.in.ff_Relative_Complement <- repFilter(fxs.in.ff_Relative_Complement, .method = "by.meta", .query = list(MHC = include("ff")))

# mutate the MHC column in the metadata.
# This is not explicitly correct but we can correct the names to be exactly as we want after the figures are plotted.
bxf.in.bb_Relative_Complement$meta <- bxf.in.bb_Relative_Complement$meta %>% mutate(MHC=recode(MHC, 'bb'='bxf'))
bxg7.in.bb_Relative_Complement$meta <- bxg7.in.bb_Relative_Complement$meta %>% mutate(MHC=recode(MHC, 'bb'='bxg7'))
bxs.in.bb_Relative_Complement$meta <- bxs.in.bb_Relative_Complement$meta %>% mutate(MHC=recode(MHC, 'bb'='bxs'))
bxf.in.ff_Relative_Complement$meta <- bxf.in.ff_Relative_Complement$meta %>% mutate(MHC=recode(MHC, 'ff'='bxf'))
fxs.in.ff_Relative_Complement$meta <- fxs.in.ff_Relative_Complement$meta %>% mutate(MHC=recode(MHC, 'ff'='fxs'))

# looks like there is an issue with the samples having the same name so they are not being compared correctly. So the metadata and the coresponding data frame names need to be modified.
rename_RC <- function(TCRlist){
  n <- TCRlist$data %>% names()
  n_mod.RC <- paste(n, "mod.RC", TCRlist$meta$MHC, sep = "_")
  names(TCRlist$data) <- n_mod.RC # new names for dataframes
  TCRlist$meta$Sample <- n_mod.RC # rename metadata columns
  return(TCRlist)
}
bxf.in.bb_Relative_Complement <- rename_RC(bxf.in.bb_Relative_Complement)
bxg7.in.bb_Relative_Complement <- rename_RC(bxg7.in.bb_Relative_Complement)
bxs.in.bb_Relative_Complement <- rename_RC(bxs.in.bb_Relative_Complement)
bxf.in.ff_Relative_Complement <- rename_RC(bxf.in.ff_Relative_Complement)
fxs.in.ff_Relative_Complement <- rename_RC(fxs.in.ff_Relative_Complement)

# Bind up data into two new data lists
# Data list 1
`TCR.CD4_bb-bxf&bxg7&bxs` <- list(
  data = c(
    bb_rep$data,
    bxf.in.bb_Relative_Complement$data,
    bxg7.in.bb_Relative_Complement$data,
    bxs.in.bb_Relative_Complement$data),
  meta = rbind(
    bb_rep$meta,
    bxf.in.bb_Relative_Complement$meta,
    bxg7.in.bb_Relative_Complement$meta,
    bxs.in.bb_Relative_Complement$meta)
)
# Data list 2
`TCR.CD4_ff-bxf&fxs` <- list(
  data = c(
    ff_rep$data,
    bxf.in.ff_Relative_Complement$data,
    fxs.in.ff_Relative_Complement$data),
  meta = rbind(
    ff_rep$meta,
    bxf.in.ff_Relative_Complement$meta,
    fxs.in.ff_Relative_Complement$meta)
)

`jaccard.CD4.bb-bxf&bxg7&bxs` <- apply_symm(`TCR.CD4_bb-bxf&bxg7&bxs`$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix
`chao_jaccard.CD4.bb-bxf&bxg7&bxs` <- apply_symm(`TCR.CD4_bb-bxf&bxg7&bxs`$data, Chao_Jaccard_Matrix)
`jaccard.CD4.ff-bxf&fxs` <- apply_symm(`TCR.CD4_ff-bxf&fxs`$data, Jaccard_Matrix)# Use the apply_symm function to generate just the value matrix
`chao_jaccard.CD4.ff-bxf&fxs` <- apply_symm(`TCR.CD4_ff-bxf&fxs`$data, Chao_Jaccard_Matrix)

`melted_jaccard.CD4.bb-bxf&bxg7&bxs` <- Value_Matrix(`jaccard.CD4.bb-bxf&bxg7&bxs`)
`melted_chao_jaccard.CD4.bb-bxf&bxg7&bxs` <- Value_Matrix(`chao_jaccard.CD4.bb-bxf&bxg7&bxs`)
`melted_jaccard.CD4.ff-bxf&fxs` <- Value_Matrix(`jaccard.CD4.ff-bxf&fxs`)
`melted_chao_jaccard.CD4.ff-bxf&fxs` <- Value_Matrix(`chao_jaccard.CD4.ff-bxf&fxs`)

`Jaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` <- Overlap.Stats(`melted_jaccard.CD4.bb-bxf&bxg7&bxs`, `TCR.CD4_bb-bxf&bxg7&bxs`)
`ChaoJaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` <- Overlap.Stats(`melted_chao_jaccard.CD4.bb-bxf&bxg7&bxs`, `TCR.CD4_bb-bxf&bxg7&bxs`)

`Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs` <- Overlap.Stats(`melted_jaccard.CD4.ff-bxf&fxs`, `TCR.CD4_ff-bxf&fxs`)
`ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs` <- Overlap.Stats(`melted_chao_jaccard.CD4.ff-bxf&fxs`, `TCR.CD4_ff-bxf&fxs`)

# self vs. 'semi self comparisons should be removed as they are not independent.'
# remove comparator and comparatee values which are similaerly named...
remove.nonindependent.samp <- function(OverlapDF){
  correctedDF <- OverlapDF %>%
    mutate(TFmatch = if_else(substr(OverlapDF$Var1, 1, 3) == substr(OverlapDF$Var2, 1, 3), 1, 0)) %>%
    filter(TFmatch == 0) %>%
    select(-TFmatch)
  return(correctedDF)
}
`Jaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` <- remove.nonindependent.samp(`Jaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs`)
`ChaoJaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` <- remove.nonindependent.samp(`ChaoJaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs`)
`Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs` <- remove.nonindependent.samp(`Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs`)
`ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs` <- remove.nonindependent.samp(`ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs`)

# DEFINE NEW COLORS FOR RELATIVE COMPLEMENT PLOTS
bxf_color_RC_bb = "#FFE8CF"
bxg7_color_RC = "#619910"
bxs_color_RC = "#F2036F"

bxf_color_RC_ff = "#21409A"
fxs_color_RC = "#996699"



# generate plots:
F.MHC_groups <- c(bxf, bxg7, bxs)
C.MHC_groups <- c(bxf_color_RC_bb, bxg7_color_RC, bxs_color_RC)
Comparator_group <- c(bxf, bxg7, bxs)
ylimit = 0.3
`Relative_Complement.bb-bxf&bxg7&bxs_CD4_Jac.plot` <- minimally_Defined.Overlap.plots(`Jaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs`, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = 0.6
`Relative_Complement.bb-bxf&bxg7&bxs_CD4_CJ.plot` <- minimally_Defined.Overlap.plots(`ChaoJaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs`, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)
# 
F.MHC_groups <- c(bxf, fxs)
C.MHC_groups <- c(bxf_color_RC_ff, fxs_color_RC)
Comparator_group <- c(bxf, fxs)
ylimit = 0.3
`Relative_Complement.ff-bxf&fxs_CD4_Jac.plot` <- minimally_Defined.Overlap.plots(`Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs`, F.MHC_groups, C.MHC_groups, Comparator_group, ylimit, Jac)
ylimit = 0.6
`Relative_Complement.ff-bxf&fxs_CD4_CJ.plot` <- minimally_Defined.Overlap.plots(`ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs`, F.MHC_groups, C.MHC_groups,Comparator_group, ylimit, Chao)

# Display Plots:
compiled_Relative_Complemention_bb_hets <- `Relative_Complement.bb-bxf&bxg7&bxs_CD4_Jac.plot`+`Relative_Complement.bb-bxf&bxg7&bxs_CD4_CJ.plot`
compiled_Relative_Complemention_ff_hets <- `Relative_Complement.ff-bxf&fxs_CD4_Jac.plot`+`Relative_Complement.ff-bxf&fxs_CD4_CJ.plot`

ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_Relative_Complemention_bb_hets_v2.pdf", compiled_Relative_Complemention_bb_hets, width=245,height=67, units = "mm")
ggsave("Figure Export/Repertoire Overlap/Bar Graphs/CD4.compiled_Relative_Complemention_ff_hets_v2.pdf", compiled_Relative_Complemention_ff_hets, width=227,height=42, units = "mm")

# bb-bxf&bxg7&bxs Relative Complement
F.MHC_groups <- c(bxf, bxg7, bxs)
Comparator_group <- c(bxf, bxg7, bxs)
`ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

`ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.bb-bxf&bxg7&bxs` %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% ANOVA.Tu()

# ff-bxf&fxs Relative Complement
F.MHC_groups <- c(bxf, fxs)
Comparator_group <- c(bxf, fxs)
# `ANOVA.Tu.Relative Complement.ff-bxf&fxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs` %>%
#   filter(Comparator %in% Comparator_group) %>%
#   filter(Comparatee %in% F.MHC_groups) %>%
#   filter(value != 1) %>% ANOVA.Tu()
# 
# `ANOVA.Tu.Relative Complement.ff-bxf&fxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs` %>%
#   filter(Comparator %in% Comparator_group) %>%
#   filter(Comparatee %in% F.MHC_groups) %>%
#   filter(value != 1) %>% ANOVA.Tu()

`Ttest.Relative Complement.ff-bxf&fxs_CD4_Jac` <- `Jaccard_Overlap_Statistics.CD4.ff-bxf&fxs` %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% Simple.T.test()

`Ttest.Relative Complement.ff-bxf&fxs_CD4_CJ` <- `ChaoJaccard_Overlap_Statistics.CD4.ff-bxf&fxs` %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1) %>% Simple.T.test()
# 
# dataset_names <- list('Tu.Jac.bb-bxf&bxg7&bxs' = `ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_Jac`, 'Tu.Jac.ff-bxf&fxs' = `ANOVA.Tu.Relative Complement.ff-bxf&fxs_CD4_Jac`,
#                       'Tu.CJ.bb-bxf&bxg7&bxs' = `ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_CJ`, 'Tu.CJ.ff-bxf&fxs' = `ANOVA.Tu.Relative Complement.ff-bxf&fxs_CD4_CJ`)
# write.xlsx(dataset_names, file = 'Figure Export/Repertoire Overlap/Statistical Tables/Relative Complement_Jaccard.and.VJ_Jaccard Tu.xlsx')
# rm(dataset_names)

dataset_names_v2 <- list('Tu.Jac.bb-bxf&bxg7&bxs' = `ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_Jac`, 'Ttest.Jac.ff-bxf&fxs' = `Ttest.Relative Complement.ff-bxf&fxs_CD4_Jac`,
                      'Tu.CJ.bb-bxf&bxg7&bxs' = `ANOVA.Tu.Relative Complement.bb-bxf&bxg7&bxs_CD4_CJ`, 'Ttest.CJ.ff-bxf&fxs' = `Ttest.Relative Complement.ff-bxf&fxs_CD4_CJ`)
write.xlsx(dataset_names_v2, file = 'Figure Export/Repertoire Overlap/Statistical Tables/Relative Complement_Jaccard.and.VJ_Jaccard Tu_v2.xlsx')
rm(dataset_names_v2)
