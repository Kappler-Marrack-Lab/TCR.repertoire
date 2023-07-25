# Author: Alex Brown
# License: MIT
# Name: TCR rarefaction plots
#----------------------------------------------------------------------------------
# Load Libraries
library(vegan)
library(iNEXT)
library(immunarch)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(tictoc)

# Set CPU Core limits. This program runs analysis in parallel to analyse the data.
# can crash your computer if too many resources are allocated.
# This will vary depending on the computer. Using the "Upper_Limit_CPU_Cores" variable should fix any system crashes.
Max_CPU_Cores = detectCores()
Upper_Limit_CPU_Cores = 2*round((Max_CPU_Cores*0.8)/2)
set.seed(42)
#----------------------------------------------------------------------------------
# Load in the repertoire data
load("DOvbeta_repData.RData")
CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'
#----------------------------------------------------------------------------------
# for TCR data you could take a single dataset and take the clone count column and randomize this, turn it into a vector. do you get the same result?
# extract clone count as a column, randmoize, convert to a vector
#----------------------------------------------------------------------------------
# Save Workspace file.
# save.image(file = "Rarefaction_Curve.RData")
# To restore your workspace, type this:
# load("Rarefaction_Curve.RData")


# Save multiple objects. These are the output file which take a long time to make. Save and load these when R is restarted.
# save(TCR.DP.iNEXT, TCR.DP.point, TCR.DP.line, TCR.CD4.iNEXT, TCR.CD4.point, TCR.CD4.line, TCR.CD8.iNEXT, TCR.CD8.point, TCR.CD8.line, file = "Rarefaction_DF.RData")
# To load the data again
# load("Rarefaction_DF.RData")
#----------------------------------------------------------------------------------
# Parallel Rareification Function
# This is a working parallelized function of iNEXT. 2-3x faster than previously
parallel_rarefaction <- function(shuffled_data){
  # if(shuffled_data %% 1 == 0) system(paste("echo 'now processing:",shuffled_data,"'"))
  out_df <- iNEXT(as.vector(shuffled_data), q=0, datatype="abundance")
  df <- fortify(out_df, type=1)
  return(df)
}
#----------------Calculate diversity curves. These take ~20 mins to run 75 reperotires -----------------------------------------
# ------- Running this with pbmclapply outside of the loops works as it should! 2-3x faster.
# Putting this first part into a function messes things up for some reason.

TCR_data_list_DP <- list()
TCR_data_list_CD4 <- list()
TCR_data_list_CD8 <- list()
#Generation of the list for mclapply
immnames <- names(TCR.DP_DoVb_DB$data)
for(i in 1:length(immnames)){
  df0 <- TCR.DP_DoVb_DB$data[[i]]$Clones
  df_shuf <- sample(df0)# shuffle rows
  tmp <- list(Data=df_shuf)
  name <- immnames[i]
  TCR_data_list_DP[name] <- tmp
}#end of for loop

immnames <- names(TCR.CD4_DoVb_DB$data)
for(i in 1:length(immnames)){
  df0 <- TCR.CD4_DoVb_DB$data[[i]]$Clones
  df_shuf <- sample(df0)# shuffle rows
  tmp <- list(Data=df_shuf)
  name <- immnames[i]
  TCR_data_list_CD4[name] <- tmp
}#end of for loop

immnames <- names(TCR.CD8_DoVb_DB$data)
for(i in 1:length(immnames)){
  df0 <- TCR.CD8_DoVb_DB$data[[i]]$Clones
  df_shuf <- sample(df0)# shuffle rows
  tmp <- list(Data=df_shuf)
  name <- immnames[i]
  TCR_data_list_CD8[name] <- tmp
}#end of for loop

list_of_rarefaction_df_DP <- pbmclapply(TCR_data_list_DP, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 10mins (1min on M1 Mac)
list_of_rarefaction_df_CD4 <- pbmclapply(TCR_data_list_CD4, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 10mins (6min on M1 Mac)
list_of_rarefaction_df_CD8 <- pbmclapply(TCR_data_list_CD8, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 10mins (2min on M1 Mac)

# DP
finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
colnames(finaldf) <- n
# Extraction of data and adding metadata information
immnames <- names(TCR.DP_DoVb_DB$data)
for(j in 1:length(immnames)){
  temp_rare <- as.data.frame(list_of_rarefaction_df_DP[j])
  temp_rare$Sample <- immnames[j]
  colnames(temp_rare) <- n
  finaldf <- rbind(finaldf, temp_rare)
}
TCR.DP.iNEXT <- merge(TCR.DP_DoVb_DB$meta,finaldf,by="Sample",all=TRUE)
TCR.DP.point <- TCR.DP.iNEXT[which(TCR.DP.iNEXT$method=="observed"),]
TCR.DP.line <- TCR.DP.iNEXT[which(TCR.DP.iNEXT$method!="observed"),]
TCR.DP.line$method <- factor(TCR.DP.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))

# CD4
immnames <- names(TCR.CD4_DoVb_DB$data)
finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
colnames(finaldf) <- n
# Extraction of data and adding metadata information
for(j in 1:length(immnames)){
  temp_rare <- as.data.frame(list_of_rarefaction_df_CD4[j])
  temp_rare$Sample <- immnames[j]
  colnames(temp_rare) <- n
  finaldf <- rbind(finaldf, temp_rare)
}
TCR.CD4.iNEXT <- merge(TCR.CD4_DoVb_DB$meta,finaldf,by="Sample",all=TRUE)
TCR.CD4.point <- TCR.CD4.iNEXT[which(TCR.CD4.iNEXT$method=="observed"),]
TCR.CD4.line <- TCR.CD4.iNEXT[which(TCR.CD4.iNEXT$method!="observed"),]
TCR.CD4.line$method <- factor(TCR.CD4.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))

# CD8
immnames <- names(TCR.CD8_DoVb_DB$data)
finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
colnames(finaldf) <- n
# Extraction of data and adding metadata information
for(j in 1:length(immnames)){
  temp_rare <- as.data.frame(list_of_rarefaction_df_CD8[j])
  temp_rare$Sample <- immnames[j]
  colnames(temp_rare) <- n
  finaldf <- rbind(finaldf, temp_rare)
}
TCR.CD8.iNEXT <- merge(TCR.CD8_DoVb_DB$meta,finaldf,by="Sample",all=TRUE)
TCR.CD8.point <- TCR.CD8.iNEXT[which(TCR.CD8.iNEXT$method=="observed"),]
TCR.CD8.line <- TCR.CD8.iNEXT[which(TCR.CD8.iNEXT$method!="observed"),]
TCR.CD8.line$method <- factor(TCR.CD8.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))

#----------------------------------------------------------------------------------
# Co-receptor
CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'

# MHC
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"
IAbHet = "b+/-"
b.f_bb = "b+f_bb"
b.f_ff = "b+f_ff"
b.g7_bb = "b+g7_bb"
b.g7_g7g7 = "b+g7_g7g7"

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
b.f_bb_color = "#5c54fb"
b.f_ff_color = "#c5584e"
b.g7_bb_color = "#54fbb5"
b.g7_g7g7_color = "#fb961a"

# Define the MHC groups and input elements for the current plot
F.MHC_groups <- c(ff, ss, bb, g7g7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color)
#----------------------Plotting Functions------------------------------------------------
# Rarefaction plotting by defined group/MHC
rarefaction_plot <- function(TCR.iNEXT, TCR.point, co_receptor, F.MHC_groups, C.MHC_groups){
  # Join columns for the line and dot data
  TCR.iNEXT[TCR.iNEXT == "N/A"] <- NA
  TCR.iNEXT <- TCR.iNEXT %>%
    unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
  TCR.iNEXT$MHC <- gsub("_NA.*", "\\1", TCR.iNEXT$MHC)
  
  TCR.point[TCR.point == "N/A"] <- NA
  TCR.point <- TCR.point %>%
    unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
  TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
  
  plot <- ggplot(subset(TCR.iNEXT, MHC %in% F.MHC_groups), aes(x=x, y=y, color=MHC, group=Sample)) + 
    geom_point(aes(shape=site), size=3, subset(TCR.point, MHC %in% F.MHC_groups), key_glyph = "rect") +
    geom_line(data = subset(TCR.iNEXT, method %in% c("interpolated") & MHC %in% F.MHC_groups), linetype = "solid", lwd=1.0, key_glyph = "rect") +
    geom_line(data = subset(TCR.iNEXT, method %in% c("extrapolated") & MHC %in% F.MHC_groups), linetype = "dashed", lwd=1.0, key_glyph = "rect") +
    geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                    fill=MHC, color=NULL), alpha=0.2) +
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, 
                      values = C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, 
                       values = C.MHC_groups) +
    
    labs(title="TCRα Species Diversity", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="Number of Sequences Examined", y="Estimated Diversity\n(unique clonotypes)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18))
  
  return(plot)
}
rare.box_plot <- function(TCR.point, co_receptor, F.MHC_groups, C.MHC_groups){
  # Join columns for the dot data
  TCR.point[TCR.point == "N/A"] <- NA
  TCR.point <- TCR.point %>%
    unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
  TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
  
  
  plot <- ggplot(subset(TCR.point, MHC %in% F.MHC_groups), aes(x=MHC, y=y, fill=MHC))+
    geom_boxplot()+
    geom_point(aes(fill = MHC), colour="black", alpha=0.2, size = 5)+
    geom_label_repel(aes(label = Sample,
                         fill = MHC), color = 'white',
                     size = 3.5, key_glyph = "rect") +
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, 
                      values=C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, 
                       values=C.MHC_groups) +
    labs(title="TCRα Species Diversity", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="MHC Background", y="Estimated Diversity\n(unique clonotypes)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18))
  return(plot)
} #currently has labels

#----------------------------------------------------------------------------------
# Define the MHC groups and input elements for the current plot
# ---------- b and f MHC ---------
F.MHC_groups <- c(bb, ff, bxf, IAbHet, b.f_bb, b.f_ff)
C.MHC_groups <- c(bb_color, ff_color, bxf_color, IAbHet_color, b.f_bb_color, b.f_ff_color)
TCR.CD4.b_f_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
TCR.CD4.b_f_rare_box_plot <- rare.box_plot(TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb.ff_CD4_rarefaction plot.pdf", width = 20, height = 10)
TCR.CD4.b_f_rare_plot + TCR.CD4.b_f_rare_box_plot
dev.off()

TCR.CD8.b_f_rare_plot <- rarefaction_plot(TCR.CD8.iNEXT, TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
TCR.CD8.b_f_rare_box_plot <- rare.box_plot(TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb.ff_CD8_rarefaction plot.pdf", width = 20, height = 10)
TCR.CD8.b_f_rare_plot + TCR.CD8.b_f_rare_box_plot
dev.off()
# ---------- bb only
bb_CD4 = "bb CD4"
bb_CD8 = "bb CD8"
bb_DP = "bb DP"
bb_color_CD4 = "#54a0fb"
bb_color_CD8 = "#427cbf"
bb_color_DP = "#203d5e"
bb_all = "CD4, CD8 and DP"
F.MHC_groups <- c(bb_CD4, bb_CD8, bb_DP)
C.MHC_groups <- c(bb_color_CD4, bb_color_CD8, bb_color_DP)
bb_only.iNEXT <- bind_rows(TCR.CD4.iNEXT, TCR.CD8.iNEXT, TCR.DP.iNEXT)
bb_only.point <- bind_rows(TCR.CD4.point, TCR.CD8.point, TCR.DP.point)
bb_only.iNEXT <- bb_only.iNEXT %>% filter(MHC == "bb") %>% unite(MHC, c("MHC", "Co-receptor"), sep = " ")
bb_only.point <- bb_only.point %>% filter(MHC == "bb") %>% unite(MHC, c("MHC", "Co-receptor"), sep = " ")
TCR.bb_only_rare_plot <- rarefaction_plot(bb_only.iNEXT, bb_only.point, bb_all, F.MHC_groups, C.MHC_groups)
TCR.bb_only_rare_box_plot <- rare.box_plot(bb_only.point, bb_all, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb only_rarefaction plot.pdf", width = 20, height = 10)
TCR.bb_only_rare_plot + TCR.bb_only_rare_box_plot
dev.off()
# ---------- b and g7 MHC ---------
F.MHC_groups <- c(bb, g7g7, bxg7, b.g7_bb, b.g7_g7g7)
C.MHC_groups <- c(bb_color, g7g7_color, bxg7_color, b.g7_bb_color, b.g7_g7g7_color)
TCR.CD4.b_g7_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
TCR.CD4.b_g7_rare_box_plot <- rare.box_plot(TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb.g7g7_CD4_rarefaction plot.pdf", width = 20, height = 10)
TCR.CD4.b_g7_rare_plot + TCR.CD4.b_g7_rare_box_plot
dev.off()

TCR.CD8.b_g7_rare_plot <- rarefaction_plot(TCR.CD8.iNEXT, TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
TCR.CD8.b_g7_rare_box_plot <- rare.box_plot(TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb.g7g7_CD8_rarefaction plot.pdf", width = 20, height = 10)
TCR.CD8.b_g7_rare_plot + TCR.CD8.b_g7_rare_box_plot
dev.off()
# ---------- CD4 het vs hom plots ---------
F.MHC_groups <- c(ff, ss, bb, g7g7, bxf, fxs, bxs, bxg7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color, bxf_color, fxs_color, bxs_color, bxg7_color)
TCR.CD4.het_VS_hom_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
TCR.CD4.het_VS_hom_rare_box_plot <- rare.box_plot(TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_het_VS_hom_CD4_rarefaction plot.pdf", width = 20, height = 10)
TCR.CD4.het_VS_hom_rare_plot + TCR.CD4.het_VS_hom_rare_box_plot
dev.off()


# ----------  CD4, CD8 bb vs b+/- # ----------
# 2021-03-06 We need to sequence more of these mice and make sure they have the right genotpye. Right now the sequences we have are very unclear.
F.MHC_groups <- c(bb, IAbHet)
C.MHC_groups <- c(bb_color, IAbHet_color)
TCR.CD4.bb_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
TCR.CD8.bb_rare_plot <- rarefaction_plot(TCR.CD8.iNEXT, TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
TCR.CD4.bb_rare_box_plot <- rare.box_plot(TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
TCR.CD8.bb_rare_box_plot <- rare.box_plot(TCR.CD8.point, CD8, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb_and_IAb_CD4_CD8_rarefaction violin_plot.pdf", width = 7, height = 10)
(TCR.CD4.bb_rare_plot + TCR.CD4.bb_rare_box_plot) / (TCR.CD8.bb_rare_plot + TCR.CD8.bb_rare_box_plot)
dev.off()



# ---------- Comparisons of DP curves from all data and singles dropped ---------
# DP cells single hit clones removed.
two.plus_TCR.DP_DoVb_DB <- TCR.DP_DoVb_DB
two.plus_TCR.DP_DoVb_DB$data <- mclapply(two.plus_TCR.DP_DoVb_DB$data, function(x) filter(x, Clones != 1), mc.cores = Upper_Limit_CPU_Cores)
#Generation of the list for mclapply
two.plus_TCR_data_list_DP <- list()
immnames <- names(two.plus_TCR.DP_DoVb_DB$data)
for(i in 1:length(immnames)){
  df0 <- two.plus_TCR.DP_DoVb_DB$data[[i]]$Clones
  df_shuf <- sample(df0)# shuffle rows
  tmp <- list(Data=df_shuf)
  name <- immnames[i]
  two.plus_TCR_data_list_DP[name] <- tmp
}#end of for loop
list_of_rarefaction_df_two.plus_DP <- pbmclapply(two.plus_TCR_data_list_DP, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 10mins
# DP two plus
finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
colnames(finaldf) <- n
# Extraction of data and adding metadata information
immnames <- names(two.plus_TCR.DP_DoVb_DB$data)
for(j in 1:length(immnames)){
  temp_rare <- as.data.frame(list_of_rarefaction_df_two.plus_DP[j])
  temp_rare$Sample <- immnames[j]
  colnames(temp_rare) <- n
  finaldf <- rbind(finaldf, temp_rare)
}
two.plus_TCR.DP.iNEXT <- merge(two.plus_TCR.DP_DoVb_DB$meta,finaldf,by="Sample",all=TRUE)
two.plus_TCR.DP.point <- two.plus_TCR.DP.iNEXT[which(two.plus_TCR.DP.iNEXT$method=="observed"),]
two.plus_TCR.DP.line <- two.plus_TCR.DP.iNEXT[which(two.plus_TCR.DP.iNEXT$method!="observed"),]
two.plus_TCR.DP.line$method <- factor(two.plus_TCR.DP.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))

F.MHC_groups <- c(bb, ff, bxf, IAbHet, b.f_bb, b.f_ff)
C.MHC_groups <- c(bb_color, ff_color, bxf_color, IAbHet_color, b.f_bb_color, b.f_ff_color)
two.plus_TCR.DP.bb_rare_plot <- rarefaction_plot(two.plus_TCR.DP.iNEXT, two.plus_TCR.DP.point, DP, F.MHC_groups, C.MHC_groups)
TCR.DP.bb_rare_plot <- rarefaction_plot(TCR.DP.iNEXT, TCR.DP.point, DP, F.MHC_groups, C.MHC_groups)
pdf("Figure Export/Rarefaction Curve/Va_bb_DP Singles removed vs all data rarefaction plot.pdf", width = 10, height = 5)
two.plus_TCR.DP.bb_rare_plot + TCR.DP.bb_rare_plot
dev.off()

