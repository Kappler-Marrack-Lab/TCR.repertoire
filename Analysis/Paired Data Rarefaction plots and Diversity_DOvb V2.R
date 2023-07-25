# Author: Alex Brown
# License: MIT
# Name: Paired data TCR rarefaction violin_plots V2
# purpose of this is to test plots which are 2x parental sequences 1x F1 sequences (minus F1 unique sequences) + Teteraparetnal sequences and calculate hills diversity coefficient
# baiscally the same as the original implimentation but changing the source data.
# this includes the 6 runs of bxg7. Two of which have a very large number of clones.
#----------------------------------------------------------------------------------
# Load Libraries
require(devtools)
install_version("iNEXT", version = "2.0.19", repos = "http://cran.us.r-project.org")

library(vegan)
library(iNEXT)
library(immunarch)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(tidyr)
library(parallel)
library(pbmcapply)
library(tictoc)
library(ggsignif)
library(Hmisc)
library(ggrepel)
Max_CPU_Cores = detectCores()
Upper_Limit_CPU_Cores = 2*round((Max_CPU_Cores*0.8)/2) # this sets the program to used 80% of CPU processing capacity. Close all other programs.
set.seed(42)
#----------------------------------------------------------------------------------
# Load in the repertoire data
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed
load(file = "Paired Parental and low frequency sequences removed.RData") # paired parental data

# Combine data from two sources:
TCR.CD4_plus.Paired <- list(
  data = c(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, TCR_DoVb_CD4_paired_AA.ones.removed.LFSR$data),
  meta = rbind(
    dplyr::select(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta, Sample, MHC, `Co-receptor`),
    dplyr::select(TCR_DoVb_CD4_paired_AA.ones.removed.LFSR$meta, Sample, MHC, `Co-receptor`)
  )
)

CD4 = 'CD4'
CD8 = 'CD8'

#----------------------------------------------------------------------------------
# Parallel Rareification Function
# This is a working parallelized function of iNEXT. 2-3x faster than previously
# this actually is plotting the number of unique clones, and is not a diversity number.
parallel_rarefaction <- function(shuffled_data){
  # if(shuffled_data %% 1 == 0) system(paste("echo 'now processing:",shuffled_data,"'"))
  out_df <- iNEXT(as.vector(shuffled_data), q=0, datatype="abundance")
  df <- fortify(out_df, type=1)
  return(df)
}
#----------------Calculate diversity curves. These take ~20 mins to run 75 reperotires -----------------------------------------
# ------- Running this with pbmclapply outside of the loops works as it should! 2-3x faster.
# Putting this first part into a function messes things up for some reason.
TCR_data_list_CD4 <- list()
# TCR_data_list_CD8 <- list()
#Generation of the list for mclapply
immnames <- names(TCR.CD4_plus.Paired$data)
for(i in 1:length(immnames)){
  df0 <- TCR.CD4_plus.Paired$data[[i]]$Clones
  df0 <- df0[! df0 %in% 1] # removes all the single clonotypes
  df_shuf <- sample(df0)# shuffle rows
  tmp <- list(Data=df_shuf)
  name <- immnames[i]
  TCR_data_list_CD4[name] <- tmp
}#end of for loop

list_of_rarefaction_df_CD4 <- pbmclapply(TCR_data_list_CD4, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 20 mins
# list_of_rarefaction_df_CD8 <- pbmclapply(TCR_data_list_CD8, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 10 mins

# CD4
immnames <- names(TCR.CD4_plus.Paired$data)
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
TCR.CD4.iNEXT <- merge(TCR.CD4_plus.Paired$meta,finaldf,by="Sample",all=TRUE)
TCR.CD4.point <- TCR.CD4.iNEXT[which(TCR.CD4.iNEXT$method=="observed"),]
TCR.CD4.line <- TCR.CD4.iNEXT[which(TCR.CD4.iNEXT$method!="observed"),]
TCR.CD4.line$method <- factor(TCR.CD4.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))

# CD8
# immnames <- names(TCR.CD8_DoVb_Paired$data)
# finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
# n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
# colnames(finaldf) <- n
# # Extraction of data and adding metadata information
# for(j in 1:length(immnames)){
#   temp_rare <- as.data.frame(list_of_rarefaction_df_CD8[j])
#   temp_rare$Sample <- immnames[j]
#   colnames(temp_rare) <- n
#   finaldf <- rbind(finaldf, temp_rare)
# }
# TCR.CD8.iNEXT <- merge(TCR.CD8_DoVb_Paired$meta,finaldf,by="Sample",all=TRUE)
# TCR.CD8.point <- TCR.CD8.iNEXT[which(TCR.CD8.iNEXT$method=="observed"),]
# TCR.CD8.line <- TCR.CD8.iNEXT[which(TCR.CD8.iNEXT$method!="observed"),]
# TCR.CD8.line$method <- factor(TCR.CD8.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
#----------------------------------------------------------------------------------

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
tetra_bb.ff = "b+f"
tetra_bb.g7g7 = "b+g7"

# Simulated Pairs
`bb+ff` = "bb+ff"
`ff+ss` = "ff+ss"
`bb+ss` = "bb+ss"
`bb+g7g7` = "bb+g7g7"


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

bb.ff_color = "#7f7fff"
ff.ss_color = "#05a0ff"
bb.ss_color = "#e501fd"
bb.g7g7_color = "#08b129"

# Define the MHC groups and input elements for the current plot
F.MHC_groups <- c(ff, ss, bb, g7g7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color)

#----------------------Plotting Functions------------------------------------------------
# rarefaction violin_plotting by defined group/MHC
# These make just the comparisons we want to make and adds t-test satistical bars
sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}
rarefaction_plot <- function(TCR.iNEXT, TCR.point, co_receptor, F.MHC_groups, C.MHC_groups){
  # # Join columns for the line and dot data
  # TCR.iNEXT[TCR.iNEXT == "N/A"] <- NA
  # TCR.iNEXT <- TCR.iNEXT %>%
  #   unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
  # TCR.iNEXT$MHC <- gsub("_NA.*", "\\1", TCR.iNEXT$MHC)
  # 
  # TCR.point[TCR.point == "N/A"] <- NA
  # TCR.point <- TCR.point %>%
  #   unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
  # TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
  
  temp <- TCR.iNEXT %>% filter(MHC %in% c("bxf", "bxg7", "bxs", "fxs", "bb+ff", "bb+g7g7", "ss+bb", "ss+ff")) %>% filter(method == "observed") #in the updated iNEXT package "Observed" is used in place of "observed"
  x_max <- max(temp$x)*1.05
  y_max <- max(temp$y)*1.05
  
  plot <- ggplot(subset(TCR.iNEXT, MHC %in% F.MHC_groups), aes(x=x, y=y, color=MHC, group=Sample)) + 
    # geom_point(aes(shape=site), size=3, subset(TCR.point, MHC %in% F.MHC_groups), key_glyph = "rect") +
    geom_line(data = subset(TCR.iNEXT, method %in% c("interpolated") & MHC %in% F.MHC_groups), linetype = "solid", lwd=1.0, key_glyph = "rect") + #in the updated iNEXT package "Rarefaction" is used in place of "interpolated"
    # geom_line(data = subset(TCR.iNEXT, method %in% c("extrapolated") & MHC %in% F.MHC_groups), linetype = "dashed", lwd=1.0, key_glyph = "rect") +
    # geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
    #                 fill=MHC, color=NULL), alpha=0.2) +
    
    # identify the max point in dataframe from geom point and make this the x upepr limit.
    scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, 
                      values = C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, 
                       values = C.MHC_groups) +
    
    labs(title="TCRα Species Diversity", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells\nwith haplotypes: ", F.MHC_groups), collapse=" "), x="Number of Sequences Examined", y="Estimated Diversity\n(unique clonotypes)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          axis.text.x = element_text(angle = 45, hjust=1),
          aspect.ratio=1)
  
  return(plot)
}
rare.box_plot <- function(TCR.point, co_receptor, F.MHC_groups, C.MHC_groups){
  MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct comparisons
  filtered.df <- TCR.point %>%
    filter(MHC %in% F.MHC_groups)
  plot <- ggplot(filtered.df, aes(x=MHC, y=y, fill=MHC))+
  # plot <- ggplot(subset(TCR.point, MHC %in% F.MHC_groups), aes(x=MHC, y=y, fill=MHC))+
    
    geom_violin()+
    geom_boxplot(width = 0.2, fill="white")+
    # stat_summary(fun.data=mean_sdl, mult=1, 
    #              geom="pointrange", color="black") +
    geom_jitter(height = 0, width = 0.1, alpha=0.2, size = 3, colour="black") +
    # geom_point(aes(fill = MHC), colour="black", alpha=0.2, size = 5)+
    geom_signif(comparisons = MHC_list,
                test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
    ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
    # geom_label_repel(aes(label = Sample,
    #                      fill = MHC), color = 'white',
    #                  size = 3.5, key_glyph = "rect") +
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, 
                      values=C.MHC_groups) +
    scale_color_manual(breaks = F.MHC_groups, 
                       values=C.MHC_groups) +
    labs(title="TCRα Species Diversity", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="MHC Background", y="Estimated Diversity\n(unique clonotypes)") +
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom", 
          legend.title=element_blank(),
          text=element_text(size=18),
          aspect.ratio=1)
  return(plot)
} #currently has labels removed
#----------------------------------------------------------------------------------
# Define the MHC groups and input elements for the current plot
# ---------- b and f MHC ---------
F.MHC_groups <- c(bxf, `bb+ff`)
C.MHC_groups <- c(bxf_color, bb.ff_color)
TCR.CD4.b_f_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
# ---------- b and g7 MHC ---------
F.MHC_groups <- c(bxg7, `bb+g7g7`)
C.MHC_groups <- c(bxg7_color, bb.g7g7_color)
TCR.CD4.b_g7_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
# ---------- b and s MHC ---------
F.MHC_groups <- c(bxs, `bb+ss`)
C.MHC_groups <- c(bxs_color, bb.ss_color)
TCR.CD4.b_s_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)
# ---------- f and s MHC ---------
F.MHC_groups <- c(fxs, `ff+ss`)
C.MHC_groups <- c(fxs_color, ff.ss_color)
TCR.CD4.f_s_rare_plot <- rarefaction_plot(TCR.CD4.iNEXT, TCR.CD4.point, CD4, F.MHC_groups, C.MHC_groups)

# --------- Plot Curves ----------
# View all the plots
pdf("Figure Export/Rarefaction Curve - Simplified/square rare plot 1.pdf", width = 10, height = 5)
(TCR.CD4.b_f_rare_plot + TCR.CD4.b_g7_rare_plot)

dev.off()

pdf("Figure Export/Rarefaction Curve - Simplified/square rare plot 2.pdf", width = 10, height = 5)
(TCR.CD4.b_s_rare_plot + TCR.CD4.f_s_rare_plot)
dev.off()
# 
# # ----------  CD4 bb vs CD4b+/- # ---------- 
# # Load in the bb and b+/- data from the CD4 total rep dataset.
# load("TCR_fixed_beta_total_rep.RData")
# `rep.CD4_bb_b.+/-` <- repFilter(TCR.CD4_DoVb_TR, .method = "by.meta", .query = list(MHC = include("bb", "b+/-")))
# 
# TCR_data_list_CD4_b.related <- list()
# #Generation of the list for mclapply
# immnames <- names(`rep.CD4_bb_b.+/-`$data)
# for(i in 1:length(immnames)){
#   df0 <- `rep.CD4_bb_b.+/-`$data[[i]]$Clones
#   df0 <- df0[! df0 %in% 1] # removes all the single clonotypes
#   df_shuf <- sample(df0)# shuffle rows
#   tmp <- list(Data=df_shuf)
#   name <- immnames[i]
#   TCR_data_list_CD4_b.related[name] <- tmp
# }#end of for loop
# 
# list_of_rarefaction_df_CD4_b.related <- pbmclapply(TCR_data_list_CD4_b.related, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores) # This is the step that takes up to 20 mins
# # CD4 samples with only bb and b+/-. These samples have not been modified.
# immnames <- names(`rep.CD4_bb_b.+/-`$data)
# finaldf <- data.frame(matrix(ncol = 10, nrow = 0))
# n <- c("datatype", "plottype", "site", "method", "order", "x", "y", "y.lwr", "y.upr", "Sample")
# colnames(finaldf) <- n
# # Extraction of data and adding metadata information
# for(j in 1:length(immnames)){
#   temp_rare <- as.data.frame(list_of_rarefaction_df_CD4_b.related[j])
#   temp_rare$Sample <- immnames[j]
#   colnames(temp_rare) <- n
#   finaldf <- rbind(finaldf, temp_rare)
# }
# TCR.CD4_b.related.iNEXT <- merge(`rep.CD4_bb_b.+/-`$meta,finaldf,by="Sample",all=TRUE)
# TCR.CD4_b.related.point <- TCR.CD4_b.related.iNEXT[which(TCR.CD4_b.related.iNEXT$method=="observed"),] #in updated iNEXT package this is now "Observed" not "observed"
# TCR.CD4_b.related.line <- TCR.CD4_b.related.iNEXT[which(TCR.CD4_b.related.iNEXT$method!="observed"),]
# TCR.CD4_b.related.line$method <- factor(TCR.CD4_b.related.line$method, c("interpolated", "extrapolated"), c("interpolation", "extrapolation"))
# 
# F.MHC_groups <- c(bb, IAbHet)
# C.MHC_groups <- c(bb_color, IAbHet_color)
# `TCR.CD4.b&b+/-_rare_plot` <- rarefaction_plot(TCR.CD4_b.related.iNEXT, TCR.CD4_b.related.point, CD4, F.MHC_groups, C.MHC_groups)
# `TCR.CD4.b&b+/-_boxplot` <- rare.box_plot(TCR.CD4_b.related.point, CD4, F.MHC_groups, C.MHC_groups)
# bb_plot <- (`TCR.CD4.b&b+/-_rare_plot`+`TCR.CD4.b&b+/-_boxplot`)
# pdf("Figure Export/Rarefaction Curve - Simplified/bb plot.pdf", width = 10, height = 5)
# bb_plot
# dev.off()
# 
# pdf("Figure Export/Rarefaction Curve - Simplified/bb plot_v2.pdf", width = 5, height = 5)
# `TCR.CD4.b&b+/-_rare_plot`
# dev.off()

