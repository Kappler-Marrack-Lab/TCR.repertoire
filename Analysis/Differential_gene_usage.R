# Author: Alex Brown
# License: MIT
#-------------------STARTING NOTES-------------------------------------------------
# 2023-01-02
# create some statistical comparisons of gene usage across haplotypes
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
#----------------------------------------------------------------------------------
# Load in the repertoire data
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed


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
DP_color = "#808080"
tetra_bb.ff_color = "#EEC0C4"
tetra_bb.g7g7_color = "#335566"
#-------------------------- Gene Usage --------------------------------------------------------
# also include the DP cells (aggrigate bb, ff, bb DObeta)
# Modify the haplotype of the DP data to be DP in general since we are just looking at V, J frequencies
TCR_DoVb_DP_AA.ones.removed.LFSR$meta$MHC <- "DP"

DP.CD4 <- list(
  data = c(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, TCR_DoVb_DP_AA.ones.removed.LFSR$data),
  meta = rbind(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta, TCR_DoVb_DP_AA.ones.removed.LFSR$meta)
)
# combine the DP and CD4 data into one.
DP.CD4 <- repFilter(DP.CD4, .method = "by.meta", .query = list(MHC = include("bb","ff","g7g7","ss","bxf","bxg7","bxs","fxs", "DP"))) #subest for het and hom reps

# rename Vs to family names
VDJ_lookup <- vroom(file = "MK_Ion_IMGT_Convert.csv", delim = ",")
V.fam.name <- function(df, VDJ_lookup){
  # Convert to IMGT format
  new <- as.data.frame(df)
  new[] <- VDJ_lookup$TRAV_Family[match(unlist(df), VDJ_lookup$MK_Ion_TRAV)]
  new.V.name <- new$V.name
  df <- df %>% mutate(V.name = new.V.name)
  return(df)
}

DP.CD4.family = DP.CD4
DP.CD4.family$data <- lapply(DP.CD4$data, V.fam.name, VDJ_lookup)
TRAV.family_TCRa.geneUsage <- geneUsage(DP.CD4.family$data, .type = c("family"), .norm = F) #quantify as numbers then average based on the TRAV family count. Set normalization to false we can deal with it later. The way it is computed here is not what we want.
TRAV.family_TCRa.geneUsage.plot <- vis(TRAV.family_TCRa.geneUsage, .by = "MHC", .meta = DP.CD4.family$meta, .plot = "box")

pdf("Figure Export/Differential Gene Usage/CD4 V usage by Family.pdf", width = 12.5, height = 4)
TRAV.family_TCRa.geneUsage.plot
dev.off()

# Since the TRAJs to do not have families they are easier to deal with
TRAJ_TCRa.geneUsage <- geneUsage(DP.CD4$data, .gene = c("musmus.TRAJ"), .norm = T) # keep the normalization since there are no families.
vis(TRAJ_TCRa.geneUsage, .by = "MHC", .meta = DP.CD4$meta, .plot = "box")

# ----------------------------
# manual generation of the plot.
TRAV.fam.gu.long <- TRAV.family_TCRa.geneUsage %>%
  gather(key = "Sample", value = "geneUsage.count", -Names, factor_key=TRUE) %>%
  left_join(DP.CD4$meta, by = "Sample") %>%
  dplyr::rename(TRAV_Family = Names) %>%
  select(TRAV_Family, Sample,  MHC, geneUsage.count)

Allele.count.df <- data.frame(TRAV_Family = c("TRAV01", "TRAV02", "TRAV03", "TRAV04", "TRAV05", "TRAV06", "TRAV07", "TRAV08", "TRAV09", "TRAV10", "TRAV11", "TRAV12", "TRAV13", "TRAV14", "TRAV15", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAV21"),
                              Allele_Count = c(1, 1, 5, 7, 3, 15, 13, 4, 6, 3, 3, 9, 13, 8, 3, 4, 1, 1, 1, 1)
)

TRAV.fam.gu.long.Allele.norm <- left_join(TRAV.fam.gu.long, Allele.count.df, by = "TRAV_Family") %>%
  replace(is.na(.), 0) %>%
  mutate(geneUsage.freq.Allele.count = geneUsage.count/Allele_Count) %>%
  group_by(Sample) %>%
  mutate(geneUsage.freq.Allele.norm = geneUsage.freq.Allele.count/sum(geneUsage.freq.Allele.count))
                                
# -----------------------------
TRAJ.gu.long <-TRAJ_TCRa.geneUsage %>%
  gather(key = "Sample", value = "geneUsage.norm", -Names, factor_key=TRUE) %>%
  left_join(DP.CD4$meta, by = "Sample") %>%
  dplyr::rename(TRAJ = Names) %>%
  select(TRAJ, Sample,  MHC, geneUsage.norm)
  

# generate plot
F.MHC_groups <- c(ff, ss, bb, g7g7, bxf, fxs, bxs, bxg7, DP)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color, bxf_color, fxs_color, bxs_color, bxg7_color, DP_color)

ggplot(TRAV.fam.gu.long.Allele.norm, aes(x = TRAV_Family, y = geneUsage.freq.Allele.norm, fill = MHC))+
  geom_bar(stat="identity", position=position_dodge())+
  # set Custom Colors
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  theme_bw()


# line plot
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

TRAV.fam.gu.long.Allele.norm.summary  <- data_summary(TRAV.fam.gu.long.Allele.norm, varname="geneUsage.freq.Allele.norm", 
             groupnames=c("MHC", "TRAV_Family"))

TRAV.gu_plot <- ggplot(TRAV.fam.gu.long.Allele.norm.summary, aes(x = TRAV_Family, y = geneUsage.freq.Allele.norm, group = MHC, color = MHC))+
  geom_errorbar(aes(ymin=geneUsage.freq.Allele.norm-sd, ymax=geneUsage.freq.Allele.norm+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  # scale_y_continuous(limits = c(0,0.08), expand = expansion(mult = c(0, 0.08))) +
  labs(title="TRAV Gene Usage normalized to number of alleles per TRAV family",
       x ="TRAV Family", y = "Normalized Gene Usage")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

 # --------------------------------------------------
 
 TRAJ.gu.long.summary  <- data_summary(TRAJ.gu.long, varname="geneUsage.norm", 
                                                       groupnames=c("MHC", "TRAJ"))
 
TRAJ.gu_plot <- ggplot(TRAJ.gu.long.summary, aes(x = TRAJ, y = geneUsage.norm, group = MHC, color = MHC))+
   geom_errorbar(aes(ymin=geneUsage.norm-sd, ymax=geneUsage.norm+sd), width=.1, 
                 position=position_dodge(0.05)) +
   geom_line() + geom_point()+
   scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
   # scale_y_continuous(limits = c(0,0.08), expand = expansion(mult = c(0, 0.08))) +
   labs(title="TRAJ Gene Usage",
        x ="TRAJ", y = "Normalized Gene Usage")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 
# -----------------------------------------------------------------------------------------------------------------------------
pdf("Figure Export/Differential Gene Usage/CD4 Normalized TRAV gene usage.pdf", width = 5.3, height = 3)
TRAV.gu_plot
dev.off()

pdf("Figure Export/Differential Gene Usage/CD4 Normalized TRAJ gene usage.pdf", width = 10, height = 3)
TRAJ.gu_plot
dev.off()


# -------------------------- Gene usage analysis by Jensen-Shannon Divergence -------------------------------------------------
TRAV.subType_TCRa.geneUsage <- geneUsage(DP.CD4$data, .norm = T) # computed on a per allele level
imm_gu_js <- geneUsageAnalysis(TRAV.subType_TCRa.geneUsage, .method = "js", .verbose = T)
imm_gu_js.PLOT <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
imm_gu_js.PLOT

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  diag(cormat) <- 1
  return(cormat)
}
Value_Matrix <- function(correlation_MX){
  lower_tri_j <- get_lower_tri(correlation_MX)# grab lower triangle of matrix
  melted_cormat <- reshape2::melt(lower_tri_j, na.rm = TRUE) # reshape to dataframe
  # This will restore the whole square and make the scale bars look the same relative to other plots.
  # melted_cormat <- reshape2::melt(correlation_MX, na.rm = FALSE) # uncomment this and comment the top to to ge the square matrix
  melted_cormat[is.na(melted_cormat)] <- 1
  melted_cormat$rnd_value <- scientific(melted_cormat$value, digits = 2)
  return(melted_cormat)
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
  Overlap_Statistics.DF <- Overlap_Statistics.DF %>%
    filter(value<1)
  return(Overlap_Statistics.DF)
}
imm_gu_js.melted <- Value_Matrix(imm_gu_js)
imm_gu_js.melted <- Overlap.Stats(imm_gu_js.melted, DP.CD4) %>% dplyr::filter(Comparator %in% c("bb","ff","g7g7", "ss")) %>% dplyr::filter(Comparatee %in% c("bb","ff","g7g7", "ss"))

# These make just the comparisons we want to make and adds t-test satistical bars
sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

F.MHC_groups <- c(ff, ss, bb, g7g7)
Comparator_group <- c(ff, ss, bb, g7g7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color)

MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct statistical  comparisons
filtered.df <- imm_gu_js.melted %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1)

imm_gu_js.ggplot <- ggplot(imm_gu_js.melted, aes(x=Comparatee, y=value, fill=Comparatee))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  geom_signif(comparisons = MHC_list,
              test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
  ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
  
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  facet_wrap(~Comparator, ncol=1, strip.position = "right") +
  ylim(0.0,0.8)+
  theme_classic()+
  labs(x ="Comparatee", y = "Jensen-Shannon Divergence")+
  theme(
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        strip.text.y.right = element_text(angle = 0, size = 14, face="bold"),
        strip.placement = "outside",
        panel.spacing.y=unit(0.0, "pt"),
        panel.border=element_rect(color="black", fill = NA, size=1)
  )
  
pdf("Figure Export/Differential Gene Usage/CD4 Jensen-Shannon Divergence.pdf", width = 4, height = 7)
imm_gu_js.ggplot
dev.off()

# -------------------------- Gene usage analysis by Coorelation  -------------------------------------------------
imm_gu_cor <- geneUsageAnalysis(TRAV.subType_TCRa.geneUsage, .method = "cor", .verbose = T)
imm_gu_cor.melted <- Value_Matrix(imm_gu_cor)
imm_gu_cor.melted <- Overlap.Stats(imm_gu_cor.melted, DP.CD4) %>% 
  dplyr::filter(Comparator %in% c("bxf","bxg7","bxs", "fxs")) %>%
  dplyr::mutate(merge = paste(Comparator, Comparatee, sep = "_")) %>%
  dplyr::filter(merge %in% c(
    "bxf_bb", "bxf_bxf", "bxf_ff",
    "bxg7_bb", "bxg7_bxg7", "bxg7_g7g7",
    "bxs_bb", "bxs_bxs", "bxs_ss",
    "fxs_ff", "fxs_fxs", "fxs_ss")
    ) %>%
  select(-merge)


F.MHC_groups <- c(bb, ff, g7g7, ss, bxf, bxg7, bxs, fxs)
Comparator_group <- c(bxf, bxg7, bxs, fxs)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color, bxf_color, bxg7_color, bxs_color, fxs_color)

MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct statistical  comparisons
filtered.df <- imm_gu_js.melted %>%
  filter(Comparator %in% Comparator_group) %>%
  filter(Comparatee %in% F.MHC_groups) %>%
  filter(value != 1)

#  ---- bxf corr ------

imm_gu_cor.melted.bxf <- filter(imm_gu_cor.melted, Comparator == "bxf")

cor.bxf.gg <- ggplot(imm_gu_cor.melted.bxf, aes(x=Comparatee, y=value, fill=Comparatee))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  # geom_signif(comparisons = MHC_list,
  #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
  # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
  # 
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  facet_wrap(~Comparator, ncol=1, strip.position = "right") +
  ylim(0.0,1.0)+
  theme_classic()+
  labs(x ="Comparatee", y = "Gene usage correlation")+
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text.y.right = element_text(angle = 0, size = 14, face="bold"),
    strip.placement = "outside",
    panel.spacing.y=unit(0.0, "pt"),
    panel.border=element_rect(color="black", fill = NA, size=1)
  )

#  ---- bxg7 corr ------
imm_gu_cor.melted.bxg7 <- filter(imm_gu_cor.melted, Comparator == "bxg7")

cor.bxg7.gg <- ggplot(imm_gu_cor.melted.bxg7, aes(x=Comparatee, y=value, fill=Comparatee))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  # geom_signif(comparisons = MHC_list,
  #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
  # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
  # 
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  facet_wrap(~Comparator, ncol=1, strip.position = "right") +
  ylim(0.0,1.0)+
  theme_classic()+
  labs(x ="Comparatee", y = "Gene usage correlation")+
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text.y.right = element_text(angle = 0, size = 14, face="bold"),
    strip.placement = "outside",
    panel.spacing.y=unit(0.0, "pt"),
    panel.border=element_rect(color="black", fill = NA, size=1)
  )
  
#  ---- bxs corr ------
imm_gu_cor.melted.bxs <- filter(imm_gu_cor.melted, Comparator == "bxs")

cor.bxs.gg <- ggplot(imm_gu_cor.melted.bxs, aes(x=Comparatee, y=value, fill=Comparatee))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  # geom_signif(comparisons = MHC_list,
  #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
  # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.
  # 
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  facet_wrap(~Comparator, ncol=1, strip.position = "right") +
  ylim(0.0,1.0)+
  theme_classic()+
  labs(x ="Comparatee", y = "Gene usage correlation")+
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text.y.right = element_text(angle = 0, size = 14, face="bold"),
    strip.placement = "outside",
    panel.spacing.y=unit(0.0, "pt"),
    panel.border=element_rect(color="black", fill = NA, size=1)
  )

#  ---- fxs corr ------
imm_gu_cor.melted.fxs <- filter(imm_gu_cor.melted, Comparator == "fxs")

cor.fxs.gg <- ggplot(imm_gu_cor.melted.fxs, aes(x=Comparatee, y=value, fill=Comparatee))+
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+
  # geom_signif(comparisons = MHC_list,
  #             test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
  # ) +  # map_signif_level=TRUE is how I had this set up before to include non-sifnificant comparisons.

  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  facet_wrap(~Comparator, ncol=1, strip.position = "right") +
  ylim(0.0,1.0)+
  theme_classic()+
  labs(x ="Comparatee", y = "Gene usage correlation")+
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    strip.text.y.right = element_text(angle = 0, size = 14, face="bold"),
    strip.placement = "outside",
    panel.spacing.y=unit(0.0, "pt"),
    panel.border=element_rect(color="black", fill = NA, size=1)
  )

pdf("Figure Export/Differential Gene Usage/CD4 Gene usage correlation.pdf", width = 7, height = 4)
(cor.bxf.gg + cor.bxg7.gg) / (cor.bxs.gg + cor.fxs.gg)
dev.off()
