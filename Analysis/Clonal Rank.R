# Clonal Rank -------------
# Purpose of this is to see if large clones appear because they are easy to generate or because they are important/selected by MHC
# this is done by comparing public high frequency clones to DP thymocytes.
#  -------------
# Detach previous packages in use:
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))

# Load Packages and requisite data
library(immunarch)
require(tidyr)
library(tidyverse)
library(gghighlight)
library(ggplot2)
library(ggrepel)
library(vroom)
library(RColorBrewer)
library(paletteer)
theme_set(theme_classic())
set.seed(42)
# ----------------------------------------
# Save and load File
# save.image(file = "Clonal Rank.RData")
# load("Clonal Rank.RData")

load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed
# -----------
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"
# -----------
ff_color ="#c5944e"
ss_color = "#fbff00"
bb_color = "#54a0fb" 
g7g7_color = "#9f104c"
# -----------


# Functions
DP_Rank_Plot_list.list <- function(ArchData, DP_list){
  # Identify top clonotypes in query dataset
  ArchData_pr <- pubRep(ArchData$data, "aa+v+j", .coding = T, .verbose = F)
  ArchData_pr <- as_tibble(ArchData_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
    mutate(rank = row_number()) %>%  select(rank, merge, sum)
  
  # Generate a proportion rank of the DP thymocytes
  DP_list_pr <- pubRep(DP_list$data, "aa+v+j", .quant = c("prop"), .coding = T, .verbose = F)
  DP_total <- as_tibble(DP_list_pr) %>%
    dplyr::mutate(prop.sum = rowSums(across(DP_list$meta %>% pull(Sample)), na.rm = TRUE)) %>%   # Sums the total clone count of the haplotypes in a public repertoire table
    dplyr::mutate(Proportion = prop.sum/Samples) %>%  arrange(desc(Proportion)) %>%# calculate the normalized proportion
    dplyr::mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
    select(merge, Proportion)
  
  # Merge Dataframes
  DP_sub <- DP_total %>%
    group_by(merge) %>%
    summarise(Proportion = sum(Proportion)) %>% # this removes and summarizes clonotpyes which show up more than once
    left_join(ArchData_pr, by = c("merge")) %>% #combine dataframes
    arrange(desc(Proportion)) %>%
    mutate(Match = ifelse(!is.na(rank), paste(merge, "\n", "Clonal Rank:", rank), NA)) #rename data for labeling
  
  Top100.nm <-deparse(substitute(ArchData))
  DP.nm <-deparse(substitute(DP_list)) #for printing the name of the DP thymocyte dataframe in title
  
  
  
  
  Rank_plot <- DP_sub %>%
    arrange(desc(Proportion)) %>%
    slice(1:2500) %>%
    ggplot(., aes(x = reorder(merge, -Proportion), y = Proportion, fill = Match)) +
    geom_bar(stat = "identity", fill = "#FF0000") +
    geom_text_repel(aes(label=Match),
                    # ylim  = 0.002,
                    force_pull   = 0, # do not pull toward data points
                    nudge_y      = 0.005,
                    direction    = "x",
                    angle        = 90,
                    hjust        = 0,
                    color = "red",
                    segment.color = "red",
                    size = 3,
                    # Add extra padding around each text label.
                    box.padding = unit(0.5, 'lines'),
                    # Add extra padding around each data point.
                    point.padding = unit(1.6, 'lines'),
                    max.iter = 1e4, max.time = 1
    ) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(expand=c(0,0.0005)) +
    scale_x_discrete(expand = c(0.01, 0)) +
    labs(title = paste("Top 100 shared clonotypes in", Top100.nm, "and their frequency in", DP.nm), x = " Top 2,500 DP Thymocyte Clonotypes (CDR3aa TRAV TRAJ)") +
    gghighlight(Match > 0, use_direct_label = FALSE, label_key = Match)
  
  return(Rank_plot)
} # input: input MHC subset list, DP thymocyte cells by name
DP_Rank_Plot_list.df <- function(ArchData, DP_df){
  # Identify top clonotypes in query dataset
  ArchData_pr <- pubRep(ArchData$data, "aa+v+j", .coding = T, .verbose = F)
  ArchData_pr <- as_tibble(ArchData_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
    mutate(rank = row_number()) %>%  select(rank, merge, sum)
  
  # Merge Dataframes
  DP_total <- DP_df%>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>% select(merge, Proportion)
  DP_sub <- DP_total %>%
    group_by(merge) %>%
    summarise(Proportion = sum(Proportion)) %>% # this removes and summarizes clonotpyes which show up more than once
    left_join(ArchData_pr, by = c("merge")) %>% #combine dataframes
    arrange(desc(Proportion)) %>%
    mutate(Match = ifelse(!is.na(rank), paste(merge, "\n", "Clonal Rank:", rank), NA)) #rename data for labeling

  Top100.nm <-deparse(substitute(ArchData))
  DP.nm <-deparse(substitute(DP_df)) #for printing the name of the DP thymocyte dataframe in title
  
  Rank_plot <- DP_sub %>%
    arrange(desc(Proportion)) %>%
    slice(1:5000) %>%
    ggplot(., aes(x = reorder(merge, -Proportion), y = Proportion, fill = Match)) +
    geom_bar(stat = "identity", fill = "#FF0000") +
    geom_text_repel(aes(label=Match),
                    # ylim  = 0.002,
                    force_pull   = 0, # do not pull toward data points
                    nudge_y      = 0.005,
                    direction    = "x",
                    angle        = 90,
                    hjust        = 0,
                    color = "red",
                    segment.color = "red",
                    size = 3,
                    # Add extra padding around each text label.
                    box.padding = unit(0.5, 'lines'),
                    # Add extra padding around each data point.
                    point.padding = unit(1.6, 'lines'),
                    max.iter = 1e4, max.time = 1
    ) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    scale_y_continuous(expand=c(0,0)) +
    labs(title = paste("Top 100 shared clonotypes in", Top100.nm, "and their frequency in", DP.nm), x = " Top 5,000 DP Thymocyte Clonotypes (CDR3aa TRAV TRAJ)") +
    gghighlight(Match > 0, use_direct_label = FALSE, label_key = Match)

  return(Rank_plot)
} # input: input MHC subset list, DP thymocyte cells by name
# ----------
CD4.bb <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb")))
CD8.bb <- repFilter(TCR_DoVb_CD8_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb")))
CD4.ff <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("ff")))
CD8.ff <- repFilter(TCR_DoVb_CD8_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("ff")))
CD4.g7g7 <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("g7g7")))
CD8.g7g7 <- repFilter(TCR_DoVb_CD8_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("g7g7")))
CD4.ss <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("ss")))
DP.bb <- repFilter(TCR_DoVb_DP_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb")))
DP.ff <- repFilter(TCR_DoVb_DP_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("ff")))
# ------------

# # Plot of top 100 clones which appear in the double positive repertoire
# bb_DP_378 <- TCR_DoVb_DP_AA.ones.removed.LFSR$data$`378 postJ`
# ff_DP_387 <- TCR_DoVb_DP_AA.ones.removed.LFSR$data$`387 postJ`
# 
# CD4.bb_rank_vsSingleDP.df <- DP_Rank_Plot_list.df(CD4.bb, bb_DP_378) # expects arch data, followed by DP thymocyte dataframe
# CD8.bb_rank_vsSingleDP.df <- DP_Rank_Plot_list.df(CD8.bb, bb_DP_378) # expects arch data, followed by DP thymocyte dataframe
# CD4.ff_rank_vsSingleDP.df <- DP_Rank_Plot_list.df(CD4.ff, ff_DP_387) # expects arch data, followed by DP thymocyte dataframe
# CD8.ff_rank_vsSingleDP.df <- DP_Rank_Plot_list.df(CD8.ff, ff_DP_387) # expects arch data, followed by DP thymocyte dataframe
# CD8.bb_ranked_with_ff.DP_vsSingleDP.df <- DP_Rank_Plot_list.df(CD8.bb, ff_DP_387) # expects arch data, followed by DP thymocyte dataframe
# CD8.ff_ranked_with_bb.DP_vsSingleDP.df <- DP_Rank_Plot_list.df(CD8.ff, bb_DP_378) # expects arch data, followed by DP thymocyte dataframe
# 
# # # Compare CD4 data against selected DP thymocytes
# all.DP.vs.CD4.bb <- DP_Rank_Plot_list.list(CD4.bb, TCR_DoVb_DP_AA.ones.removed.LFSR)
# all.DP.vs.CD4.ff <- DP_Rank_Plot_list.list(CD4.ff, TCR_DoVb_DP_AA.ones.removed.LFSR)
# all.DP.vs.CD4.g7g7 <- DP_Rank_Plot_list.list(CD4.g7g7, TCR_DoVb_DP_AA.ones.removed.LFSR)
# all.DP.vs.CD4.ss <- DP_Rank_Plot_list.list(CD4.ss, TCR_DoVb_DP_AA.ones.removed.LFSR)
# 
# # Compare CD4 data against all preselection thymocytes
# ps.vs.CD4.bb <- DP_Rank_Plot_list.list(CD4.bb, TCR_DoVb_DP_AA.ones.removed.LFSR)
# ps.vs.CD4.ff <- DP_Rank_Plot_list.list(CD4.ff, TCR_DoVb_DP_AA.ones.removed.LFSR)
# ps.vs.CD4.g7g7 <- DP_Rank_Plot_list.list(CD4.g7g7, TCR_DoVb_DP_AA.ones.removed.LFSR)
# ps.vs.CD4.ss <- DP_Rank_Plot_list.list(CD4.ss, TCR_DoVb_DP_AA.ones.removed.LFSR)
# 
# pdf("Figure Export/DP Thymocyte Rank/Rank of top 2.5K clones in pooled DP thymocyte populations.pdf", width = 8, height = 8)
# (all.DP.vs.CD4.bb + all.DP.vs.CD4.ff)/
# (all.DP.vs.CD4.g7g7 + all.DP.vs.CD4.ss)
# dev.off()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Map all top 100 shared clones in CD4 to a single preselection plot:
# Identify top clonotypes in query dataset
# bb pr
CD4.bb_pr <- pubRep(CD4.bb$data, "aa+v+j", .coding = T, .verbose = F)
CD4.bb_pr <- as_tibble(CD4.bb_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
  mutate(rank = row_number()) %>% mutate(MHC = "bb") %>%  select(rank, merge, sum, MHC)

# ff pr
CD4.ff_pr <- pubRep(CD4.ff$data, "aa+v+j", .coding = T, .verbose = F)
CD4.ff_pr <- as_tibble(CD4.ff_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
  mutate(rank = row_number()) %>% mutate(MHC = "ff") %>%  select(rank, merge, sum, MHC)

# g7g7 pr
CD4.g7g7_pr <- pubRep(CD4.g7g7$data, "aa+v+j", .coding = T, .verbose = F)
CD4.g7g7_pr <- as_tibble(CD4.g7g7_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
  mutate(rank = row_number()) %>% mutate(MHC = "g7g7") %>%  select(rank, merge, sum, MHC)

# ss pr
CD4.ss_pr <- pubRep(CD4.ss$data, "aa+v+j", .coding = T, .verbose = F)
CD4.ss_pr <- as_tibble(CD4.ss_pr) %>% mutate(sum = rowSums(across(where(is.numeric)))) %>% arrange(desc(sum)) %>% slice_max(sum, n = 100) %>% mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
  mutate(rank = row_number()) %>% mutate(MHC = "ss") %>%  select(rank, merge, sum, MHC)

ArchData_pr <- rbind(CD4.bb_pr, CD4.ff_pr, CD4.g7g7_pr, CD4.ss_pr)

# Generate a proportion rank of the DP thymocytes
TCR_DoVb_DP_AA.ones.removed.LFSR_pr <- pubRep(TCR_DoVb_DP_AA.ones.removed.LFSR$data, "aa+v+j", .quant = c("prop"), .coding = T, .verbose = F)
DP_total <- as_tibble(TCR_DoVb_DP_AA.ones.removed.LFSR_pr) %>%
  dplyr::mutate(prop.sum = rowSums(across(TCR_DoVb_DP_AA.ones.removed.LFSR$meta %>% pull(Sample)), na.rm = TRUE)) %>%   # Sums the total clone count of the haplotypes in a public repertoire table
  dplyr::mutate(Proportion = prop.sum/Samples) %>%  arrange(desc(Proportion)) %>%# calculate the normalized proportion
  dplyr::mutate(merge = paste(V.name, CDR3.aa, J.name)) %>%
  select(merge, Proportion)

# Merge Dataframes
DP_sub <- DP_total %>%
  group_by(merge) %>%
  summarise(Proportion = sum(Proportion)) %>% # this removes and summarizes clonotpyes which show up more than once
  left_join(ArchData_pr, by = c("merge")) %>% #combine dataframes
  arrange(desc(Proportion)) %>%
  mutate(Match = ifelse(!is.na(rank), paste(merge, "\n", "CD4", MHC, "Clonal Rank:", rank), NA)) #rename data for labeling


# scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
# scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  
F.MHC_groups <- c(bb, ff, g7g7, ss)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color)
  
Rank_plot <- DP_sub %>%
  arrange(desc(Proportion)) %>%
  slice(1:3000) %>%
  ggplot(., aes(x = reorder(merge, -Proportion), y = Proportion, fill = MHC)) +
  geom_bar(stat = "identity") +
  geom_text_repel(aes(label=Match, color=MHC),
                  # ylim  = 0.002,
                  force_pull   = 0, # do not pull toward data points
                  nudge_y      = 0.005,
                  direction    = "x",
                  angle        = 90,
                  hjust        = 0,
                  # color = "red",
                  # segment.color = "red",
                  size = 3,
                  # Add extra padding around each text label.
                  box.padding = unit(0.5, 'lines'),
                  # Add extra padding around each data point.
                  point.padding = unit(1.6, 'lines'),
                  max.iter = 1e4, max.time = 1
  ) +
  
  # scale_fill_discrete(
  #   name = "MHC",
  #   # The same color scale will apply to both of these aesthetics.
  #   aesthetics = c("fill", "segment.color")
  # )+
  # 
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(expand=c(0,0.0005)) +
  scale_x_discrete(expand = c(0.01, 0)) +
  labs(title = paste("Top 100 shared clonotypes per MHC-II haplotype and their frequency in the preselection DP thymocyte repertoire"), x = "3,000 most frequent preselection DP thymocyte sequences\n(TRAV-CDR3α-TRAJ)")

  # gghighlight(Match > 0, use_direct_label = FALSE, label_key = Match)


pdf("Figure Export/DP Thymocyte Rank/Rank of top 3K clones in the preselection DP thymocyte repertoire.pdf", width = 10, height = 6)
Rank_plot
dev.off()