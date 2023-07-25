library(immunarch)
library(tidyverse)

load(file = "Paired Parental and low frequency sequences removed.RData") # paired parental file
load(file = "TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData") # only ones removed
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed

load(file = "counts tables filtered.LFSR.ones.removed.RData")
# primary files:
# TCR_DoVb_CD4_AA.ones.removed
# TCR_DoVb_CD4_AA.ones.removed.LFSR


# TCR_DoVb_CD4_paired_AA.ones.removed.LFSR

# pr_pr_CD4s_filtered.LFSR.ones.removed.LFSR.ones.removed


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

CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'

# MHC homozygous
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"

# MHC heterozygous
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"

# MHC hemizygous
IAbHet = "b+/-"

# Tetraparetnal
tetra_bb.ff = "bb&ff"
tetra_bb.g7g7 = "bb&g7g7"

# Simulated Pairs
`bb+ff` = "bb+ff"
`ff+ss` = "ff+ss"
`bb+ss` = "bb+ss"
`bb+g7g7` = "bb+g7g7"
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Calculate ave number of unique sequences (CDR3AA), and functional reads per mouse (2 functions)
CDR3AA.count <- function(df){
  df.out <- df %>%
    mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    distinct(merge) %>%
    nrow()
  return(df.out)
}

Clone.Calc <- function(df_filtered, df_unfiltered){
  
  # ave number of unique sequences (CDR3AA) per mouse after applying a filter cutoff
  clonotypes_df <- sapply(df_filtered$data, nrow)
  CD4.df <- df_filtered$meta
  CD4.df$V.J.CDR3AA_Filtered <- clonotypes_df
  
  # calculated the number of unique sequeneces (V-CDR3AA-J) per mouse prior to applying a filter cutoff.
  clonotypes.AA_df <- sapply(df_unfiltered$data, CDR3AA.count)
  CD4.df$V.J.CDR3AA_Raw <- clonotypes.AA_df

  # ave number of functional reads per mouse (sum of the clones column)
  Total.seqs_unfiltered <- sapply(df_unfiltered$data, function(i) sum(i$Clones))
  Total.seqs_filtered <- sapply(df_filtered$data, function(i) sum(i$Clones))
  CD4.df$total.seq.unfiltered <- Total.seqs_unfiltered
  CD4.df$total.seq.filtered <- Total.seqs_filtered
  CD4.df$MHC <- factor(CD4.df$MHC)
  return(CD4.df)
}

clones_count_df <- Clone.Calc(TCR_DoVb_CD4_AA.ones.removed.LFSR, TCR_DoVb_CD4_AA.ones.removed) # this has additional columns than the table used later
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Use the counts table to calculate the number of unique sequences when two parents are added together (use the data which has the filtering cutoff applied)
parental_combinations <- function(df, MHC_A, MHC_B){
  # identify parental runs which were a single mouse run at once
  group_A <- df$meta %>% filter(Mouse_Count  == 1, MHC == MHC_A) %>% pull(Sample)
  group_B <- df$meta %>% filter(Mouse_Count  == 1, MHC == MHC_B) %>% pull(Sample)
  Unique_Combination <- expand.grid(group_A, group_B)# Apply expand.grid function
  # names(Unique_Combination) <- c(MHC_A, MHC_B) 
  Unique_Combination <- Unique_Combination %>%
    mutate(MHC = paste(MHC_A, MHC_B, sep = "+")) %>%
    mutate(Sample  = paste(Unique_Combination[[1]], Unique_Combination[[2]], sep = "+"))
  return(Unique_Combination)
}

# Create a master table of all the parental combinations I can make:
PC.df <- rbind(
  parental_combinations(TCR_DoVb_CD4_AA.ones.removed.LFSR, "bb", "ff"),
  parental_combinations(TCR_DoVb_CD4_AA.ones.removed.LFSR, "bb", "g7g7"),
  parental_combinations(TCR_DoVb_CD4_AA.ones.removed.LFSR, "bb", "ss"),
  parental_combinations(TCR_DoVb_CD4_AA.ones.removed.LFSR, "ff", "ss")
)

# Pull out the paired parental vectors
vec_A <- PC.df %>% pull(Var1)
vec_B <- PC.df %>% pull(Var2)

# Need to unfactor these vectors
vec_A <- as.character(vec_A)
vec_B <- as.character(vec_B)

# Set up a for loop to compare the shared observations in the desired two parental MHC runs.
my_vec <- double()
for(i in 1:length(vec_A)){
  my_out <- pr_CD4s_filtered.LFSR.ones.removed %>% select(vec_A[i], vec_B[i]) %>% filter(.[[1]] > 0 | .[[2]] > 0)  %>% nrow()
  my_vec <- c(my_vec, my_out) # my_vec is the returned object
}

# Add this data onto the paired parental data
PC.df$V.J.CDR3AA_Filtered <- my_vec
PC.df$Mouse_Count <- "1+1"

# combine the dataframes and select the columns of interest.
clones_count_df_v2 <- rbind(
  select(clones_count_df, Sample, MHC, Mouse_Count, V.J.CDR3AA_Filtered),
  select(PC.df, Sample, MHC, Mouse_Count, V.J.CDR3AA_Filtered)
) # this contains the normalized unique counts information for all the comparisons we are interested in.


# Calculate totals of each group from the counts table which has the filtering cutoff applied
# Totals need to be calculated with equal numbers of samples per group:
bb_totals <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`264 postJ` > 0 | `267 postJ` > 0 | `270 postJ` > 0)  %>% nrow() # bb: 264, 267, 270 (total of 3 mice)
ff_totals <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`222 postj` > 0 | `225 postJ` > 0 | `228 postJ` > 0)  %>% nrow() # ff: 222, 225, 228 (total of 3 mice)
g7g7_totals <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`493_postJ` > 0 | `495 postj` > 0 | `503 postj` > 0)  %>% nrow() # g7g7: 493, 495, 503 (total of 3 mice)
# ---
bxf_totals <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`388 postJ` > 0 | `394_postJ` > 0 | `421 postJ` > 0)  %>% nrow() # bxf: 388, 394, 421 (total of 3 mice)
bxg7_totals <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`455 postj` > 0 | `631 postJ` > 0 | `643_postJ` > 0)  %>% nrow() # bxg7: 455, 631, 643 (total of 3 mice)
# ---
`bb&ff_totals` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`528_postJ` > 0 | `529_postJ` > 0 | `532_postJ` > 0 | `533_postJ` > 0 | `537 postJ` > 0)  %>% nrow() # bb&ff: 528, 529, 532, 533, 537 (total of 3 mice)
`bb&g7g7_totals` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`570A and B postj` > 0 | `574_postJ` > 0 | `578_postJ` > 0 | `579 postj` > 0)  %>% nrow() # bb&g7g7: 570, 574, 578, 579 (total of 3 mice)
# ---
`bb+ff_totals` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`264 postJ` > 0 | `267 postJ` > 0 | `270 postJ` > 0 | `222 postj` > 0 | `225 postJ` > 0 | `228 postJ` > 0)  %>% nrow() #(3x bb + 3x ff)
`bb+g7g7_totals` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`264 postJ` > 0 | `267 postJ` > 0 | `270 postJ` > 0 | `493_postJ` > 0 | `495 postj` > 0 | `503 postj` > 0)  %>% nrow() #(3x bb + 3x g7g7)

rep_totals <- data.frame(
  MHC = c(bb, ff, g7g7,
          bxf, bxg7,
          "bb&ff", "bb&g7g7",
          "bb+ff", "bb+g7g7"),
  
  V.J.CDR3AA_Filtered_totals = c(bb_totals, ff_totals, g7g7_totals,
                                 bxf_totals, bxg7_totals,
                                 `bb&ff_totals`, `bb&g7g7_totals`,
                                 `bb+ff_totals`, `bb+g7g7_totals`)
)


# Recalculate totals for the bb and b+/- runs only (figure 5):
`bb_totals_x5` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`264 postJ` > 0 | `267 postJ` > 0 | `270 postJ` > 0 | `505_postJ_two_mice` > 0)  %>% nrow() # bb: 264, 267, 270, 505 (total of 5 mice)
`b+/-_totals_x5` <- pr_CD4s_filtered.LFSR.ones.removed %>% filter(`603VJ_postPy` > 0 | `607 postJ` > 0 | `619VJ_postPy` > 0 | `624 postJ` > 0 | `655 postJ` > 0)  %>% nrow() # b+/-: 603, 607, 619, 624, 655 (total of 5 mice)

bb_vs_Iab_totals <- data.frame(
  MHC = c(bb, "b+/-"),
  V.J.CDR3AA_Filtered_totals = c(`bb_totals_x5`, `b+/-_totals_x5`)
)



# ------------ PLOT FUNCTIONS ------------------------
ggAve.counts <- function(df, y.Lim, F.MHC_groups, C.MHC_groups){
  df <- df %>% filter(MHC %in% F.MHC_groups)
  df$MHC <- factor(df$MHC, levels = rev(F.MHC_groups))
  
  ggout <- ggplot(subset(df, MHC %in% F.MHC_groups), aes(x=MHC, y=V.J.CDR3AA_Filtered, fill=MHC))+
    geom_bar(
      stat="summary",
      fun=mean
    ) +
    
    stat_summary(fun.data=mean_se, geom = "linerange")+
    geom_jitter(shape = 1, width = 0.1, size = 0.5)+
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    coord_flip()+
    theme_classic()+
    scale_y_continuous(expand = c(0,0), limits = c(0, y.Lim)) +
    theme(legend.position = "none",
          axis.text.y = element_text(face="bold", size = 12, color = "black"),
          axis.text.x = element_text(color = "black"),
    ) +
    labs(x="MHC")
  return(ggout)
}
# plot total V-CDR3AA-J counts per MHC
ggAve.totals <- function(df, y.Lim, F.MHC_groups, C.MHC_groups){
  df <- df %>% filter(MHC %in% F.MHC_groups)
  df$MHC <- factor(df$MHC, levels = rev(F.MHC_groups))
  ggout <- ggplot(subset(df, MHC %in% F.MHC_groups), aes(x=MHC, y=V.J.CDR3AA_Filtered_totals, fill=MHC))+
    geom_col() +
    # geom_text(position=position_dodge(width=1), hjust=-1, fontface="italic", color="grey")+
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    coord_flip()+
    theme_classic()+
    scale_y_continuous(expand = c(0,0), limits = c(0, y.Lim)) +
    theme(legend.position = "none",
          axis.text.y = element_text(face="bold", size = 12, color = "black"),
          axis.text.x = element_text(color = "black"),
    ) +
    
    labs(x="MHC")
  return(ggout)
}

ggAve.reads <- function(df, y.Lim, F.MHC_groups, C.MHC_groups){
  df <- df %>% filter(MHC %in% F.MHC_groups)
  df$MHC <- factor(df$MHC, levels = rev(F.MHC_groups))
  
  ggout <- ggplot(subset(df, MHC %in% F.MHC_groups), aes(x=MHC, y=total.seq.filtered, fill=MHC))+
    geom_bar(
      stat="summary",
      fun=mean
    ) +
    
    stat_summary(fun.data=mean_se, geom = "linerange")+
    geom_jitter(shape = 1, width = 0.1, size = 0.5)+
    # set Custom Colors
    scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
    coord_flip()+
    theme_classic()+
    scale_y_continuous(expand = c(0,0), limits = c(0, y.Lim)) +
    theme(legend.position = "none",
          axis.text.y = element_text(face="bold", size = 12, color = "black"),
          axis.text.x = element_text(color = "black"),
    ) +
    labs(x="MHC")
  return(ggout)
}


# Define a consistent Y- axis for the plots
clonotype_count.yLim_V.J.CDR3AA <- max(clones_count_df_v2$V.J.CDR3AA_Filtered)*1.05
bb_IAb_yLim <- max(bb_vs_Iab_totals$V.J.CDR3AA_Filtered_totals)*1.05 
tetra_totals_yLim <- max(rep_totals$V.J.CDR3AA_Filtered_totals)*1.05 
# -----------------------------------------------------
# Average number of total functional sequences:
# ggAve.reads(clones_count_df, max(clones_count_df$total.seq.filtered)*1.05, F.MHC_groups, C.MHC_groups)

# -----------------------------------------------------
# Plot groups: bf, bg7, sb, sf
# Group bf: bb, bxf, ff, bb+ff
F.MHC_groups <- c(bb, bxf, ff, `bb+ff`)
C.MHC_groups <- c(bb_color, bxf_color, ff_color, bb.ff_color)
bf_p1 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")

F.MHC_groups <- c(bb, bxg7, g7g7, `bb+g7g7`)
C.MHC_groups <- c(bb_color, bxg7_color, g7g7_color, bb.g7g7_color)
bg7_p1 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")

F.MHC_groups <- c(bb, bxs, ss, `bb+ss`)
C.MHC_groups <- c(bb_color, bxs_color, ss_color, bb.ss_color)
bs_p1 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")

F.MHC_groups <- c(ff, fxs, ss, `ff+ss`)
C.MHC_groups <- c(ff_color, fxs_color, ss_color, ff.ss_color)
fs_p1 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")

averave_number_of_unique_seq_corrected <- bf_p1/bg7_p1/bs_p1/fs_p1
ggsave("Figure Export/Clonotype and Seq average 2/ave unique sequences Corrected Fig3.pdf", plot = averave_number_of_unique_seq_corrected,  width = 8.6, height = 15, units = "cm")
# -----------------------------------------------------
# bb vs b+/- samples
F.MHC_groups <- c(bb, IAbHet)
C.MHC_groups <- c(bb_color, IAbHet_color)
b.Iab_p1 <- ggAve.counts(clones_count_df_v2, 16000, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")
b.Iab_totals <- ggAve.totals(bb_vs_Iab_totals, bb_IAb_yLim, F.MHC_groups, C.MHC_groups) + labs(y="Total number of\n unique sequences\n(V-CDR3-J)")
ggsave("Figure Export/Clonotype and Seq average 2/ave.bb and IAb unique sequences Corrected.pdf", plot = b.Iab_p1,  width = 8.0, height = 3.2, units = "cm")
ggsave("Figure Export/Clonotype and Seq average 2/total.bb and IAb unique sequences Corrected.pdf", plot = b.Iab_totals,  width = 8.0, height = 3.2, units = "cm")

# -----------------------------------------------------
# Plot groups: bf, bg7, sb, sf
# Group bf: bb, bxf, ff, bb+ff
F.MHC_groups <- c(bb, bxf, ff, tetra_bb.ff, `bb+ff`)
C.MHC_groups <- c(bb_color, bxf_color, ff_color, tetra_bb.ff_color, bb.ff_color)
bf_p2 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")
bf_totals <- ggAve.totals(rep_totals, tetra_totals_yLim, F.MHC_groups, C.MHC_groups) + labs(y="Total number of\n unique sequences\n(V-CDR3-J)")

F.MHC_groups <- c(bb, bxg7, g7g7, tetra_bb.g7g7, `bb+g7g7`)
C.MHC_groups <- c(bb_color, bxg7_color, g7g7_color, tetra_bb.g7g7_color, bb.g7g7_color)
bg7_p2 <- ggAve.counts(clones_count_df_v2, clonotype_count.yLim_V.J.CDR3AA, F.MHC_groups, C.MHC_groups) + labs(y="Average number of\n unique sequences\n(V-CDR3-J)")
bg7_totals <- ggAve.totals(rep_totals, tetra_totals_yLim, F.MHC_groups, C.MHC_groups) + labs(y="Total number of\n unique sequences\n(V-CDR3-J)")


averave_number_of_unique_seq_corrected_tetra <- bf_p2/bg7_p2
total_number_of_unique_seq_corrected_tetra <- bf_totals/bg7_totals
ggsave("Figure Export/Clonotype and Seq average 2/ave.tetra unique sequences Corrected Fig5.pdf", plot = averave_number_of_unique_seq_corrected_tetra,  width = 10, height = 10, units = "cm")
ggsave("Figure Export/Clonotype and Seq average 2/total.tetra unique sequences Corrected Fig5.pdf", plot = total_number_of_unique_seq_corrected_tetra,  width = 10, height = 10, units = "cm")


# --------------- Generate Tables --------------------
require(rstatix)
bf_UniqueSeq.A <- clones_count_df_v2 %>% filter(MHC %in% c(bb, bxf, ff, tetra_bb.ff, `bb+ff`)) %>% tukey_hsd(V.J.CDR3AA_Filtered ~ MHC)
bg7_UniqueSeq.A <- clones_count_df_v2 %>% filter(MHC %in% c(bb, bxg7, g7g7, tetra_bb.g7g7, `bb+g7g7`)) %>% tukey_hsd(V.J.CDR3AA_Filtered ~ MHC)
sb_UniqueSeq.A <- clones_count_df_v2 %>% filter(MHC %in% c(bb, bxs, ss, `bb+ss`)) %>% tukey_hsd(V.J.CDR3AA_Filtered ~ MHC)
sf_UniqueSeq.A <- clones_count_df_v2 %>% filter(MHC %in% c(ff, fxs, ss, `ff+ss`)) %>% tukey_hsd(V.J.CDR3AA_Filtered ~ MHC)

bb.Iab <- clones_count_df_v2 %>% filter(MHC %in% c(bb, IAbHet)) %>% t_test(V.J.CDR3AA_Filtered ~ MHC) %>% add_significance("p")

l <- list(bf_UniqueSeq.A = bf_UniqueSeq.A, bg7_UniqueSeq.A = bg7_UniqueSeq.A, sb_UniqueSeq.A = sb_UniqueSeq.A, sf_UniqueSeq.A = sf_UniqueSeq.A, bb.Iab = bb.Iab)
openxlsx::write.xlsx(l, file = "Figure Export/Clonotype and Seq average 2/Unique sequences by group filtering cutoff applied.xlsx")




