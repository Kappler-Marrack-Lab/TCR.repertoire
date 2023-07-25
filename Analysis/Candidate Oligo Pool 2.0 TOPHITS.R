library(immunarch)
library(ggrepel)
library(vroom)
library(ggplot2)
library(svglite)
require(tidyr)
library(viridis)
library(RColorBrewer)
library(scales)
library(paletteer)
set.seed(42)

# Shared Clonotypes:
# interpreting these plots:
#  Size =  freq <- log10(sqrt(as.numeric(df[, 1]) * df[, 2])) / 2
#  Color = pnt.cols <- log(df[, 1] / df[, 2])
# Visualise correlation of public clonotype frequencies in pairs of repertoires.
# X and Y axis are frequency counts of V subfamily and CDR3 AA public clonotypes, with log10 transformation.
# high frequency in X rep is = right hand
# high frequency in Y rep is top
# way to interpret these figures:
# Axis: higher on one side means the public clone is high frequency in that reperotire.
# Size: size is a relative scale of how big each clone is in both reperotires. Also tracks with a white color.
# color is how much more a clone from one reperotire is reperesented in the sharing vs. another clone.
# R^2 is the gemetric mean of the plot. high R2 should correlate to equal contribution between two reptoerires.
# an x or y shifted plot means that the one reperotire is contributing more to the public sharing than the other.
# High presence of values on the direct center, especially large values suggests that these public clones are important for both reperotires.
# May still be useful to document in the top left corner what the frequency of sharing is in the two compared reperotires, and total Vsub gene and CDR3AA shared clones are.
# can also do this with V, J and CDR3 mates.


# -----------------------------------------------------
# To load the modified public clonotypes function
load("vis_public_clonotypes_labels.RData")
load("Canidate Oligo Pool 2.0 TOPHITS.RData")
load("position_bunch.RData")

# save.image(file = "Canidate Oligo Pool 2.0 TOPHITS.RData")
# -----------------------------------------------------

load(file = "TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData") # filtered and ones removed
load(file = "counts tables filtered.ones.removed.RData") # counts table with only the ones removed

# -----
# Import anticipated CDR3 sequences
MK_convert <- vroom("MK_Ion_IMGT_Convert.csv", delim = ",")
CDR.add <- function(VDJ_lookup, df){
  # Convert to IMGT format
  new <- as.data.frame(df)
  new[] <- VDJ_lookup$CDR1[match(unlist(df), VDJ_lookup$MK_Ion_TRAV)]
  new.CDR1 <- new$V.name
  new[] <- VDJ_lookup$CDR2[match(unlist(df), VDJ_lookup$MK_Ion_TRAV)]
  new.CDR2 <- new$V.name
  df <- df %>% mutate(CDR1 = new.CDR1) %>% mutate(CDR2 = new.CDR2)
  return(df)
}
# ----

pr_CD4_merge_prop <- pubRep(TCR_DoVb_CD4_AA.ones.removed$data, "aa+v+j", .quant = c("prop"), .verbose = T)

pr_CD4_merge <- pubRep(TCR_DoVb_CD4_AA.ones.removed$data, "aa+v+j", .verbose = T)

pr_full_CD4_sample.cnt <- as_tibble(pr_CD4_merge) %>%
  # Sums the total clone count of the haplotypes in a public reperotire table
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(ff.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(g7g7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(bxf.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(bxg7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(`b+f.sum` = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(`b+g7.sum` = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample)), na.rm = TRUE)) %>%
  
  # Sums the number of mice (samples) which a given clone appears
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample))))) %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample))))) %>%
  dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample))))) %>%
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample))))) %>%
  dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample))))) %>%
  dplyr::mutate(bf_tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample))))) %>%
  dplyr::mutate(bg7_tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample))))) %>%
  
  # Filter the relevant columns
  select(CDR3.aa, V.name, J.name,
         bb.sum, bb.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample),
         ff.sum, ff.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample),
         g7g7.sum, g7g7.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample),
         bxf.sum, bxf.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample),
         bxg7.sum, bxg7.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample),
         `b+f.sum`, bf_tetra.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample),
         `b+g7.sum`, bg7_tetra.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample)
         )


pr_full_CD4_sample.prop <- as_tibble(pr_CD4_merge_prop) %>%
  # Sums the total clone count of the haplotypes in a public reperotire table
  dplyr::mutate(bb = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(ff = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(g7g7 = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(bxf = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(bxg7 = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(`b+f` = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample)), na.rm = TRUE)) %>%
  dplyr::mutate(`b+g7` = rowSums(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample)), na.rm = TRUE)) %>%
  
  # Sums the number of mice (samples) which a given clone appears
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample))))) %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample))))) %>%
  dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample))))) %>%
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample))))) %>%
  dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample))))) %>%
  dplyr::mutate(bf_tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample))))) %>%
  dplyr::mutate(bg7_tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample))))) %>%
  
  # calculate the average proportion
  dplyr::mutate(bb.ave.prop = bb/bb.Samples) %>%
  dplyr::mutate(ff.ave.prop = ff/ff.Samples) %>%
  dplyr::mutate(g7g7.ave.prop = g7g7/g7g7.Samples) %>%
  dplyr::mutate(bxf.ave.prop = bxf/bxf.Samples) %>%
  dplyr::mutate(bxg7.ave.prop = bxg7/bxg7.Samples) %>%
  dplyr::mutate(`b+f.ave.prop` = `b+f`/bf_tetra.Samples) %>%
  dplyr::mutate(`b+g7.ave.prop` = `b+g7`/bg7_tetra.Samples) %>%
  
  # Filter the relevant columns
  select(CDR3.aa, V.name, J.name,
         bb.ave.prop, bb.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bb') %>% pull(Sample),
         ff.ave.prop, ff.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'ff') %>% pull(Sample),
         g7g7.ave.prop, g7g7.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'g7g7') %>% pull(Sample),
         bxf.ave.prop, bxf.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxf') %>% pull(Sample),
         bxg7.ave.prop, bxg7.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'bxg7') %>% pull(Sample),
         `b+f.ave.prop`, bf_tetra.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+f') %>% pull(Sample),
         `b+g7.ave.prop`, bg7_tetra.Samples, TCR_DoVb_CD4_AA.ones.removed$meta %>% filter(MHC == 'b+g7') %>% pull(Sample)
  )

# This has the samples stored based on average frequency values but also pulls in the sum of counts form another dataframe
pr_full_CD4_mnml <- pr_full_CD4_sample.prop %>%
  mutate(bb.sum = pr_full_CD4_sample.cnt$bb.sum) %>%
  mutate(ff.sum = pr_full_CD4_sample.cnt$ff.sum) %>%
  mutate(g7g7.sum = pr_full_CD4_sample.cnt$g7g7.sum) %>%
  mutate(bxf.sum = pr_full_CD4_sample.cnt$bxf.sum) %>%
  mutate(bxg7.sum = pr_full_CD4_sample.cnt$bxg7.sum) %>%
  mutate(`b+f.sum` = pr_full_CD4_sample.cnt$`b+f.sum`) %>%
  mutate(`b+g7.sum` = pr_full_CD4_sample.cnt$`b+g7.sum`) %>%
  
  # Select out the columns which are relevant.
  select(CDR3.aa, V.name, J.name,
         bb.ave.prop, bb.sum, bb.Samples,
         ff.ave.prop, ff.sum, ff.Samples,
         g7g7.ave.prop, g7g7.sum, g7g7.Samples,
         bxf.ave.prop, bxf.sum, bxf.Samples,
         bxg7.ave.prop, bxg7.sum, bxg7.Samples,
         `b+f.ave.prop`, `b+f.sum`, bf_tetra.Samples,
         `b+g7.ave.prop`, `b+g7.sum`, bg7_tetra.Samples
  )
  
pr_full_CD4_mnml <- CDR.add(MK_convert, pr_full_CD4_mnml)
# Correct public repertoire data sets which should be used for the oligo pool
in.bb.NOT.ff.or.bxf <- pr_full_CD4_mnml %>%
  arrange(desc(bb.ave.prop))%>%
  filter(bb.Samples >= 5, ff.Samples == 0, bxf.Samples == 0) %>%
  select(CDR3.aa, V.name, J.name, CDR1, CDR2, bb.Samples, ff.Samples, bxf.Samples, bf_tetra.Samples, bb.ave.prop, bb.sum) %>%
  mutate(V.CDR3 = paste(V.name, CDR3.aa))

in.ff.NOT.bb.or.bxf <- pr_full_CD4_mnml %>%
  arrange(desc(ff.ave.prop))%>%
  filter(ff.Samples >= 5, bb.Samples == 0, bxf.Samples == 0) %>%
  select(CDR3.aa, V.name, J.name, CDR1, CDR2, ff.Samples, bb.Samples, bxf.Samples, bf_tetra.Samples, ff.ave.prop, ff.sum) %>%
  mutate(V.CDR3 = paste(V.name, CDR3.aa))

in.bb.NOT.g7g7.or.bxg7 <- pr_full_CD4_mnml %>%
  arrange(desc(bb.ave.prop))%>%
  filter(bb.Samples >= 4, g7g7.Samples == 0, bxg7.Samples == 0) %>%
  select(CDR3.aa, V.name, J.name, CDR1, CDR2, bb.Samples, g7g7.Samples, bxg7.Samples, bg7_tetra.Samples, bb.ave.prop, bb.sum) %>%
  mutate(V.CDR3 = paste(V.name, CDR3.aa))

in.g7g7.NOT.bb.or.bxg7 <- pr_full_CD4_mnml %>%
  arrange(desc(g7g7.ave.prop))%>%
  filter(g7g7.Samples >= 3, bb.Samples == 0, bxg7.Samples == 0) %>%
  select(CDR3.aa, V.name, J.name, CDR1, CDR2, g7g7.Samples, bb.Samples, bxg7.Samples, bg7_tetra.Samples, g7g7.ave.prop, g7g7.sum) %>%
  mutate(V.CDR3 = paste(V.name, CDR3.aa))

# NA values are all pretty much appearing because the particular clonotypes are showing up in the F1 or rarely in the other parent and the F1. Can serve as controls but should be minimized.
# Re-Curate list of candidate TCRs for a new Oligo pool. Remove F1 TCRs which appear in the list of parental sequences which happen to share the same TRAV family and CDR3 as a clonotype which appears in the 
subset_ImmuArchList <- function(TCR_DoVb_TR, condition){
  C_meta <- filter(TCR_DoVb_TR$meta, `Co-receptor` == condition | MHC == condition)
  C_samples <- C_meta$Sample
  subset.data <- TCR_DoVb_TR$data %>%
    magrittr::extract(C_samples)
  TCR.NEW_DoVb_TR <- list(
    data = subset.data,
    meta = C_meta
  )
  return(TCR.NEW_DoVb_TR)
}
bxf = "bxf"
bxg7 = "bxg7"

TCR_DoVb_CD4_AA.ones.removed_bxf <- subset_ImmuArchList(TCR_DoVb_CD4_AA.ones.removed, bxf)
TCR_DoVb_CD4_AA.ones.removed_bxg7 <- subset_ImmuArchList(TCR_DoVb_CD4_AA.ones.removed, bxg7)

# Generate a list of the public repertoires of the F1
pr.aav.bxf <- pubRep(TCR_DoVb_CD4_AA.ones.removed_bxf$data, "aa+v+j", .verbose = F)
pr.aav.bxg7 <- pubRep(TCR_DoVb_CD4_AA.ones.removed_bxg7$data, "aa+v+j", .verbose = F)
# 
pr.aav.bxf <- as_tibble(pr.aav.bxf) %>% dplyr::rename(F1_Count = Samples) %>% select(CDR3.aa, V.name, J.name, F1_Count)
pr.aav.bxg7 <- as_tibble(pr.aav.bxg7) %>% dplyr::rename(F1_Count = Samples) %>% select(CDR3.aa, V.name, J.name, F1_Count)
pr.aav.bxf <- CDR.add(MK_convert, pr.aav.bxf)
pr.aav.bxg7 <- CDR.add(MK_convert, pr.aav.bxg7)

# Working lists of in one parent not the other parent and not the F1
# in.bb.NOT.ff.or.bxf
# in.ff.NOT.bb.or.bxf
# in.bb.NOT.g7g7.or.bxg7
# in.g7g7.NOT.bb.or.bxg7 

# TRAVs in use:
# TRAV03-3DN
# TRAV05-4A
# TRAV06-5A 
# TRAV07-6A
# TRAV10-1A 
# TRAV14-3A 

# Plots to see what are good choices for TRAV family members
PR <- function(df, f1.pr){
  out.df <- right_join(f1.pr, df, by = c("CDR3.aa" = "CDR3.aa")) %>%
    dplyr::rename(V.name.F1 = V.name.x) %>%
    filter(paste(CDR1.y, CDR2.y) != paste(CDR1.x, CDR2.x)) %>%
    distinct(V.CDR3, .keep_all = TRUE)
  out.df$bool = 1
  return(out.df)
}

in.bb.NOT.ff.or.bxf_plot <- PR(in.bb.NOT.ff.or.bxf, pr.aav.bxf)
in.ff.NOT.bb.or.bxf_plot <- PR(in.ff.NOT.bb.or.bxf, pr.aav.bxf)
in.bb.NOT.g7g7.or.bxg7_plot <- PR(in.bb.NOT.g7g7.or.bxg7, pr.aav.bxg7)
in.g7g7.NOT.bb.or.bxg7_plot <- PR(in.g7g7.NOT.bb.or.bxg7, pr.aav.bxg7)

Blues.pal <- brewer.pal(6, "Blues")
# Certain Blues are commented out based on replicate occurance which was use in the oligo Pool 2.0 selection criteria.
Plot1.in.bb.NOT.ff.or.bxf <- ggplot(in.bb.NOT.ff.or.bxf_plot[order(in.bb.NOT.ff.or.bxf_plot$bb.Samples), ], aes(y=V.name.y, x=bool, fill = as.factor(bb.Samples))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x ="# clonotypes found in bb \n and fully absent in ff and bxf", y = "TRAV", fill = "Replicate \n Occurance") +
  theme_classic() +
  scale_fill_manual(values = c(
                              # "3" = Blues.pal[3],
                               # "4" = Blues.pal[4],
                               "5" = Blues.pal[5],
                               "6" = Blues.pal[6])) 

Plot2.in.ff.NOT.bb.or.bxf <- ggplot(in.ff.NOT.bb.or.bxf_plot[order(in.ff.NOT.bb.or.bxf_plot$ff.Samples), ], aes(y=V.name.y, x=bool, fill = as.factor(ff.Samples))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x ="# clonotypes found in ff \n and fully absent in bb and bxf", y = "TRAV", fill = "Replicate \n Occurance") +
  theme_classic() +
  scale_fill_manual(values = c(
                               # "3" = Blues.pal[3],
                               # "4" = Blues.pal[4],
                               "5" = Blues.pal[5],
                               "6" = Blues.pal[6])) 

Plot3.in.bb.NOT.g7g7.or.bxg7  <- ggplot(in.bb.NOT.g7g7.or.bxg7_plot[order(in.bb.NOT.g7g7.or.bxg7_plot$bb.Samples), ], aes(y=V.name.y, x= bool, fill = as.factor(bb.Samples))) + 
  geom_bar(stat = "identity") +
  labs(x ="# clonotypes found in bb \n and fully absent in g7g7 and bxg7", y = "TRAV", fill = "Replicate \n Occurance") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = c(
                              # "3" = Blues.pal[3],
                               "4" = Blues.pal[4],
                               "5" = Blues.pal[5],
                               "6" = Blues.pal[6])) 

Plot4.in.g7g7.NOT.bb.or.bxg7 <- ggplot(in.g7g7.NOT.bb.or.bxg7_plot[order(in.g7g7.NOT.bb.or.bxg7_plot$g7g7.Samples), ], aes(y=V.name.y, x=bool, fill = as.factor(g7g7.Samples))) +
  geom_bar(stat = "identity") +
  labs(x ="# clonotypes found in g7g7 \n and fully absent in bb and bxg7", y = "TRAV", fill = "Replicate \n Occurance") +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() + 
  scale_fill_manual(values = c("3" = Blues.pal[3] #,
                               # "4" = Blues.pal[4],
                               # "5" = Blues.pal[5],
                               # "6" = Blues.pal[6]
                               )) 

groupTRAV.plot <- (Plot1.in.bb.NOT.ff.or.bxf + Plot2.in.ff.NOT.bb.or.bxf) / (Plot3.in.bb.NOT.g7g7.or.bxg7 + Plot4.in.g7g7.NOT.bb.or.bxg7)
ggsave('Figure Export/oligo pool 2.0/Canidate oPool TRAVs with no matching CDR1 or CDR2.pdf', groupTRAV.plot,width=7,height=13)


# one way of curating the list is to see what the total clone count is. However, because certain runs have differetn numbers of reads using a frequency based formula would be better.
# for counts the cut off is 50. Determine a frequency cut off:

# Don't remove the clones which have the same V Family member and CDR3 from the F1 sequences... These are probbaly usefull still. Still look and see if there are some which are really simmilar though.
# Remove clones which have the same anticipated CDR1 and CDR2 sequence
TRAV.filter <- function(df, f1.pr, TRAV){
  out.df <- right_join(f1.pr, df, by = c("CDR3.aa" = "CDR3.aa")) %>%
    dplyr::rename(V.name.F1 = V.name.x) %>%
    filter(V.name.y == TRAV) %>%
    filter(paste(CDR1.y, CDR2.y) != paste(CDR1.x, CDR2.x)) %>%
    distinct(V.CDR3, .keep_all = TRUE)
  return(out.df)
}

TRAV03 = 'TRAV03-3DN'
TRAV05 = 'TRAV05-4DN'
TRAV06 = 'TRAV06-5A'
TRAV07 = 'TRAV07-6A'
TRAV10 = 'TRAV10-1A'
TRAV14 = 'TRAV14-3A'
# --------
TRAV3_in.bb.NOT.ff.or.bxf_prop <- TRAV.filter(in.bb.NOT.ff.or.bxf, pr.aav.bxf, TRAV03) %>% arrange(desc(bb.ave.prop)) %>% filter(bb.ave.prop > 1.5e-05) #35
TRAV6_in.bb.NOT.ff.or.bxf_prop <- TRAV.filter(in.bb.NOT.ff.or.bxf, pr.aav.bxf, TRAV06) %>% arrange(desc(bb.ave.prop)) %>% filter(bb.ave.prop > 1.5e-05) #37
TRAV10_in.bb.NOT.ff.or.bxf_prop <- TRAV.filter(in.bb.NOT.ff.or.bxf, pr.aav.bxf, TRAV10) %>% arrange(desc(bb.ave.prop)) %>% filter(bb.ave.prop > 1.5e-05) #42
# --------
# TRAV3_in.ff.NOT.bb.or.bxf_prop <- TRAV.filter(in.ff.NOT.bb.or.bxf, pr.aav.bxf, TRAV03) %>% arrange(desc(ff.ave.prop)) %>% filter(ff.ave.prop > 1.5e-05) #98
TRAV7_in.ff.NOT.bb.or.bxf_prop <- TRAV.filter(in.ff.NOT.bb.or.bxf, pr.aav.bxf, TRAV07) %>% arrange(desc(ff.ave.prop)) %>% filter(ff.ave.prop > 1.5e-05) #116
TRAV14_in.ff.NOT.bb.or.bxf_prop <- TRAV.filter(in.ff.NOT.bb.or.bxf, pr.aav.bxf, TRAV14) %>% arrange(desc(ff.ave.prop)) %>% filter(ff.ave.prop >1.5e-05) #104
# --------
TRAV5_in.g7g7.NOT.bb.or.bxg7_prop <- TRAV.filter(in.g7g7.NOT.bb.or.bxg7, pr.aav.bxg7, TRAV05) %>% arrange(desc(g7g7.ave.prop)) %>% filter(g7g7.ave.prop > 1.5e-05) #49
TRAV6_in.g7g7.NOT.bb.or.bxg7_prop <- TRAV.filter(in.g7g7.NOT.bb.or.bxg7, pr.aav.bxg7, TRAV06) %>% arrange(desc(g7g7.ave.prop)) %>% filter(g7g7.ave.prop > 1.5e-05) #112
# TRAV7_in.g7g7.NOT.bb.or.bxg7_prop <- TRAV.filter(in.g7g7.NOT.bb.or.bxg7, pr.aav.bxg7, TRAV07) %>% arrange(desc(g7g7.ave.prop)) %>% filter(g7g7.ave.prop > 1.5e-05) #24
# --------
TRAV5_in.bb.NOT.g7g7.or.bxg7_prop <- TRAV.filter(in.bb.NOT.g7g7.or.bxg7, pr.aav.bxg7, TRAV05) %>% arrange(desc(bb.ave.prop)) #25
TRAV6_in.bb.NOT.g7g7.or.bxg7_prop <- TRAV.filter(in.bb.NOT.g7g7.or.bxg7, pr.aav.bxg7, TRAV06) %>% arrange(desc(bb.ave.prop)) #53
TRAV10_in.bb.NOT.g7g7.or.bxg7_prop <- TRAV.filter(in.bb.NOT.g7g7.or.bxg7, pr.aav.bxg7, TRAV10) %>% arrange(desc(bb.ave.prop)) #30
# --------
TOPHITS_in.bb.NOT.ff.or.bxf <- rbind(TRAV3_in.bb.NOT.ff.or.bxf_prop, TRAV6_in.bb.NOT.ff.or.bxf_prop, TRAV10_in.bb.NOT.ff.or.bxf_prop) %>% arrange(desc(bb.ave.prop)) %>% head(100)
TOPHITS_in.ff.NOT.bb.or.bxf <- rbind(TRAV7_in.ff.NOT.bb.or.bxf_prop, TRAV14_in.ff.NOT.bb.or.bxf_prop) %>% arrange(desc(ff.ave.prop)) %>% head(100)
TOPHITS_in.bb.NOT.g7g7.or.bxg7 <- rbind(TRAV5_in.bb.NOT.g7g7.or.bxg7_prop, TRAV6_in.bb.NOT.g7g7.or.bxg7_prop, TRAV10_in.bb.NOT.g7g7.or.bxg7_prop) %>% arrange(desc(bb.ave.prop)) %>% head(100)
TOPHITS_in.g7g7.NOT.bb.or.bxg7 <- rbind(TRAV5_in.g7g7.NOT.bb.or.bxg7_prop, TRAV6_in.g7g7.NOT.bb.or.bxg7_prop) %>% arrange(desc(g7g7.ave.prop)) %>% head(100)
# --------
nrow(TOPHITS_in.bb.NOT.ff.or.bxf) + nrow(TOPHITS_in.ff.NOT.bb.or.bxf) + nrow(TOPHITS_in.bb.NOT.g7g7.or.bxg7) + nrow(TOPHITS_in.g7g7.NOT.bb.or.bxg7)
# --------
vroom_write(TOPHITS_in.bb.NOT.ff.or.bxf, "/Volumes/SSD/Analysis/stitchr/Re-checked files/TOPHITS_in.bb.NOT.ff.or.bxf.csv", delim = ",")
vroom_write(TOPHITS_in.ff.NOT.bb.or.bxf, "/Volumes/SSD/Analysis/stitchr/Re-checked files/TOPHITS_in.ff.NOT.bb.or.bxf.csv", delim = ",")
vroom_write(TOPHITS_in.bb.NOT.g7g7.or.bxg7, "/Volumes/SSD/Analysis/stitchr/Re-checked files/TOPHITS_in.bb.NOT.g7g7.or.bxg7.csv", delim = ",")
vroom_write(TOPHITS_in.g7g7.NOT.bb.or.bxg7, "/Volumes/SSD/Analysis/stitchr/Re-checked files/TOPHITS_in.g7g7.NOT.bb.or.bxg7.csv", delim = ",")
# -----------------------------------------
#  VJ position bunch
z.df  <- in.bb.NOT.ff.or.bxf %>% filter(V.name %in% c(TRAV03, TRAV06, TRAV10))

df.pt.1 <- TOPHITS_in.bb.NOT.ff.or.bxf %>% select(V.name.y, J.name.y, bb.ave.prop) %>% rename(TRAV = V.name.y, TRAJ = J.name.y, ave.prop = bb.ave.prop) %>% mutate(Class = "in bb NOT\n ff or bxf") 
df.pt.2 <- TOPHITS_in.ff.NOT.bb.or.bxf %>% select(V.name.y, J.name.y, ff.ave.prop) %>% rename(TRAV = V.name.y, TRAJ = J.name.y, ave.prop = ff.ave.prop) %>% mutate(Class = "in ff NOT\n bb or bxf")
df.pt.3 <- TOPHITS_in.bb.NOT.g7g7.or.bxg7 %>% select(V.name.y, J.name.y, bb.ave.prop) %>% rename(TRAV = V.name.y, TRAJ = J.name.y, ave.prop = bb.ave.prop) %>% mutate(Class = "in bb NOT\n g7g7 or bxg7")
df.pt.4 <- TOPHITS_in.g7g7.NOT.bb.or.bxg7 %>% select(V.name.y, J.name.y, g7g7.ave.prop) %>% rename(TRAV = V.name.y, TRAJ = J.name.y, ave.prop = g7g7.ave.prop) %>% mutate(Class = "in g7g7\n NOT bb or bxg7")
TopHits.400x <- rbind(df.pt.1, df.pt.2, df.pt.3, df.pt.4)

# Set color scale
my_breaks = c(5.0e-05, 1.00e-05, 5.0e-04, 1.00e-04, 0.50e-03, 1.00e-03, 0.50e-02, 1.00e-02)

TopHits.VJ.position_bunch <- ggplot(TopHits.400x, aes(x = TRAJ, y = TRAV, color = ave.prop)) +
  geom_point(position = position_bunch(shape = "hex", width = 1.5), size=0.5, alpha = 0.9) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral")),
                        name = "Clonotype Ave. Freq.",
                        trans = "log",
                        breaks = my_breaks,
                        labels = my_breaks
  ) +
  theme_minimal()+
  facet_grid(rows = vars(Class), scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position = "bottom",
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        plot.title=element_text(hjust=0),
        strip.text=element_text(hjust=0),
        panel.spacing.x=unit(-1.0, "pt"),
        panel.spacing.y=unit(-1.0, "pt"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=7),
        legend.key.size=unit(0.2, "cm"),
        legend.key.width=unit(2.2, "cm"),
        legend.title.align=0.5,
        legend.box.just = "center",
        strip.text.y.right = element_text(angle = 0),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1),
         colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5)
  ) +
  labs(x="TRAJ", y="TRAV")

ggsave("Figure Export/oligo pool 2.0/Canidate oPool TRAVs VJ Freq.pdf", TopHits.VJ.position_bunch, width=170,height=75, units = "mm")
# -----------------------------------------



# Find why Sample 513 and 517 don't have the clone the  say they do...
target <- tibble(CDR3.aa = "CAVSALSNNRIFF", V.name  = "TRAV03-3DN")
tc <- trackClonotypes(TCR_DoVb_CD4_AA.ones.removed$data, target, .col = "aa+v")
vis(tc)
# THEY ARE NOT SHOWING UP BECAUSE THEY ONLY OCCUR ONCE IN THAT GIVEN MOUSE AND WE THROW OUT ALL THE SINGLES FROM ANALYSYS
# Pippa's data shows that these sequences show up 2x and 3x respectively check the data before immune arch...
# ------------------------------------------
# Find sequences which are in both parents, but not in the F1 from the CD4 counts table

# Clones which appear consistently in bb and ff but never in bxf
bb.ff.NOT.bxf_counts <- pr_full_CD4_mnml %>%
  filter(bxf.sum == 0) %>%
  filter(bb.Samples > 2) %>%
  filter(ff.Samples > 2) %>%
  mutate(bb.ff.sum = bb.sum + ff.sum) %>%
  mutate(bb.ff.Samples = bb.Samples + ff.Samples) %>%
  arrange(desc(bb.ff.sum)) %>%
  arrange(desc(bb.ff.Samples)) %>%
  select(CDR3.aa, V.name, J.name, bb.ff.sum, bb.ff.Samples,
         bb.ave.prop, bb.sum, bb.Samples,
         ff.ave.prop, ff.sum, ff.Samples,
         bxf.ave.prop, bxf.sum, bxf.Samples)

# This is a check to make sure the same CDR3 doesn't appear in the bxf groups.
bxf_CDR3 <- pr_full_CD4_mnml %>% filter(bxf.sum > 1) %>% select(CDR3.aa, bxf.Samples)
right_join(bxf_CDR3, bb.ff.NOT.bxf_counts) %>% filter(bxf.Samples > 1) 

# Clones which appear consistently in bb and g7g7 but never in bxg7
bb.g7g7.NOT.bxg7_counts <- pr_full_CD4_mnml %>%
  filter(bxg7.sum == 0) %>%
  filter(bb.Samples > 2) %>%
  filter(g7g7.Samples > 2) %>%
  mutate(bb.g7g7.sum = bb.sum + g7g7.sum) %>%
  mutate(bb.g7g7.Samples = bb.Samples + g7g7.Samples) %>%
  arrange(desc(bb.g7g7.sum)) %>%
  arrange(desc(bb.g7g7.Samples)) %>%
  select(CDR3.aa, V.name, J.name, bb.g7g7.sum, bb.g7g7.Samples,
         bb.ave.prop, bb.sum, bb.Samples,
         g7g7.ave.prop, g7g7.sum, g7g7.Samples,
         bxg7.ave.prop, bxg7.sum, bxg7.Samples)

# This is a check to make sure the same CDR3 doesn't appear in the bxg7 groups.
bxg7_CDR3 <- pr_full_CD4_mnml %>% filter(bxg7.sum > 1) %>% select(CDR3.aa, bxg7.Samples)
right_join(bxg7_CDR3, bb.g7g7.NOT.bxg7_counts) %>% filter(bxg7.Samples > 1)

# Check to make sure the iNKT (TRAV11-TRAJ18) or MAIT (TRAV01-TRAJ33) associated TCRalpha is not in the list
# Manually checked that this is not the case.


# Export
vroom_write(bb.ff.NOT.bxf_counts, "stitchr/sequences not in F1/bb.ff.NOT.bxf_counts.csv", delim = ",")
vroom_write(bb.g7g7.NOT.bxg7_counts, "stitchr/sequences not in F1/bb.g7g7.NOT.bxg7_counts.csv", delim = ",")


# Old code to recheck the oligo pools
# ------------------
# Read in the oligo pool dataset

# Read in the precursor data to the oligo pool dataset. Effectively the same sequences
# in.bb_NOT.bxf_CD4_CDR3.match <- vroom("/Volumes/SSD/Analysis/stitchr/Stitched_VDJ/in.bb_NOT.bxf_CD4_CDR3.match.csv")
# in.ff_NOT.bxf_CD4_CDR3.match <- vroom("/Volumes/SSD/Analysis/stitchr/Stitched_VDJ/in.ff_NOT.bxf_CD4_CDR3.match.csv")
# in.bb_NOT.bxg7_CD4_CDR3.match <- vroom("/Volumes/SSD/Analysis/stitchr/Stitched_VDJ/in.bb_NOT.bxg7_CD4_CDR3.match.csv")
# in.g7g7_NOT.bxg7_CD4_CDR3.match <- vroom("/Volumes/SSD/Analysis/stitchr/Stitched_VDJ/in.g7g7_NOT.bxg7_CD4_CDR3.match.csv")
# # --
# in.bb_NOT.bxf_CD4_CDR3.match_select <- in.bb_NOT.bxf_CD4_CDR3.match %>% mutate(V.CDR3 = paste(V.name.y, CDR3.aa)) %>% select(TCR_name, V.CDR3)
# in.ff_NOT.bxf_CD4_CDR3.match_select <- in.ff_NOT.bxf_CD4_CDR3.match %>% mutate(V.CDR3 = paste(V.name.y, CDR3.aa)) %>% select(TCR_name, V.CDR3)
# in.bb_NOT.bxg7_CD4_CDR3.match_select <- in.bb_NOT.bxg7_CD4_CDR3.match %>% mutate(V.CDR3 = paste(V.name.y, CDR3.aa)) %>% select(TCR_name, V.CDR3)
# in.g7g7_NOT.bxg7_CD4_CDR3.match_select <- in.g7g7_NOT.bxg7_CD4_CDR3.match %>% mutate(V.CDR3 = paste(V.name.y, CDR3.aa)) %>% select(TCR_name, V.CDR3)
# # --
# RECHECKED.in.bb.NOT.ff.or.bxf <- right_join(in.bb.NOT.ff.or.bxf, in.bb_NOT.bxf_CD4_CDR3.match_select, by = c("V.CDR3" = "V.CDR3"))
# RECHECKED.in.ff.NOT.bb.or.bxf <- right_join(in.ff.NOT.bb.or.bxf, in.ff_NOT.bxf_CD4_CDR3.match_select, by = c("V.CDR3" = "V.CDR3"))
# RECHECKED.in.bb.NOT.g7g7.or.bxg7 <- right_join(in.bb.NOT.g7g7.or.bxg7, in.bb_NOT.bxg7_CD4_CDR3.match_select, by = c("V.CDR3" = "V.CDR3"))
# RECHECKED.in.g7g7.NOT.bb.or.bxg7 <- right_join(in.g7g7.NOT.bb.or.bxg7, in.g7g7_NOT.bxg7_CD4_CDR3.match_select, by = c("V.CDR3" = "V.CDR3"))

# vroom_write(RECHECKED.in.bb.NOT.ff.or.bxf, "/Volumes/SSD/Analysis/stitchr/Re-checked files/RECHECKED.in.bb.NOT.ff.or.bxf.csv", delim = ",")
# vroom_write(RECHECKED.in.ff.NOT.bb.or.bxf, "/Volumes/SSD/Analysis/stitchr/Re-checked files/RECHECKED.in.ff.NOT.bb.or.bxf.csv", delim = ",")
# vroom_write(RECHECKED.in.bb.NOT.g7g7.or.bxg7, "/Volumes/SSD/Analysis/stitchr/Re-checked files/RECHECKED.in.bb.NOT.g7g7.or.bxg7.csv", delim = ",")
# vroom_write(RECHECKED.in.g7g7.NOT.bb.or.bxg7, "/Volumes/SSD/Analysis/stitchr/Re-checked files/RECHECKED.in.g7g7.NOT.bb.or.bxg7.csv", delim = ",")

# ---
# Find out why NA values are present

# in.bb_NOT.bxf_NA <- vroom("/Volumes/SSD/Analysis/stitchr/Re-checked files/in.bb_NOT.bxf_NA.csv")
# in.bb_NOT.bxf_NA <- in.bb_NOT.bxf_NA %>% select(CDR3.aa, V.name)

# pr_full_CD4_mnml_tst <- pr_full_CD4_mnml %>%
#   select(CDR3.aa, V.name, J.name, bb.Samples, ff.Samples, g7g7.Samples, bxf.Samples, bxg7.Samples) %>%
#   mutate(V.CDR3 = paste(V.name, CDR3.aa))
# 
#  
# NA.in.bb.NOT.ff.or.bxf <- RECHECKED.in.bb.NOT.ff.or.bxf %>% filter(is.na(CDR3.aa)) %>% select(V.CDR3, TCR_name)
# NA.in.ff.NOT.bb.or.bxf <- RECHECKED.in.ff.NOT.bb.or.bxf %>% filter(is.na(CDR3.aa)) %>% select(V.CDR3, TCR_name)
# NA.in.bb.NOT.g7g7.or.bxg7 <- RECHECKED.in.bb.NOT.g7g7.or.bxg7 %>% filter(is.na(CDR3.aa)) %>% select(V.CDR3, TCR_name)
# NA.in.g7g7.NOT.bb.or.bxg7 <- RECHECKED.in.g7g7.NOT.bb.or.bxg7 %>% filter(is.na(CDR3.aa)) %>% select(V.CDR3, TCR_name)
# 
# 
# NA.in.bb.NOT.ff.or.bxf.2 <- right_join(pr_full_CD4_mnml_tst, NA.in.bb.NOT.ff.or.bxf, by = c("V.CDR3" = "V.CDR3"))
# NA.in.ff.NOT.bb.or.bxf.2 <- right_join(pr_full_CD4_mnml_tst, NA.in.ff.NOT.bb.or.bxf, by = c("V.CDR3" = "V.CDR3"))
# NA.in.bb.NOT.g7g7.or.bxg7.2 <- right_join(pr_full_CD4_mnml_tst, NA.in.bb.NOT.g7g7.or.bxg7, by = c("V.CDR3" = "V.CDR3"))
# NA.in.g7g7.NOT.bb.or.bxg7.2 <- right_join(pr_full_CD4_mnml_tst, NA.in.g7g7.NOT.bb.or.bxg7, by = c("V.CDR3" = "V.CDR3"))

