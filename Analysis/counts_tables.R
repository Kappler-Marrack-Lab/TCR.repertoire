# Purpose:
# Create a new counts table.
# Table 1: Contains original data minus sequences which occur once in TRAV-CDR3AA-TRAJ
# Table 2: CDR3NT sequences are aggrigated, frequency is recalcualted after agrigation. a frequency cutoff is applied. ones are then removed.
# For both the frequency values reported are based on the total productive reads in the run (ie, not filtered and ones retained).

# library loading:
library(immunarch)
library(tidyverse)

# Load source data:
load("TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData")
load("TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData")

# Functions:
# A function to output a counts table of the unique sequences per group and include group statistics
pr_filter <- function(rep.list){
  
  df.q <- pubRep(rep.list$data, "v+j+aa", .quant =c("prop"), .verbose = F) %>% as_tibble() # frequency based calculation
  df.c <- pubRep(rep.list$data, "v+j+aa", .quant =c("count"), .verbose = F) %>% as_tibble() # frequency based calculation
  
  
  
  # Calculate the group total productive reads for each MHC type in use:
  bb_totals <- rep.list$meta %>% filter(MHC == "bb") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  ff_totals <- rep.list$meta %>% filter(MHC == "ff") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  g7g7_totals <- rep.list$meta %>% filter(MHC == "g7g7") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  ss_totals <- rep.list$meta %>% filter(MHC == "ss") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  bxf_totals <- rep.list$meta %>% filter(MHC == "bxf") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  bxg7_totals <- rep.list$meta %>% filter(MHC == "bxg7") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  bxs_totals <- rep.list$meta %>% filter(MHC == "bxs") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  fxs_totals <- rep.list$meta %>% filter(MHC == "fxs") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  `bb&ff_totals` <- rep.list$meta %>% filter(MHC == "bb&ff") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  `bb&g7g7_totals` <- rep.list$meta %>% filter(MHC == "bb&g7g7") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  `b+/-_totals` <- rep.list$meta %>% filter(MHC == "b+/-") %>% summarise(Total = sum(Total_Productive_Reads))  %>% pull(Total)
  
  
  
  df.c.out <- df.c %>%
    dplyr::mutate(bb.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb') %>% pull(Sample)))))%>%
    dplyr::mutate(ff.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'ff') %>% pull(Sample)))))%>%
    dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)))))%>%
    dplyr::mutate(ss.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'ss') %>% pull(Sample))))) %>%
    dplyr::mutate(bxf.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxf') %>% pull(Sample)))))%>%
    dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)))))%>%  
    dplyr::mutate(fxs.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'fxs') %>% pull(Sample))))) %>%
    dplyr::mutate(bxs.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxs') %>% pull(Sample))))) %>%
    dplyr::mutate(bf_tetra.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample))))) %>%
    dplyr::mutate(bg7_tetra.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample))))) %>%
    dplyr::mutate(`b+/-.Samples` = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'b+/-') %>% pull(Sample))))) %>%
    
    dplyr::mutate(bb.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bb') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(ff.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'ff') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(g7g7.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(ss.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'ss') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(bxf.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxf') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(bxg7.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(fxs.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'fxs') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(bxs.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxs') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(`bb&ff.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(`bb&g7g7.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
    dplyr::mutate(`b+/-.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'b+/-') %>% pull(Sample)), na.rm = TRUE)) %>%
    
    # mutate(bb.ave.prop = (bb.sum/max(bb.Samples))/sum(bb.sum/max(bb.Samples))) %>%
    # mutate(ff.ave.prop = (ff.sum/max(ff.Samples))/sum(ff.sum/max(ff.Samples))) %>%
    # mutate(g7g7.ave.prop = (g7g7.sum/max(g7g7.Samples))/sum(g7g7.sum/max(g7g7.Samples))) %>%
    # mutate(ss.ave.prop = (ss.sum/max(ss.Samples))/sum(ss.sum/max(ss.Samples))) %>%
    # mutate(bxf.ave.prop = (bxf.sum/max(bxf.Samples))/sum(bxf.sum/max(bxf.Samples))) %>%
    # mutate(bxg7.ave.prop = (bxg7.sum/max(bxg7.Samples))/sum(bxg7.sum/max(bxg7.Samples))) %>%
    # mutate(fxs.ave.prop = (fxs.sum/max(fxs.Samples))/sum(fxs.sum/max(fxs.Samples))) %>%
    # mutate(bxs.ave.prop = (bxs.sum/max(bxs.Samples))/sum(bxs.sum/max(bxs.Samples))) %>%
    # mutate(`bb&ff.ave.prop` = (`bb&ff.sum`/max(bf_tetra.Samples))/sum(`bb&ff.sum`/max(bf_tetra.Samples))) %>%
    # mutate(`bb&g7g7.ave.prop` = (`bb&g7g7.sum`/max(bg7_tetra.Samples))/sum(`bb&g7g7.sum`/max(bg7_tetra.Samples))) %>%
    # mutate(`b+/-.ave.prop` = (`b+/-.sum`/max(`b+/-.Samples`))/sum(`b+/-.sum`/max(`b+/-.Samples`))) %>%

    mutate(bb.ave.prop = bb.sum/bb_totals) %>%
    mutate(ff.ave.prop = ff.sum/ff_totals) %>%
    mutate(g7g7.ave.prop = g7g7.sum/g7g7_totals) %>%
    mutate(ss.ave.prop = ss.sum/ss_totals) %>%
    mutate(bxf.ave.prop = bxf.sum/bxf_totals) %>%
    mutate(bxg7.ave.prop = bxg7.sum/bxg7_totals) %>%
    mutate(bxs.ave.prop = bxs.sum/bxs_totals) %>%
    mutate(fxs.ave.prop = fxs.sum/fxs_totals) %>%
    mutate(`bb&ff.ave.prop` = `bb&ff.sum`/`bb&ff_totals`) %>%
    mutate(`bb&g7g7.ave.prop` = `bb&g7g7.sum`/`bb&g7g7_totals`) %>%
    mutate(`b+/-.ave.prop` = `b+/-.sum`/`b+/-_totals`) %>%
    
    
    dplyr::select(
      V.name, J.name, CDR3.aa,
      bb.ave.prop,
      ff.ave.prop,
      g7g7.ave.prop,
      ss.ave.prop,
      bxf.ave.prop,
      bxg7.ave.prop,
      fxs.ave.prop,
      bxs.ave.prop,
      `bb&ff.ave.prop`,
      `bb&g7g7.ave.prop`,
      `b+/-.ave.prop`
    )
  
  
  
  
  
  
  df_out <- df.q %>%
    # Sums the number of mice (samples) which a given clone appears
    dplyr::mutate(bb.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb') %>% pull(Sample)))))%>%
    dplyr::mutate(ff.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'ff') %>% pull(Sample)))))%>%
    dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)))))%>%
    dplyr::mutate(ss.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'ss') %>% pull(Sample))))) %>%
    dplyr::mutate(bxf.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxf') %>% pull(Sample)))))%>%
    dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)))))%>%
    dplyr::mutate(fxs.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'fxs') %>% pull(Sample))))) %>%
    dplyr::mutate(bxs.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bxs') %>% pull(Sample))))) %>%
    dplyr::mutate(bf_tetra.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample))))) %>%
    dplyr::mutate(bg7_tetra.Samples = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample))))) %>%
    dplyr::mutate(`b+/-.Samples` = rowSums(!is.na(across(rep.list$meta %>% filter(MHC == 'b+/-') %>% pull(Sample))))) %>%
    
    
    # dplyr::mutate(bb.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bb') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(ff.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'ff') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(g7g7.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(ss.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'ss') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(bxf.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxf') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(bxg7.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(fxs.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'fxs') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(bxs.sum = rowSums(across(rep.list$meta %>% filter(MHC == 'bxs') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(`bb&ff.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(`bb&g7g7.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample)), na.rm = TRUE)) %>%
    # dplyr::mutate(`b+/-.sum` = rowSums(across(rep.list$meta %>% filter(MHC == 'b+/-') %>% pull(Sample)), na.rm = TRUE)) %>%
    # 
    # mutate(bb.ave.prop = (bb.sum/max(bb.Samples))/sum(bb.sum/max(bb.Samples))) %>%
    # mutate(ff.ave.prop = (ff.sum/max(ff.Samples))/sum(ff.sum/max(ff.Samples))) %>%
    # mutate(g7g7.ave.prop = (g7g7.sum/max(g7g7.Samples))/sum(g7g7.sum/max(g7g7.Samples))) %>%
    # mutate(ss.ave.prop = (ss.sum/max(ss.Samples))/sum(ss.sum/max(ss.Samples))) %>%
    # mutate(bxf.ave.prop = (bxf.sum/max(bxf.Samples))/sum(bxf.sum/max(bxf.Samples))) %>%
    # mutate(bxg7.ave.prop = (bxg7.sum/max(bxg7.Samples))/sum(bxg7.sum/max(bxg7.Samples))) %>%
    # mutate(fxs.ave.prop = (fxs.sum/max(fxs.Samples))/sum(fxs.sum/max(fxs.Samples))) %>%
    # mutate(bxs.ave.prop = (bxs.sum/max(bxs.Samples))/sum(bxs.sum/max(bxs.Samples))) %>%
    # mutate(`bb&ff.ave.prop` = (`bb&ff.sum`/max(bf_tetra.Samples))/sum(`bb&ff.sum`/max(bf_tetra.Samples))) %>%
    # mutate(`bb&g7g7.ave.prop` = (`bb&g7g7.sum`/max(bg7_tetra.Samples))/sum(`bb&g7g7.sum`/max(bg7_tetra.Samples))) %>%
  # mutate(`b+/-.ave.prop` = (`b+/-.sum`/max(`b+/-.Samples`))/sum(`b+/-.sum`/max(`b+/-.Samples`))) %>%
  # 
  
  # mutate(bb.ave.prop = (bb.sum/max(bb.Samples))) %>%
  # mutate(ff.ave.prop = (ff.sum/max(ff.Samples))) %>%
  # mutate(g7g7.ave.prop = (g7g7.sum/max(g7g7.Samples))) %>%
  # mutate(ss.ave.prop = (ss.sum/max(ss.Samples))) %>%
  # mutate(bxf.ave.prop = (bxf.sum/max(bxf.Samples))) %>%
  # mutate(bxg7.ave.prop = (bxg7.sum/max(bxg7.Samples))) %>%
  # mutate(fxs.ave.prop = (fxs.sum/max(fxs.Samples))) %>%
  # mutate(bxs.ave.prop = (bxs.sum/max(bxs.Samples))) %>%
  # mutate(`bb&ff.ave.prop` = (`bb&ff.sum`/max(bf_tetra.Samples))) %>%
  # mutate(`bb&g7g7.ave.prop` = (`bb&g7g7.sum`/max(bg7_tetra.Samples))) %>%
  # mutate(`b+/-.ave.prop` = (`b+/-.sum`/max(`b+/-.Samples`))) %>%
  
  # Reference the average proportion from the numerical counts table:
  mutate(bb.ave.prop = df.c.out$bb.ave.prop) %>%
    mutate(ff.ave.prop = df.c.out$ff.ave.prop) %>%
    mutate(g7g7.ave.prop = df.c.out$g7g7.ave.prop) %>%
    mutate(ss.ave.prop = df.c.out$ss.ave.prop) %>%
    mutate(bxf.ave.prop = df.c.out$bxf.ave.prop) %>%
    mutate(bxg7.ave.prop = df.c.out$bxg7.ave.prop) %>%
    mutate(fxs.ave.prop = df.c.out$fxs.ave.prop) %>%
    mutate(bxs.ave.prop = df.c.out$bxs.ave.prop) %>%
    mutate(`bb&ff.ave.prop` = df.c.out$`bb&ff.ave.prop`) %>%
    mutate(`bb&g7g7.ave.prop` = df.c.out$`bb&g7g7.ave.prop`) %>%
    mutate(`b+/-.ave.prop` = df.c.out$`b+/-.ave.prop`) %>%
    
    
    
    dplyr::select(
      V.name, J.name, CDR3.aa,
      bb.Samples, bb.ave.prop, c(rep.list$meta %>% filter(MHC == 'bb') %>% pull(Sample)),
      ff.Samples, ff.ave.prop, c(rep.list$meta %>% filter(MHC == 'ff') %>% pull(Sample)),
      g7g7.Samples, g7g7.ave.prop, c(rep.list$meta %>% filter(MHC == 'g7g7') %>% pull(Sample)),
      ss.Samples, ss.ave.prop, c(rep.list$meta %>% filter(MHC == 'ss') %>% pull(Sample)),
      bxf.Samples, bxf.ave.prop, c(rep.list$meta %>% filter(MHC == 'bxf') %>% pull(Sample)),
      bxg7.Samples, bxg7.ave.prop, c(rep.list$meta %>% filter(MHC == 'bxg7') %>% pull(Sample)),
      fxs.Samples, fxs.ave.prop, c(rep.list$meta %>% filter(MHC == 'fxs') %>% pull(Sample)),
      bxs.Samples, bxs.ave.prop, c(rep.list$meta %>% filter(MHC == 'bxs') %>% pull(Sample)),
      bf_tetra.Samples, `bb&ff.ave.prop`, c(rep.list$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample)),
      bg7_tetra.Samples, `bb&g7g7.ave.prop`, c(rep.list$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample)),
      `b+/-.Samples`, `b+/-.ave.prop`, c(rep.list$meta %>% filter(MHC == 'b+/-') %>% pull(Sample))
    )
  
  df_out <- df_out %>% mutate_all(~replace(., is.na(.), 0))
  
  return(df_out)
}


# Generate the organized counts table from the rep data:
pr_CD4s_filtered.ones.removed <- pr_filter(TCR_DoVb_CD4_AA.ones.removed)
pr_CD8s_filtered.ones.removed <- pr_filter(TCR_DoVb_CD8_AA.ones.removed)

pr_CD4s_filtered.LFSR.ones.removed <- pr_filter(TCR_DoVb_CD4_AA.ones.removed.LFSR)
pr_CD8s_filtered.LFSR.ones.removed <- pr_filter(TCR_DoVb_CD8_AA.ones.removed.LFSR)

# Save objects
save(pr_CD4s_filtered.ones.removed, pr_CD8s_filtered.ones.removed,
     file = "counts tables filtered.ones.removed.RData")

# Save objects
save(pr_CD4s_filtered.LFSR.ones.removed, pr_CD8s_filtered.LFSR.ones.removed,
     file = "counts tables filtered.LFSR.ones.removed.RData")



print("task finished.")
print("Data stored in:")
print("counts tables filtered.ones.removed.RData")
print("counts tables filtered.LFSR.ones.removed.RData")


library(vroom)

# Export tables to CSV files
vroom_write(pr_CD4s_filtered.ones.removed, file = "Figure Export/Counts Tables/pr_CD4s_filtered.ones.removed.csv", delim = ",")
vroom_write(pr_CD8s_filtered.ones.removed, file = "Figure Export/Counts Tables/pr_CD8s_filtered.ones.removed.csv", delim = ",")

vroom_write(pr_CD4s_filtered.LFSR.ones.removed, file = "Figure Export/Counts Tables/pr_CD4s_filtered.LFSR.ones.removed.csv", delim = ",")
vroom_write(pr_CD8s_filtered.LFSR.ones.removed, file = "Figure Export/Counts Tables/pr_CD8s_filtered.LFSR.ones.removed.csv", delim = ",")
