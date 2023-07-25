library(immunarch)
library(tidyverse)
library(tidyr)
library(ggrepel)
require(tidyr)
library(emdbook)
library(grid)
library(eulerr)
library(viridis)
library(vroom)
library(stringr)
set.seed(42)

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


load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed
load(file = "counts tables filtered.LFSR.ones.removed.RData")

# Make sure the number of samples in each circle is the same
pr_CD4_merge <- pubRep(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, "v+j+aa", .verbose = T)

# bb_samples = c("264 postJ", "267 postJ", "505 postj.one_mouse")
# bxf_samples = c("394 postJ" , "421 postJ" , "461LNC.SC_postPy.one_mouse")
# ff_samples = c("222 postj", "228 postJ" , "546 postJ.one_mouse")
MHC.counts <- TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% count(MHC)
bb.ff_MHC.counts <- MHC.counts %>% filter(MHC %in% c("bb", "bxf", "ff")) #this will allow you to dynamically identify the maximum number of mice per group to use in the euler diagrams
bb.g7g7_MHC.counts <- MHC.counts %>% filter(MHC %in% c("bb", "bxg7", "g7g7"))  #this will allow you to dynamically identify the maximum number of mice per group to use in the euler diagrams
bb.ss_MHC.counts <- MHC.counts %>% filter(MHC %in% c("bb", "bxs", "ss"))
ff.ss_MHC.counts <- MHC.counts %>% filter(MHC %in% c("ff", "fxs", "ss"))



# Select minimum sample count in the group to equalize the samples being evaluated
bb.ff_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n)))))) %>%
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n)))))) %>%
  dplyr::mutate(bxf.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n)))))) %>%
  dplyr::mutate(ff.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(bb.ff_MHC.counts$n))), na.rm = TRUE)) %>%
  dplyr::filter(bb.sum + bxf.sum + ff.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         bb.Samples, bxf.Samples, ff.Samples,
         bb.sum, bxf.sum, ff.sum
         )
  
bb.g7g7_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n)))))) %>%
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n)))))) %>%
  dplyr::mutate(bxg7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n)))))) %>%
  dplyr::mutate(g7g7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample) %>% head(n=min(bb.g7g7_MHC.counts$n))), na.rm = TRUE)) %>%
  dplyr::filter(bb.sum + bxg7.sum + g7g7.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         bb.Samples, bxg7.Samples, g7g7.Samples,
         bb.sum, bxg7.sum, g7g7.sum
         )

bb.ss_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(bxs.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxs') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(bxs.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxs') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(ss.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(ss.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample) %>% head(n=min(bb.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  dplyr::filter(bb.sum + bxs.sum + ss.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         bb.Samples, bxs.Samples, ss.Samples,
         bb.sum, bxs.sum, ss.sum
  )

ff.ss_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(ff.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(fxs.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'fxs') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(fxs.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'fxs') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(ss.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n)))))) %>%
  dplyr::mutate(ss.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample) %>% head(n=min(ff.ss_MHC.counts$n))), na.rm = TRUE)) %>%
  dplyr::filter(ff.sum + fxs.sum + ss.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         ff.Samples, fxs.Samples, ss.Samples,
         ff.sum, fxs.sum, ss.sum
  )
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# substitute F1 for tetraparentals
`bb&ff_MHC.counts` <- MHC.counts %>% filter(MHC %in% c("bb", "bb&ff", "ff", "bxf")) #this will allow you to dynamically identify the maximum number of mice per group to use in the euler diagrams
`bb&g7g7_MHC.counts` <- MHC.counts %>% filter(MHC %in% c("bb", "bb&g7g7", "g7g7", "bxg7"))  #this will allow you to dynamically identify the maximum number of mice per group to use in the euler diagrams

bb.ff.tetra_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n)))))) %>%
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n)))))) %>%
  dplyr::mutate(tetra.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb&ff') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n)))))) %>%
  dplyr::mutate(ff.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n)))))) %>%
  dplyr::mutate(bxf.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample) %>% head(n=min(`bb&ff_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::filter(bb.sum + tetra.sum + ff.sum + bxf.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         bb.Samples, tetra.Samples, ff.Samples, bxf.Samples,
         bb.sum, tetra.sum, ff.sum, bxf.sum
  )

bb.g7g7.tetra_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n)))))) %>%
  dplyr::mutate(bb.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(tetra.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n)))))) %>%
  dplyr::mutate(tetra.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb&g7g7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n)))))) %>%
  dplyr::mutate(g7g7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n)))))) %>%
  dplyr::mutate(bxg7.sum = rowSums(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample) %>% head(n=min(`bb&g7g7_MHC.counts`$n))), na.rm = TRUE)) %>%
  
  dplyr::filter(bb.sum + tetra.sum + g7g7.sum + bxg7.sum > 0) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa,
         bb.Samples, tetra.Samples, g7g7.Samples, bxg7.Samples,
         bb.sum, tetra.sum, g7g7.sum, bxg7.sum
  )
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Euler plots don't pay attention to the number of clones present per haplotype just are they present or not.
# euler function
euler_clonotypes <- function(clone.counts, Par1_color, F1_color, Par2_color) {
  # par 1
  Par1 <-  clone.counts %>%
    filter(.[[5]] > 0) %>% select(clonotype)
  # F1 
  F1 <- clone.counts %>%
    filter(.[[6]] > 0) %>% select(clonotype)
  # par 2
  Par2 <- clone.counts %>%
    filter(.[[7]] > 0) %>% select(clonotype)

  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    F1  = F1$clonotype,
    Par2 = Par2$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c(names(clone.counts)[5], names(clone.counts)[6], names(clone.counts)[7])
  map.name <- str_sub(map.name,1,nchar(map.name)-8)
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x)
  
  # generate the euler diagrams with desired colors
  euler_plot <- plot(eu.x,
       quantities = TRUE,
       fills = c(Par1_color, F1_color, Par2_color)
  )
  
 return(euler_plot) 
}


euler_clonotypes(bb.ff_clone.counts, bb_color, bxf_color, ff_color)
euler_clonotypes(bb.g7g7_clone.counts, bb_color, bxg7_color, g7g7_color)
euler_clonotypes(bb.ss_clone.counts, bb_color, bxs_color, ss_color)
euler_clonotypes(ff.ss_clone.counts, ff_color, fxs_color, ss_color)

ggsave('Figure Export/Euler Diagram/CD4 bb vs ff vs bxf.pdf', euler_clonotypes(bb.ff_clone.counts, bb_color, bxf_color, ff_color) ,width=3,height=3)
ggsave('Figure Export/Euler Diagram/CD4 bb vs g7g7 vs bxg7.pdf', euler_clonotypes(bb.g7g7_clone.counts, bb_color, bxg7_color, g7g7_color) ,width=3,height=3)
ggsave('Figure Export/Euler Diagram/CD4 bb vs ss vs bxs.pdf', euler_clonotypes(bb.ss_clone.counts, bb_color, bxs_color, ss_color) ,width=3,height=3)
ggsave('Figure Export/Euler Diagram/CD4 ff vs ss vs fxs.pdf', euler_clonotypes(ff.ss_clone.counts, ff_color, fxs_color, ss_color) ,width=3,height=3)

# -----------------------------------------------------------------------------------------------------------------
# tetra plots:
euler_clonotypes(bb.ff.tetra_clone.counts, bb_color, tetra_bb.ff_color, ff_color)
euler_clonotypes(bb.g7g7.tetra_clone.counts, bb_color, tetra_bb.g7g7_color, g7g7_color)

# -----------------------------------------------------------------------------------------------------------------
# maually select out the samples to use in the plots:
# bb, ff, bxf 5x pergroup
b.f.bxf_dat <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(
  MHC = include("bb", "bxf", "ff"),
  Sample = exclude("394_postJ", "463_465LNC.SC_postPy.two_mice")
))
pr_b.f.bxf_dat <- pubRep(b.f.bxf_dat$data, "v+j+aa", .verbose = F) %>% as_tibble() %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(b.f.bxf_dat$meta %>% filter(MHC == 'bb') %>% pull(Sample)))))%>%
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(b.f.bxf_dat$meta %>% filter(MHC == 'bxf') %>% pull(Sample)))))%>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(b.f.bxf_dat$meta %>% filter(MHC == 'ff') %>% pull(Sample)))))%>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa, bb.Samples, bxf.Samples, ff.Samples)




euler_clonotypes(pr_b.f.bxf_dat, bb_color, bxf_color, ff_color)

# test that the number of unique clones in the F1 are correct
pr_b.f.bxf_dat %>%
  mutate(sum = rowSums(across(c(bb.Samples, ff.Samples)))) %>%
  filter(sum == 0)


# -----------------------------------------------------------------------------------------------------------------
# maually select out the samples to use in the plots:
# ss, ss, fxs 5x pergroup
f.s.fxs_dat <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(
  MHC = include("ff", "fxs", "ss"),
  Sample = include(
    "222 postj", "546 postJ.two_mice", "550_postJ_two_mice", # ff
    "255 postj", "306 postJ", "438 postJ", # ss
    "397 postJ",	"400 postJ", "403 postJ" # fxs
  )
))
pr_f.s.fxs_dat <- pubRep(f.s.fxs_dat$data, "v+j+aa", .verbose = F) %>% as_tibble() %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(f.s.fxs_dat$meta %>% filter(MHC == 'ff') %>% pull(Sample)))))%>%
  dplyr::mutate(fxs.Samples = rowSums(!is.na(across(f.s.fxs_dat$meta %>% filter(MHC == 'fxs') %>% pull(Sample)))))%>%
  dplyr::mutate(ss.Samples = rowSums(!is.na(across(f.s.fxs_dat$meta %>% filter(MHC == 'ss') %>% pull(Sample)))))%>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa, ff.Samples, fxs.Samples, ss.Samples)

euler_clonotypes(pr_f.s.fxs_dat, ff_color, fxs_color, ss_color)
# -----------------------------------------------------------------------------------------------------------------
b.s.bxs_dat <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(
  MHC = include("bb", "bxs", "ss"),
  Sample = include(
    "270 postJ", "505_postJ_two_mice", "264 postJ", # bb
    "255 postj", "306 postJ", "438 postJ", # ss
    "406 postJ", "409 postJ", "412 postJ" # bxs
  )
))
pr_b.s.bxs_dat <- pubRep(b.s.bxs_dat$data, "v+j+aa", .verbose = F) %>% as_tibble() %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(b.s.bxs_dat$meta %>% filter(MHC == 'bb') %>% pull(Sample)))))%>%
  dplyr::mutate(bxs.Samples = rowSums(!is.na(across(b.s.bxs_dat$meta %>% filter(MHC == 'bxs') %>% pull(Sample)))))%>%
  dplyr::mutate(ss.Samples = rowSums(!is.na(across(b.s.bxs_dat$meta %>% filter(MHC == 'ss') %>% pull(Sample)))))%>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, V.name, J.name, CDR3.aa, bb.Samples, bxs.Samples, ss.Samples)

euler_clonotypes(pr_b.s.bxs_dat, bb_color, bxs_color, ss_color)
# -----------------------------------------------------------------------------------------------------------------

ggsave('Figure Export/Euler Diagram/CD4 bb vs ss vs bxs_v2.pdf', euler_clonotypes(pr_b.s.bxs_dat, bb_color, bxs_color, ss_color) ,width=3,height=3)
ggsave('Figure Export/Euler Diagram/CD4 ff vs ss vs fxs_v2.pdf', euler_clonotypes(pr_f.s.fxs_dat, ff_color, fxs_color, ss_color) ,width=3,height=3)
# -----------------------------------------------------------------------------------------------------------------







# -----------------------------------------------------------------------------------------------------------------
# Filter only bb samples:
bb.5x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=5)) %>%
  replace(is.na(.), 0)
names(bb.5x_clone.counts) <- c("clonotype", "bb Mouse #1", "bb Mouse #2", "bb Mouse #3", "bb Mouse #4", "bb Mouse #5")
bb.1 <-  bb.5x_clone.counts %>% filter(.[[2]] > 0)
bb.2 <-  bb.5x_clone.counts %>% filter(.[[3]] > 0)
bb.3 <-  bb.5x_clone.counts %>% filter(.[[4]] > 0)
bb.4 <-  bb.5x_clone.counts %>% filter(.[[5]] > 0)
bb.5 <-  bb.5x_clone.counts %>% filter(.[[6]] > 0)

bb.list <- list(
  n1 = bb.1$clonotype,
  n2 = bb.2$clonotype,
  n3 = bb.3$clonotype,
  n4 = bb.4$clonotype,
  n5 = bb.5$clonotype
)

names(bb.list) <- c("bb Mouse #1", "bb Mouse #2", "bb Mouse #3", "bb Mouse #4", "bb Mouse #5")
eu.bb <- euler(bb.list)
eu.bb_plot <- plot(eu.bb, fill_alpha=0.7,
                   quantities = TRUE,
                   fills = c("#66cdaa", "#ffa500", "#00ff00", "#0000ff", "#1e90ff")
)
ggsave('Figure Export/Euler Diagram/CD4 bb compiled euler diagram.pdf', eu.bb_plot,width=5,height=5)
# -------------------------------------------------------------------------------------------------------------------------------


# make a concentric euler diagram to illistrate how many clones are shared in the different animals:
ff.6x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample) %>% head(n=6)) %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ff') %>% pull(Sample))))) %>%
  filter(ff.Samples > 0) %>%
  replace(is.na(.), 0)
  
ff.1 <-  ff.6x_clone.counts %>% filter(ff.Samples > 0)
ff.2 <-  ff.6x_clone.counts %>% filter(ff.Samples > 1)
ff.3 <-  ff.6x_clone.counts %>% filter(ff.Samples > 2)
ff.4 <-  ff.6x_clone.counts %>% filter(ff.Samples > 3)
ff.5 <-  ff.6x_clone.counts %>% filter(ff.Samples > 4)
ff.6 <-  ff.6x_clone.counts %>% filter(ff.Samples > 5)

ff.list <- list(
  n1 = ff.1$clonotype,
  n2 = ff.2$clonotype,
  n3 = ff.3$clonotype,
  n4 = ff.4$clonotype,
  n5 = ff.5$clonotype,
  n6 = ff.6$clonotype
)
names(ff.list) <- c("ff: 1", "ff: 2", "ff: 3", "ff: 4", "ff: 5", "ff: 6")
eu.ff <- euler(ff.list)
eu.ff_plot <- plot(eu.ff, fill_alpha=0.7,
                   quantities = TRUE,
                   fills = c("#f3e9db","#e7d4b8", "#dcbe94", "#d0a971", "#c5944e", "#9d763e")
)
ggsave('Figure Export/Euler Diagram/CD4 ff concentric euler diagram.pdf', eu.ff_plot,width=5,height=5)
# ----------------------------------------------------------------------------------------------------------------------------------------

# make a concentric euler diagram to illistrate how many clones are shared in the dibberent animals:
bb.5x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample) %>% head(n=5)) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bb') %>% pull(Sample))))) %>%
  filter(bb.Samples > 0) %>%
  replace(is.na(.), 0)

bb.1 <-  bb.5x_clone.counts %>% filter(bb.Samples > 0)
bb.2 <-  bb.5x_clone.counts %>% filter(bb.Samples > 1)
bb.3 <-  bb.5x_clone.counts %>% filter(bb.Samples > 2)
bb.4 <-  bb.5x_clone.counts %>% filter(bb.Samples > 3)
bb.5 <-  bb.5x_clone.counts %>% filter(bb.Samples > 4)

bb.list <- list(
  n1 = bb.1$clonotype,
  n2 = bb.2$clonotype,
  n3 = bb.3$clonotype,
  n4 = bb.4$clonotype,
  n5 = bb.5$clonotype
)
names(bb.list) <- c("bb: 1", "bb: 2", "bb: 3", "bb: 4", "bb: 5")
eu.bb <- euler(bb.list)
eu.bb_plot <- plot(eu.bb, fill_alpha=0.7,
                   quantities = TRUE,
                   fills = c("#dcecfe","#bad9fd", "#98c6fc", "#76b3fb", "#54a0fb")
)
ggsave('Figure Export/Euler Diagram/CD4 bb concentric euler diagram.pdf', eu.bb_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------

# make a concentric euler diagram to illistrate how many clones are shared in the dig7g7erent animals:
g7g7.6x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample) %>% head(n=6)) %>%
  dplyr::mutate(g7g7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'g7g7') %>% pull(Sample))))) %>%
  filter(g7g7.Samples > 0) %>%
  replace(is.na(.), 0)

g7g7.1 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 0)
g7g7.2 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 1)
g7g7.3 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 2)
g7g7.4 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 3)
g7g7.5 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 4)
g7g7.6 <-  g7g7.6x_clone.counts %>% filter(g7g7.Samples > 5)

g7g7.list <- list(
  n1 = g7g7.1$clonotype,
  n2 = g7g7.2$clonotype,
  n3 = g7g7.3$clonotype,
  n4 = g7g7.4$clonotype,
  n5 = g7g7.5$clonotype,
  n6 = g7g7.6$clonotype
)
names(g7g7.list) <- c("g7g7: 1", "g7g7: 2", "g7g7: 3", "g7g7: 4", "g7g7: 5", "g7g7: 6")
eu.g7g7 <- euler(g7g7.list)
eu.g7g7_plot <- plot(eu.g7g7, fill_alpha=0.7,
                     quantities = TRUE,
                     fills = c("#ebcfdb","#d89fb7", "#c56f93", "#b23f6f", "#9f104c", "#8f0e44")
)
ggsave('Figure Export/Euler Diagram/CD4 g7g7 concentric euler diagram.pdf', eu.g7g7_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------
# make a concentric euler diagram to illistrate how many clones are shared in the dibxferent animals:
bxf.6x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample) %>% head(n=6)) %>%
  dplyr::mutate(bxf.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxf') %>% pull(Sample))))) %>%
  filter(bxf.Samples > 0) %>%
  replace(is.na(.), 0)

bxf.1 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 0)
bxf.2 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 1)
bxf.3 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 2)
bxf.4 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 3)
bxf.5 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 4)
bxf.6 <-  bxf.6x_clone.counts %>% filter(bxf.Samples > 5)

bxf.list <- list(
  n1 = bxf.1$clonotype,
  n2 = bxf.2$clonotype,
  n3 = bxf.3$clonotype,
  n4 = bxf.4$clonotype,
  n5 = bxf.5$clonotype,
  n6 = bxf.6$clonotype
)
names(bxf.list) <- c("bxf: 1", "bxf: 2", "bxf: 3", "bxf: 4", "bxf: 5", "bxf: 6")
eu.bxf <- euler(bxf.list)
eu.bxf_plot <- plot(eu.bxf, fill_alpha=0.7,
                    quantities = TRUE,
                    fills = c("#e5e5cc","#cccc99", "#b2b266", "#999932", "#808000", "#737300")
)

ggsave('Figure Export/Euler Diagram/CD4 bxf concentric euler diagram.pdf', eu.bxf_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------

# make a concentric euler diagram to illistrate how many clones are shared in the dibxg7erent animals:
bxg7.6x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample) %>% head(n=5)) %>%
  dplyr::mutate(bxg7.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxg7') %>% pull(Sample))))) %>%
  filter(bxg7.Samples > 0) %>%
  replace(is.na(.), 0)

bxg7.1 <-  bxg7.6x_clone.counts %>% filter(bxg7.Samples > 0)
bxg7.2 <-  bxg7.6x_clone.counts %>% filter(bxg7.Samples > 1)
bxg7.3 <-  bxg7.6x_clone.counts %>% filter(bxg7.Samples > 2)
bxg7.4 <-  bxg7.6x_clone.counts %>% filter(bxg7.Samples > 3)
bxg7.5 <-  bxg7.6x_clone.counts %>% filter(bxg7.Samples > 4)

bxg7.list <- list(
  n1 = bxg7.1$clonotype,
  n2 = bxg7.2$clonotype,
  n3 = bxg7.3$clonotype,
  n4 = bxg7.4$clonotype,
  n5 = bxg7.5$clonotype
)
names(bxg7.list) <- c("bxg7: 1", "bxg7: 2", "bxg7: 3", "bxg7: 4", "bxg7: 5")
eu.bxg7 <- euler(bxg7.list)
eu.bxg7_plot <- plot(eu.bxg7, fill_alpha=0.7,
                     quantities = TRUE,
                     fills = c("#fddbf6","#fbb8ee", "#fa94e6", "#f871de", "#f74ed6")
)
ggsave('Figure Export/Euler Diagram/CD4 bxg7 concentric euler diagram.pdf', eu.bxg7_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------
# make a concentric euler diagram to illistrate how many clones are shared in the difxserent animals:
fxs.3x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'fxs') %>% pull(Sample) %>% head(n=3)) %>%
  dplyr::mutate(fxs.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'fxs') %>% pull(Sample))))) %>%
  filter(fxs.Samples > 0) %>%
  replace(is.na(.), 0)

fxs.1 <-  fxs.3x_clone.counts %>% filter(fxs.Samples > 0)
fxs.2 <-  fxs.3x_clone.counts %>% filter(fxs.Samples > 1)
fxs.3 <-  fxs.3x_clone.counts %>% filter(fxs.Samples > 2)

fxs.list <- list(
  n1 = fxs.1$clonotype,
  n2 = fxs.2$clonotype,
  n3 = fxs.3$clonotype
)
names(fxs.list) <- c("fxs: 1", "fxs: 2", "fxs: 3")
eu.fxs <- euler(fxs.list)
eu.fxs_plot <- plot(eu.fxs, fill_alpha=0.7,
                    quantities = TRUE,
                    fills = c("#fdbf99", "#fb8f4c", "#fa5f00")
)

ggsave('Figure Export/Euler Diagram/CD4 fxs concentric euler diagram.pdf', eu.fxs_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------
# make a concentric euler diagram to illistrate how many clones are shared in the dibxserent animals:
bxs.3x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxs') %>% pull(Sample) %>% head(n=3)) %>%
  dplyr::mutate(bxs.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'bxs') %>% pull(Sample))))) %>%
  filter(bxs.Samples > 0) %>%
  replace(is.na(.), 0)

bxs.1 <-  bxs.3x_clone.counts %>% filter(bxs.Samples > 0)
bxs.2 <-  bxs.3x_clone.counts %>% filter(bxs.Samples > 1)
bxs.3 <-  bxs.3x_clone.counts %>% filter(bxs.Samples > 2)

bxs.list <- list(
  n1 = bxs.1$clonotype,
  n2 = bxs.2$clonotype,
  n3 = bxs.3$clonotype
)
names(bxs.list) <- c("bxs: 1", "bxs: 2", "bxs: 3")
eu.bxs <- euler(bxs.list)
eu.bxs_plot <- plot(eu.bxs, fill_alpha=0.7,
                    quantities = TRUE,
                    fills = c("#a3fe99","#5efe4d", "#1afe02")
)

ggsave('Figure Export/Euler Diagram/CD4 bxs concentric euler diagram.pdf', eu.bxs_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------
# make a concentric euler diagram to illistrate how many clones are shared in the disserent animals:
ss.4x_clone.counts <- as_tibble(pr_CD4_merge) %>%
  dplyr::mutate(clonotype = paste0(CDR3.aa, "_", V.name, "_", J.name)) %>%
  select(clonotype, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample) %>% head(n=4)) %>%
  dplyr::mutate(ss.Samples = rowSums(!is.na(across(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == 'ss') %>% pull(Sample))))) %>%
  filter(ss.Samples > 0) %>%
  replace(is.na(.), 0)

ss.1 <-  ss.4x_clone.counts %>% filter(ss.Samples > 0)
ss.2 <-  ss.4x_clone.counts %>% filter(ss.Samples > 1)
ss.3 <-  ss.4x_clone.counts %>% filter(ss.Samples > 2)
ss.4 <-  ss.4x_clone.counts %>% filter(ss.Samples > 3)


ss.list <- list(
  n1 = ss.1$clonotype,
  n2 = ss.2$clonotype,
  n3 = ss.3$clonotype,
  n4 = ss.4$clonotype
)
names(ss.list) <- c("ss: 1", "ss: 2", "ss: 3", "ss: 4")
eu.ss <- euler(ss.list)
eu.ss_plot <- plot(eu.ss, fill_alpha=0.7,
                   quantities = TRUE,
                   fills = c("#feffe5", "#fdff99", "#fcff4c", "#fbff00")
)

ggsave('Figure Export/Euler Diagram/CD4 ss concentric euler diagram.pdf', eu.ss_plot,width=5,height=5)

# ----------------------------------------------------------------------------------------------------------------------------------------

# Select just the bb and ff parents for a simple diagram 
  # par 1
  Par1 <-  bb.ff_clone.counts %>%
    filter(.[[5]] > 0) %>% select(clonotype)
  # par 2
  Par2 <- bb.ff_clone.counts %>%
    filter(.[[7]] > 0) %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    Par2 = Par2$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c(names(bb.ff_clone.counts)[5], names(bb.ff_clone.counts)[7])
  map.name <- str_sub(map.name,1,nchar(map.name)-8)
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x)
  
  # generate the euler diagrams with desired colors
  euler_plot_bb.ff_parent <- plot(eu.x,
                     quantities = TRUE,
                     fills = c(bb_color, ff_color)
  )
  ggsave('Figure Export/Euler Diagram/CD4 parents bb vs ff.pdf', euler_plot_bb.ff_parent,width=3,height=3)
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  
  # Select just the bxf and tetra. bb&ff for a 2-way overlap diagram 
  # tetra. bb&ff
  tetra <-  bb.ff.tetra_clone.counts %>%
    filter(.[[6]] > 0) %>% select(clonotype)
  # bxf
  F1 <- bb.ff.tetra_clone.counts %>%
    filter(.[[8]] > 0) %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    tetra = tetra$clonotype,
    F1 = F1$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c(names(bb.ff.tetra_clone.counts)[6], names(bb.ff.tetra_clone.counts)[8])
  map.name <- str_sub(map.name,1,nchar(map.name)-8)
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x)
  
  # generate the euler diagrams with desired colors
  euler_plot_tetra.vs.bxf <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(tetra_bb.ff_color, bxf_color)
  )
  
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  
  # Select just the bxf and tetra. bb&ff for a 2-way overlap diagram 
  Par1 <-  bb.ff.tetra_clone.counts %>%
    filter(.[[5]] > 0) %>% select(clonotype)
  # tetra. bb&ff
  tetra <-  bb.ff.tetra_clone.counts %>%
    filter(.[[6]] > 0) %>% select(clonotype)
  # par 2
  Par2 <- bb.ff.tetra_clone.counts %>%
    filter(.[[7]] > 0) %>% select(clonotype)
  # bxf
  F1 <- bb.ff.tetra_clone.counts %>%
    filter(.[[8]] > 0) %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    tetra = tetra$clonotype,
    Par2 = Par2$clonotype,
    F1 = F1$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c(names(bb.ff.tetra_clone.counts)[5], names(bb.ff.tetra_clone.counts)[6], names(bb.ff.tetra_clone.counts)[7], names(bb.ff.tetra_clone.counts)[8])
  map.name <- str_sub(map.name,1,nchar(map.name)-8)
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  euler_plot_tetra.vs.bxf <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(bb_color,tetra_bb.ff_color,ff_color, bxf_color)
  )
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  
  # Select just the bxf and tetra. bb&ff for a 2-way overlap diagram 
  Par1 <-  bb.g7g7.tetra_clone.counts %>%
    filter(.[[5]] > 0) %>% select(clonotype)
  # tetra. bb&g7g7
  tetra <-  bb.g7g7.tetra_clone.counts %>%
    filter(.[[6]] > 0) %>% select(clonotype)
  # par 2
  Par2 <- bb.g7g7.tetra_clone.counts %>%
    filter(.[[7]] > 0) %>% select(clonotype)
  # bxg7
  F1 <- bb.g7g7.tetra_clone.counts %>%
    filter(.[[8]] > 0) %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    tetra = tetra$clonotype,
    Par2 = Par2$clonotype,
    F1 = F1$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c(names(bb.g7g7.tetra_clone.counts)[5], names(bb.g7g7.tetra_clone.counts)[6], names(bb.g7g7.tetra_clone.counts)[7], names(bb.g7g7.tetra_clone.counts)[8])
  map.name <- str_sub(map.name,1,nchar(map.name)-8)
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  euler_plot_tetra.vs.bxg7 <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(bb_color,tetra_bb.g7g7_color,g7g7_color, bxg7_color)
  )  
  
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  # Subset only the most consistent clones and see how the plots appear:
  
  # Select just the bxf and tetra. bb&ff for a 2-way overlap diagram 
  Par1 <-  bb.5 %>% select(clonotype)
  # par 2
  Par2 <- ff.5 %>% select(clonotype)
  #  All of bxf
  F1.All <- bxf.1 %>% select(clonotype)


  # top hits of bxf
  F1.top <- bxf.5 %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    Par2 = Par2$clonotype,
    F1.All = F1.All$clonotype,
    F1.top = F1.top$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c("bb 5x", "ff 5-6x", "bxf 1-6x", "bxf 5-6x")
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  parental.consistent_bxf <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(bb_color,ff_color, "#e5e5cc", bxf_color)
  )
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  # Subset only the most consistent clones and see how the plots appear:
  Par1 <-  bb.5 %>% select(clonotype)
  # par 2
  Par2 <- g7g7.5 %>% select(clonotype)
  #  All of bxg7
  F1.All <- bxg7.1 %>% select(clonotype)

  # top hits of bxg7
  F1.top <- bxg7.5 %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    Par2 = Par2$clonotype,
    F1.All = F1.All$clonotype,
    F1.top = F1.top$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c("bb 5x", "g7g7 5-6x", "bxg7 1-6x", "bxg7 5-6x")
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  parental.consistent_bxg7 <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(bb_color,g7g7_color, "#fddbf6", bxg7_color)
  )
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  # Subset only the most consistent clones and see how the plots appear:
  Par1 <-  bb.5 %>% select(clonotype)
  # par 2
  Par2 <- ss.4 %>% select(clonotype)
  #  All of bxs
  F1.All <- bxs.1 %>% select(clonotype)
  
  # top hits of bxs
  F1.top <- bxs.3 %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    Par2 = Par2$clonotype,
    F1.All = F1.All$clonotype,
    F1.top = F1.top$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c("bb 5x", "ss 4x", "bxs 1-3x", "bxs 3x")
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  parental.consistent_bxs <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(bb_color, ss_color, "#a3fe99", bxs_color)
  )
  
  # ----------------------------------------------------------------------------------------------------------------------------------------
  # Subset only the most consistent clones and see how the plots appear:
  Par1 <-  ff.5 %>% select(clonotype)
  # par 2
  Par2 <- ss.4 %>% select(clonotype)
  #  All of fxs
  F1.All <- fxs.1 %>% select(clonotype)
  
  # top hits of fxs
  F1.top <- fxs.3 %>% select(clonotype)
  
  # make a list of the clonotypes present in each particular haplotype
  x <- list(
    Par1 = Par1$clonotype,
    Par2 = Par2$clonotype,
    F1.All = F1.All$clonotype,
    F1.top = F1.top$clonotype
  )
  
  # Map the actual names of the clonotypes to the plots
  map.name = c("ff 5-6x", "ss 4x", "fxs 1-3x", "fxs 3x")
  names(x) <- map.name
  
  # Generate euler dataframe
  eu.x <- euler(x, shape = "ellipse")
  
  # generate the euler diagrams with desired colors
  parental.consistent_fxs <- plot(eu.x,
                                  quantities = TRUE,
                                  fills = c(ff_color, ss_color, "#fdbf99", fxs_color)
  )
  
ggsave('Figure Export/Euler Diagram/CD4 parental.consistent_bxf.pdf', parental.consistent_bxf,width=4,height=4)
ggsave('Figure Export/Euler Diagram/CD4 parental.consistent_bxg7.pdf', parental.consistent_bxg7,width=4,height=4)
ggsave('Figure Export/Euler Diagram/CD4 parental.consistent_bxs.pdf', parental.consistent_bxs,width=4,height=4)
ggsave('Figure Export/Euler Diagram/CD4 parental.consistent_fxs.pdf', parental.consistent_fxs,width=4,height=4)
