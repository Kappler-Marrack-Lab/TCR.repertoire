# Combine Parental Runs to generate a paired parental dataset which has been filtered.
library(immunarch)
library(tidyverse)


# load the filtered data objects
load("TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # loads objects: TCR_DoVb_CD4_AA.ones.removed.LFSR, CD4s_filtered



Combine_DF <- function(x){
  x2 <- x %>% 
    dplyr::select(Clones, Proportion, CDR3.aa, V.name, J.name) %>%
    mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
    group_by(CDR3.aa, V.name, J.name) %>%
    dplyr::summarise(Clones = sum(Clones)) %>% 
    as.data.frame()
  x3  <- x2 %>%
    dplyr::mutate(Proportion = Clones/sum(Clones)) %>%
    # filter(Proportion >= 5e-06) %>%
    add_column(CDR3.nt = NA) %>% add_column(D.name = NA) %>% add_column(V.end = NA) %>%
    add_column(D.start = NA) %>% add_column(D.end = NA) %>% add_column(J.start = NA) %>%
    add_column(VJ.ins = NA) %>% add_column(VD.ins = NA) %>% add_column(DJ.ins = NA) %>% add_column(Sequence = NA) %>%
    select(Clones, Proportion, CDR3.nt, CDR3.aa, V.name, D.name, J.name, V.end,  D.start, D.end,  J.start, VJ.ins, VD.ins, DJ.ins, Sequence)
  return(x3)
}
# ----------------------

# calculate the number of unique sequences when two parents are added together (use the data which has the filtering cutoff applied)
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

# Create a master table of all the relavant parental combinations I can make:
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

# Set up a for loop to combine the unique sequences in the desired two parental MHC runs.
my_list <- list() # define new empty list
for(i in 1:length(vec_A)){
  joint.DF <- Combine_DF(
    rbind(
      # Access list by name in vector index:
      TCR_DoVb_CD4_AA.ones.removed.LFSR[[1]][[vec_A[i]]],
      TCR_DoVb_CD4_AA.ones.removed.LFSR[[1]][[vec_B[i]]]
    )
  )
  
  my_list <- append(my_list, list(as.tibble(joint.DF))) # add this dataframe to the end of the list, etc.
}
names(my_list) <- PC.df$Sample # rename the list elements based on dataframe of names

# create a new metadata file for the list
paired.meta <- PC.df %>%
  select(Sample, MHC) %>%
  mutate(Tetraparental_Cell_Source = NA) %>%
  mutate(Cell_Number = NA) %>%
  mutate(`Co-receptor` = "CD4") %>%
  mutate(Mouse_Count = 2) %>%
  mutate(Concern = "Paired Parental Runs")

# Combine data and meta data in new list
TCR_DoVb_CD4_paired_AA.ones.removed.LFSR <- list(
  data = my_list,
  meta = paired.meta
)
# ------------------------------------------------------------------------------------------
# Save objects
save(TCR_DoVb_CD4_paired_AA.ones.removed.LFSR, 
     file = "Paired Parental and low frequency sequences removed.RData")


print("task finished.")
print("Data stored in:")
print("Paired Parental and low frequency sequences removed.RData")