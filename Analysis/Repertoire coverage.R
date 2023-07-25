# Recreate plots based on traditional ecological survey sampling techniques

# Note plots depicting the total number of sequences in figure 5 are only using 3x mice each and these plots use all mice per haplotype,
# sorted by the samples with the greates to least number of unique sequences.
# http://www.countrysideinfo.co.uk/3howto.htm#Line%20Transect%20Method
library(tidyverse)
library(ggsignif)
library(patchwork)
library(ggeasy)
library(reshape2)
library(immunarch)

# --------------------------------------------------------------------------------------------------
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed


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
`ss+ff` = "ss+ff"
`ss+bb` = "ss+bb"
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
ss.ff_color = "#05a0ff"
ss.bb_color = "#e501fd"
bb.g7g7_color = "#08b129"

DP = "DP"
DP_color = "#808080"
# --------------------------------------------------------------------------------------------------
# Make a function of the Repertoire coverage plot for each haplotype 

gg.systematic.sampling.plot <- function(rep.dat, query, MHC_color){
  
  MHC.Rep <- repFilter(rep.dat, .method = "by.meta", .query = list(MHC = include(query)))
  # MHC.Rep$meta
  
  # calculate the number of unique sequences obtained per run:
  CDR3AA.count <- function(df){
    df.out <- df %>%
      mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
      distinct(merge) %>%
      nrow()
    return(df.out)
  }
  clonotypes.AA_df <- sapply(MHC.Rep$data, CDR3AA.count)
  
  
  # Apply vector
  merge.VDJ <- function(df){
    df.out <- df %>%
      dplyr::mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
      select(merge)
    return(df.out)
  } # Creates a column of V-CDR3-J information
  MHC.Rep$data <- lapply(MHC.Rep$data, merge.VDJ) # apply this new column to all of the data.
  
  
  
  # obtain the order I want from the meta data:
  Sample.Order <- MHC.Rep$meta %>%
    mutate(count = clonotypes.AA_df) %>%
    arrange(desc(count)) %>%
    # arrange(desc(Mouse_Count)) %>%
    pull(Sample)
  
  
  
  # aggregate the data and change sample names to tallys
  df.depth <- bind_rows(MHC.Rep$data, .id = "column_label") %>% 
    as_tibble()
  df.depth <-  left_join(data.frame(column_label = Sample.Order),    # Reorder data frame
                         df.depth,
                         by = "column_label")
  
  df.depth$column_label <- factor(df.depth$column_label, levels=c(Sample.Order)) #set the speciifc order I want
  df.depth <- df.depth %>%  mutate(depth = unclass(column_label))
  
  # Cumulative counts of unique vales
  MHC.depth <- df.depth %>%
    #group_by(group)%>% # if you have a third variable and you want to achieve the same results for each group
    mutate(cum_unique_entries = cumsum(!duplicated(merge))) %>%
    group_by(depth) %>% # add group variable for more layers
    summarise(cum_unique_entries = last(cum_unique_entries))
  
  zeros <- data.frame(depth = 0, cum_unique_entries = 0)
  new_dat <- rbind(zeros, MHC.depth)
  # --------------------------------------------------------------------
  
  # # aggregate the data and change sample names to tallys
  # df.depth <- bind_rows(MHC.Rep$data, .id = "column_label") %>% 
  #   as_tibble() %>%
  #   dplyr::mutate(column_label = as.factor(column_label)) %>%
  #   mutate(depth = unclass(column_label))
  # 
  # # Cumulative counts of unique vales
  # MHC.depth <- df.depth %>%
  #   #group_by(group)%>% # if you have a third variable and you want to achieve the same results for each group
  #   mutate(cum_unique_entries = cumsum(!duplicated(merge))) %>%
  #   group_by(depth) %>% # add group variable for more layers
  #   summarise(cum_unique_entries = last(cum_unique_entries))
  # 
  # zeros <- data.frame(depth = 0, cum_unique_entries = 0)
  # new_dat <- rbind(zeros, MHC.depth)
  # 
  
  # plot the data
  gg.out <- ggplot(data = new_dat, aes(x = depth, y = cum_unique_entries))+
    geom_point(color = MHC_color, size = 2.5)+
      geom_line(color = MHC_color, linewidth = 1.25)+
    theme_bw()+
    scale_x_continuous(labels = as.character(new_dat$depth), breaks = new_dat$depth)+
    ylim(NA, 32000)+
    labs(x = "# Unique Sequencing\nRuns", y = "Cumulative Number of\nUnique Sequences")+
    ggtitle(paste(query, "Samples")) +
    theme(
      plot.title = element_text(size=12, face = "bold", color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      
      panel.background = element_rect(fill = "#e6e7e8",
                                      colour = "#e6e7e8",
                                      linewidth = 0.25, linetype = "dashed"),
      panel.grid.major = element_line(linewidth = 0.25, linetype = 'dashed',
                                      colour = "white"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
      
    )
    
  
  return(gg.out)
}

plot.bb <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, bb, bb_color)
plot.ff <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, ff, ff_color)
plot.g7g7 <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, g7g7, g7g7_color)
plot.bxf <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxf, bxf_color)
plot.bxg7 <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxg7, bxg7_color)

plot.IAb <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, IAbHet, IAbHet_color)
plot.ss <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, ss, ss_color)

plot.bxs <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxs, bxs_color)
plot.fxs <- gg.systematic.sampling.plot  (TCR_DoVb_CD4_AA.ones.removed.LFSR, fxs, fxs_color)

`plot.bb&ff` <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, tetra_bb.ff, tetra_bb.ff_color)
`plot.bb&g7g7` <- gg.systematic.sampling.plot(TCR_DoVb_CD4_AA.ones.removed.LFSR, tetra_bb.g7g7, tetra_bb.g7g7_color)


# # Also plot the DP samples:
# TCR_DoVb_DP_AA.ones.removed.LFSR$meta$MHC <- "DP"
# # remove sample 376 because it has a very low number of sequences
# TCR_DoVb_DP_AA.ones.removed.LFSR <- repFilter(TCR_DoVb_DP_AA.ones.removed.LFSR, .method = "by.meta",
#           .query = list(Sample = exclude("376 postJ")),
#           .match="substring")
# 

ggsave("Figure Export/Repertoire coverage/Repertoire coverage.bb.pdf", plot = plot.bb,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.ff.pdf", plot = plot.ff,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.g7g7.pdf", plot = plot.g7g7,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.bxf.pdf", plot = plot.bxf,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.bxg7.pdf", plot = plot.bxg7,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.IAb het.pdf", plot = plot.IAb,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.ss.pdf", plot = plot.ss,  width = 6, height = 6, units = "cm")
# ggsave("Figure Export/Repertoire coverage/Repertoire coverage.DP.pdf", plot = plot.DP,  width = 6, height = 6, units = "cm")


ggsave("Figure Export/Repertoire coverage/Repertoire coverage.bxs.pdf", plot = plot.bxs,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage.fxs.pdf", plot = plot.fxs,  width = 6, height = 6, units = "cm")

ggsave("Figure Export/Repertoire coverage/Repertoire coverage. tetra bb&ff.pdf", plot = `plot.bb&ff`,  width = 6, height = 6, units = "cm")
ggsave("Figure Export/Repertoire coverage/Repertoire coverage. tetra bb&g7g7.pdf", plot = `plot.bb&g7g7`,  width = 6, height = 6, units = "cm")

# # -------------------------------------------------------------------------------------------------------------------------------------------
# 
# gg.systematic.sampling.df <- function(rep.dat, query){
#   
#   MHC.Rep <- repFilter(rep.dat, .method = "by.meta", .query = list(MHC = include(query)))
#   MHC.Rep$meta
#   
#   # Apply vector
#   merge.VDJ <- function(df){
#     df.out <- df %>%
#       dplyr::mutate(merge = paste(V.name, CDR3.aa, J.name, sep = "_")) %>%
#       select(merge)
#     return(df.out)
#   } # Creates a column of V-CDR3-J information
#   MHC.Rep$data <- lapply(MHC.Rep$data, merge.VDJ) # apply this new column to all of the data.
#   
#   # aggregate the data and change sample names to tallys
#   df.depth <- bind_rows(MHC.Rep$data, .id = "column_label") %>% 
#     as_tibble() %>%
#     dplyr::mutate(column_label = as.factor(column_label)) %>%
#     mutate(depth = unclass(column_label))
#   
#   # Cumulative counts of unique vales
#   MHC.depth <- df.depth %>%
#     #group_by(group)%>% # if you have a third variable and you want to achieve the same results for each group
#     mutate(cum_unique_entries = cumsum(!duplicated(merge))) %>%
#     group_by(depth) %>% # add group variable for more layers
#     summarise(cum_unique_entries = last(cum_unique_entries))
#   
#   zeros <- data.frame(depth = 0, cum_unique_entries = 0)
#   new_dat <- rbind(zeros, MHC.depth)
#   new_dat$MHC <- query
#   return(new_dat)
# }
# ss.pooled <- rbind(
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, bb),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, ff),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, g7g7),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxf),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxg7),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, ss),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, bxs),
#   gg.systematic.sampling.df(TCR_DoVb_CD4_AA.ones.removed.LFSR, fxs)
# )
# 
# 
# 
# F.MHC_groups <- c(bb, ff, g7g7, ss, bxf, bxg7, bxs, fxs)
# C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color, bxf_color, bxg7_color, bxs_color, fxs_color)
# 
# ss.pooled.plot <- ggplot(data = ss.pooled, aes(x = depth, y = cum_unique_entries, color = MHC))+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#   scale_color_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
#   # scale_x_continuous(labels = as.character(new_dat$depth), breaks = new_dat$depth)+
#   labs(x = "# Unique Sequencing\nRuns", y = "Cumulative Number of\nUnique Sequences")
# 
# ggsave("Figure Export/Repertoire coverage/Repertoire coverage.pooled groups.pdf", plot = ss.pooled.plot,  width = 9, height = 6, units = "cm")
# 

  