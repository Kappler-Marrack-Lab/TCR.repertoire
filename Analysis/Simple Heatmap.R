library(stringr)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(scales)
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(data.table)  # faster fread() and better weekdays()
library(lubridate)   # date manipulation
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(ggplot2)
library(patchwork)
library(dendsort)
library(seriation)
library(circlize)
library(ComplexHeatmap)
library(corrplot)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(circlize)
library(colorspace)
library(GetoptLong)
library(devtools)
library(superheat)
library(stringr)
library(immunarch)
library(vroom)
library(parallel)
library(pbmcapply)
library(emdbook)
library(ggsignif)
set.seed(42)

#=================================
# Read in data
#=================================

load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed


# 
# TCR_DoVb_TR_AA.ones.removed.LFSR <- list(
#   data = c(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, TCR_DoVb_CD4_AA.ones.removed.LFSR$data, TCR_DoVb_DP_AA.ones.removed.LFSR$data),
#   meta = rbind(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta, TCR_DoVb_CD4_AA.ones.removed.LFSR$meta, TCR_DoVb_DP_AA.ones.removed.LFSR$meta)
# )
# 

# Remove the following TRAVs.
# They do no appear to applify anything:
# TRAV05-1A, TRAV09-2N, TRAV09-4N, TRAV15-2DN
filter_nonfunctional <- function(x){
  x <- x %>%
    filter((V.name %in% c("TRAV05-1A", "TRAV09-2N", "TRAV09-4N", "TRAV15-2DN")) == FALSE) %>%
    filter(Clones > 1)
  return(x)
}

TCR_DoVb_TR_AA.ones.removed.LFSR$data <- lapply(TCR_DoVb_TR_AA.ones.removed.LFSR$data, filter_nonfunctional)
TCR_DoVb_CD4_AA.ones.removed.LFSR$data <- lapply(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, filter_nonfunctional)
TCR_DoVb_CD8_AA.ones.removed.LFSR$data <- lapply(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, filter_nonfunctional)
TCR_DoVb_DP_AA.ones.removed.LFSR$data <- lapply(TCR_DoVb_DP_AA.ones.removed.LFSR$data, filter_nonfunctional)


# Reference Data.
TCR_DoVb_TR_AA.ones.removed.LFSR$meta
TCR_DoVb_CD4_AA.ones.removed.LFSR$meta
TCR_DoVb_CD8_AA.ones.removed.LFSR$meta
TCR_DoVb_DP_AA.ones.removed.LFSR$meta
# -----------------------------------------
# This re-assigns the databases lazily to add the tetraprental sequences I am interested in
# TCR_DoVb_CD4_AA.ones.removed.LFSR <- TCR_DoVb_CD4_AA.ones.removed.LFSR_TetraComb
# TCR_DoVb_CD8_AA.ones.removed.LFSR <- TCR_DoVb_CD8_AA.ones.removed.LFSR_TetraComb
# unique(TCR_DoVb_CD4_AA.ones.removed.LFSR$meta$MHC)
# TCR_DoVb_CD4_AA.ones.removed.LFSR$meta$Tetraparental_Cell_Source <- NA
# TCR_DoVb_CD8_AA.ones.removed.LFSR$meta$Tetraparental_Cell_Source <- NA
  
#Test to make sure clones are greater than 1
min(TCR_DoVb_CD4_AA.ones.removed.LFSR$data$`228 postJ`$Clones)
min(TCR_DoVb_TR_AA.ones.removed.LFSR$data$`228 postJ`$Clones)
# ---------------------------------------------------------------------------------------------------------
# Identify the size of the VJ matrix which needs to be made.
# See what is the number of unique TRAV and TRAJs which appear Save these as dimensions and use for a later matrix
V.vec <- mclapply(TCR_DoVb_TR_AA.ones.removed.LFSR$data, function(V) unique(V[["V.name"]]), mc.cores = 4)
J.vec <- mclapply(TCR_DoVb_TR_AA.ones.removed.LFSR$data, function(J) unique(J[["J.name"]]), mc.cores = 4)

V.vec <- as_tibble(unique(unlist(V.vec)))
names(V.vec) <- "V.name"
V.vec <- V.vec %>% arrange(V.name) %>% filter(V.name != "Vfam-Vsub")
J.vec <- as_tibble(unique(unlist(J.vec)))
names(J.vec) <- "J.name"
J.vec <- J.vec %>% arrange(J.name) %>% filter(J.name != "J")

V = as.character(V.vec$V.name)
J = as.character(J.vec$J.name)
# ---------------------------------------------------------------------------------------------------------
# Define Functions:
VJheat <- function(immunarch.list, V.vec, J.vec){
  immnames <- names(immunarch.list$data) # First section makes the sample dataframe.
  nc <- ncol(immunarch.list$data[[1]])
  spectraCDR3 <- data.frame(matrix(ncol = nc, nrow = 0))
  n <- names(immunarch.list$data[[1]])
  colnames(spectraCDR3) <- n
  
  for(i in 1:length(immnames)){
    print(paste0("Extracting sample and metadata from: ", immnames[i]))
    spectraCDR3_temp <- immunarch.list$data[[i]]
    spectraCDR3_temp$Sample <- immnames[i]
    spectraCDR3 <- rbind(spectraCDR3, spectraCDR3_temp)
  }
  spectraCDR3 <- merge(immunarch.list$meta, spectraCDR3, by="Sample",all=TRUE)
  
  # MHC Choice
  VJheat_min <- as_tibble(data.frame(
    "V.name" = spectraCDR3$V.name, 
    "J.name" = spectraCDR3$J.name,
    "MHC" = spectraCDR3$MHC, #This is equal to saying "MHC" = spectraCDR3$MHC, #"MHC" = spectraCDR3[Meta.Factor]
    "Clones" = spectraCDR3$Clones,
    stringsAsFactors = TRUE))
  #  I think the clones with matching TRAV, TRAJ need to be collapsed to remove duplicate TRAV / TRAJ pairs
  VJheat_min <- VJheat_min  %>%
    dplyr::group_by(V.name, J.name, MHC)  %>%
    dplyr::mutate(Clones = base::sum(Clones)) %>%
    unique() %>% ungroup() %>%
    dplyr::arrange(MHC, V.name, J.name)
  VJheat_expand <- VJheat_min %>% expand(MHC, V.name, J.name)
  VJheat_full <- VJheat_min %>% dplyr::right_join(VJheat_expand)
  VJheat_full$Clones <- VJheat_full$Clones %>% replace_na(0)
  
  # This expands the heatmap data to include missing VJ values and ensures all combinations are included.
  # Requires V and J vectors be included as input variables as these are dependent on the assembly software.
  # Pulls in the specific MHC variables from the input data
  MHC_vec = as.character(unique(VJheat_full$MHC))
  V = as.character(V.vec$V.name)
  J = as.character(J.vec$J.name)
  VJdim <- as_tibble(expand.grid(V.name = V, J.name = J, MHC = MHC_vec))
  
  VJheat_full <- merge(VJheat_full, VJdim, by = c("V.name", "J.name", "MHC"), all = TRUE) #Generates NA values where there was previously not a column
  VJheat_full$Clones <- VJheat_full$Clones %>% replace_na(0) # convert NA to 0 so it can be plotted
  # Remove artifact rows which somehow show up in some heatmaps
  VJheat_full <- VJheat_full %>% arrange(V.name) %>% filter(V.name != "Vfam-Vsub")
  VJheat_full <- VJheat_full %>% arrange(J.name) %>% filter(J.name != "J")
  
  # set ordering
  #set vector of levels you want
  Vlevels <- V.vec$V.name
  Jlevels <- J.vec$J.name
  #reorder factors
  VJheat_full$V.name <- factor(VJheat_full$V.name,levels=Vlevels)
  VJheat_full$J.name <- factor(VJheat_full$J.name, levels=Jlevels)
  
  return(VJheat_full) # this is the end of the function and what the function returns to me.
} #generates the source data for the heatmap
VJheat.norm <- function(VJheatmap){
  VJheatmap_norm <- VJheatmap %>%
    dplyr::group_by(MHC) %>%
    mutate(clones_norm = Clones/sum(Clones)) %>%
    select(-Clones) %>%
    # rename(Clones = clones_norm) # this was working but now its not...
    mutate(Clones = clones_norm)  %>% select(-clones_norm)
  return(VJheatmap_norm)
} # Option to pass a normalized heatmap
makeHeatMap <- function(VJheat_full){
  heatscale.min <- min(VJheat_full$Clones) # 0 is min
  heatscale.max <- max(VJheat_full$Clones) # 463643 is max.
  
  #this is for automated equal Log-value range
  Max_val <- sort(unique(VJheat_full$Clones), TRUE)[1] # first highest value
  Min_val <- sort(unique(VJheat_full$Clones), FALSE)[1] # first lowest value
  second.lowest <- sort(unique(VJheat_full$Clones), FALSE)[2] # second lowest value
  value_range <- lseq(from=second.lowest, to=Max_val, length.out=10) # change lseq() to seq() to make this a linear scale.
  value_range <- append(0, value_range)
  value_range <- signif(value_range, 3)
  
  hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
  # make the heatmap for the specified groups in a given dataset. 
  HM <- ggplot(VJheat_full, aes(x = J.name, y = V.name, fill = Clones)) +
    geom_tile(color="white", size=0.1) +
    scale_fill_gradientn(colours = hm.palette(100),
                         na.value = "black",
                         name = "Freq. TCRα\n Sequences",
                         trans = "log",
                         breaks = value_range,
                         labels = value_range,
                         
                         # guide="legend"                       
    ) +
    coord_equal() +
    facet_wrap(~VJheat_full[[3]], ncol=4) + #change this back to ncol=7
    labs(x=NULL, y=NULL, title="TRAV / TRAJ Heatmap")  +
    theme_tufte(base_family="Helvetica") +
    theme(axis.ticks=element_blank()) +
    theme(axis.text=element_text(size=5)) +
    theme(panel.border=element_blank()) +
    theme(plot.title=element_text(hjust=0)) +
    theme(strip.text=element_text(hjust=0)) +
    theme(panel.spacing.x=unit(0.5, "cm")) +
    theme(panel.spacing.y=unit(0.5, "cm")) +
    theme(legend.title=element_text(size=6)) +
    theme(legend.title.align=1) +
    theme(legend.text=element_text(size=6)) +
    theme(legend.position="bottom") +
    guides(fill = guide_colourbar(barwidth = 20, barheight = 1)) +
    theme(legend.key.size=unit(0.2, "cm")) +
    theme(legend.key.width=unit(1, "cm")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(HM)
} # Makes ggplot based heatmaps
# ---------------------------------------------------------------------------------------------------------
# Generate Heatmaps
VJheat_CD4 <- VJheat(TCR_DoVb_CD4_AA.ones.removed.LFSR, V.vec, J.vec)
VJheat_CD8 <- VJheat(TCR_DoVb_CD8_AA.ones.removed.LFSR, V.vec, J.vec)
VJheat_DP <- VJheat(TCR_DoVb_DP_AA.ones.removed.LFSR, V.vec, J.vec)

# Modify the names of the DP MHC so that we colapse all the DP runs together.
TCR.DP.combine <- TCR_DoVb_DP_AA.ones.removed.LFSR
TCR.DP.combine$meta$MHC <- "DP"
VJheat_DP.combine <- VJheat(TCR.DP.combine, V.vec, J.vec)
# ---------------------------------------------------------------------------------------------------------
# Generate ggplots
VJheat_CD4_ggplot <- makeHeatMap(VJheat_CD4)
VJheat_CD8_ggplot <- makeHeatMap(VJheat_CD8)
VJheat_DP_ggplot <- makeHeatMap(VJheat_DP)
VJheat_DP_ggplot <- makeHeatMap(VJheat_DP)
makeHeatMap(VJheat_DP.combine)

pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC DP.pdf", width = 8, height = 7)
VJheat_DP_ggplot
dev.off()
pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC CD4.pdf", width = 20, height = 15)
VJheat_CD4_ggplot
dev.off()
pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC CD8.pdf", width = 20, height = 15)
VJheat_CD8_ggplot
dev.off()
# ---------------------------------------------------------------------------------------------------------
# Merge of all data, discriminated by co-receptor. This is for making a single plot where everything can be compaired together.
VJheat_DP$Co_receptor = "DP"
VJheat_DP.combine$Co_receptor = "DP"
VJheat_CD4$Co_receptor = "CD4"
VJheat_CD8$Co_receptor = "CD8"

VJheat_DP <- within(VJheat_DP,  MHC <- paste(MHC, Co_receptor, sep = " "))
VJheat_CD4 <- within(VJheat_CD4,  MHC <- paste(MHC, Co_receptor, sep = " "))
VJheat_CD8 <- within(VJheat_CD8,  MHC <- paste(MHC, Co_receptor, sep = " "))
VJheat_all_data <- do.call("rbind", list(VJheat_DP.combine, VJheat_CD4, VJheat_CD8))
# ---------------------------------------------------------------------------------------------------------
# Normalize the merged heatmap data and plot it.
VJheat_all_data.norm <- VJheat.norm(VJheat_all_data)
VJheat_all_data.norm_ggplot <- makeHeatMap(VJheat_all_data.norm)
pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC All Data - Normalzied.pdf", width = 20, height = 20)
VJheat_all_data.norm_ggplot
dev.off()

# pull in relevant CD4 and DP data and have it organized
VJheat_DP_CD4_data <- do.call("rbind", list(VJheat_DP.combine, VJheat_CD4))
VJheat_DP_CD4_data.norm <- VJheat.norm(VJheat_DP_CD4_data)
unique(VJheat_DP_CD4_data$MHC)
MHC.order <-  c("bb CD4",  "ff CD4", "g7g7 CD4",  "ss CD4",
"bxf CD4","fxs CD4", "bxg7 CD4", "bxs CD4",
"DP")
VJheat_DP_CD4_data.norm <- VJheat_DP_CD4_data.norm %>% filter(MHC %in% MHC.order) %>% arrange(match(MHC, MHC.order))

# Find second lowest value in the dataset:
VJfreqs <- VJheat_DP_CD4_data.norm$Clones
min(VJfreqs[VJfreqs != min(VJfreqs)])
max(VJfreqs)


pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC relevant CD4 DP - Normalzied.pdf", width = 20, height = 20)
makeHeatMap(VJheat_DP_CD4_data.norm) + facet_wrap(~factor(MHC, levels=MHC.order), ncol = 4)
dev.off()
# # ---------------------------------------------------------------------------------------------------------
# # Test the DP thymus data using all previous runs
# # Load in this data.
# # file_path = "/Volumes/SSD/Analysis/DP only" # This is the specific location of the data on the lab computer.
# # TCR.DP_DoVb_TR_all <- repLoad(file_path)
# # heatmap_thymus <- VJheat(TCR.DP_DoVb_TR_all, V.vec, J.vec)
# # heatmap_thymus.norm <- VJheat.norm(heatmap_thymus)
# # heatmap_thymus.norm.ggplot <- makeHeatMap(heatmap_thymus.norm)
# # pdf("Figure Export/VJ Usage Simple Heatmap/Va_Heatmap by MHC All Thymus  Data - Normalzied.pdf", width = 10, height = 10)
# # heatmap_thymus.norm.ggplot
# # dev.off()
# 
# # ---------------------------------------------------------------------------------------------------------
# # Create a scale for heatmaps which is log based and automatic
# # WORKING
# # Max_val <- sort(unique(VJheat_all_data.norm$Clones), TRUE)[1] # first highest value
# # Min_val <- sort(unique(VJheat_all_data.norm$Clones), FALSE)[1] # first lowest value
# # second.lowest <- sort(unique(VJheat_all_data.norm$Clones), FALSE)[2] # second lowest value
# # value_range <- lseq(from=second.lowest, to=Max_val, length.out=10)
# # value_range <- append(0, value_range)
# # value_range <- signif(value_range, 3)
# # ---------------------------------------------------------------------------------------------------------
# ff ="ff"
# ss = "ss"
# bb = "bb"
# g7g7 = "g7g7"
# bxf = "bxf"
# fxs = "fxs"
# bxs = "bxs"
# bxg7 = "bxg7"
# IAbHet = "b+/-"
# b.f = "b+f"
# b.f_Simulated = "b+fSimulated"
# b.g7 = "b+g7"
# b.g7_Simulated = "b+g7Simulated"
# tetra_bb.ff = "tetraparental bb&ff"
# tetra_bb.g7g7 = "tetraparental bb&g7g7"
# 
# # Define MHC group colors:
# ff_color ="#c5944e"
# ss_color = "#fbff00"
# bb_color = "#54a0fb" 
# g7g7_color = "#9f104c"
# bxf_color = "#808000"
# fxs_color = "#fa5f00"
# bxs_color = "#1afe02"
# bxg7_color = "#f74ed6"
# IAbHet_color = "#bddbff"
# b.f_color = "#EEC0C4"
# b.f_Simulated_color = "#ad2b37"
# b.g7_color = "#335566"
# b.g7_Simulated_color = "#17d2e3"
# tetra_bb.ff_color = "#EEC0C4"
# tetra_bb.g7g7_color = "#335566"
# 
# CD4 = 'CD4'
# CD8 = 'CD8'
# DP = 'DP'
# # ---------------------------------------------------------------------------------------------------------
# 
# # ------- immuneArch testing--------------------------------------------------------------------------------------------------
# # Find the number of zeros. Holes in the repertoire.
# 
# # dimensions on the VJ plots = 81 Vs * 43 Js = 3,483
# Tens.and.Zeros <- function(VJheat){
#   zeros.A <- VJheat[which(VJheat$Clones == 0),]$Clones
#   ones.A <- VJheat[which(VJheat$Clones == 1),]$Clones
#   tens.A <- VJheat[which(VJheat$Clones < 11),]$Clones
#   hundred.A <- VJheat[which(VJheat$Clones < 101),]$Clones
#   thousand.A <- VJheat[which(VJheat$Clones < 1001),]$Clones
#   tenthousand.A <- VJheat[which(VJheat$Clones < 10001),]$Clones
# 
#   print(paste0("Zeros:  ", length(zeros.A) ))
#   print(paste0("Ones:  ", length(ones.A) ))
#   print(paste0("10 and Lower:  ", length(tens.A) ))
#   print(paste0("100 and Lower:  ", length(hundred.A) ))
#   print(paste0("1,000 and Lower:  ", length(thousand.A) ))
#   print(paste0("10,000 and Lower:  ", length(tenthousand.A) ))
#   
# }
# 
# # ------------------
# 
# bb_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bb CD4")
# bb_CD4_data_gg <- makeHeatMap(bb_CD4_data)
# bb_CD4_data_gg
# 
# ff_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "ff CD4")
# bxf_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bxf CD4")
# b.f_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "b+f CD4")
# g7g7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "g7g7 CD4")
# bxg7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bxg7 CD4")
# b.g7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "b+g7 CD4")
# 
# bb_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bb CD8")
# ff_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "ff CD8")
# bxf_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bxf CD8")
# b.f_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "b+f CD8")
# g7g7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "g7g7 CD8")
# bxg7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bxg7 CD8")
# b.g7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "b+g7 CD8")
# # ----------------------
# 
# 
# 
# Tens.and.Zeros(bb_CD4_data)
# Tens.and.Zeros(ff_CD4_data)
# Tens.and.Zeros(b.f_CD4_data)
# Tens.and.Zeros(bxf_CD4_data)
# 
# imm_rare <- repClonality(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, .method = "rare")
# imm_top <- repClonality(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
# vis(imm_rare) + vis(imm_rare, .by = "MHC", .meta = TCR_DoVb_CD8_AA.ones.removed.LFSR$meta)
# vis(imm_top) + vis(imm_top, .by = "MHC", .meta = TCR_DoVb_CD8_AA.ones.removed.LFSR$meta)
# 
# imm_gu <- geneUsage(TCR_DoVb_CD8_AA.ones.removed.LFSR$data, "musmus.trav", .norm = T)
# vis(imm_gu, .by = "MHC", .meta = TCR_DoVb_CD8_AA.ones.removed.LFSR$meta, .plot = "box")
# # ---------------------------------------------------------------------------------------------------------
# # Generates box plots of the frequecy of zeros, tens and hundreds extracted from VJ pairings in the heatmaps.
# # Assess for each sample.
# # require(dplyr)
# sigFunc = function(x){
#   if(x < 0.001){"***"} 
#   else if(x < 0.01){"**"}
#   else if(x < 0.05){"*"}
#   else{NA}
# }
# 
# VJheat_zeros <- function(immunarch.list, V.vec, J.vec){
#   immnames <- names(immunarch.list$data) # First section makes the sample dataframe.
#   nc <- ncol(immunarch.list$data[[1]])
#   spectraCDR3 <- data.frame(matrix(ncol = nc, nrow = 0))
#   n <- names(immunarch.list$data[[1]])
#   colnames(spectraCDR3) <- n
#   
#   for(i in 1:length(immnames)){
#     print(paste0("Extracting sample and metadata from: ", immnames[i]))
#     spectraCDR3_temp <- immunarch.list$data[[i]]
#     spectraCDR3_temp$Sample <- immnames[i]
#     spectraCDR3 <- rbind(spectraCDR3, spectraCDR3_temp)
#   }
#   spectraCDR3 <- merge(immunarch.list$meta, spectraCDR3, by="Sample",all=TRUE)
#   
#   # MHC Choice
#   VJheat_min <- as_tibble(data.frame(
#     "V.name" = spectraCDR3$V.name, 
#     "J.name" = spectraCDR3$J.name,
#     "MHC" = spectraCDR3$MHC, #This is equal to saying "MHC" = spectraCDR3$MHC, #"MHC" = spectraCDR3[Meta.Factor]
#     "Clones" = spectraCDR3$Clones,
#     "Sample" = spectraCDR3$Sample,
#     stringsAsFactors = TRUE))
#   # VJheat_min <- spread(VJheat_min, MHC, Sample)
#   
#   #  I think the clones with matching TRAV, TRAJ need to be collapsed to remove duplicate TRAV / TRAJ pairs
#   VJheat_min <- VJheat_min  %>%
#     dplyr::group_by(V.name, J.name, MHC, Sample)  %>%
#     dplyr::mutate(Clones = base::sum(Clones)) %>%
#     unique() %>% ungroup() %>%
#     dplyr::arrange(MHC, V.name, J.name, Sample)
#   
#   col_order <- c("V.name", "J.name", "MHC",
#                  "Sample", "Clones")
#   VJheat_min <- VJheat_min[, col_order]
#   VJheat_min <- VJheat_min %>% arrange(V.name, J.name, MHC, Sample)
#   
#   VJheat_min$MHC = paste(VJheat_min$MHC,"|",VJheat_min$Sample)
#   VJheat_min <- VJheat_min %>%  select(-Sample) # this is a trick to eliminate duplicates
#   
#   MHC_vec = as.character(unique(VJheat_min$MHC))
#   V = as.character(V.vec$V.name)
#   J = as.character(J.vec$J.name)
#   VJdim <- as_tibble(expand.grid(V.name = V, J.name = J, MHC = MHC_vec))
#   VJdim <- VJdim %>% arrange(V.name, J.name, MHC)
#   
#   VJheat_max <- merge(VJheat_min, VJdim, by = c("V.name", "J.name", "MHC"), all = TRUE) #Generates NA values where there was previously not a column
#   VJheat_max$Clones <- VJheat_max$Clones %>% replace_na(0) # convert NA to 0 so it can be plotted
#   VJheat_max <- VJheat_max %>% separate(MHC, c("MHC", "Sample") , sep = "([|])")
#   # Remove artifact rows which somehow show up in some heatmaps
#   VJheat_max <- VJheat_max %>% arrange(V.name, J.name, MHC, Sample) %>% filter(V.name != "Vfam-Vsub")
#   VJheat_max <- VJheat_max %>% arrange(V.name, J.name, MHC, Sample) %>% filter(J.name != "J")
# 
#   VJheat_max_wide <- VJheat_max %>%
#     dplyr::group_by(MHC, Sample) %>%
#     dplyr::summarise(Zero = sum(Clones == 0),
#               Ten = sum(Clones < 11),
#               Hundred = sum(Clones < 101),
#               Thousand = sum(Clones < 1001))
# 
#   VJheat_max_long <- gather(VJheat_max_wide, Clone_Appearance, Count, Zero:Thousand, factor_key=TRUE)
#   VJheat_max_long$MHC <- gsub(" ", "", VJheat_max_long$MHC)
#   return(VJheat_max_long) # this is the end of the function and what the function returns to me.
# } #generates the source data for the heatmap
# zeros.box_plot <- function(VJheat_Z, co_receptor, F.MHC_groups, C.MHC_groups){
#   # # Join columns for the dot data
#   # TCR.point[TCR.point == "N/A"] <- NA
#   # TCR.point <- TCR.point %>%
#   #   unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
#   # TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
#   # 
#   
#   plot <- ggplot(subset(VJheat_Z, MHC %in% F.MHC_groups), aes(x=Clone_Appearance, y=Count, fill=MHC))+
#     geom_boxplot()+
#     # set Custom Colors
#     scale_fill_manual(breaks = F.MHC_groups, 
#                       values=C.MHC_groups) +
#     scale_color_manual(breaks = F.MHC_groups, 
#                        values=C.MHC_groups) +
#     labs(title="TCRα Relative VJ pairing abundance", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="Group", y="Abundance of VJ Pairs") +
#     theme_bw(base_size = 18) +
#     theme(legend.position = "bottom", 
#           legend.title=element_blank(),
#           text=element_text(size=18))
#   return(plot)
# }
# 
# zeros.box_plot_stats <- function(VJheat_Z, co_receptor, F.MHC_groups, C.MHC_groups){
#   # # Join columns for the dot data
#   # TCR.point[TCR.point == "N/A"] <- NA
#   # TCR.point <- TCR.point %>%
#   #   unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
#   # TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
#   # 
#   MHC_list <- as.list(as.data.frame(combn(F.MHC_groups,2))) # converts MHC groups into the correct comparisons
#   filtered.df <- VJheat_Z %>%
#     filter(MHC %in% F.MHC_groups)
#   
#   plot <- ggplot(data=filtered.df, aes(x=MHC, y=Count, fill=MHC))+
#     # geom_boxplot()+
#     geom_bar(stat="summary")+
#     geom_point(alpha = 0.1 )+
#     stat_summary(fun.data="mean_se")+
#     geom_signif(comparisons = MHC_list,
#                 test = 't.test', map_signif_level=sigFunc, step_increase = 0.1, tip_length = 0
#     ) +
#     facet_wrap(~Clone_Appearance) +
#     # set Custom Colors
#     scale_fill_manual(breaks = F.MHC_groups, 
#                       values=C.MHC_groups) +
#     scale_color_manual(breaks = F.MHC_groups, 
#                        values=C.MHC_groups) +
#     scale_y_continuous(limits = c(0,5000), expand = expansion(mult = c(0, .1))) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none"
#     ) +
#     labs(title="TCRα Relative VJ pairing abundance", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="Group", y="Abundance of Unselected\n TRAV / TRAJ Pairs")
# 
# 
#   return(plot)
# }
# 
# zeros.box_plot_stats_drop <- function(VJheat_Z, co_receptor, F.MHC_groups, C.MHC_groups, F.MHC_groups_drop){
#   # # Join columns for the dot data
#   # TCR.point[TCR.point == "N/A"] <- NA
#   # TCR.point <- TCR.point %>%
#   #   unite(MHC, c("MHC", "Tetraparental_Cell_Source"), sep = "_")
#   # TCR.point$MHC <- gsub("_NA.*", "\\1", TCR.point$MHC)
#   # 
#   MHC_list <- as.list(as.data.frame(combn(F.MHC_groups_drop,2))) # converts MHC groups into the correct comparisons
#   filtered.df <- VJheat_Z %>%
#     filter(MHC %in% F.MHC_groups)
#   
#   plot <- ggplot(data=filtered.df, aes(x=MHC, y=Count, fill=MHC))+
#     # geom_boxplot()+
#     geom_bar(stat="summary")+
#     geom_point(alpha = 0.1 )+
#     stat_summary(fun.data="mean_se")+
#     geom_signif(comparisons = MHC_list,
#                 test = 't.test', map_signif_level=TRUE, step_increase = 0.1
#     ) +
#     facet_wrap(~Clone_Appearance) +
#     # set Custom Colors
#     scale_fill_manual(breaks = F.MHC_groups, 
#                       values=C.MHC_groups) +
#     scale_color_manual(breaks = F.MHC_groups, 
#                        values=C.MHC_groups) +
#     labs(title="TCRα Relative VJ pairing abundance", subtitle=paste0(c(co_receptor, " Fixed DOβ T cells with haplotypes: ", F.MHC_groups), collapse=" "), x="Group", y="Abundance of VJ Pairs") +
#     theme_bw(base_size = 18) +
#     theme(legend.position = "bottom", 
#           legend.title=element_blank(),
#           text=element_text(size=18))
#   return(plot)
# }
# 
# CD4_z <- VJheat_zeros(TCR_DoVb_CD4_AA.ones.removed.LFSR, V.vec, J.vec)
# CD8_z <- VJheat_zeros(TCR_DoVb_CD8_AA.ones.removed.LFSR, V.vec, J.vec)
# 
# # load in the merged tetraparental files and the simulated tetraparental files.
# # These are in the "DOvbeta_repData.RData"
# TCR.CD4_Tet_merge_DoVb_TR$meta$Simulated_Status[TCR.CD4_Tet_merge_DoVb_TR$meta$Simulated_Status == "N/A"] <- NA
# TCR.CD8_Tet_merge_DoVb_TR$meta$Simulated_Status[TCR.CD8_Tet_merge_DoVb_TR$meta$Simulated_Status == "N/A"] <- NA
# TCR.CD4_Tet_merge_DoVb_TR$meta$MHC = paste0(TCR.CD4_Tet_merge_DoVb_TR$meta$MHC, " ", TCR.CD4_Tet_merge_DoVb_TR$meta$Simulated_Status)
# TCR.CD8_Tet_merge_DoVb_TR$meta$MHC = paste0(TCR.CD8_Tet_merge_DoVb_TR$meta$MHC, " ", TCR.CD8_Tet_merge_DoVb_TR$meta$Simulated_Status)
# TCR.CD4_Tet_merge_DoVb_TR$meta$MHC <- str_replace_all(TCR.CD4_Tet_merge_DoVb_TR$meta$MHC, " NA", "")
# TCR.CD8_Tet_merge_DoVb_TR$meta$MHC <- str_replace_all(TCR.CD8_Tet_merge_DoVb_TR$meta$MHC, " NA", "")
# 
# CD4_z_tetra <- VJheat_zeros(TCR.CD4_Tet_merge_DoVb_TR, V.vec, J.vec)
# CD8_z_tetra <- VJheat_zeros(TCR.CD8_Tet_merge_DoVb_TR, V.vec, J.vec)
# 
# 
# # Remove the b+f and b+g7 from the original data
# CD4_z_merged_tet_removed <- CD4_z %>% dplyr::group_by(Clone_Appearance) %>% 
#   filter((MHC %in% c("b+f", "b+g7")) == FALSE)
# CD8_z_merged_tet_removed <- CD8_z %>% dplyr::group_by(Clone_Appearance) %>% 
#   filter((MHC %in% c("b+f", "b+g7")) == FALSE)
# # bind the original data with the "merged tetraparental + simulated tetraparental groups"
# CD4_z_merged <- rbind(CD4_z_tetra, CD4_z_merged_tet_removed)
# CD4_z_merged <- CD4_z_merged %>% dplyr::group_by(Clone_Appearance)
# 
# CD8_z_merged <- rbind(CD8_z_tetra, CD8_z_merged_tet_removed)
# CD8_z_merged <- CD8_z_merged %>% dplyr::group_by(Clone_Appearance)
# 
# # Gather only the Zeros data (V1)
# CD4_z_merged_gaps <- CD4_z_merged %>% filter(Clone_Appearance == "Zero")
# CD8_z_merged_gaps <- CD8_z_merged %>% filter(Clone_Appearance == "Zero")
# 
# # generate plots (V1)
# F.MHC_groups <- c(bb, ff, bxf, b.f, b.f_Simulated)
# C.MHC_groups <- c(bb_color, ff_color, bxf_color, b.f_color, b.f_Simulated_color)
# CD4_z_plot_b.f_tetra.sim <- zeros.box_plot_stats(CD4_z_merged_gaps, CD4, F.MHC_groups, C.MHC_groups) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# CD8_z_plot_b.f_tetra.sim <- zeros.box_plot_stats(CD8_z_merged_gaps, CD8, F.MHC_groups, C.MHC_groups) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#  
# F.MHC_groups <- c(bb, g7g7, bxg7, b.g7, b.g7_Simulated)
# C.MHC_groups <- c(bb_color, g7g7_color, bxg7_color, b.g7_color, b.g7_Simulated_color)
# F.MHC_groups_drop <- c(bb, g7g7, bxg7, b.g7_Simulated) # use this function variable for when the error bars cannot be made for the comparisons you want because of too few data points.
# CD4_z_plot_b.g7_tetra.sim <- zeros.box_plot_stats_drop(CD4_z_merged_gaps, CD4, F.MHC_groups, C.MHC_groups, F.MHC_groups_drop) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# CD8_z_plot_b.g7_tetra.sim <- zeros.box_plot_stats_drop(CD8_z_merged_gaps, CD8, F.MHC_groups, C.MHC_groups, F.MHC_groups_drop) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # Gather only the Zeros data (V2)
# CD4_z_merged_gaps <- CD4_z %>% filter(Clone_Appearance == "Zero")
# CD8_z_merged_gaps <- CD8_z %>% filter(Clone_Appearance == "Zero")
# CD4_z_merged_gaps$MHC <- gsub('tetraparental','tetraparental ', CD4_z_merged_gaps$MHC)
# CD8_z_merged_gaps$MHC <- gsub('tetraparental','tetraparental ', CD8_z_merged_gaps$MHC)
# 
# # generate plots (V2)
# F.MHC_groups <- c(bb, ff, bxf, tetra_bb.ff)
# C.MHC_groups <- c(bb_color, ff_color, bxf_color, tetra_bb.ff_color)
# CD4_z_plot_b.f_tetra.sim <- zeros.box_plot_stats(CD4_z_merged_gaps, CD4, F.MHC_groups, C.MHC_groups) 
# CD8_z_plot_b.f_tetra.sim <- zeros.box_plot_stats(CD8_z_merged_gaps, CD8, F.MHC_groups, C.MHC_groups)
# 
# F.MHC_groups <- c(bb, g7g7, bxg7, tetra_bb.g7g7)
# C.MHC_groups <- c(bb_color, g7g7_color, bxg7_color, tetra_bb.g7g7_color)
# CD4_z_plot_b.g7_tetra.sim <- zeros.box_plot_stats(CD4_z_merged_gaps, CD4, F.MHC_groups, C.MHC_groups)
# CD8_z_plot_b.g7_tetra.sim <- zeros.box_plot_stats(CD8_z_merged_gaps, CD8, F.MHC_groups, C.MHC_groups)
# 
# pdf("Figure Export/VJ Usage Simple Heatmap/VaJa utilization by count + stats + Simulated parental sum zeros only v2.pdf", width = 5, height = 15)
# (CD4_z_plot_b.f_tetra.sim + CD8_z_plot_b.f_tetra.sim) / (CD4_z_plot_b.g7_tetra.sim + CD8_z_plot_b.g7_tetra.sim)
# dev.off()
# 
# # find VJ elements with zeros and make a table of them for the F1 and tetraparental
# # VJheat_CD4_norm_subset <- VJheat_all_data.norm %>% filter(Co_receptor == CD4)
# # makeHeatMap(VJheat_CD4_norm_subset)
# 
# F1_tet_zeros_dif <- function(VJheat_all_data.norm, F.MHC_groups){
#   VJheat_CD4_norm_subset <- VJheat_all_data.norm %>%
#     filter(MHC %in% F.MHC_groups) %>%
#     filter(Clones == 0) %>% mutate(merge = paste(V.name, J.name))
#     bxg7_set <- VJheat_CD4_norm_subset %>% filter(MHC == F.MHC_groups[1])
#     tetra.bb.ff_set <- VJheat_CD4_norm_subset %>% filter(MHC == F.MHC_groups[2])
#     zeros_unique_to_in_bxg7_set <- setdiff(bxg7_set$merge, tetra.bb.ff_set$merge)
#     bxg7_tet_zeros_subset.HM.Norm <- VJheat_all_data.norm %>%
#       filter(MHC %in% F.MHC_groups) %>%
#       mutate(merge = paste(V.name, J.name)) %>%
#       filter(merge %in% zeros_unique_to_in_bxg7_set)
#     bxg7_tet_zeros_subset.HM.Norm$MHC <- factor(bxg7_tet_zeros_subset.HM.Norm$MHC, levels = F.MHC_groups)
#     # ggHeat <- makeHeatMap(bxg7_tet_zeros_subset.HM.Norm)
#   return(bxg7_tet_zeros_subset.HM.Norm)
# }
# 
# F.MHC_groups <- c("bxf CD4", "tetraparental bb&ff CD4", "bb CD4", "ff CD4", "bb DP")
# plot1 <- F1_tet_zeros_dif(VJheat_all_data.norm, F.MHC_groups)
# bxf_Zeros_from_F1 <- plot1 %>% filter(MHC == "bxf CD4")
# 
# F.MHC_groups <- c("bxg7 CD4", "tetraparental bb&g7g7 CD4", "bb CD4", "g7g7 CD4")
# plot2 <- F1_tet_zeros_dif(VJheat_all_data.norm, F.MHC_groups)
# bxg7_Zeros_from_F1 <- plot2 %>% filter(MHC == "bxg7 CD4")
# 
# 
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ Missing in F1 present in tetraparental.pdf", width = 13, height = 13)
# makeHeatMap(plot1) / makeHeatMap(plot2)
# dev.off()
# 
# # Save multiple objects. These are the output file which take a long time to make. Save and load these when R is restarted.
# save(bxf_Zeros_from_F1, bxg7_Zeros_from_F1, file = "CD4_Zeros_from_F1.RData") # use this to scan reperotire diversity once the VJ pairs which we have excluded are accounted for
# 
# 
# # --------------------------- Find VJ elements with largest difference ------------------------------------------------------------------------------
# # Use the VJ pairs with the greatest difference or highest values and pull these out for comparisons.
# 
# # Find Diffentailly expressed VJ pairs
# # Function requires 3 inputs. Data frame 1, Data frame 2 and the number of values you want to plot from each dataframe
# # From the complex Heat map identify the values which are differntially expressed as VJ pairs between the groupings
# VJheat_uncombined <- function(immunarch.list, V.vec, J.vec){
#   immnames <- names(immunarch.list$data) # First section makes the sample dataframe.
#   nc <- ncol(immunarch.list$data[[1]])
#   spectraCDR3 <- data.frame(matrix(ncol = nc, nrow = 0))
#   n <- names(immunarch.list$data[[1]])
#   colnames(spectraCDR3) <- n
#   
#   for(i in 1:length(immnames)){
#     print(paste0("Extracting sample and metadata from: ", immnames[i]))
#     spectraCDR3_temp <- immunarch.list$data[[i]]
#     spectraCDR3_temp$Sample <- immnames[i]
#     spectraCDR3 <- rbind(spectraCDR3, spectraCDR3_temp)
#   }
#   spectraCDR3 <- merge(immunarch.list$meta, spectraCDR3, by="Sample",all=TRUE)
#   
#   # MHC Choice
#   VJheat_min <- as_tibble(data.frame(
#     "V.name" = spectraCDR3$V.name, 
#     "J.name" = spectraCDR3$J.name,
#     "MHC" = spectraCDR3$MHC, #This is equal to saying "MHC" = spectraCDR3$MHC, #"MHC" = spectraCDR3[Meta.Factor]
#     "Clones" = spectraCDR3$Clones,
#     "Sample" = spectraCDR3$Sample,
#     stringsAsFactors = TRUE))
#   # VJheat_min <- spread(VJheat_min, MHC, Sample)
#   
#   #  I think the clones with matching TRAV, TRAJ need to be collapsed to remove duplicate TRAV / TRAJ pairs
#   VJheat_min <- VJheat_min  %>%
#     dplyr::group_by(V.name, J.name, MHC, Sample)  %>%
#     dplyr::mutate(Clones = base::sum(Clones)) %>%
#     unique() %>% ungroup() %>%
#     dplyr::arrange(MHC, V.name, J.name, Sample)
#   
#   col_order <- c("V.name", "J.name", "MHC",
#                  "Sample", "Clones")
#   VJheat_min <- VJheat_min[, col_order]
#   VJheat_min <- VJheat_min %>% arrange(V.name, J.name, MHC, Sample)
#   
#   VJheat_min$MHC = paste(VJheat_min$MHC,"|",VJheat_min$Sample)
#   VJheat_min <- VJheat_min %>%  select(-Sample) # this is a trick to eliminate duplicates
#   
#   MHC_vec = as.character(unique(VJheat_min$MHC))
#   V = as.character(V.vec$V.name)
#   J = as.character(J.vec$J.name)
#   VJdim <- as_tibble(expand.grid(V.name = V, J.name = J, MHC = MHC_vec))
#   VJdim <- VJdim %>% arrange(V.name, J.name, MHC)
#   
#   VJheat_max <- merge(VJheat_min, VJdim, by = c("V.name", "J.name", "MHC"), all = TRUE) #Generates NA values where there was previously not a column
#   VJheat_max$Clones <- VJheat_max$Clones %>% replace_na(0) # convert NA to 0 so it can be plotted
#   VJheat_max <- VJheat_max %>% separate(MHC, c("MHC", "Sample") , sep = "([|])")
#   # Remove artifact rows which somehow show up in some heatmaps
#   VJheat_max <- VJheat_max %>% arrange(V.name, J.name, MHC, Sample) %>% filter(V.name != "Vfam-Vsub")
#   VJheat_max <- VJheat_max %>% arrange(V.name, J.name, MHC, Sample) %>% filter(J.name != "J")
#   
#   VJheat_max$MHC <- gsub(" ", "", VJheat_max$MHC)
#   
#   return(VJheat_max) # this is the end of the function and what the function returns to me.
# } #generates the source data for the heatmap
# CD4_Uncombined <- VJheat_uncombined(TCR_DoVb_CD4_AA.ones.removed.LFSR, V.vec, J.vec)
# CD8_Uncombined <- VJheat_uncombined(TCR_DoVb_CD8_AA.ones.removed.LFSR, V.vec, J.vec)
# 
# 
# head(VJheat_all_data)
# head(VJheat_all_data.norm)
# N.num = 10
# F.MHC_groups <- c(bb, ff, bxf, b.f)
# C.MHC_groups <- c(bb_color, ff_color, bxf_color, b.f_color)
# # input: VJheat_all_data.norm, CD4_Uncombined, N.num, F.MHC_groups, C.MHC_groups, co_receptor
# 
# # top_VJ_pairs <- function(VJheat_all_data.norm, CD4_Uncombined, N.num, F.MHC_groups, C.MHC_groups, co_receptor){
# #   for(i in 1:length(F.MHC_groups)){
# #   bb_CD4_data <- VJheat_all_data.norm %>% filter(MHC == paste(F.MHC_groups[1], co_receptor))
# #   
# #   
# #   }
# # }
# # 
# # union_list <- Reduce(union, list(top_bb_CD4, top_ff_CD4, top_b.f_CD4, top_bxf_CD4))
# # union_list <- union_list %>% select(V.name, J.name) %>% unique() %>%
# # --------------------------
# save(VJheat_all_data.norm, file = "VJheat_all_data.norm.RData")
# b.g7_CD4_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "bxg7 CD4", "g7g7 CD4", "b+g7 CD4"))
# b.g7_CD8_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD8", "bxg7 CD8", "g7g7 CD8", "b+g7 CD8"))
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ-b.g7_CD4andCD8.pdf", width = 15, height = 15)
# makeHeatMap(b.g7_CD4_VJ_gg) / makeHeatMap(b.g7_CD8_VJ_gg)
# dev.off()
# 
# b.f_CD4_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "bxf CD4", "ff CD4", "b+f CD4"))
# b.f_CD8_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD8", "bxf CD8", "ff CD8", "b+f CD8"))
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ-b.f_CD4andCD8.pdf", width = 15, height = 15)
# makeHeatMap(b.f_CD4_VJ_gg) / makeHeatMap(b.f_CD8_VJ_gg)
# dev.off()
# 
# b.s_CD4_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "bxs CD4", "ss CD4"))
# f.s_CD4_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("ff CD4", "fxs CD4", "ss CD4"))
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ-b.s_CD4and_f.sCD8.pdf", width = 15, height = 15)
# makeHeatMap(b.s_CD4_VJ_gg) / makeHeatMap(f.s_CD4_VJ_gg)
# dev.off()
# 
# # MHChom comparisons
# MHC.hom_CD4_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "ff CD4", "g7g7 CD4", "ss CD4"))
# MHC.hom_CD8_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD8", "ff CD8", "g7g7 CD8"))
# MHC.hom_CD4.CD8_VJ_gg <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "ff CD4", "g7g7 CD4", "ss CD4", "bb CD8", "ff CD8", "g7g7 CD8"))
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ-MHC_homCD4CD8.pdf", width = 15, height = 15)
# # makeHeatMap(MHC.hom_CD4_VJ_gg) / makeHeatMap(MHC.hom_CD8_VJ_gg)
# makeHeatMap(MHC.hom_CD4.CD8_VJ_gg)
# dev.off()
# 
# unique(VJheat_all_data.norm$MHC)
# MHC.hom_CD4_VJ_bxg7.tetra.b.g7_CD4 <- VJheat_all_data.norm %>% filter(MHC %in% c("bxg7 CD4", "tetraparental bb&g7g7 CD4"))
# pdf("Figure Export/VJ Usage Simple Heatmap/MHC.hom_CD4_VJ_bxg7.tetra.b.g7_CD4.pdf", width = 15, height = 15)
# # makeHeatMap(MHC.hom_CD4_VJ_gg) / makeHeatMap(MHC.hom_CD8_VJ_gg)
# makeHeatMap(MHC.hom_CD4_VJ_bxg7.tetra.b.g7_CD4)
# dev.off()
# 
# # Combine CD4 and CD8 samples to compare with thymus data.
# bb_CD4_CD8_DP <- VJheat_all_data.norm %>% mutate(MHC=recode(MHC, 
#                          `bb CD4`="bb CD4+CD8",
#                          `bb CD8`="bb CD4+CD8")
#                          ) %>%
#                 filter(MHC %in% c("bb CD4+CD8", "bb DP")) %>%
#                 group_by(MHC, V.name, J.name) %>%
#                 mutate(Clones = sum(Clones))
# temp_bb <- VJheat_all_data.norm %>% filter(MHC %in% c("bb CD4", "bb CD8"))
# bb_CD4_CD8_DP <- rbind(bb_CD4_CD8_DP, temp_bb)
# bb_CD4_CD8_DP <- bb_CD4_CD8_DP %>% mutate(MHC=recode(MHC,
#                           `bb CD4`="1 bb CD4",
#                           `bb CD8`="2 bb CD8",
#                           `bb CD4+CD8`="3 bb CD4+CD8",
#                           `bb DP`="4 bb DP"))
# ff_CD4_CD8_DP <- VJheat_all_data.norm %>% mutate(MHC=recode(MHC, 
#                           `ff CD4`="ff CD4+CD8",
#                           `ff CD8`="ff CD4+CD8")
#                           ) %>%
#                 filter(MHC %in% c("ff CD4+CD8", "ff DP"))
#                 group_by(MHC, V.name, J.name) %>%
#                 mutate(Clones = sum(Clones))
# temp_ff <- VJheat_all_data.norm %>% filter(MHC %in% c("ff CD4", "ff CD8"))
# ff_CD4_CD8_DP <- rbind(ff_CD4_CD8_DP, temp_ff)
# ff_CD4_CD8_DP <- ff_CD4_CD8_DP %>% mutate(MHC=recode(MHC,
#                                                      `ff CD4`="1 ff CD4",
#                                                      `ff CD8`="2 ff CD8",
#                                                      `ff CD4+CD8`="3 ff CD4+CD8",
#                                                      `ff DP`="4 ff DP"))
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ-bb_CD4_CD8_DP and ff_CD4_CD8_DP.pdf", width = 15, height = 15)
# makeHeatMap(bb_CD4_CD8_DP) / makeHeatMap(ff_CD4_CD8_DP)
# dev.off()
# # unique(df$MHC)
# # --------------------------
# 
# 
# bb_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bb CD4")
# bb_CD4_data_gg <- makeHeatMap(bb_CD4_data)
# bb_CD4_data_gg
# 
# ff_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "ff CD4")
# bxf_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bxf CD4")
# b.f_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "b+f CD4")
# g7g7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "g7g7 CD4")
# bxg7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "bxg7 CD4")
# b.g7_CD4_data <- VJheat_all_data.norm %>% filter(MHC == "b+g7 CD4")
# 
# bb_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bb CD8")
# ff_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "ff CD8")
# bxf_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bxf CD8")
# b.f_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "b+f CD8")
# g7g7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "g7g7 CD8")
# bxg7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "bxg7 CD8")
# b.g7_CD8_data <- VJheat_all_data.norm %>% filter(MHC == "b+g7 CD8")
# 
# # go through each separated CD4 data and pull out the top 10 clones from each
# top_bb_CD4 <- ungroup(bb_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_ff_CD4 <- ungroup(ff_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_bxf_CD4 <- ungroup(bxf_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_b.f_CD4 <- ungroup(b.f_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_g7g7_CD4 <- ungroup(g7g7_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_bxg7_CD4 <- ungroup(bxg7_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_b.g7_CD4 <- ungroup(b.g7_CD4_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# 
# # go through each separated CD8 data and pull out the top 10 clones from each
# top_bb_CD8 <- ungroup(bb_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_ff_CD8 <- ungroup(ff_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_bxf_CD8 <- ungroup(bxf_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_b.f_CD8 <- ungroup(b.f_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_g7g7_CD8 <- ungroup(g7g7_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_bxg7_CD8 <- ungroup(bxg7_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# top_b.g7_CD8 <- ungroup(b.g7_CD8_data) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# 
# colsToUse_b.f_CD4 <- as_tibble(dplyr::union(top_bb_CD4[,c("V.name", "J.name")], top_ff_CD4[,c("V.name", "J.name")], top_b.f_CD4[,c("V.name", "J.name")], top_bxf_CD4[,c("V.name", "J.name")])) 
# colsToUse_b.f_CD4 <-  colsToUse_b.f_CD4 %>% slice(1:N.num) # Only collect the top N columns
# 
# colsToUse_b.g7_CD4 <- as_tibble(dplyr::union(top_bb_CD4[,c("V.name", "J.name")], top_g7g7_CD4[,c("V.name", "J.name")], top_b.g7_CD4[,c("V.name", "J.name")], top_bxg7_CD4[,c("V.name", "J.name")])) 
# colsToUse_b.g7_CD4 <-  colsToUse_b.g7_CD4 %>% slice(1:N.num) # Only collect the top N columns
# 
# colsToUse_b.f_CD8 <- as_tibble(dplyr::union(top_bb_CD8[,c("V.name", "J.name")], top_ff_CD8[,c("V.name", "J.name")], top_b.f_CD8[,c("V.name", "J.name")], top_bxf_CD8[,c("V.name", "J.name")])) 
# colsToUse_b.f_CD8 <-  colsToUse_b.f_CD8 %>% slice(1:N.num) # Only collect the top N columns
# 
# colsToUse_b.g7_CD8 <- as_tibble(dplyr::union(top_bb_CD8[,c("V.name", "J.name")], top_g7g7_CD8[,c("V.name", "J.name")], top_b.g7_CD8[,c("V.name", "J.name")], top_bxg7_CD8[,c("V.name", "J.name")])) 
# colsToUse_b.g7_CD8 <-  colsToUse_b.g7_CD8 %>% slice(1:N.num) # Only collect the top N columns
# 
# 
# 
# # Use the search column to pull out indexed data needed from the MHC matching samples from the CD4_Uncombined or CD8_Uncombined data.
# 
# # This is not the best function. Should be made more dynamic and less hardcoded.
# top_VJ <- function(colsToUse, data_Uncombined, F.MHC_groups, C.MHC_groups){
#   colsToUse <- colsToUse %>%
#     mutate(VJ = paste0(V.name, "|", J.name))
#   VJ.df1.M <- data_Uncombined %>%
#     mutate(VJ = paste0(V.name, "|", J.name)) %>% 
#     subset(MHC %in% F.MHC_groups) %>%
#     subset(VJ %in% colsToUse$VJ) %>%
#     filter(!(VJ == "TRAV11|TRAJ18")) # Remove the invariant TRAV11, TRAJ18 Pair from the dataset. This is an artifact of the data.
#   
#   plot <- ggplot(data = VJ.df1.M, aes(x = reorder(VJ, -Clones), y = Clones, fill = MHC)) +
#     geom_boxplot(width = 0.75, position = position_dodge()) +
#     scale_y_continuous(expand = expansion(mult = c(0, .1))) +
#     coord_flip() +
#     theme_classic() +
#     theme(legend.position="bottom", legend.box="vertical") +
#     scale_fill_manual(breaks = F.MHC_groups, 
#                       values=C.MHC_groups) +
#     scale_color_manual(breaks = F.MHC_groups, 
#                        values=C.MHC_groups) +
#     labs(
#       x ="TRAV|TRAJ Pair",
#       y = "Frequency",
#       fill = "MHC \n Group"
#     )
#   return(plot)
# }
# F.MHC_groups <- c(bb, ff, bxf, b.f)
# C.MHC_groups <- c(bb_color, ff_color, bxf_color, b.f_color)
# top_VJ_CD4_b.f <- top_VJ(colsToUse_b.f_CD4, CD4_Uncombined, F.MHC_groups, C.MHC_groups)
# top_VJ_CD8_b.f <- top_VJ(colsToUse_b.f_CD8, CD8_Uncombined, F.MHC_groups, C.MHC_groups)
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ top pairs b and f groups.pdf", width = 9, height = 5)
# top_VJ_CD4_b.f + top_VJ_CD8_b.f
# dev.off()
# 
# F.MHC_groups <- c(bb, g7g7, bxg7, b.g7)
# C.MHC_groups <- c(bb_color, g7g7_color, bxg7_color, b.g7_color)
# top_VJ_CD4_b.g7 <- top_VJ(colsToUse_b.g7_CD4, CD4_Uncombined, F.MHC_groups, C.MHC_groups)
# top_VJ_CD8_b.g7 <- top_VJ(colsToUse_b.g7_CD8, CD8_Uncombined, F.MHC_groups, C.MHC_groups)
# 
# pdf("Figure Export/VJ Usage Simple Heatmap/VJ top pairs b and g7 groups.pdf", width = 9, height = 5)
# top_VJ_CD4_b.g7 + top_VJ_CD8_b.g7
# dev.off()
# 
# 
# 
# # -----------------------------------------------------------------------------------------------------------------------------------------------------
# # most different genes between two groups
# dif <- ungroup(bb_CD4_data) %>% select(c("V.name", "J.name")) %>% mutate(Clones = abs(bb_CD4$Clones - bxf_CD4_data$Clones)) %>% arrange(desc(Clones)) %>% slice(1:N.num)
# # highest expressed genes among various groups comapred to each other.

