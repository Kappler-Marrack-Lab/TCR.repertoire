# Query for alpha chains which have the WT DOalpha from DO-11.10 T cell hyrbidoma from the  samples sets
# First focus on CD4, CD8 from bb, g7g7, and bxg7
library(immunarch)
library(ggplot2)
library(parallel)
library(tidyverse)
library(viridis)
library(scales)
# set.seed(42)
#=================================
# Read in data
#=================================
# ---------------------------------------------------------------------------------------------------------
# Save Workspace file.
# save.image(file = "DOalpha Query.RData")
# # To restore your workspace, type this:
# load("DOalpha Query.RData")
# ---------------------------------------------------------------------------------------------------------
load("TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData")


# Reference Data.
# ---------------------------------------------------------------------------------------------------------
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"
IAbHet = "b+/-"
`bb&ff` = "bb&ff"
`bb&g7g7` = "bb&g7g7"

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
`bb&ff_color` = "#EEC0C4"
`bb&g7g7_color` = "#335566"


# -------------------------------------------------------------------------------------------------------
target_CAAS <- tibble(CDR3.aa = "CAASPNYNVLYF", V.name = "TRAV05-4DN", J.name = "TRAJ21")
target_CAAA <- tibble(CDR3.aa = "CAAAPNYNVLYF", V.name = "TRAV05-4DN", J.name = "TRAJ21") # this is the JSB autoreactive mutation
# target <- tibble(V.name = "TRAV05-4DN", J.name = "TRAJ21")

extract_clonotype_CD4 <- function(list, target){
  tc <- trackClonotypes(list$data, target)
  tc_df <- as_tibble(tc)
  tc_df <- mutate(tc_df, Clonotype = paste(tc_df$V.name, tc_df$CDR3.aa, tc_df$J.name))
  clone <- tc_df$Clonotype
  
  tc_df_2 <- select(tc_df, -c(Clonotype, CDR3.aa, V.name, J.name))
  tc_df_2 <- as.data.frame(t(as.matrix(tc_df_2)))
  tc_df_2 <- tibble::rownames_to_column(tc_df_2, "Sample")
  tc_df_2$Clonotype <- clone
  names(tc_df_2) <- c("Sample", "Proportion", "Clonotype")
  tc_df <- merge(list$meta, tc_df_2, by="Sample",all=TRUE)
  return(tc_df)
}

extract_clonotype_CD8 <- function(list, target){
  tc <- trackClonotypes(list$data, target)
  tc_df <- as_tibble(tc)
  tc_df <- mutate(tc_df, Clonotype = paste(tc_df$V.name, tc_df$CDR3.aa, tc_df$J.name))
  clone <- tc_df$Clonotype
  
  tc_df_2 <- select(tc_df, -c(Clonotype, CDR3.aa, V.name, J.name))
  tc_df_2 <- as.data.frame(t(as.matrix(tc_df_2)))
  tc_df_2 <- tibble::rownames_to_column(tc_df_2, "Sample")
  tc_df_2$Clonotype <- clone
  names(tc_df_2) <- c("Sample", "Proportion", "Clonotype")
  tc_df <- merge(list$meta, tc_df_2, by="Sample",all=TRUE)
  return(tc_df)
}

CD4_CAAS <- extract_clonotype_CD4(TCR_DoVb_CD4_AA.ones.removed, target_CAAS)
CD8_CAAS <- extract_clonotype_CD8(TCR_DoVb_CD8_AA.ones.removed, target_CAAS)
CD4_CAAA <- extract_clonotype_CD4(TCR_DoVb_CD4_AA.ones.removed, target_CAAA)
CD8_CAAA <- extract_clonotype_CD8(TCR_DoVb_CD8_AA.ones.removed, target_CAAA)

F.MHC_groups <- c(bb, ff, g7g7, ss,
                  bxf, bxg7, bxs, fxs,
                  `bb&ff`, `bb&g7g7`, IAbHet)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color,
                  bxf_color, bxg7_color, bxs_color, fxs_color,
                  `bb&ff_color`, `bb&g7g7_color`, IAbHet_color)  
 
CD4_CAAS$MHC <- factor(CD4_CAAS$MHC, levels = F.MHC_groups)
CD4_CAAA$MHC <- factor(CD4_CAAA$MHC, levels = F.MHC_groups)

DOalpha_WT_query_CD4_plot <- ggplot(CD4_CAAS, aes(fill=MHC, y=Proportion, x=MHC)) + 
  # geom_bar(position="dodge", stat="summary", fun = mean)+
  stat_summary(fun.data=mean_se, geom = "linerange", position=position_dodge(width=1))+
  geom_jitter(size = 2, position = position_jitterdodge(3), shape = 21, color = "black")+
  theme_classic()+
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(paste0("Frequency of Clonotype:", paste(target_CAAS[1,], collapse = " "), "\nacross CD4 MHC haplotypes"))

DOalpha_S93A_query_CD4_plot <- ggplot(CD4_CAAA, aes(fill=MHC, y=Proportion, x=MHC)) + 
  # geom_bar(position="dodge", stat="summary", fun = mean)+
  stat_summary(fun.data=mean_se, geom = "linerange", position=position_dodge(width=1))+
  geom_jitter(size = 2, position = position_jitterdodge(3), shape = 21, color = "black")+
  theme_classic()+
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(paste0("Frequency of Clonotype:", paste(target_CAAA[1,], collapse = " "), "\nacross CD4 MHC haplotypes"))

DOalpha_WT_query_CD4_plot / DOalpha_S93A_query_CD4_plot

pdf("Figure Export/DObeta clonotypes/DOb11.10 alpha clonotypes.pdf", width = 4.5, height = 5)
DOalpha_WT_query_CD4_plot / DOalpha_S93A_query_CD4_plot
dev.off()
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# recreate plot above with fewer groups to simplify the plot:
F.MHC_groups <- c(bb,
                  bxf, bxg7, bxs,
                  `bb&ff`, `bb&g7g7`, IAbHet)
C.MHC_groups <- c(bb_color,
                  bxf_color, bxg7_color, bxs_color,
                  `bb&ff_color`, `bb&g7g7_color`, IAbHet_color)  

CD4_CAAS$MHC <- factor(CD4_CAAS$MHC, levels = F.MHC_groups)
CD4_CAAA$MHC <- factor(CD4_CAAA$MHC, levels = F.MHC_groups)

CD4_CAAS.sub <- CD4_CAAS %>% filter(MHC %in% F.MHC_groups)
CD4_CAAA.sub <- CD4_CAAA %>% filter(MHC %in% F.MHC_groups)

DOalpha_WT_query_CD4_plot_minimal <- ggplot(CD4_CAAS.sub, aes(fill=MHC, y=Proportion, x=MHC)) + 
  # geom_bar(position="dodge", stat="summary", fun = mean)+
  stat_summary(fun.data=mean_se, geom = "linerange", position=position_dodge(width=1))+
  geom_jitter(size = 2, position = position_jitterdodge(3), shape = 21, color = "black")+
  theme_classic()+
  
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-6, 0.0003))+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(paste0("Frequency of Clonotype:", paste(target_CAAS[1,], collapse = " "), "\nacross CD4 MHC haplotypes"))

DOalpha_S93A_query_CD4_plot_minimal <- ggplot(CD4_CAAA.sub, aes(fill=MHC, y=Proportion, x=MHC)) + 
  # geom_bar(position="dodge", stat="summary", fun = mean)+
  stat_summary(fun.data=mean_se, geom = "linerange", position=position_dodge(width=1))+
  geom_jitter(size = 2, position = position_jitterdodge(3), shape = 21, color = "black")+
  theme_classic()+
  scale_fill_manual(breaks = F.MHC_groups, values=C.MHC_groups) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e-6, 0.0003))+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(paste0("Frequency of Clonotype:", paste(target_CAAA[1,], collapse = " "), "\nacross CD4 MHC haplotypes"))

DOalpha_WT_query_CD4_plot_minimal / DOalpha_S93A_query_CD4_plot_minimal

pdf("Figure Export/DObeta clonotypes/DOb11.10 alpha clonotypes_minimal.pdf", width = 3.75, height = 5)
DOalpha_WT_query_CD4_plot_minimal / DOalpha_S93A_query_CD4_plot_minimal
dev.off()


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Average frequency per group:
CD4_CAAS %>% group_by(MHC) %>% 
  summarize(avg = scientific(mean(Proportion), digits = 3), n = n(), 
            sd = sd(Proportion), se = sd/sqrt(n))

tc_CD4 <- trackClonotypes(TCR.CD4_DoVb_DB$data, target)
irp_tc_CD4 <- vis(tc_CD4, .plot = "smooth")
tc_CD8 <- trackClonotypes(TCR.CD8_DoVb_DB$data, target)
irp_tc_CD8 <- vis(tc_CD8, .plot = "smooth")
irp_tc_CD4 + irp_tc_CD8

# Top CD4 = 0.0002602783 = ff
# Top CD8 = 0.0001672742 = g7g7 

# ----------------------------
# Filter the dataset for the V and J the WT alpha resides in.
# then look at the frequency of all the CDR3 in this space and see what the frequencey is per haplotype.
VJ_filter <- function(x){
  y <- x %>% filter(V.name == "TRAV05-4DN") %>% filter(J.name == "TRAJ21") %>% filter(Clones > 1) %>%
    select(Clones, Proportion, CDR3.aa, CDR3.nt)
  return(y)
}
DOalpha_CDR3_Query <- mclapply(TCR.CD4_DoVb_DB$data, VJ_filter, mc.cores = 6)

nc <- ncol(DOalpha_CDR3_Query[[1]])
DOalphaCDR3 <- data.frame(matrix(ncol = nc, nrow = 0))
n <- names(DOalpha_CDR3_Query[[1]])
colnames(DOalphaCDR3) <- n

# Extract the list names
immnames <- names(DOalpha_CDR3_Query) # First section makes the sample dataframe.
for(i in 1:length(immnames)){
  print(paste0("Extracting sample and metadata from: ", immnames[i]))
  temp <- DOalpha_CDR3_Query[[i]]
  temp$Sample <- immnames[i]
  DOalphaCDR3 <- rbind(DOalphaCDR3, temp)
}
DOalphaCDR3 <- merge(TCR.CD4_DoVb_DB$meta, DOalphaCDR3, by="Sample",all=TRUE)
g <- unique(DOalphaCDR3$CDR3.aa)
length(g)

VJ_expand <- DOalphaCDR3 %>%  expand(CDR3.aa, MHC)
VJ_full <- DOalphaCDR3 %>% select(Clones, CDR3.aa, MHC)  %>% dplyr::right_join(VJ_expand)
VJ_full$Clones <- VJ_full$Clones %>% replace_na(0)
VJ_full <- VJ_full %>% group_by(MHC, CDR3.aa) %>% dplyr::summarise(Clones = sum(Clones))
ggplot(VJ_full, aes(x = CDR3.aa, y = MHC, fill = Clones)) +
  geom_tile(color="white", size=0.1) +
  scale_fill_viridis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# ---------------------------------------------------------------------------------------------------------
# test the presence of a wobble in the usage of CDR3 codons.
wobble_test <- DOalphaCDR3 %>% filter(CDR3.aa == "CAAAPNYNVLYF")
g <- unique(wobble_test$CDR3.nt)
length(g)
# boxplot
# wobble plot
wobble_plot <- ggplot(wobble_test, aes(fill=CDR3.nt, y=Proportion, x=MHC)) + 
  geom_boxplot() +
  # geom_point(fill=alpha("red", 0.3),alpha=0.7) +
  ggtitle(paste0("Frequency of CAAAPNYNVLYF  Clonotypes\n ", " from TRAV05-4DN + TRAJ21 pairs", "across CD4 MHC haplotypes"))
# ff samples ahve a bias to the known CDR3 dna clone, but this could be due to germline DNA sequence. However other CDR3.nt variants still exist which is promising.
pdf("CAAAPNYNVLYF Dobalpha clonotypes.pdf", width = 14, height = 8)
wobble_plot
dev.off()

# boxplot
ggplot(DOalphaCDR3, aes(fill=MHC, y=Proportion, x=CDR3.aa)) + 
  geom_boxplot() +
  # geom_point(fill=alpha("red", 0.3),alpha=0.7) +
  ggtitle(paste0("Frequency of Clonotypes ", " from TRAV05-4DN + TRAJ21 pairs", "across CD4 MHC haplotypes"))

# heatmap
ggplot(DOalphaCDR3, aes(fill=MHC, y=Proportion, x=CDR3.aa)) + 
  geom_tile(color="white", size=0.1) +
  ggtitle(paste0("Frequency of Clonotypes ", " from TRAV05-4DN + TRAJ21 pairs", "across CD4 MHC haplotypes"))

