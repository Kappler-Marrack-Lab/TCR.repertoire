# Plot flow data into a heatmap
# TCR_clones_flow_ANOVA
library(tidyverse)
library(ggpubr)
library(vroom)
# library(rstatix) # don't need this for the heatmaps
library(readxl)
library(multcomp)
library(ggprism)
library(immunarch)
library(RColorBrewer)

# Colors
# ff ="ff SC, IL-4 prestim (Fresh SC)"
# bb = "bb SC, IL-4 prestim (Fresh SC)"
# g7g7 = "g7g7 SC, IL-4 prestim (Fresh SC)"
# bxf = "bxf SC, IL-4 prestim (Fresh SC)"
# bxg7 = "bxg7 SC, IL-4 prestim (Fresh SC)"
ff ="ff"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
bxg7 = "bxg7"
tetra_bb.ff = "tetra. bb&ff"
tetra_bb.g7g7 = "tetra. bb&g7g7"

# Define MHC group colors:
ff_color ="#c5944e"
bb_color = "#54a0fb" 
g7g7_color = "#9f104c"
bxf_color = "#808000"
bxg7_color = "#f74ed6"
tetra_bb.ff_color = "#EEC0C4"
tetra_bb.g7g7_color = "#335566"


load(file = "TCR_fixed_beta_total_rep.CDR3AA.ones.removed.RData") # filtered and ones removed
load(file = "counts tables filtered.ones.removed.RData") # counts table with only the ones removed



TCR.response.df <- read_xlsx("/Users/LappyPro/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Kappler Lab/_Lab Databases/TCR Clones.xlsx",
                             sheet = "Master Response Table")

# VDJ_lookup <- vroom(file = "MK_Ion_IMGT_Convert.csv", delim = ",")
# IMGT_to_arden <- function(VDJ_lookup, df){
#   # Convert to MK format
#   new <- as.data.frame(df)
#   new[] <- VDJ_lookup$MK_Ion_TRAV[match(unlist(df), VDJ_lookup$IMGT_TRAV)]
#   new.V.name <- new$TRAV
#   new[] <- VDJ_lookup$MK_Ion_TRAJ[match(unlist(df), VDJ_lookup$IMGT_TRAJ)]
#   new.J.name <- new$TRAJ
#   df <- df %>% mutate(Arden_TRAV = new.V.name) %>% mutate(Arden_TRAJ = new.J.name)
#   new <- as.data.frame(df)
#   new[] <- VDJ_lookup$CDR1[match(unlist(df), VDJ_lookup$MK_Ion_TRAV)]
#   new.CDR1 <- new$Arden_TRAV
#   new[] <- VDJ_lookup$CDR2[match(unlist(df), VDJ_lookup$MK_Ion_TRAV)]
#   new.CDR2 <- new$Arden_TRAV
#   df <- df %>% mutate(CDR1 = new.CDR1) %>% mutate(CDR2 = new.CDR2) %>%
#     relocate(c("Arden_TRAV", "CDR1", "CDR2"), .after = TRAV) %>%
#     relocate(Arden_TRAJ, .after = TRAJ)
#   return(df)
# }

# TCR.response.df <- IMGT_to_arden(VDJ_lookup, TCR.response.df)

# Create a new counts table of the following groups:
# bb, ff, g7g7, bxf, bxg7, bb&ff tetra, bb&g7g7 tetra


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Apply the pr_CD4s_filtered.ones.removed matrix instead of what we had before. This is contains the missing INS sequences.
# Select out the columns which are relevant.
pr.CD4.mat.P2M.DMR <- pr_CD4s_filtered.ones.removed %>%
  dplyr::select(V.name, J.name, CDR3.aa,
              bb.ave.prop, bb.Samples,
              ff.ave.prop, ff.Samples,
              g7g7.ave.prop, g7g7.Samples,
              bxf.ave.prop, bxf.Samples,
              bxg7.ave.prop, bxg7.Samples,
              `bb&ff.ave.prop`, bf_tetra.Samples,
              `bb&g7g7.ave.prop`, bg7_tetra.Samples
)

# organize a table with the relevant clone counts information and identify NA values.
pr.CD4.mat.P2M.DMR <- pr.CD4.mat.P2M.DMR %>% mutate(V.J.CDR3 = paste(V.name, J.name, CDR3.aa, sep = "_"))
TCR.response.df <- TCR.response.df %>% mutate(V.J.CDR3 = paste(Arden_TRAV, Arden_TRAJ, TRA_CDR3, sep = "_"))
joint.df <- TCR.response.df %>% left_join(pr.CD4.mat.P2M.DMR, by = "V.J.CDR3") %>%
  dplyr::select(-V.name, -J.name, -CDR3.aa, -V.J.CDR3) %>%
  relocate(bb.ave.prop:bg7_tetra.Samples, .after = TRA_CDR3) %>%
  # mutate(across(bb.ave.prop:bg7_tetra.Samples, ~na_if(., "NaN"))) %>%
  mutate(across(c(`Matched oPool Seq #`, `Notes`), ~na_if(., "NA"))) %>%
  mutate(across(where(is.character), ~na_if(., "n/a"))) %>%
  # mutate(across(c(`TCRbeta MFI (Run 1)`:`ZsGreen_NFAT+ MFI (Run 3)`), ~na_if(., "n/a"))) %>% # 2023.07.01 Commented out. No n/a values here and causes an issue with dplyr
  dplyr::mutate(APC = factor(`Antigen Presenting Cell (APC)`)) %>%
  dplyr::mutate(APC = sub(" .*", "", APC)) %>%
  rowwise() %>% 
  dplyr::mutate(mean.NFATZsGreen.pos=mean(c(`ZsGreen_NFAT+ (Run 1)`, `ZsGreen_NFAT+ (Run 2)`, `ZsGreen_NFAT+ (Run 3)`))) %>%
  ungroup()%>%
  dplyr::mutate(Group = recode(Group, "in bb not g7g7 or bxg7" = "in bb NOT g7g7 or bxg7"))
# 
# # Read in the excel file which has notations of what are the clones incompletely negatively selected
# Backchecked_TCR_responses <- read_excel("Figure Export/5KC NFAT Reporter/Data Files/TCR_clones.xlsx", sheet = "TCR_clones") %>% na_if("NA")
# 
# 
# Known_oPool_seq <- read_excel("/Users/LappyPro/Downloads/Match.xlsx")
# 
# TCRdf_joint <- left_join(Backchecked_TCR_responses, Known_oPool_seq, by = "Matched oPool Seq #") %>%
#   dplyr::select(Group, `Sequence ID`, `Matched oPool Seq #`, TRAV, Arden_TRAV, TRAJ, Arden_TRAJ, TRA_CDR3,
#                 TRAV_Opool, TRAJ_Opool, TRA_CDR3_Opool, V.name_Opool, J.name_Opool, CDR3.aa_Opool)
# vroom_write(TCRdf_joint, file = "/Users/LappyPro/Downloads/TCRdf_joint.csv", delim = ",")
# 
# 
# 
# joint.df$Notes <- Backchecked_TCR_responses$Notes
# 


# Add in variables to APC group to view repertoire frequency in different groupings:
# Easy way to do this is to create duplicated rows and assign new variable name. These will only be used to plot
# repertoire frequency and not NFAT response.
duplicated_rows.tetra.bbff <- joint.df %>%
  dplyr::filter(Group %in% c("in bb NOT ff or bxf", "in ff NOT bb or bxf")) %>%
  dplyr::filter(APC == "bb") %>%
  dplyr::mutate(APC = recode(APC, "bb" = "tetra. bb&ff")) %>%
  dplyr::mutate(mean.NFATZsGreen.pos = NA)

duplicated_rows.tetra.bbg7g7 <- joint.df %>%
  dplyr::filter(Group %in% c("in bb NOT g7g7 or bxg7", "in g7g7 NOT bb or bxg7")) %>%
  dplyr::filter(APC == "bb") %>%
  dplyr::mutate(APC = recode(APC, "bb" = "tetra. bb&g7g7")) %>%
  dplyr::mutate(mean.NFATZsGreen.pos = NA)

joint.df <- rbind(joint.df, duplicated_rows.tetra.bbff, duplicated_rows.tetra.bbg7g7)


joint.df.long_NFAT.pos <- joint.df %>%
  # gather(key = Replicate, value = Percent_NFAT_Pos., c(`ZsGreen_NFAT+ (Run 1)`, `ZsGreen_NFAT+ (Run 2)`, `ZsGreen_NFAT+ (Run 3)`), factor_key=TRUE) %>%
  
  dplyr::mutate(Sequence.ID = factor(`Sequence ID`)) %>%
  dplyr::mutate(bf_tetra.Samples = if_else(bf_tetra.Samples > 0, "Present in bb&ff\nTetraparental", "Absent in bb&ff\nTetraparental", "NA")) %>%
  dplyr::mutate(bg7_tetra.Samples = if_else(bg7_tetra.Samples > 0, "Present in bb&g7g7\nTetraparental", "Absent in bb&g7g7\nTetraparental", "NA")) %>%
  dplyr::mutate(Incomplete_Negative_Selection = if_else(Notes == "INS", "INS", "NA")) %>%
  mutate(Incomplete_Negative_Selection = na_if(Incomplete_Negative_Selection, "NA")) %>%
  dplyr::mutate(V.CDR3.J = paste(Arden_TRAV, TRA_CDR3, Arden_TRAJ, sep = "_")) %>%
  
  # Create a column based on the average frequency in the repertoire for the specific TCR sequences we have tested.
  dplyr::mutate(Frequency = case_when(APC == "bb" ~ bb.ave.prop,
                                      APC == "ff" ~ ff.ave.prop,
                                      APC == "g7g7" ~ g7g7.ave.prop,
                                      APC == "bxf" ~ bxf.ave.prop,
                                      APC == "bxg7" ~ bxg7.ave.prop,
                                      APC == "tetra. bb&ff" ~ `bb&ff.ave.prop`,
                                      APC == "tetra. bb&g7g7" ~ `bb&g7g7.ave.prop`
                                      )
                ) %>%
  
  dplyr::mutate(Frequency = na_if(Frequency, 0)) %>%
  dplyr::select(Group, Incomplete_Negative_Selection, Sequence.ID, V.CDR3.J, APC,
                bf_tetra.Samples, bg7_tetra.Samples,  mean.NFATZsGreen.pos, Frequency,
                bb.ave.prop, ff.ave.prop, g7g7.ave.prop,
                bxf.ave.prop, bxg7.ave.prop, `bb&ff.ave.prop`, `bb&g7g7.ave.prop`) %>%
  
  # Need to remove TCR: 129, 178, 298, 286, 300, 294 because they were cloned onto the wrong TRAV gene backbone...
  dplyr::filter(Sequence.ID %in% c("TCR_129", "TCR_178", "TCR_298", "TCR_286", "TCR_300", "TCR_294", "TCR_287", "TCR_288", "TCR_026", "gBlock TCR_1047", "gBlock TCR_1037", "gBlock TCR_1034") == FALSE) %>%
  dplyr::mutate(Sequence.ID = sub(".*TCR_", "", Sequence.ID))

# Separate out the incomplete negative selection values:
ISN_sub <- joint.df.long_NFAT.pos %>%
  filter(Incomplete_Negative_Selection == "INS")

F1.absent_sub <- joint.df.long_NFAT.pos %>%
  filter(is.na(Incomplete_Negative_Selection)) #remove the incompleately negatively selected samples. We can plot these separatly


# ------------------- HEATMAP FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# singular functions for making heatmaps if needed to do so separately.
Avatar.heatmap.rf <- function(df, Groupings, F.MHC_groups, C.MHC_groups){
  palette_purples <- colorRampPalette(colors = c("white", "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d"))(12)
  
  plot.1 <- ggplot(subset(df, Group %in% Groupings), aes(x = APC, y = Sequence.ID, fill = Frequency))+
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) +
    coord_equal()+
    theme_classic()+
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position="left",
      plot.margin = unit(c(t=5.5, r=0, b=5.5, l=5.5), "pt")
    )+
    scale_fill_gradientn(colours = palette_purples, guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    # theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")+
    # geom_vline(xintercept = seq(0.5, length(df$Sequence.ID), by = 1), color="gray", size=.5, alpha=.5) + # set vertical lines between x groups
    labs(x = "MHC II", y = "TCR Avatar #", fill = "Repertoire\nFrequency")+
    ggtitle(Groupings)
  return(plot.1)
}
Avatar.heatmap.ZsGreen <- function(df, Groupings, F.MHC_groups, C.MHC_groups){
  palette_greens <- colorRampPalette(colors = c("white", "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"))(12)
  
  plot.1 <- ggplot(subset(df, Group %in% Groupings), aes(x = APC, y = Sequence.ID, fill = mean.NFATZsGreen.pos))+
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) +
    coord_equal()+
    theme_classic()+
    theme(
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = unit(c(t=5.5, r=5.5, b=5.5, l=0), "pt")
    )+
    scale_fill_gradientn(colours = palette_greens, limits = c(0, 100), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    # theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")+
    # geom_vline(xintercept = seq(0.5, length(df$Sequence.ID), by = 1), color="gray", size=.5, alpha=.5) + # set vertical lines between x groups
    labs(x = "APC", fill = "% NFAT-ZsGreen\nPositive")+
    ggtitle(Groupings)
  # scale_y_continuous(expand = c(0,0), limits = c(0, 120), breaks=c(0,25,50,75, 100))
  return(plot.1)
}

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

Upper.lim <- my.max(joint.df.long_NFAT.pos$Frequency)
Lower.lim <- my.min(joint.df.long_NFAT.pos$Frequency)

# Single function which combines both heatmap plots for simplicity.
Avatar.heatmap <- function(df1, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim){
  # df1 for Repertoire Frequency
  # df2 for %NFAT-ZsGreen Positive
  df2 <- df1 %>% dplyr::filter(APC %in% c("tetra. bb&ff", "tetra. bb&g7g7") == FALSE)
  palette_purples <- colorRampPalette(colors = c("white", "#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d"))(12)
  palette_greens <- colorRampPalette(colors = c("white", "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"))(12)
  # Set plot 1 limits:

  # -----------------------------------------------------------------------------------------------------------------------
  plot.1 <- ggplot(subset(df1, Group %in% Groupings), aes(x = APC, y = Sequence.ID, fill = Frequency))+
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) +
    coord_equal()+
    theme_classic()+
    theme(
      legend.text=element_text(size=11),
      legend.title=element_text(size=13),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color="black", size = 11),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title=element_text(size=13),
      legend.position="left",
      plot.margin = unit(c(t=5.5, r=0, b=5.5, l=5.5), "pt")
    )+
    scale_fill_gradientn(
      colours = palette_purples,
      limits = c(Lower.lim, Upper.lim),
      trans = "log10",
      labels = scales::label_log(base = 10, digits = 3),
      breaks = c(1e-06, 1e-05, 1e-04, 1e-03),
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      )+
    labs(x = "Repertoire\nFreq. in:", y = "TCR Avatar #", fill = "Repertoire\nFrequency")

  # -----------------------------------------------------------------------------------------------------------------------
  
  plot.2 <- ggplot(subset(df2, Group %in% Groupings), aes(x = APC, y = Sequence.ID, fill = mean.NFATZsGreen.pos))+
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) +
    coord_equal()+
    theme_classic()+
    theme(
      legend.text=element_text(size=11),
      legend.title=element_text(size=13),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black", size = 11),
      axis.title.x = element_text(size=13),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      plot.title=element_text(size=14, face = "bold"),
      plot.margin = unit(c(t=5.5, r=5.5, b=5.5, l=0), "pt")
    )+
    scale_fill_gradientn(colours = palette_greens, limits = c(0, 100), guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    labs(x = "APC", fill = "% NFAT-ZsGreen\nPositive")
  
  heatmap.set <- plot.1 + plot.2
  return(heatmap.set)
}
F.MHC_groups <- c(bb, ff, g7g7, bxf, bxg7, tetra_bb.ff, tetra_bb.g7g7)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, bxf_color, bxg7_color, tetra_bb.ff_color, tetra_bb.g7g7_color)    
# ------------------- PLOTTING HEATMAPS -----------------------------------------------------------------------------------------------------------------------------------------------------------------
F.MHC_groups <- c(bb, ff, g7g7, bxf, bxg7)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, bxf_color, bxg7_color)    
Groupings <- "Control"
p1.Control <- Avatar.heatmap(F1.absent_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim)
p1.Control

F.MHC_groups <- c(bb, ff, g7g7, bxf, bxg7, tetra_bb.ff, tetra_bb.g7g7)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, bxf_color, bxg7_color, tetra_bb.ff_color, tetra_bb.g7g7_color)  
Groupings <- "in bb NOT ff or bxf"
p2.Tetra_present <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bf_tetra.Samples == "Present in bb&ff\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb, bb&ff NOT ff or bxf")
p2.Tetra_absent <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bf_tetra.Samples == "Absent in bb&ff\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb NOT bb&ff, ff or bxf")
p2.INS <- Avatar.heatmap(ISN_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb NOT ff and INS in bxf")

Groupings <- "in ff NOT bb or bxf"
p3.Tetra_present <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bf_tetra.Samples == "Present in bb&ff\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in ff, bb&ff NOT bb or bxf")
p3.Tetra_absent <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bf_tetra.Samples == "Absent in bb&ff\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in ff NOT bb&ff, bb or bxf")
p3.INS <- Avatar.heatmap(ISN_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in ff NOT bb and INS in bxf")

Groupings <- "in bb NOT g7g7 or bxg7"
p4.Tetra_present <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bg7_tetra.Samples == "Present in bb&g7g7\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb, bb&g7g7 NOT g7g7 or bxg7")
p4.Tetra_absent <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bg7_tetra.Samples == "Absent in bb&g7g7\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb NOT bb&g7g7, g7g7 or bxg7")
p4.INS <- Avatar.heatmap(ISN_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb NOT g7g7 and INS in bxg7")

Groupings <- "in g7g7 NOT bb or bxg7"
p5.Tetra_present <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bg7_tetra.Samples == "Present in bb&g7g7\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in g7g7, bb&g7g7 NOT bb or bxg7")
p5.Tetra_absent <- Avatar.heatmap(F1.absent_sub %>% dplyr::filter(bg7_tetra.Samples == "Absent in bb&g7g7\nTetraparental"), Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in g7g7 NOT bb&g7g7, bb or bxg7")
p5.INS <- Avatar.heatmap(ISN_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in g7g7 NOT bb and INS in bxg7")

Groupings <- "bb& ff NOT bxf"
p6 <- Avatar.heatmap(F1.absent_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb & ff NOT bxf") # + facet_grid(.~bf_tetra.Samples, scales = "free", space = "free")

Groupings <- "bb & g7g7 NOT bxg7"
p7 <- Avatar.heatmap(F1.absent_sub, Groupings, F.MHC_groups, C.MHC_groups, Upper.lim, Lower.lim) + ggtitle("in bb & g7g7 NOT bxg7") # + facet_grid(.~bg7_tetra.Samples, scales = "free", space = "free")

# -------------------------------- Save Plots ----------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/Control TCRs.pdf", width=7, height=1.5, useDingbats=FALSE)
p1.Control
dev.off()
# ----------------------------- Set the width of sets to 7.2mm for each plot (measure manually) ------------------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb, bb&ff NOT ff or bxf TCRs.pdf", width=8.0, height=5.2, useDingbats=FALSE)
p2.Tetra_present
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb NOT bb&ff, ff or bxf TCRs.pdf", width=7.0, height=2.7, useDingbats=FALSE)
p2.Tetra_absent
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb NOT ff and INS in bxf TCRs.pdf", width=5.2, height=4.48, useDingbats=FALSE)
p2.INS
dev.off()
# -----------------------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in ff, bb&ff NOT bb or bxf TCRs.pdf", width=6, height=3.6, useDingbats=FALSE)
p3.Tetra_present
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in ff NOT bb&ff, bb or bxf TCRs.pdf", width=7.6, height=6.25, useDingbats=FALSE)
p3.Tetra_absent
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in ff NOT bb and INS in bxf TCRs.pdf", width=7.6, height=2.52, useDingbats=FALSE)
p3.INS
dev.off()
# -----------------------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb, bb&g7g7 NOT g7g7 or bxg7 TCRs.pdf", width=6, height=3.3, useDingbats=FALSE)
p4.Tetra_present
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb NOT bb&g7g7, g7g7 or bxg7 TCRs.pdf", width=6, height=5.0, useDingbats=FALSE)
p4.Tetra_absent
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in bb NOT g7g7 and INS in bxg7 TCRs.pdf", width=7, height=2.78, useDingbats=FALSE)
p4.INS
dev.off()
# -----------------------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in g7g7, bb&g7g7 NOT bb or bxg7 TCRs.pdf", width=6, height=4.25, useDingbats=FALSE)
p5.Tetra_present
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in g7g7 NOT bb&g7g7 bb or bxg7 TCRs.pdf", width=6, height=4.35, useDingbats=FALSE)
p5.Tetra_absent
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/in g7g7 NOT bb and INS in bxg7 TCRs.pdf", width=7, height=2.8, useDingbats=FALSE)
p5.INS
dev.off()
# -----------------------------
pdf("Figure Export/5KC NFAT Reporter/Heatmaps/bb& ff NOT bxf TCRs.pdf", width=6, height=2.9, useDingbats=FALSE)
p6
dev.off()

pdf("Figure Export/5KC NFAT Reporter/Heatmaps/bb & g7g7 NOT bxg7 TCRs.pdf", width=6, height=1.57, useDingbats=FALSE)
p7
dev.off()
