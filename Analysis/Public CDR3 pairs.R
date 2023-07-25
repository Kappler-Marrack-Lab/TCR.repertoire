library(immunarch)
library(ggrepel)
library(vroom)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(GGally)
library(svglite)
require(tidyr)
library(R.utils)
library(viridis)
library(RColorBrewer)
library(hablar)
library(latex2exp)
library(scico)
library(scales) # to access break formatting functions
require(scales)

theme_set(theme_classic())
set.seed(42)


load("TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData")
load("counts tables filtered.LFSR.ones.removed.RData")
# ______

CDR3_pairs <- pr_CD4s_filtered.LFSR.ones.removed %>% select(
  bb.Samples, 
  ff.Samples,
  g7g7.Samples,
  ss.Samples,
  bxf.Samples,
  bxg7.Samples,
  fxs.Samples,
  bxs.Samples
)

# # --------- Pairs plot from base R NOT IN USE -----------------
# CDR3_pairs_log <- CDR3_pairs %>% replace(is.na(.), 0) %>% log10() %>% rationalize() %>% replace(is.na(.), 0) # log10 transform data and then remove NA and -Inf values and replace with 0
# 
# # Defined upper pannel to simplify plot
# upper.panel<-function(x, y){
#   points(x,y, pch = ifelse(x == 0 | y == 0, 4, 19), col = ifelse(x == 0 | y == 0, rgb(1, 0, 0, 0.05), rgb(0, 0, 0, 0.05))) #values with 0 (x,y) are set to red, can also change the shape
#   R <- round(cor(x, y), digits=2) # True R squared which includes zero values on periphery
#   x[x==0] <- min(x)
#   y[y==0] <- min(y)
#   Rom <- round(cor(x, y), digits=2) # Rom  = 'R^2 ommited', which excludes zero values highlighted in red
#   txt <- paste0("R = ", R, "\n",
#                 "Rom = ", Rom) # can't figure out the math type
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   text(0.6, 0.8, txt)
# }
# 
# 
# CDR3.pairs_plot <- pairs(CDR3_pairs_log[,1:8],
#                        lower.panel = NULL,
#                        upper.panel = upper.panel
# )
# 
# 
# pdf("Figure Export/Public Repertoire/V.J.CDR3 pairs.pdf", width = 10, height = 10) # Open a pdf file #
# # 2. Create a plot
# CDR3.pairs_plot <- pairs(CDR3_pairs_log[,1:8],
#                        lower.panel = NULL,
#                        upper.panel = upper.panel
# )
# dev.off() # Close the pdf file
# 
# png("Figure Export/Public Repertoire/V.J.CDR3 pairs.png", width = 10, height = 10) # Open a pdf file #
# # 2. Create a plot
# CDR3.pairs_plot <- pairs(CDR3_pairs_log[,1:8],
#                          lower.panel = NULL,
#                          upper.panel = upper.panel
# )
# dev.off() # Close the png file


# ---------  ggPairs plots (2022-06-16) -----------------------------------------
# # # Remove values which are zeros for both the haplotypes
# CDR3_pairs_select <- CDR3_pairs %>% replace(is.na(.), 0) %>% select(bb.ave.prop, bxf.ave.prop) %>% filter(!(bb.ave.prop == 0 & bxf.ave.prop == 0))
# 
# # create separate data where zeros are filtered out
# CDR3_pairs_select_zeros_dropped <- CDR3_pairs_select %>% filter(!(bb.ave.prop == 0)) %>% filter(!(bxf.ave.prop == 0))
# CDR3_pairs_select_zeros_only <- CDR3_pairs_select %>% filter(bb.ave.prop == 0 | bxf.ave.prop == 0)
# 
# # create a plot for the % of clones in each part of the diagram
# set_total <- nrow(CDR3_pairs_select)
# shared <- nrow(CDR3_pairs_select_zeros_dropped) #clones intersection
# x_only <- sum(CDR3_pairs_select_zeros_only$bb.ave.prop > 0) #clones only on X axis
# y_only <- sum(CDR3_pairs_select_zeros_only$bxf.ave.prop > 0) #clones only on Y axis
# 
# set_sharing <- data.frame(group = c("shared","bb unique", "bxf unique"),
#                                  counts = c(shared, x_only, y_only),
#                                  percentage = c(signif((shared/set_total)*100, 2), signif((x_only/set_total)*100, 2), signif((y_only/set_total)*100, 2)))
# 
# 
# C.MHC_groups <- c(bb_color, bxf_color, "grey")
# 
# sets_plot <- ggplot(set_sharing, aes(x = "", y = counts, fill = group)) +     
#   geom_bar(position = position_stack(), stat = "identity", width = .7) +
#   scale_fill_manual(values=C.MHC_groups) +
#   geom_text(aes(label = paste0(percentage,"%")), position = position_stack(vjust = 0.5)) +
#   coord_flip() +
#   scale_x_discrete(expand=c(0,0)) +
#   # scale_y_discrete(expand=c(0, 0)) +
#   theme(aspect.ratio=0.075) +
#   labs(
#     x = "",
#     y = "# Unique Sequences",
#     fill = "Group")
# 
# # Set zero values to some value lower than everything else in the data.
# Xbreak <- min(CDR3_pairs_select_zeros_dropped$bb.ave.prop) - 0.000002 # set x-break here
# Ybreak <- min(CDR3_pairs_select_zeros_dropped$bxf.ave.prop) - 0.000002 #set y-break here
# 
# # Re-assign this new minimal value to all 'zeros'.
# CDR3_pairs_select_zeros_only['bb.ave.prop'][CDR3_pairs_select_zeros_only['bb.ave.prop'] == 0] <- Xbreak
# CDR3_pairs_select_zeros_only['bxf.ave.prop'][CDR3_pairs_select_zeros_only['bxf.ave.prop'] == 0] <- Ybreak
# 
# # Set a facet variable (used for R^2 value)
# CDR3_pairs_select_combine <- rbind(CDR3_pairs_select_zeros_dropped, CDR3_pairs_select_zeros_only)
# CDR3_pairs_select_combine$facet <- ifelse(CDR3_pairs_select_combine$bb.ave.prop == min(CDR3_pairs_select_combine$bb.ave.prop), 1,
#                                           ifelse(CDR3_pairs_select_combine$bxf.ave.prop == min(CDR3_pairs_select_combine$bxf.ave.prop), 3, 2))
# # facet 2 =  main data (shared clones), 3 = Zeros on X axis (clones only in group X), 2 = Zeros on Y axis (clones only in group Y)
# 
# # 2D Histogram
# my_breaks = c(0, 0.1e2, 0.5e2, 0.1e3, 0.5e3, 0.1e4, 0.5e4, 0.1e5, 0.5e5, 1.0e5, 5.0e5)
# getPalette = colorRampPalette(scico(length(my_breaks), palette = 'lajolla'))
# getPalette = colorRampPalette(brewer.pal(9, "BuPu"))
# 
# 
# Density_plot <- ggplot(CDR3_pairs_select_combine, aes(x = bb.ave.prop, y = bxf.ave.prop)) +
#   # Plotting
#   geom_density_2d_filled(
#     breaks = my_breaks,
#     contour_var = "count") +
#   
#   # Set Colors
#   # scale_fill_viridis_d(option = "rocket") +
#   # scico::scale_fill_scico(palette = "tofino") +
#   # scale_fill_brewer(palette = "Spectral") +
#   scale_fill_manual(values = getPalette(length(my_breaks))) +
# 
# 
#   # Axis Modification
#   scale_x_log10(breaks = c(Xbreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1)),
#                 expand = c(Xbreak, Xbreak)) +
#   scale_y_log10(breaks = c(Ybreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1)),
#                 expand = c(Ybreak, Ybreak)) +
#   annotation_logticks() +
#   
#   # bounding lines
#   geom_vline(xintercept = Xbreak+0.000001, linetype = "solid", size = 1.5, color = "white") +
#   geom_hline(yintercept = Ybreak+0.000001, linetype = "solid", size = 1.5, color = "white") +
#   geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
#   
#   # Stat line for just the shared clones
#   stat_poly_line(data = subset(CDR3_pairs_select_combine, facet == 2), color = "Red") +
#   stat_poly_eq(data = subset(CDR3_pairs_select_combine, facet == 2), label.x = "right", label.y = "top", color = "Red") +
#   
#   
#   # Theme
#   theme_bw() +
#   theme(aspect.ratio=1)+
#   labs(
#     x = "bb Sequence Frequency",
#     y = "bxf Sequence Frequency")
# 
# 
# sets_plot/Density_plot
#----------------------------------------------------------------------------------
# Co-receptor
CD4 = 'CD4'
CD8 = 'CD8'
DP = 'DP'

# MHC
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"
IAbHet = "b+/-"
b.f_bb = "b+f_bb"
b.f_ff = "b+f_ff"
b.g7_bb = "b+g7_bb"
b.g7_g7g7 = "b+g7_g7g7"

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
b.f_bb_color = "#5c54fb"
b.f_ff_color = "#c5584e"
b.g7_bb_color = "#54fbb5"
b.g7_g7g7_color = "#fb961a"

# Define the MHC groups and input elements for the current plot
F.MHC_groups <- c(ff, ss, bb, g7g7)
C.MHC_groups <- c(ff_color, ss_color, bb_color, g7g7_color)

CDR3_pairs_copy <- pr_CD4s_filtered.LFSR.ones.removed %>%
  select(
    bb.ave.prop, ff.ave.prop, g7g7.ave.prop, ss.ave.prop, bxf.ave.prop, bxg7.ave.prop, fxs.ave.prop, bxs.ave.prop
  )

#----------------------------------------------------------------------------------
# Create a function to replicate the above ggplot for desired pairs
names(CDR3_pairs_copy) <- c("bb", "ff", "g7g7", "ss", "bxf", "bxg7", "fxs", "bxs")
# input to the function: CDR3_pairs

F.MHC_groups <- c(bb, bxf)
subset(CDR3_pairs_copy, select=F.MHC_groups)
# subset(CDR3_pairs_copy, MHC %in% F.MHC_groups)

ggCloneCor <- function(counts.table, group_1, group_2, C.MHC_groups){
  # Remove values which are zeros for both the haplotypes

  F.MHC_groups = c(group_1, group_2)
  CDR3_pairs_select <- counts.table %>% replace(is.na(.), 0) %>% subset(select=F.MHC_groups) %>% filter(!(.[[1]] == 0 & .[[2]] == 0))
  CDR3_pairs_select_zeros_dropped <- CDR3_pairs_select %>% filter(!(.[[1]] == 0)) %>% filter(!(.[[2]] == 0))
  CDR3_pairs_select_zeros_only <- CDR3_pairs_select %>% filter(.[[1]] == 0 | .[[2]] == 0)
  
  
  # create a plot for the % of clones in each part of the diagram
  set_total <- nrow(CDR3_pairs_select)
  shared <- nrow(CDR3_pairs_select_zeros_dropped) #clones intersection
  x_only <- sum(CDR3_pairs_select_zeros_only[[1]] > 0) #clones only on X axis
  y_only <- sum(CDR3_pairs_select_zeros_only[[2]] > 0) #clones only on Y axis

  x_name <- paste(deparse(substitute(group_1)), "unique")
  y_name <- paste(deparse(substitute(group_2)), "unique")
  set_sharing <- data.frame(group = factor(c("shared", x_name, y_name),
                                           levels = c(x_name, y_name, "shared")),
                            counts = c(shared, x_only, y_only),
                            percentage = c(signif((shared/set_total)*100, 2), signif((x_only/set_total)*100, 2), signif((y_only/set_total)*100, 2)),
                            color = c("grey", C.MHC_groups)) %>% arrange(group)
  
  sets_plot <- ggplot(set_sharing, aes(x = group, y = counts, fill = group)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=set_sharing$color) +
    # geom_text(aes(label = paste0(percentage,"%")), position = position_stack(vjust = 0.5)) +
    # coord_flip() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0, 0)) + #this is causing the y values to not be displayed
    # guides(fill = guide_legend(reverse = TRUE))+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          legend.position="none", legend.direction="vertical", aspect.ratio=6)+
    labs(
      x = "",
      y = "# Unique Sequences",
      fill = "Group")
  
  
  # Original sets plot
  # sets_plot <- ggplot(set_sharing, aes(x = "", y = counts, fill = group)) +
  #   geom_bar(position = position_stack(), stat = "identity", width = .7) +
  #   scale_fill_manual(values=set_sharing$color) +
  #   geom_text(aes(label = paste0(percentage,"%")), position = position_stack(vjust = 0.5)) +
  #   coord_flip() +
  #   scale_x_discrete(expand=c(0,0)) +
  #   # scale_y_discrete(expand=c(0, 0)) +
  #   theme(aspect.ratio=0.075) +
  #   theme(legend.position="bottom") +
  #   guides(fill = guide_legend(reverse = TRUE))+
  #   labs(
  #     x = "",
  #     y = "# Unique Sequences",
  #     fill = "Group")
  
  

  # # This works for most of the plots.
  # Xbreak <- 3e-07
  # Ybreak <- 3e-07
  # split_line <- 1e-06
  
  # Reassign to fix some plots which look a little off:
  Xbreak <- 2e-07
  Ybreak <- 2e-07
  split_line <- 7e-07
  
  
  # Re-assign this new minimal value to all 'zeros'.
  CDR3_pairs_select_zeros_only[[1]][CDR3_pairs_select_zeros_only[[1]] == 0] <- Xbreak
  CDR3_pairs_select_zeros_only[[2]][CDR3_pairs_select_zeros_only[[2]] == 0] <- Ybreak
  
  # Set a facet variable (used for R^2 value)
  CDR3_pairs_select_combine <- rbind(CDR3_pairs_select_zeros_dropped, CDR3_pairs_select_zeros_only)
  CDR3_pairs_select_combine$facet <- ifelse(CDR3_pairs_select_combine[[1]] == min(CDR3_pairs_select_combine[[1]]), 1,
                                            ifelse(CDR3_pairs_select_combine[[2]] == min(CDR3_pairs_select_combine[[2]]), 3, 2))
  # facet 2 =  main data (shared clones), 3 = Zeros on X axis (clones only in group X), 2 = Zeros on Y axis (clones only in group Y)

  # 2D Histogram
  # my_breaks = c(0, 0.1e2, 0.5e2, 0.1e3, 0.5e3, 0.1e4, 0.5e4, 0.1e5, 0.5e5, 1.0e5, 5.0e5) # breaks for discrete Sequence count
  my_breaks = c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1) #breaks for a 0-1 density 
  # getPalette = colorRampPalette(scico(length(my_breaks), palette = 'lajolla'))
  getPalette = colorRampPalette(brewer.pal(9, "BuPu"))
  

  Density_plot <- ggplot(CDR3_pairs_select_combine, aes_string(x = names(CDR3_pairs_select_combine)[1], y = names(CDR3_pairs_select_combine)[2])) +
    # Plotting
    
    # # THis gives the actual counts which appears to not be tottaly accurate or something is wrong with it...
    # geom_density_2d_filled(
    #   breaks = my_breaks,
    #   contour_var = "count") +
    
    geom_density_2d_filled(contour_var = "ndensity", breaks = my_breaks)+
    

    # Set Colors
    # scale_fill_viridis_d(option = "rocket") +
    # scico::scale_fill_scico(palette = "tofino") +
    # scale_fill_brewer(palette = "Spectral") +
    scale_fill_manual(values = getPalette(length(my_breaks))) +


    # Axis Modification
    scale_x_log10(breaks = c(Xbreak, 10^(-6:-1)),
                  labels = c(0, math_format()(-6:-1)),
                  expand = c(Xbreak, Xbreak)) +
    scale_y_log10(breaks = c(Ybreak, 10^(-6:-1)),
                  labels = c(0, math_format()(-6:-1)),
                  expand = c(Ybreak, Ybreak)) +
    annotation_logticks() +

    
    geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
    
    # Stat line for just the shared clones
    stat_poly_line(data = subset(CDR3_pairs_select_combine, facet == 2), color = "Red") +
    stat_poly_eq(data = subset(CDR3_pairs_select_combine, facet == 2), label.x = "right", label.y = "top", color = "Red", size = 6) +
    
    # Add boxes which defined a color:
    annotate("rect", xmin=0, xmax=split_line, ymin=0, ymax=split_line, fill="white")+ #bottom corner
    
    # identification of groups by color lines
    annotate(geom = 'segment', x = Inf, xend = Inf, y = 0, yend = split_line,  color = C.MHC_groups[1], size = 7) + #x axis non-shared
    annotate(geom = 'segment', x = 0, xend = split_line, y = Inf, yend = Inf,  color = C.MHC_groups[2], size = 7) + #y axis non-shared
    annotate(geom = 'segment', x = split_line, xend = Inf, y = Inf, yend = Inf,  color = "grey", size = 7) + #shared
    annotate(geom = 'segment', x = Inf, xend = Inf, y = split_line, yend = Inf,  color = "grey", size = 7) + #shared
    
    # bounding lines
    geom_vline(xintercept = split_line, linetype = "solid", size = 4, color = "white") +
    geom_hline(yintercept = split_line, linetype = "solid", size = 4, color = "white") +
    
    # Theme
    theme_bw() +
    theme(aspect.ratio=1) +
    theme(legend.position="none") + #this will hide the levels color coding.
    labs(
      x = paste(deparse(substitute(group_1)), "Sequence Frequency"),
      y = paste(deparse(substitute(group_2)), "Sequence Frequency"))
  
  # Organize Layout
  # legend <- cowplot::get_legend(sets_plot)
  # sets_plot <- sets_plot + theme(legend.position="none")
  # pairs_plot <- Density_plot | (legend / sets_plot)
  # 
  # layout <- "
  # ##BBBB
  # AACCDD
  # ##CCDD
  # "
  # 
  # plot_spacer()
  # p1 + p2 + p3 + p4 + 
  #   plot_layout(design = layout)
  # 
  pairs_plot <- Density_plot+sets_plot
  
  # return(suppressWarnings(print(pairs_plot)))
  return(pairs_plot)
}

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs bxf_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, bxf, C.MHC_groups <- c(bb_color, bxf_color))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ff vs bxf_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, ff, bxf, C.MHC_groups <- c(ff_color, bxf_color))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs ff_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, ff, C.MHC_groups <- c(bb_color, ff_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs bxg7_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, bxg7, C.MHC_groups <- c(bb_color, bxg7_color)) # issues here something is breaking in the function. Something is causing it to not assign the zeros to a new value...
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J g7g7 vs bxg7_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, g7g7, bxg7, C.MHC_groups <- c(g7g7_color, bxg7_color))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs g7g7_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, g7g7, C.MHC_groups <- c(bb_color, g7g7_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs bxs_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, bxs, C.MHC_groups <- c(bb_color, bxs_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ss vs bxs_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, ss, bxs, C.MHC_groups <- c(ss_color, bxs_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb vs ss_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, bb, ss, C.MHC_groups <- c(bb_color, ss_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ff vs fxs_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, ff, fxs, C.MHC_groups <- c(ff_color, fxs_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ss vs fxs_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, ss, fxs, C.MHC_groups <- c(ss_color, fxs_color))
dev.off()

pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ff vs ss_V3.pdf", width = 4, height = 4)
ggCloneCor(CDR3_pairs_copy, ff, ss, C.MHC_groups <- c(ff_color, ss_color))
dev.off()


# saveRDS(ggCloneCor, "ggCloneCor.rds")
# Identify V.CDR3.J pairs based on the same haplotype
# also!
# Combine average frequency of multiple repertoires of the same haplotype together and compare the repertoires then.
# Since the other comparisons are of multiple repertoires you might expect multiple animals to increase the shared number of clones.

# n = 6 for CD4 bb and CD4 ff unique repertoire runs
multi_title_CD4 <- function(df, hap, n) {
  rename(df, "{hap}.{n}" := TCR_DoVb_CD4_AA.ones.removed.LFSR$meta %>% filter(MHC == hap) %>% pull(Sample) %>% nth(n))
}
pr_CD4_merge_prop <- pubRep(TCR_DoVb_CD4_AA.ones.removed.LFSR$data, "v+j+aa", .quant = c("prop"), .verbose = T)


V.CDR3.J_sample_replicates_dummy.names <- pr_CD4_merge_prop %>%
  multi_title_CD4(bb, 1) %>% multi_title_CD4(bb, 2) %>% multi_title_CD4(bb, 3) %>% multi_title_CD4(bb, 4) %>% multi_title_CD4(bb, 5) %>%
  multi_title_CD4(ff, 1) %>% multi_title_CD4(ff, 2) %>% multi_title_CD4(ff, 3) %>% multi_title_CD4(ff, 4) %>% multi_title_CD4(ff, 5) %>% multi_title_CD4(ff, 6) %>%
  multi_title_CD4(g7g7, 1) %>% multi_title_CD4(g7g7, 2) %>% multi_title_CD4(g7g7, 3) %>% multi_title_CD4(g7g7, 4) %>% multi_title_CD4(g7g7, 5) %>% multi_title_CD4(g7g7, 6) %>%
  multi_title_CD4(ss, 1) %>% multi_title_CD4(ss, 2) %>% multi_title_CD4(ss, 3) %>% multi_title_CD4(ss, 4) %>%
  
  select(
         bb.1, bb.2, bb.3, bb.4, bb.5,
         ff.1, ff.2, ff.3, ff.4, ff.5, ff.6,
         g7g7.1, g7g7.2, g7g7.3, g7g7.4, g7g7.5, g7g7.6,
         ss.1, ss.2, ss.3, ss.4
         ) %>%
  as_tibble() %>%
  dplyr::mutate(bb.25 = rowSums(across(c(bb.2, bb.5)), na.rm = TRUE)/2) %>%
  dplyr::mutate(bb.134 = rowSums(across(c(bb.1, bb.3, bb.4)), na.rm = TRUE)/3) %>%
  dplyr::mutate(ff.135 = rowSums(across(c(ff.1, ff.3, ff.5)), na.rm = TRUE)/3) %>%
  dplyr::mutate(ff.246 = rowSums(across(c(ff.2, ff.4, ff.6)), na.rm = TRUE)/3) %>%
  dplyr::mutate(g7g7.135 = rowSums(across(c(g7g7.1, g7g7.3, g7g7.5)), na.rm = TRUE)/3) %>%
  dplyr::mutate(g7g7.246 = rowSums(across(c(g7g7.2, g7g7.4, g7g7.6)), na.rm = TRUE)/3) %>%
  dplyr::mutate(ss.12 = rowSums(across(c(ss.1, ss.2)), na.rm = TRUE)/2) %>%
  dplyr::mutate(ss.34 = rowSums(across(c(ss.3, ss.4)), na.rm = TRUE)/2)


V.CDR3.J_sample_replicates_dummy.names[V.CDR3.J_sample_replicates_dummy.names == 0] <- NA


pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J bb.25 vs bb.134_V2.pdf", width = 4, height = 4)
ggCloneCor(V.CDR3.J_sample_replicates_dummy.names, group_1 = "bb.25", group_2 = "bb.134", C.MHC_groups <- c("#5b54fb", "#54f4fb"))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ff.135 vs ff.246_V3.pdf", width = 4, height = 4)
ggCloneCor(V.CDR3.J_sample_replicates_dummy.names, group_1 = "ff.135", group_2 = "ff.246", C.MHC_groups <- c("#c08657", "#b6a161"))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J g7g7.135 vs g7g7.246_V3.pdf", width = 4, height = 4)
ggCloneCor(V.CDR3.J_sample_replicates_dummy.names, group_1 = "g7g7.135", group_2 = "g7g7.246", C.MHC_groups <- c("#6f0b67", "#BB5781"))
dev.off()
pdf(file = "Figure Export/Public Repertoire/2D Density/V.CDR3.J/gg V.CDR3.J ss.12 vs ss.34_V2.pdf", width = 4, height = 4)
ggCloneCor(V.CDR3.J_sample_replicates_dummy.names, group_1 = "ss.12", group_2 = "ss.34", C.MHC_groups <- c("#B0B300", "#FDF193"))
dev.off()





# # EXTRA TEST PLOTS
# # PLOT OF GEOM POINT 2 with zero values present
# ggplot(CDR3_pairs_select_combine, aes(bb.ave.prop, bxf.ave.prop, alpha = 0.05)) +
#   
#   # Plotting
#   geom_point(shape = 16, size = 5, show.legend = FALSE) +
#   
#   # Axis Modification
#   scale_x_log10(breaks = c(Xbreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1))) +
#   scale_y_log10(breaks = c(Ybreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1))) +
#   annotation_logticks() +
#   
#   # bounding lines
#   geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
#   geom_vline(xintercept = Xbreak+0.000001, linetype = "dashed", size = 0.5) +
#   geom_hline(yintercept = Ybreak+0.000001, linetype = "dashed", size = 0.5) +
#   
#   # Theme
#   theme_bw() +
#   theme(aspect.ratio=1)
# 
# 
# # 2D Density plot
# ggplot(CDR3_pairs_select_combine, aes(bb.ave.prop, bxf.ave.prop)) +
#   
#   # Plotting
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#   geom_point(shape = 16, size = 5, alpha = 0.01, show.legend = FALSE) +
#   scico::scale_fill_scico(palette = "bilbao") +
#   # scico::scale_fill_scico(palette = "lajolla") +
#   # scale_fill_viridis_c(option = "magma") +
#   # scale_fill_distiller(palette= "Spectral", direction=1) +
#   
#   # Axis Modification
#   scale_x_log10(breaks = c(Xbreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1))) +
#   scale_y_log10(breaks = c(Ybreak, 10^(-6:-1)),
#                 labels = c(0, math_format()(-6:-1))) +
#   annotation_logticks() +
#   
#   # bounding lines
#   geom_abline(intercept = 0, slope = 1, size = 0.5, linetype = "dashed") +
#   geom_vline(xintercept = Xbreak+0.000001, linetype = "dashed", size = 0.5) +
#   geom_hline(yintercept = Ybreak+0.000001, linetype = "dashed", size = 0.5) +
#   
#   # Theme
#   theme_bw() +
#   theme(aspect.ratio=1)

