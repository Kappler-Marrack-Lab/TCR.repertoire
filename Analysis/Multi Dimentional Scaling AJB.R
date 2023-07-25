library(tidyverse)
library(readxl)
library(ggpubr)
library(ggrepel)
library(Polychrome)
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
library(ggprism)
library(umap)
set.seed(42)
# -----------------------------------------------------
load(file = "TCR_fixed_beta_total_rep.CDR3AA.LFSR.ones.removed.RData") # filtered and ones removed

# exclude the tetraparental, b+/- sequences we are not interested in
genotypes <- repFilter(TCR_DoVb_CD4_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb","ff","g7g7","ss","bxf","bxg7","bxs","fxs"))) #subest for het and hom reps
# 
# genotypes <- list(
#   data = repFilter(TCR.CD4_DoVb_TR, .method = "by.meta", .query = list(Tetraparental_Cell_Source = exclude("bb", "ff", "g7g7")))$data,
#   meta = repFilter(TCR.CD4_DoVb_TR, .method = "by.meta", .query = list(Tetraparental_Cell_Source = exclude("bb", "ff", "g7g7")))$meta
# )
# Generate counts table
Rep_Freq <- pubRep(genotypes$data, "v+j+aa", .quant = c("prop"), .verbose = T) # frequency table
Rep_Counts <-  pubRep(genotypes$data, "v+j+aa", .verbose = T)
Rep_Counts[is.na(Rep_Counts)] = 0 #replace NA with zeros

# #fix data that have incorrect normalization
# counts_joined$count643 <- counts_joined$count643 * 100
# counts_joined$count619 <- counts_joined$count619 * 100
# counts_joined$count264 <- counts_joined$count264 /402034 * 10000

# subset to clones found at least 1 in 100,000 in at least 1 rep
min.cutoff <- Rep_Counts[] %>% apply(.,1,max) > 0.00001
# min.cutoff <- Rep_Counts[] %>% apply(.,1,max) > 0.1
Rep_Counts.filtered <- Rep_Counts[min.cutoff,]

# Inverse hyperbolic sine transformation (asinh) used here as an alternative to a log transform.
# The data has both zeros and a long right tail, since log isn't defined at zero and asinh is.
countsDists <- dist( t( asinh(Rep_Counts.filtered[,5:ncol(Rep_Counts.filtered)] ) ) )

# compute multidimentional scaling (mds)
mds <- countsDists %>% cmdscale() %>% as_tibble(rownames = "Sample", .name_repair = ~ c("Dim1", "Dim2"))
mds <- left_join(mds,genotypes$meta,by="Sample")
# ----------------------------------------------------------------------------------------------------------------------------------
# exclude the tetraparental, b+/- sequences we are not interested in
genotypes.CD8 <- repFilter(TCR_DoVb_CD8_AA.ones.removed.LFSR, .method = "by.meta", .query = list(MHC = include("bb","ff","g7g7","bxf","bxg7"))) #subest for het and hom reps

# Generate counts table
Rep_Freq.CD8 <- pubRep(genotypes.CD8$data, "v+j+aa", .quant = c("prop"), .verbose = T) # frequency table
Rep_Counts.CD8 <-  pubRep(genotypes.CD8$data, "v+j+aa", .verbose = T)
Rep_Counts.CD8[is.na(Rep_Counts.CD8)] = 0 #replace NA with zeros

# subset to clones found at least 1 in 100,000 in at least 1 rep
min.cutoff.CD8 <- Rep_Counts.CD8[] %>% apply(.,1,max) > 0.00001
# min.cutoff <- Rep_Counts[] %>% apply(.,1,max) > 0.1
Rep_Counts.filtered.CD8 <- Rep_Counts.CD8[min.cutoff.CD8,]

# Inverse hyperbolic sine transformation (asinh) used here as an alternative to a log transform.
# The data has both zeros and a long right tail, since log isn't defined at zero and asinh is.
countsDists.CD8 <- dist( t( asinh(Rep_Counts.filtered.CD8[,5:ncol(Rep_Counts.filtered.CD8)] ) ) )

# compute multidimentional scaling (mds)
mds.CD8 <- countsDists.CD8 %>% cmdscale() %>% as_tibble(rownames = "Sample", .name_repair = ~ c("Dim1", "Dim2"))
mds.CD8 <- left_join(mds.CD8, genotypes.CD8$meta,by="Sample")
# ----------------------------------------------------------------------------------------------------------------------------------

#  compute UMAP
# library(palmerpenguins)
# penguins <- penguins %>% 
#   drop_na() %>%
#   select(-year)%>%
#   mutate(ID=row_number())
# 
# # create a dataframe with all categorical variables with the unique row ID.
# penguins_meta <- penguins %>%
#   select(ID, species, island, sex)
# 
# # ============== EXAMPLE INPUT DATA ================
# penguins %>%
#   select(where(is.numeric)) %>%
#   column_to_rownames("ID") %>% scale() %>% head()

# ============== MY INPUT DATA ================
# rep.transpose <- t( asinh(Rep_Counts.filtered[,5:ncol(Rep_Counts.filtered)] ) ) #this makes a bad projection for a umap
rep.transpose <- t( scale(Rep_Counts.filtered[,5:ncol(Rep_Counts.filtered)] ) )
rep.transpose %>% as_tibble()
umap.rep <- rep.transpose %>% umap()

umap.rep_df <- umap.rep$layout %>%
  as_tibble(rownames = "Sample", .name_repair = ~ c("Dim1", "Dim2")) %>%
  left_join(genotypes$meta,by="Sample")

# ----------------------------------------------------------------------------------------------------------------------------------
# MHC
ff ="ff"
ss = "ss"
bb = "bb"
g7g7 = "g7g7"
bxf = "bxf"
fxs = "fxs"
bxs = "bxs"
bxg7 = "bxg7"

# Define MHC group colors:
ff_color ="#c5944e"
ss_color = "#fbff00"
bb_color = "#54a0fb" 
g7g7_color = "#9f104c"
bxf_color = "#808000"
fxs_color = "#fa5f00"
bxs_color = "#1afe02"
bxg7_color = "#f74ed6"


# Define the MHC groups and input elements for the current plot
F.MHC_groups <- c(bb, ff, g7g7, ss, bxf, bxg7, bxs, fxs)
C.MHC_groups <- c(bb_color, ff_color, g7g7_color, ss_color, bxf_color, bxg7_color, bxs_color, fxs_color)
# plot
p1 <- ggplot(mds, aes(x = Dim1, y = Dim2,label = MHC)) + 
  geom_point(aes(fill = factor(MHC)),
             shape=21,
             size = 4,
             alpha = 0.8,
             color = "black",
             position=position_jitter(h=5,w=5)) + 
  scale_fill_manual(breaks = F.MHC_groups, 
                     values = C.MHC_groups) +
  coord_cartesian(clip = "off")+
  theme_prism(border = TRUE)+
  labs(fill = "MHC")+   
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title=element_text(face="bold", size=16),
        legend.text=element_text(size=14),
        panel.background = element_rect(colour = "black", linewidth=1),
        aspect.ratio=1
        )

# plot UMAP
p2 <- ggplot(umap.rep_df, aes(x = Dim1, y = Dim2,label = MHC)) + 
  geom_point(aes(fill = factor(MHC)),
             shape=21,
             size = 4,
             alpha = 0.8,
             color = "black",
             position=position_jitter(h=5,w=5)) + 
  scale_fill_manual(breaks = F.MHC_groups, 
                    values = C.MHC_groups) +
  coord_cartesian(clip = "off")+
  theme_prism(border = TRUE)+
  labs(fill = "MHC")+   
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title=element_text(face="bold", size=16),
        legend.text=element_text(size=14),
        panel.background = element_rect(colour = "black", linewidth=1),
        aspect.ratio=1
  )


# write pdf
pdf("Figure Export/Clustering/Rep_Counts_mds.pdf", width = 5, height = 5)
p1
dev.off()

p3 <- ggplot(mds.CD8, aes(x = Dim1, y = Dim2,label = MHC)) + 
  geom_point(aes(fill = factor(MHC)),
             shape=21,
             size = 4,
             alpha = 0.8,
             color = "black",
             position=position_jitter(h=5,w=5)) + 
  scale_fill_manual(breaks = F.MHC_groups, 
                    values = C.MHC_groups) +
  coord_cartesian(clip = "off")+
  theme_prism(border = TRUE)+
  labs(fill = "MHC")+   
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title=element_text(face="bold", size=16),
        legend.text=element_text(size=14),
        panel.background = element_rect(colour = "black", linewidth=1),
        aspect.ratio=1
  )

pdf("Figure Export/Clustering/Rep_Counts_mds_CD8.pdf", width = 4, height = 4)
p3
dev.off()

# # ----------------------------------------------------------------------------------------------------------------------------------
# # Application of the Immune arch mds, tsne, clustering approaches to repertoire overlap data looks terrible.
# # tottaly unhelpful...
# genotypes.OV.MOR <- repOverlap(genotypes$data, .method = "morisita", .verbose = T)
# genotypes.OV.JAC <- repOverlap(genotypes$data, .method = "jaccard", .verbose = T)
# # Apply different analysis algorithms to the matrix of public clonotypes:
# # "mds" - Multi-dimensional Scaling
# genotypes.OV.MOR.mds <- repOverlapAnalysis(genotypes.OV.MOR, "mds")
# genotypes.OV.MOR.mds <- genotypes.OV.MOR.mds$points %>%  as_tibble(rownames = "Sample") %>% left_join(genotypes$meta,by="Sample")
# genotypes.OV.MOR.mds <- genotypes.OV.MOR.mds %>% rename(Dim1 = V1) %>% rename(Dim2 = V2)
# 
# genotypes.OV.MOR.tsne <- repOverlapAnalysis(genotypes.OV.MOR, "tsne")
# genotypes.OV.MOR.tsne <- genotypes.OV.MOR.tsne %>%  as_tibble(rownames = "Sample") %>% left_join(genotypes$meta,by="Sample")
# genotypes.OV.MOR.tsne <- genotypes.OV.MOR.tsne %>% rename(Dim1 = DimI) %>% rename(Dim2 = DimII)
# 
# ggplot(genotypes.OV.MOR.mds, aes(x = Dim1, y = Dim2,label = MHC, colour = MHC, fill = MHC, shape = Tetraparental_Cell_Source)) + 
#   geom_point(size = 4, alpha = 0.9) + 
#   # geom_text_repel() +
#   scale_fill_manual(breaks = F.MHC_groups, 
#                     values = C.MHC_groups) +
#   scale_color_manual(breaks = F.MHC_groups, 
#                      values = C.MHC_groups) +
#   scale_shape_manual(values = c(15,17,19)) + 
#   theme_bw() + 
#   theme(axis.line = element_line(color='black'),
#         plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank())
# 
# repOverlapAnalysis(genotypes.OV.MOR, "tsne") %>% vis()
# ggplot(genotypes.OV.MOR.tsne, aes(x = Dim1, y = Dim2,label = MHC, colour = MHC, fill = MHC, shape = Tetraparental_Cell_Source)) + 
#   geom_point(size = 4, alpha = 0.9) + 
#   # geom_text_repel() +
#   scale_fill_manual(breaks = F.MHC_groups, 
#                     values = C.MHC_groups) +
#   scale_color_manual(breaks = F.MHC_groups, 
#                      values = C.MHC_groups) +
#   scale_shape_manual(values = c(15,17,19)) + 
#   theme_bw() + 
#   theme(axis.line = element_line(color='black'),
#         plot.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank())
# 
# 
# # Attempt at a UMAP projection of the bb composite repertoire
# library(readr)
# library(umap)
# custom.config = umap.defaults
# custom.config$random_state = 123
# CD4_Counts_table <- read_csv("Clones Counts Table CD4.csv")
# 
# # pull in the CD4 bb and bxf counts table and pull out only  values which are here.
# bb_bxf_select <- CD4_Counts_table %>%
#   mutate(Clonotype = paste(V.name, CDR3.aa, J.name, sep = "|")) %>%
#   select(Clonotype, bb.ave.prop, bb.sum, bb.Samples, bxf.ave.prop, bxf.sum, bxf.Samples) %>%
#   replace(is.na(.), 0) %>%
#   filter(!(bb.ave.prop == 0 & bxf.ave.prop == 0))
# # Label the values which are unique to bb, bxf and common to both.
# bb.bxf.shared <- bb_bxf_select %>% filter(!(bb.sum == 0)) %>% filter(!(bxf.sum == 0)) %>% dplyr::mutate(group = "shared")
# 
# # find the values which contain zero hits in either and pull these out to not confuse them with the values which are shared between repertoires
# CDR3_pairs_select_zeros_only <- bb_bxf_select %>% filter(bb.sum == 0 | bxf.sum == 0)
# bb.unique <- CDR3_pairs_select_zeros_only %>% filter(bb.sum > 0) %>% dplyr::mutate(group = "bb unique")
# bxf.unique <- CDR3_pairs_select_zeros_only %>% filter(bxf.sum > 0) %>% dplyr::mutate(group = "bxf unique")
# bb_ff.bound <- rbind(bb.bxf.shared, bb.unique, bxf.unique) %>% dplyr::mutate(group = as.factor(group))
# Rep_Counts_data <- bb_ff.bound %>% select(bb.ave.prop, bb.sum, bb.Samples, bxf.ave.prop, bxf.sum, bxf.Samples)
# Rep_Counts_data_min1 <- bb_ff.bound %>% select(bb.sum, bb.Samples, bxf.sum, bxf.Samples)
# Rep_Counts_labels <- bb_ff.bound[, "group"]
# Rep_Counts.umap = umap(Rep_Counts_data)
# Rep_Counts.min1.umap = umap(Rep_Counts_data_min1)
# # Rep_Counts.umap_all = Rep_Counts.umap
# umap.df <- data.frame(x.umap = Rep_Counts.min1.umap$layout[,1],
#                  y.umap = Rep_Counts.min1.umap$layout[,2],
#                  group = Rep_Counts_labels)
# bb_ff.bound_new <- bb_ff.bound %>% mutate(x.umap = umap.df$x.umap) %>% mutate(y.umap = umap.df$y.umap)
# ggplot(umap.df, aes(x.umap, y.umap, colour = group)) +
#   geom_point() +
#   theme_minimal()
# 
# # Try to filter out values which are a total count less than 50 in the sum column.
# # These seem to to be the ones which separate oddly.
# bb_ff.bound.HF <- bb_ff.bound %>% filter(!(bb.sum > 50)) %>% filter(!(bxf.sum > 50))
# bb_ff.bound.HF.data <- bb_ff.bound.HF %>% select(bb.sum, bb.Samples, bxf.sum, bxf.Samples)
# bb_ff.bound.HF.labels <- bb_ff.bound.HF[, "group"]
# bb_ff.bound.HF.umap = umap(bb_ff.bound.HF.data)
# 
# bb_ff.bound.HF <- bb_ff.bound.HF %>%
#   mutate(x.umap = bb_ff.bound.HF.umap$layout[,1]) %>%
#   mutate(y.umap = bb_ff.bound.HF.umap$layout[,2])
# ggplot(bb_ff.bound.HF, aes(x.umap, y.umap, colour = group)) +
#   geom_point() +
#   theme_minimal()
# # This looks better but its still looks really wierd... Basically not as useful as an euler plot.
# 
# 
# # plot.iris(Rep_Counts.umap, Rep_Counts_labels)
# ## This is not working. Inroduces catagorical variables which are not allowed. 
# # bb.select <- bb_bxf_select %>%
# #   select(Clonotype, bb.ave.prop, bb.sum, bb.Samples) %>%
# #   mutate(MHC = "bb") %>%
# #   rename(ave.prop = bb.ave.prop) %>% rename(sum = bb.sum) %>% rename(samples = bb.Samples)
# # bxf.select <- bb_bxf_select %>%
# #   select(Clonotype, bxf.ave.prop, bxf.sum, bxf.Samples) %>%
# #   mutate(MHC = "bxf") %>%
# #   rename(ave.prop = bxf.ave.prop) %>% rename(sum = bxf.sum) %>% rename(samples = bxf.Samples)
# # bb_ff.bound <- rbind(bb.select, bxf.select)
# # bb_ff.bound$MHC <- as.factor(bb_ff.bound$MHC)
# # # Select out the data and labels columns
# # Rep_Counts_data <- bb_ff.bound %>% select(Clonotype, ave.prop, sum, samples)
# # Rep_Counts_labels <- bb_ff.bound[, "MHC"]
# 
# 
# 
# 
# 
# 
# # Try to pivot each set of values (prop, sum, samples) and have just bb and bxf become the names for labels.
# # The Clonotype information will not really be use but we can keep it until the end.
# 
# 
# #   rename(bb = bb.ave.prop) %>% rename(bxf = bxf.ave.prop) %>%
# #   pivot_longer(cols = c("bb", "bxf"), names_to = "MHC", values_to = "ave.prop")
# # bb_bxf_select$MHC = as.factor(bb_bxf_select$MHC)
# # Rep_Counts_data = bb_bxf_select[, "ave.prop"]
# # Rep_Counts_labels = bb_bxf_select$MHC
# 
# # # Attempt to make a UMAP here...
# # Rep.umap = umap(Rep_Counts_data)
# # 
# # # Test UMAPS
# # install.packages("remotes")
# # remotes::install_github("jlmelville/snedata")
# # library(snedata)
# # fashion <- download_fashion_mnist()
# # mnist <- download_mnist()
# # # install.packages("devtools")
# # devtools::install_github("jlmelville/coil20")
# # library(coil20)
# # coil20 <- download_coil20(verbose = TRUE)
# 
# 
# iris.data = iris[, grep("Sepal|Petal", colnames(iris))]
# iris.labels = iris[, "Species"]