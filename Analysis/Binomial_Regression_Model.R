library(tidyverse)
library(rstatix)
library(immunarch)
library(emmeans)
library(ggrastr)
options(ggrastr.default.dpi=900)

load("V.CDR3.J Counts Matrix.rData")

pr_full_CD4_mnml%>%head()

# Evaluate binomial regression seprately for each parrent.
# We want to first identify the fraction of the parrents which bare a given sequence. Evaluate this for each parental haplotype.
# filter for sequences which are NOT zero in both haplotypes of interest.
pr_full_CD4_mnml %>% names()
bb.bxf_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, bb.Samples, bxf.Samples) %>% mutate(merge = paste0(bb.Samples, bxf.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(bb.frac = bb.Samples/max(bb.Samples)) %>% mutate(bxf.Samples.successes = bxf.Samples) %>% mutate(bxf.Samples.failures = max(bxf.Samples)-bxf.Samples) %>% mutate(bxf.frac = bxf.Samples/max(bxf.Samples))
ff.bxf_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, ff.Samples, bxf.Samples) %>% mutate(merge = paste0(ff.Samples, bxf.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(ff.frac = ff.Samples/max(ff.Samples)) %>% mutate(bxf.Samples.successes = bxf.Samples) %>% mutate(bxf.Samples.failures = max(bxf.Samples)-bxf.Samples) %>% mutate(bxf.frac = bxf.Samples/max(bxf.Samples))

bb.bxg7_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, bb.Samples, bxg7.Samples) %>% mutate(merge = paste0(bb.Samples, bxg7.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(bb.frac = bb.Samples/max(bb.Samples)) %>% mutate(bxg7.Samples.successes = bxg7.Samples) %>% mutate(bxg7.Samples.failures = max(bxg7.Samples)-bxg7.Samples) %>% mutate(bxg7.frac = bxg7.Samples/max(bxg7.Samples))
g7g7.bxg7_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, g7g7.Samples, bxg7.Samples) %>% mutate(merge = paste0(g7g7.Samples, bxg7.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(g7g7.frac = g7g7.Samples/max(g7g7.Samples)) %>% mutate(bxg7.Samples.successes = bxg7.Samples) %>% mutate(bxg7.Samples.failures = max(bxg7.Samples)-bxg7.Samples) %>% mutate(bxg7.frac = bxg7.Samples/max(bxg7.Samples))

bb.bxs_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, bb.Samples, bxs.Samples) %>% mutate(merge = paste0(bb.Samples, bxs.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(bb.frac = bb.Samples/max(bb.Samples)) %>% mutate(bxs.Samples.successes = bxs.Samples) %>% mutate(bxs.Samples.failures = max(bxs.Samples)-bxs.Samples) %>% mutate(bxs.frac = bxs.Samples/max(bxs.Samples))
ss.bxs_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, ss.Samples, bxs.Samples) %>% mutate(merge = paste0(ss.Samples, bxs.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(ss.frac = ss.Samples/max(ss.Samples)) %>% mutate(bxs.Samples.successes = bxs.Samples) %>% mutate(bxs.Samples.failures = max(bxs.Samples)-bxs.Samples) %>% mutate(bxs.frac = bxs.Samples/max(bxs.Samples))

ff.fxs_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, ff.Samples, fxs.Samples) %>% mutate(merge = paste0(ff.Samples, fxs.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(ff.frac = ff.Samples/max(ff.Samples)) %>% mutate(fxs.Samples.successes = fxs.Samples) %>% mutate(fxs.Samples.failures = max(fxs.Samples)-fxs.Samples) %>% mutate(fxs.frac = fxs.Samples/max(fxs.Samples))
ss.fxs_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, ss.Samples, fxs.Samples) %>% mutate(merge = paste0(ss.Samples, fxs.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(ss.frac = ss.Samples/max(ss.Samples)) %>% mutate(fxs.Samples.successes = fxs.Samples) %>% mutate(fxs.Samples.failures = max(fxs.Samples)-fxs.Samples) %>% mutate(fxs.frac = fxs.Samples/max(fxs.Samples))

bb.IAbhemi_data <- pr_full_CD4_mnml %>% select(V.name, CDR3.aa, J.name, `b+/-.Samples`, bb.Samples) %>% mutate(merge = paste0(`b+/-.Samples`, bb.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(IAbhemi.frac = `b+/-.Samples`/max(`b+/-.Samples`)) %>% mutate(bb.frac = bb.Samples/max(bb.Samples)) %>% mutate(`b+/-.Samples.successes` = `b+/-.Samples`) %>% mutate(`b+/-.Samples.failures` = max(`b+/-.Samples`)-`b+/-.Samples`)

head(bb.bxf_data) # look at the organized data.
# "frac" columns are the proportion of isolates which bear a given sequence. 
# linear model requires a response variableto have a matrix of "sucesses" and "failures" column. Warning messages appear if I comute this on the fly. so having it defined as a column is useful.
# creating a function  wrapper around this is an issue... So just do all of these individually.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.bb.bxf <- glm(cbind(bxf.Samples.successes, bxf.Samples.failures) ~ bb.frac, family = 'quasibinomial', data = bb.bxf_data)
yhat.df <- emmeans(model.bb.bxf, ~ bb.frac, at=list(bb.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bxf <- ggplot(bb.bxf_data, aes(x = bb.frac))+
  rasterise(geom_jitter(aes(y=bxf.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxf sequence\n presence in a bb sample")+
  labs(x ="Fraction of bb isolates\n with sequence", y = "Probability of b×f in bb")+
  theme_classic()


# new_plot = probability_plot.bb.bxf + rasterise(rasterise(geom_jitter(aes(y=bxf.frac), alpha = 0.01), dpi = 900)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.ff.bxf <- glm(cbind(bxf.Samples.successes, bxf.Samples.failures) ~ ff.frac, family = 'quasibinomial', data = ff.bxf_data)
yhat.df <- emmeans(model.ff.bxf, ~ ff.frac, at=list(ff.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ff.bxf <- ggplot(ff.bxf_data, aes(x = ff.frac))+
  rasterise(geom_jitter(aes(y=bxf.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxf sequence\n presence in a ff sample")+
  labs(x ="Fraction of ff isolates\n with sequence", y = "Probability of b×f in ff")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.bb.bxg7 <- glm(cbind(bxg7.Samples.successes, bxg7.Samples.failures) ~ bb.frac, family = 'quasibinomial', data = bb.bxg7_data)
yhat.df <- emmeans(model.bb.bxg7, ~ bb.frac, at=list(bb.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bxg7 <- ggplot(bb.bxg7_data, aes(x = bb.frac))+
  rasterise(geom_jitter(aes(y=bxg7.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxg7 sequence\n presence in a bb sample")+
  labs(x ="Fraction of bb isolates\n with sequence", y = "Probability of b×g7 in bb")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.g7g7.bxg7 <- glm(cbind(bxg7.Samples.successes, bxg7.Samples.failures) ~ g7g7.frac, family = 'quasibinomial', data = g7g7.bxg7_data)
yhat.df <- emmeans(model.g7g7.bxg7, ~ g7g7.frac, at=list(g7g7.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.g7g7.bxg7 <- ggplot(g7g7.bxg7_data, aes(x = g7g7.frac))+
  rasterise(geom_jitter(aes(y=bxg7.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxg7 sequence\n presence in a g7g7 sample")+
  labs(x ="Fraction of g7g7 isolates\n with sequence", y = "Probability of b×g7 in g7g7")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.bb.bxs <- glm(cbind(bxs.Samples.successes, bxs.Samples.failures) ~ bb.frac, family = 'quasibinomial', data = bb.bxs_data)
yhat.df <- emmeans(model.bb.bxs, ~ bb.frac, at=list(bb.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bxs <- ggplot(bb.bxs_data, aes(x = bb.frac))+
  rasterise(geom_jitter(aes(y=bxs.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of b×s sequence\n presence in a bb sample")+
  labs(x ="Fraction of bb isolates\n with sequence", y = "Probability of b×s in bb")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.ss.bxs <- glm(cbind(bxs.Samples.successes, bxs.Samples.failures) ~ ss.frac, family = 'quasibinomial', data = ss.bxs_data)
yhat.df <- emmeans(model.ss.bxs, ~ ss.frac, at=list(ss.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ss.bxs <- ggplot(ss.bxs_data, aes(x = ss.frac))+
  rasterise(geom_jitter(aes(y=bxs.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxs sequence\n presence in a ss sample")+
  labs(x ="Fraction of ss isolates\n with sequence", y = "Probability of b×s in ss")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.ff.fxs <- glm(cbind(fxs.Samples.successes, fxs.Samples.failures) ~ ff.frac, family = 'quasibinomial', data = ff.fxs_data)
yhat.df <- emmeans(model.ff.fxs, ~ ff.frac, at=list(ff.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ff.fxs <- ggplot(ff.fxs_data, aes(x = ff.frac))+
  rasterise(geom_jitter(aes(y=fxs.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of fxs sequence\n presence in a ff sample")+
  labs(x ="Fraction of ff isolates\n with sequence", y = "Probability of f×s in ff")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.ss.fxs <- glm(cbind(fxs.Samples.successes, fxs.Samples.failures) ~ ss.frac, family = 'quasibinomial', data = ss.fxs_data)
yhat.df <- emmeans(model.ss.fxs, ~ ss.frac, at=list(ss.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ss.fxs <- ggplot(ss.fxs_data, aes(x = ss.frac))+
  rasterise(geom_jitter(aes(y=fxs.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of fxs sequence\n presence in a ss sample")+
  labs(x ="Fraction of ss isolates\n with sequence", y = "Probability of f×s in ss")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Test Probability of b+/- in bb sample
model.bb.IAbhemi <- glm(cbind(`b+/-.Samples.successes`, `b+/-.Samples.failures`) ~ bb.frac, family = 'quasibinomial', data = bb.IAbhemi_data)
yhat.df <- emmeans(model.bb.IAbhemi, ~ bb.frac, at=list(bb.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame()
probability_plot.bb.IAbhemi <- ggplot(bb.IAbhemi_data, aes(x = bb.frac))+
  rasterise(geom_jitter(aes(y=IAbhemi.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of b+/- sequence\n presence in a bb sample")+
  labs(x ="Fraction of bb isolates\n with sequence", y = "Probability of b+/- in bb")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(file = "Figure Export/Binomial Regression/probability_plots.pdf", width = 2.75, height = 2.75)
  probability_plot.bb.bxf
  probability_plot.ff.bxf
  probability_plot.bb.bxg7
  probability_plot.g7g7.bxg7
  probability_plot.bb.bxs
  probability_plot.ss.bxs
  probability_plot.ff.fxs
  probability_plot.ss.fxs
  probability_plot.bb.IAbhemi
dev.off()

ggsave(filename = "Figure Export/Binomial Regression/probability_plots_1.pdf", plot = probability_plot.bb.bxf, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_2.pdf", plot = probability_plot.ff.bxf, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_3.pdf", plot = probability_plot.bb.bxg7, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_4.pdf", plot = probability_plot.g7g7.bxg7, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_5.pdf", plot = probability_plot.bb.bxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_6.pdf", plot = probability_plot.ss.bxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_7.pdf", plot = probability_plot.ff.fxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_8.pdf", plot = probability_plot.ss.fxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_9.pdf", plot = probability_plot.bb.IAbhemi, width = 2.75, height = 2.75, device='pdf')

# ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
# there would be no sequences which fall into that category. But we can use this to identify how incompleate our sequencing data is.

# Function to do all the filtering:
tetra.filter <- function(pr.df, Par.Samples, F1.Samples, tetra.Samples){
  out_data <- pr.df %>%
    select(V.name, CDR3.aa, J.name, {{Par.Samples}}, {{F1.Samples}}, {{tetra.Samples}}) %>%
    mutate(merge = paste0({{Par.Samples}}, {{F1.Samples}}, {{tetra.Samples}})) %>%
    filter(!merge == "000") %>% select(!merge) #  limit selection to sequences in the set
  # filter out the sequences we don't want to evaluate
  DF1 <- out_data %>% filter({{Par.Samples}} > 0 & {{F1.Samples}} == 0 & {{tetra.Samples}} == 0)
  DF2 <- out_data %>% filter({{F1.Samples}} > 0 & {{Par.Samples}} == 0 & {{tetra.Samples}} == 0)
  DF3 <- out_data %>% filter({{Par.Samples}} > 0 & {{F1.Samples}} == 0 & {{tetra.Samples}} > 0)
  DF4 <- out_data %>% filter({{Par.Samples}} == 0 & {{F1.Samples}} > 0 & {{tetra.Samples}} > 0)
  DF5 <- out_data %>% filter({{Par.Samples}} == 0 & {{F1.Samples}} == 0 & {{tetra.Samples}} > 0) #sequences which are in the tetraparental but not the others
  DF6 <- rbind(DF1, DF2, DF3, DF4, DF5) %>% mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name))
  
  out_ <- out_data %>%
    mutate(V.CDR3.J = paste(V.name, CDR3.aa, J.name)) %>%
    filter(!V.CDR3.J %in% c(DF6$V.CDR3.J)) %>%
    select(-V.CDR3.J)
  
  # out_ <- rbind(out_, DF5)
  
  out_data <- out_ %>%
    mutate(Par.F1.Samples = {{Par.Samples}} + {{F1.Samples}}) %>%
    mutate(Par.F1.frac = Par.F1.Samples/max(Par.F1.Samples)) %>%
    mutate(tetra.successes = {{tetra.Samples}}) %>%
    mutate(tetra.failures = max({{tetra.Samples}})-{{tetra.Samples}}) %>%
    mutate(tetra.frac = {{tetra.Samples}}/max({{tetra.Samples}}))
  
  return(out_data)
}
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bb_bxf.vs.bfTetra_data <- tetra.filter(pr_full_CD4_mnml, bb.Samples, bxf.Samples, bf_tetra.Samples)

model.bb_bxf.vs.bfTetra <- glm(cbind(tetra.successes, tetra.failures) ~ Par.F1.frac, family = 'quasibinomial', data = bb_bxf.vs.bfTetra_data)
yhat.df <- emmeans(model.bb_bxf.vs.bfTetra, ~ Par.F1.frac, at=list(Par.F1.frac=seq(0.18, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bxf.bfTetra <- ggplot(bb_bxf.vs.bfTetra_data, aes(x = Par.F1.frac))+
  rasterise(geom_jitter(aes(y=tetra.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxf sequence\n presence in a bb sample")+
  labs(x ="Proportion of bb and b×f\n samples with sequence", y = "Probability of bb&ff tetra.\nshared with bb and b×f")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ff_bxf.vs.bfTetra_data <- tetra.filter(pr_full_CD4_mnml, ff.Samples, bxf.Samples, bf_tetra.Samples)

model.ff_bxf.vs.bfTetra <- glm(cbind(tetra.successes, tetra.failures) ~ Par.F1.frac, family = 'quasibinomial', data = ff_bxf.vs.bfTetra_data)
yhat.df <- emmeans(model.ff_bxf.vs.bfTetra, ~ Par.F1.frac, at=list(Par.F1.frac=seq(0.16, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ff.bxf.bfTetra <- ggplot(ff_bxf.vs.bfTetra_data, aes(x = Par.F1.frac))+
  rasterise(geom_jitter(aes(y=tetra.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxf sequence\n presence in a bb sample")+
  labs(x ="Proportion of ff and b×f\n samples with sequence", y = "Probability of bb&ff tetra.\nshared with ff and b×f")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bb_bxg7.vs.bg7Tetra_data <- tetra.filter(pr_full_CD4_mnml, bb.Samples, bxg7.Samples, bg7_tetra.Samples)

model.bb_bxg7.vs.bg7Tetra <- glm(cbind(tetra.successes, tetra.failures) ~ Par.F1.frac, family = 'quasibinomial', data = bb_bxg7.vs.bg7Tetra_data)
yhat.df <- emmeans(model.bb_bxg7.vs.bg7Tetra, ~ Par.F1.frac, at=list(Par.F1.frac=seq(0.18, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bxg7.bg7Tetra <- ggplot(bb_bxg7.vs.bg7Tetra_data, aes(x = Par.F1.frac))+
  rasterise(geom_jitter(aes(y=tetra.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxg7 sequence\n presence in a bb sample")+
  labs(x ="Proportion of bb and b×g7\n samples with sequence", y = "Probability of bb&g7g7 tetra.\nshared with bb and b×g7")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
g7g7_bxg7.vs.bg7Tetra_data <- tetra.filter(pr_full_CD4_mnml, g7g7.Samples, bxg7.Samples, bg7_tetra.Samples)

model.g7g7_bxg7.vs.bg7Tetra <- glm(cbind(tetra.successes, tetra.failures) ~ Par.F1.frac, family = 'quasibinomial', data = g7g7_bxg7.vs.bg7Tetra_data)
yhat.df <- emmeans(model.g7g7_bxg7.vs.bg7Tetra, ~ Par.F1.frac, at=list(Par.F1.frac=seq(0.22, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.g7g7.bxg7.bg7Tetra <- ggplot(g7g7_bxg7.vs.bg7Tetra_data, aes(x = Par.F1.frac))+
  rasterise(geom_jitter(aes(y=tetra.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bxg7 sequence\n presence in a g7g7 sample")+
  labs(x ="Proportion of g7g7 and b×g7\n samples with sequence", y = "Probability of bb&g7g7 tetra.\nshared with g7g7 and b×g7")+
  theme_classic()

# -------------
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_10.pdf", plot = probability_plot.bb.bxf.bfTetra, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_11.pdf", plot = probability_plot.ff.bxf.bfTetra, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_12.pdf", plot = probability_plot.bb.bxg7.bg7Tetra, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_13.pdf", plot = probability_plot.g7g7.bxg7.bg7Tetra, width = 2.75, height = 2.75, device='pdf')


# -------------------------------------
# use Binomial regression modeling to correlate sequences derived from one MHC homozygous animal and shared with F1 sequences from different haplotypes.
# ie: are the sequnces whihc are derived from bb from one F1 the same in other haplotypes?

load("Intersection.CD4.rep.RData")
# bxf, bxg7, bxs
TCR.CD4_bb.F1.intersections.pr <- pubRep(TCR.CD4_bb.F1.intersections$data, "v+j+aa", .verbose = T)
TCR.CD4_bb.F1.intersections.pr <- as_tibble(TCR.CD4_bb.F1.intersections.pr) %>%
  dplyr::mutate(bb.interse.bxf = rowSums(!is.na(across(TCR.CD4_bb.F1.intersections$meta %>% filter(MHC == 'bxf') %>% pull(Sample))))) %>%
  dplyr::mutate(bb.interse.bxg7 = rowSums(!is.na(across(TCR.CD4_bb.F1.intersections$meta %>% filter(MHC == 'bxg7') %>% pull(Sample))))) %>%
  dplyr::mutate(bb.interse.bxs = rowSums(!is.na(across(TCR.CD4_bb.F1.intersections$meta %>% filter(MHC == 'bxs') %>% pull(Sample))))) %>%
  dplyr::mutate(merge = paste0(bb.interse.bxf, bb.interse.bxg7, bb.interse.bxs)) %>% filter(!merge == "000") %>%
  select(V.name, J.name, CDR3.aa, bb.interse.bxf, bb.interse.bxg7, bb.interse.bxs)

TCR.CD4_ff.F1.intersections.pr <- pubRep(TCR.CD4_ff.F1.intersections$data, "v+j+aa", .verbose = T)
TCR.CD4_ff.F1.intersections.pr <- as_tibble(TCR.CD4_ff.F1.intersections.pr) %>%
  dplyr::mutate(ff.interse.bxf = rowSums(!is.na(across(TCR.CD4_ff.F1.intersections$meta %>% filter(MHC == 'bxf') %>% pull(Sample))))) %>%
  dplyr::mutate(ff.interse.fxs = rowSums(!is.na(across(TCR.CD4_ff.F1.intersections$meta %>% filter(MHC == 'fxs') %>% pull(Sample))))) %>%
  dplyr::mutate(merge = paste0(ff.interse.bxf, ff.interse.fxs)) %>% filter(!merge == "00") %>%
  select(V.name, J.name, CDR3.aa, ff.interse.bxf, ff.interse.fxs)

intersection.filter.plot <- function(pr.df, F1.x, F1.y){
  out_data <- pr.df %>%
    select(V.name, CDR3.aa, J.name, {{F1.x}}, {{F1.y}}) %>%
    mutate(merge = paste0({{F1.x}}, {{F1.y}})) %>%
    filter(!merge == "00") %>% select(!merge) %>% #  limit selection to sequences in the set
    mutate(F1.x.frac = {{F1.x}}/max({{F1.x}})) %>%
    mutate(F1.y.successes = {{F1.y}}) %>%
    mutate(F1.y.failures = max({{F1.y}})-{{F1.y}}) %>%
    mutate(F1.y.frac = {{F1.y}}/max({{F1.y}}))
  # # Build regression model.
  model.qb <- glm(cbind(F1.y.successes, F1.y.failures) ~ F1.x.frac, family = 'quasibinomial', data = out_data)
  yhat.df <- emmeans(model.qb, ~ F1.x.frac, at=list(F1.x.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.

  probability_plot <- ggplot(out_data, aes(x = F1.x.frac))+
    rasterise(geom_jitter(aes(y=F1.y.frac), alpha = 0.01))+
    geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
    geom_line( data=yhat.df, aes(y=prob), color='red') +
    xlim(0,1)+
    ylim(0,1)+
    coord_fixed()+
    theme_classic()
  
  return(probability_plot)
}
# -----------------------------------------------------------------------------------------------------------
probability_plot.bxf.bxg7 <- intersection.filter.plot(TCR.CD4_bb.F1.intersections.pr, bb.interse.bxf, bb.interse.bxg7) + labs(x ="Fraction of bb∩b×f isolates\n with sequence", y = "Probability of bb∩b×g7\n in bb∩b×f")
probability_plot.bxf.bxs <- intersection.filter.plot(TCR.CD4_bb.F1.intersections.pr, bb.interse.bxf, bb.interse.bxs) + labs(x ="Fraction of bb∩b×f isolates\n with sequence", y = "Probability of bb∩b×s\n in bb∩b×f")
probability_plot.bxg7.bxs <- intersection.filter.plot(TCR.CD4_bb.F1.intersections.pr, bb.interse.bxg7, bb.interse.bxs) + labs(x ="Fraction of bb∩b×g7 isolates\n with sequence", y = "Probability of bb∩b×s\n in bb∩b×g7")
probability_plot.bxf.fxs <- intersection.filter.plot(TCR.CD4_ff.F1.intersections.pr, ff.interse.bxf, ff.interse.fxs) + labs(x ="Fraction of ff∩b×f isolates\n with sequence", y = "Probability of ff∩f×s\n in ff∩b×f")
 #--------------------------::::::::::::::::::::::::::::::::::----------------------------------------------- 
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_14.pdf", plot = probability_plot.bxf.bxg7, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_15.pdf", plot = probability_plot.bxf.bxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_16.pdf", plot = probability_plot.bxg7.bxs, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plots_17.pdf", plot = probability_plot.bxf.fxs, width = 2.75, height = 2.75, device='pdf')
# ▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦▦
# generate a probability model for 1/2 the bb+ff simulated union + remaining 1/2 of bb
# generate a probability model for 1/2 the bb+ff simulated union + remaining 1/2 of ff
load("TCR_fixed_beta_total_rep.RData")
load("DOvbeta_paired_repData.RData")
bb.ff_data_separate <- repFilter(TCR.CD4_DoVb_TR, .method = "by.meta", .query = list(MHC = include("bb","ff"), Sample = include("222 postj", "463_465LNC.SC_postPy", "550 postJ", "264 postJ", "270 postJ"))) # select out the bb and ff data
bb.ff_data_Paired <- repFilter(TCR.CD4_DoVb_Paired, .method = "by.meta", .query = list(MHC = include("bb+ff"), Sample = include("267 postJ_225 postJ.bb.ff.CD4", "505 postj.one_mouse_228 postJ.bb.ff.CD4", "509 postj.one_mouse_546 postJ.one_mouse.bb.ff.CD4"))) # select out bb+ff paired data

bb.ff_data <- list(
  data = c(bb.ff_data_separate$data, bb.ff_data_Paired$data),
  meta = rbind(select(bb.ff_data_separate$meta, Sample, MHC, `Co-receptor`), bb.ff_data_Paired$meta)
)

bb.ff_data.pr <- pubRep(bb.ff_data$data, "v+j+aa", .verbose = T)
bb.ff_data.pr <- as_tibble(bb.ff_data.pr) %>%
  dplyr::mutate(bb.Samples = rowSums(!is.na(across(bb.ff_data$meta %>% filter(MHC == 'bb') %>% pull(Sample))))) %>%
  dplyr::mutate(ff.Samples = rowSums(!is.na(across(bb.ff_data$meta %>% filter(MHC == 'ff') %>% pull(Sample))))) %>%
  dplyr::mutate(bbff.Samples = rowSums(!is.na(across(bb.ff_data$meta %>% filter(MHC == 'bb+ff') %>% pull(Sample))))) %>%
  select(V.name, J.name, CDR3.aa, bb.Samples, ff.Samples, bbff.Samples)

bb.bbff_data <- bb.ff_data.pr %>% select(V.name, CDR3.aa, J.name, bb.Samples, bbff.Samples) %>% mutate(merge = paste0(bb.Samples, bbff.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(bb.frac = bb.Samples/max(bb.Samples)) %>% mutate(bbff.Samples.successes = bbff.Samples) %>% mutate(bbff.Samples.failures = max(bbff.Samples)-bbff.Samples) %>% mutate(bbff.frac = bbff.Samples/max(bbff.Samples))
ff.bbff_data <- bb.ff_data.pr %>% select(V.name, CDR3.aa, J.name, ff.Samples, bbff.Samples) %>% mutate(merge = paste0(ff.Samples, bbff.Samples)) %>% filter(!merge == "00") %>% select(!merge) %>% mutate(ff.frac = ff.Samples/max(ff.Samples)) %>% mutate(bbff.Samples.successes = bbff.Samples) %>% mutate(bbff.Samples.failures = max(bbff.Samples)-bbff.Samples) %>% mutate(bbff.frac = bbff.Samples/max(bbff.Samples))

model.bb.bbff <- glm(cbind(bbff.Samples.successes, bbff.Samples.failures) ~ bb.frac, family = 'quasibinomial', data = bb.bbff_data)
yhat.df <- emmeans(model.bb.bbff, ~ bb.frac, at=list(bb.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.bb.bbff <- ggplot(bb.bbff_data, aes(x = bb.frac))+
  rasterise(geom_jitter(aes(y=bbff.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bbff sequence\n presence in a bb sample")+
  labs(x ="Fraction of bb isolates\n with sequence", y = "Probability of bb+ff in bb")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model.ff.bbff <- glm(cbind(bbff.Samples.successes, bbff.Samples.failures) ~ ff.frac, family = 'quasibinomial', data = ff.bbff_data)
yhat.df <- emmeans(model.ff.bbff, ~ ff.frac, at=list(ff.frac=seq(0, 1, by = 0.01)), type='response') %>% as.data.frame() # this is like the predict function but it has 95% confidence intervals.
probability_plot.ff.bbff <- ggplot(ff.bbff_data, aes(x = ff.frac))+
  rasterise(geom_jitter(aes(y=bbff.frac), alpha = 0.01))+
  geom_ribbon(data=yhat.df, aes(ymin=asymp.LCL, ymax=asymp.UCL), fill='salmon', alpha=.4) +
  geom_line( data=yhat.df, aes(y=prob), color='red') +
  xlim(0,1)+
  ylim(0,1)+
  coord_fixed()+
  # labs(title="Probability of bbff sequence\n presence in a ff sample")+
  labs(x ="Fraction of ff isolates\n with sequence", y = "Probability of bb+ff in ff")+
  theme_classic()
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ggsave(filename = "Figure Export/Binomial Regression/probability_plot.bb_bbff.pdf", plot = probability_plot.bb.bbff, width = 2.75, height = 2.75, device='pdf')
ggsave(filename = "Figure Export/Binomial Regression/probability_plot.ff_bbff.pdf", plot = probability_plot.ff.bbff, width = 2.75, height = 2.75, device='pdf')
