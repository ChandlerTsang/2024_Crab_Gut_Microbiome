⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
####  Plot mantel phylosymbiosis plots with statistics  ####
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

## Load library
require(ggplot2)
require(plyr)
require(ggpmisc)
library(broom)
library(dplyr)


## Define file path
dist_fp <- '../phylosymbiosis_LC3/phylo.distance_Jaccard_table.csv'
dist_fp_Bray <- '../phylosymbiosis/phylo.distance_Bray_table.csv'
dist_fp_uwu <- '../phylosymbiosis/phylo.distance_Unweighted_UniFrac_table.csv'
dist_fp_wu <- '../phylosymbiosis/phylo.distance_Weighted_UniFrac_table.csv'


## Create plot directory
plot_dir <- '../phylosymbiosis_LC3/plots'
dir.create(file.path(plot_dir))


## Read DMs
dist_df <- read.csv(dist_fp, header = TRUE)
dist_df_Bray <- read.csv(dist_fp_Bray, header = TRUE)
dist_df_uwu <- read.csv(dist_fp_uwu, header = TRUE)
dist_df_wu <- read.csv(dist_fp_wu, header = TRUE)


## Divergence time calculation
dist_df$divergence <- dist_df$phylo / 2
dist_df_Bray$divergence <- dist_df_Bray$phylo / 2
dist_df_uwu$divergence <- dist_df_uwu$phylo / 2
dist_df_wu$divergence <- dist_df_wu$phylo / 2


## Remove same species comparison
dist_df <- subset(dist_df, host_1 != host_2)
dist_df_Bray <- subset(dist_df_Bray, host_1 != host_2)
dist_df_uwu <- subset(dist_df_uwu, host_1 != host_2)
dist_df_wu <- subset(dist_df_wu, host_1 != host_2)




⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
#########################################
##  Figure 4A - Jaccard vs Divergence  ##
#########################################
phylo_jaccard_all <- ggplot(dist_df,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_jaccard_all,
       filename = file.path(plot_dir, 'phylo_jaccard_all_stat_1.png'), 
       width=4,
       height=4)


##########################################################
##  Supplementary Figure 6 - Bray-Curtis vs Divergence  ##
##########################################################
phylo_Bray_all <- ggplot(dist_df_Bray,
                        aes(x = divergence, y = Bray)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Bray-Curtis)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Bray_all,
       filename = file.path(plot_dir, 'phylo_Bray_all_stat.png'), 
       width=4,
       height=4)


#################################################################
##  Supplementary Figure 6 - Unweighted UniFrac vs Divergence  ##
#################################################################
phylo_uwu_all <- ggplot(dist_df_uwu,
                        aes(x = divergence, y = Unweighted_UniFrac)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Unweighted UniFrac)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_uwu_all,
       filename = file.path(plot_dir, 'phylo_uwu_all_stat.png'), 
       width=4,
       height=4)


###############################################################
##  Supplementary Figure 6 - Weighted UniFrac vs Divergence  ##
###############################################################
phylo_wu_all <- ggplot(dist_df_wu,
                        aes(x = divergence, y = Weighted_UniFrac)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Weighted UniFrac)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_wu_all,
       filename = file.path(plot_dir, 'phylo_wu_all_stat.png'), 
       width=4,
       height=4)
	   

######################################################
##  Supplementary Figure 5 - Jaccard vs Divergence  ##
######################################################
dist_df_FMH = subset(dist_df, (host_1_gutregion == host_2_gutregion))

phylo_jaccard_FMH <- ggplot(dist_df_FMH,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() +
  facet_wrap(~ host_1_gutregion) +   
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_jaccard_FMH, 
       filename = file.path(plot_dir, 'phylo_Jaccard_FMH_stat.png'), 
       width=8,
       height=3)




⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
#####################################################################
##  Supplementary Figure 7 - Jaccard vs Divergence in each Family  ##
#####################################################################

## Jaccard - Dotillidae
data_Dotillidae = subset(dist_df, (host_1_family == 'Dotillidae') & (host_2_family == 'Dotillidae'))

phylo_Dotillidae_jaccard_all <- ggplot(data_Dotillidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Dotillidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Dotillidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


## Jaccard - Macrophthalmidae
data_Macrophthalmidae = subset(dist_df, (host_1_family == 'Macrophthalmidae') & (host_2_family == 'Macrophthalmidae'))

phylo_Macrophthalmidae_jaccard_all <- ggplot(data_Macrophthalmidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Macrophthalmidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Macrophthalmidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


## Jaccard - Portunidae
data_Portunidae = subset(dist_df, (host_1_family == 'Portunidae') & (host_2_family == 'Portunidae'))

phylo_Portunidae_jaccard_all <- ggplot(data_Portunidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Portunidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Portunidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


## Jaccard - Ocypodidae
data_Ocypodidae = subset(dist_df, (host_1_family == 'Ocypodidae') & (host_2_family == 'Ocypodidae'))

phylo_Ocypodidae_jaccard_all <- ggplot(data_Ocypodidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Ocypodidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Ocypodidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


## Jaccard - Sesarmidae
data_Sesarmidae = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae'))

phylo_Sesarmidae_jaccard_all <- ggplot(data_Sesarmidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Sesarmidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Sesarmidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


## Jaccard - Varunidae
data_Varunidae = subset(dist_df, (host_1_family == 'Varunidae') & (host_2_family == 'Varunidae'))

phylo_Varunidae_jaccard_all <- ggplot(data_Varunidae,
                        aes(x = divergence, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  coord_cartesian(ylim=c(0.6, 1)) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_Varunidae_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_Varunidae_jaccard_all_stat.png'), 
       width=3,
       height=3)


################################################################################
##  ALTERNATIVES - plot the previous 6 plots in one command using facet_wrap  ##
################################################################################

##define families with more than two species sampled
crab_families <- c('Dotillidae', 'Macrophthalmidae', 'Ocypodidae', 'Portunidae', 'Varunidae', 'Sesarmidae')


##exclude Mictyridae and Grapsidae since there are less than 3 species
data_full_jaccard = subset(dist_df, (within_family == 'True') & (host_1_family %in% crab_families))

phylo_full_jaccard_perfamily <- ggplot(data_full_jaccard,
                        aes(x = divergence,
                            y = Jaccard)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  facet_wrap(~ host_1_family) +  
  coord_cartesian(ylim=c(0.5, 1)) + 
  labs(x = "Divergence time (Ma)",
       y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_full_jaccard_perfamily, 
       filename = file.path(plot_dir, 'phylo_full_jaccard_perfamily.png'), 
       width=8,
       height=6)




⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
########################
##  Figure 4B and 4C  ##
########################

## plot Jaccard distance against microhabitat preference distance
phylo_jaccard_all <- ggplot(dist_df,
                        aes(x = microhabitat_jaccard, y = Jaccard)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Microhabitat Distance (Jaccard)", y = "Microbiome Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_jaccard_all_microhabitat_Jaccard_microbiome_stat.png'), 
       width=4,
       height=4)


## plot microhabitat preference distance against phylogenetic distance
phylo_jaccard_all <- ggplot(dist_df,
                        aes(x = divergence, y = microhabitat_jaccard)) +
  geom_point(alpha=0.05) + 
  stat_smooth(method='lm', formula = y ~ x, color = "#E69F00") + 
  theme_bw() + 
  labs(x = "Divergence time (Ma)", y = "Microhabitat Distance (Jaccard)") +
  stat_fit_glance(method = "cor.test",
                    label.y = "bottom",
                    method.args = list(formula = ~ x + y),
                    mapping = aes(label = sprintf('r[Pearson]~"="~%.3f~~italic(P)~"="~%.2g',
                                  after_stat(estimate), after_stat(p.value))),
                    parse = TRUE) +
  theme(legend.position = "none")

ggsave(phylo_jaccard_all, 
       filename = file.path(plot_dir, 'phylo_jaccard_all_microhabitat_Jaccard_divergence_stat.png'), 
       width=4,
       height=4)
