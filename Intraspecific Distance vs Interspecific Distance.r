**********************************************************
**** Intraspecific Distance vs Interspecific Distance ****
**********************************************************
####################
##  load library  ##
####################

library(ggplot2)
library(plyr)
library(ggpmisc)
library(broom)
library(ggpubr)



######################
##  data wrangling  ##
######################

## path to the 4 beta-diversity distances
dist_fp_Jaccard <- '../phylosymbiosis_LC3/phylo.distance_Jaccard_table.csv'
dist_fp_Bray <- '../phylosymbiosis_LC3/phylo.distance_Bray_table.csv'
dist_fp_Unweighted_UniFrac <- '../phylosymbiosis_LC3/phylo.distance_Unweighted_UniFrac_table.csv'
dist_fp_Weighted_UniFrac <- '../phylosymbiosis_LC3/phylo.distance_Weighted_UniFrac_table.csv'

## make plot directory
plot_dir <- '../phylosymbiosis_LC3/inter-intra-plots'
dir.create(file.path(plot_dir))

## read table
dist_df_Jaccard <- read.csv(dist_fp_Jaccard, header = TRUE)
dist_df_Bray <- read.csv(dist_fp_Bray, header = TRUE)
dist_df_Unweighted_UniFrac <- read.csv(dist_fp_Unweighted_UniFrac, header = TRUE)
dist_df_Weighted_UniFrac <- read.csv(dist_fp_Weighted_UniFrac, header = TRUE)

## Interspecific distance - filter data comparison from the same family
## Note: comparing within species distance to among species distance can easily inflate the comparison
##       due to the large distances between families, so we perform the comparison in each family separately
dist_df_inter_intra_specific <- subset(dist_df_Jaccard, host_1_family == host_2_family)
dist_df_inter_intra_specific <- subset(dist_df_Bray, host_1_family == host_2_family)
dist_df_inter_intra_specific <- subset(dist_df_Unweighted_UniFrac, host_1_family == host_2_family)
dist_df_inter_intra_specific <- subset(dist_df_Weighted_UniFrac, host_1_family == host_2_family)

## Use arcsin-transformed distance for better visualization
dist_df_inter_intra_specific['Jaccard'] <- asin(sqrt(dist_df_inter_intra_specific['Jaccard']))
dist_df_inter_intra_specific['Bray'] <- asin(sqrt(dist_df_inter_intra_specific['Bray']))
dist_df_inter_intra_specific['Unweighted_UniFrac'] <- asin(sqrt(dist_df_inter_intra_specific['Unweighted_UniFrac']))
dist_df_inter_intra_specific['Weighted_UniFrac'] <- asin(sqrt(dist_df_inter_intra_specific['Weighted_UniFrac']))



#####################
##  visualization  ##
#####################

##Jaccard
beta_boxplot <- ggplot(dist_df_inter_intra_specific,
                        aes(x = host_1_family, y = Jaccard, fill = within_species)) +
						geom_boxplot(outlier.colour = "grey", outlier.fill = 0.5) +
						scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + 
  labs(x = "Family", y = "Microbiome Distance (arcsin-transformed Jaccard)")

ggsave(beta_boxplot, 
       filename = file.path(plot_dir, 'beta_boxplot_Jaccard_arcsine.png'), 
       width=10,
       height=8)

##Bray
beta_boxplot <- ggplot(dist_df_inter_intra_specific,
                        aes(x = host_1_family, y = Bray, fill = within_species)) +
						geom_boxplot(outlier.colour = "grey", outlier.fill = 0.5) +
						scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + 
  labs(x = "Family", y = "Microbiome Distance (arcsin-transformed Bray)")


ggsave(beta_boxplot, 
       filename = file.path(plot_dir, 'beta_boxplot_Bray_arcsine.png'), 
       width=10,
       height=8)

##Unweighted_UniFrac
beta_boxplot <- ggplot(dist_df_inter_intra_specific,
                        aes(x = host_1_family, y = Unweighted_UniFrac, fill = within_species)) +
						geom_boxplot(outlier.colour = "grey", outlier.fill = 0.5) +
						scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + 
  labs(x = "Family", y = "Microbiome Distance (arcsin-transformed Unweighted_UniFrac)")

ggsave(beta_boxplot, 
       filename = file.path(plot_dir, 'beta_boxplot_Unweighted_UniFrac_arcsine.png'), 
       width=10,
       height=8)

##Weighted_UniFrac
beta_boxplot <- ggplot(dist_df_inter_intra_specific,
                        aes(x = host_1_family, y = Weighted_UniFrac, fill = within_species)) +
						geom_boxplot(outlier.colour = "grey", outlier.fill = 0.5) +
						scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + 
  labs(x = "Family", y = "Microbiome Distance (arcsin-transformed Weighted_UniFrac)")

ggsave(beta_boxplot, 
       filename = file.path(plot_dir, 'beta_boxplot_Weighted_UniFrac_arcsine.png'), 
       width=10,
       height=8)



#######################################################################
##  The following plots are for statistics - Wilcoxon Rank Sum Test  ##
#######################################################################

## Jaccard
beta_boxplot_Jaccard_stat <- ggboxplot(dist_df_inter_intra_specific,
                x = "within_species", y = "Jaccard",
                color = "host_1_family", short.panel.labs = FALSE) +
				facet_wrap(.~ host_1_family) +
      stat_compare_means(comparisons = list(c("True", "False")),
      aes(label = paste0("p = ", after_stat(p.format))),
	  method = "wilcox")

ggsave(beta_boxplot_Jaccard_stat, 
       filename = file.path(plot_dir, 'beta_boxplot_Jaccard_arcsine_stat.png'), 
       width=10,
       height=12)


## Bray
beta_boxplot_Bray_stat <- ggboxplot(dist_df_inter_intra_specific,
                x = "within_species", y = "Bray",
                color = "host_1_family", short.panel.labs = FALSE) +
				facet_wrap(.~ host_1_family) +
      stat_compare_means(comparisons = list(c("True", "False")),
      aes(label = paste0("p = ", after_stat(p.format))),
	  method = "wilcox")

ggsave(beta_boxplot_Bray_stat, 
       filename = file.path(plot_dir, 'beta_boxplot_Bray_arcsine_stat.png'), 
       width=10,
       height=12)


## Unweighted_UniFrac
beta_boxplot_Unweighted_UniFrac_stat <- ggboxplot(dist_df_inter_intra_specific,
                x = "within_species", y = "Unweighted_UniFrac",
                color = "host_1_family", short.panel.labs = FALSE) +
				facet_wrap(.~ host_1_family) +
      stat_compare_means(comparisons = list(c("True", "False")),
      aes(label = paste0("p = ", after_stat(p.format))),
	  method = "wilcox")

ggsave(beta_boxplot_Unweighted_UniFrac_stat, 
       filename = file.path(plot_dir, 'beta_boxplot_Unweighted_UniFrac_arcsine_stat.png'), 
       width=10,
       height=12)


## Weighted_UniFrac
beta_boxplot_Weighted_UniFrac_stat <- ggboxplot(dist_df_inter_intra_specific,
                x = "within_species", y = "Weighted_UniFrac",
                color = "host_1_family", short.panel.labs = FALSE) +
				facet_wrap(.~ host_1_family) +
      stat_compare_means(comparisons = list(c("True", "False")),
      aes(label = paste0("p = ", after_stat(p.format))),
	  method = "wilcox")

ggsave(beta_boxplot_Weighted_UniFrac_stat, 
       filename = file.path(plot_dir, 'beta_boxplot_Weighted_UniFrac_arcsine_stat.png'), 
       width=10,
       height=12)
