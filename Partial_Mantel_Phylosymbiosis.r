⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
####  Phylosymbiosis with Partial mantel statistics  ####
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

********************
**  load library  **
********************
library(plyr)
library(broom)
library(dplyr)
library(psych)




****************************
**  Get Jaccard Distance  **
****************************
dist_fp <- '../phylosymbiosis_LC3/phylo.distance_Jaccard_table.csv'

plot_dir <- '../phylosymbiosis_LC3/plots'
dir.create(file.path(plot_dir))

dist_df <- read.csv(dist_fp, header = TRUE)

dist_df$divergence <- dist_df$phylo / 2

dist_df <- subset(dist_df, host_1 != host_2)




**********************************************************************
****  Multiple Regression on Phylosymbiosis - Jaccard All vs All  ****
**********************************************************************
## this section can be done recursively with a loop for large trees

##############
##  node84  ##
##############
## full dataset
node84 = dist_df

## partial correlation
pcor_node84 <- partial.r(data=node84, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node84



##############
##  node83  ##
##############
## adjust dataset
node83 = subset(dist_df, (host_1_family == 'Portunidae') & (host_2_family == 'Portunidae'))

## partial correlation
pcor_node83 <- partial.r(data=node83, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node83



##############
##  node82  ##
##############
## adjust dataset
node82 = subset(dist_df, (host_1_family != 'Portunidae') & (host_2_family != 'Portunidae'))

## partial correlation
pcor_node82 <- partial.r(data=node82, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node82



##############
##  node77  ##
##############
## adjust dataset - Dotillidae
node77 = subset(dist_df, (host_1_family == 'Dotillidae') & (host_2_family == 'Dotillidae'))

## partial correlation
pcor_node77 <- partial.r(data=node77, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node77



##############
##  node76  ##
##############
## adjust dataset
node76 = subset(node82, (host_1_family != 'Dotillidae') & (host_2_family != 'Dotillidae'))

## partial correlation
pcor_node76 <- partial.r(data=node76, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node76



##############
##  node69  ##
##############
## adjust dataset - Ocypodidae
node69 = subset(dist_df, (host_1_family == 'Ocypodidae') & (host_2_family == 'Ocypodidae'))

## partial correlation
pcor_node69 <- partial.r(data=node69, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node69



##############
##  node66  ##
##############
## adjust dataset - Ocypodidae
node66 = subset(node69, (host_1_family == 'Ocypodidae') & (host_2_family == 'Ocypodidae')
                     & (host_1_genus != 'Tubuca') & (host_2_genus != 'Tubuca'))

## partial correlation
pcor_node66 <- partial.r(data=node66, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node66



##############
##  node64  ##
##############
## adjust dataset
node64 = subset(node66, (host_1_family == 'Ocypodidae') & (host_2_family == 'Ocypodidae')
                     & (host_1_genus != 'Austruca') & (host_2_genus != 'Austruca'))


## partial correlation
pcor_node64 <- partial.r(data=node64, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node64



##############
##  node68  ##
##############
## adjust dataset
node68 = subset(node76, (host_1_family != 'Ocypodidae') & (host_2_family != 'Ocypodidae'))

## partial correlation
pcor_node68 <- partial.r(data=node68, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node68



##############
##  node58  ##
##############
## adjust dataset
node58 = subset(node68, (host_1_family != 'Grapsidae') & (host_2_family != 'Grapsidae'))

## partial correlation
pcor_node58 <- partial.r(data=node58, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node58



##############
##  node54  ##
##############
## adjust dataset
node54 = subset(node58, (host_1_family != 'Sesarmidae') & (host_2_family != 'Sesarmidae'))

## partial correlation
pcor_node54 <- partial.r(data=node54, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node54



##############
##  node16  ##
##############
## adjust dataset
node16 = subset(node54, (host_1_family != 'Mictyridae') & (host_2_family != 'Mictyridae'))

## partial correlation
pcor_node16 <- partial.r(data=node16, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node16



##############
##  node14  ##
##############
## adjust dataset - Macrophthalmidae
node14 = subset(dist_df, (host_1_family == 'Macrophthalmidae') & (host_2_family == 'Macrophthalmidae'))

## partial correlation
pcor_node14 <- partial.r(data=node14, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node14



#############
##  node4  ##
#############
## adjust dataset
node4 = subset(node14, (host_1_family == 'Macrophthalmidae') & (host_2_family == 'Macrophthalmidae')
                     & (host_1 != 'Macrophthalmus_convexus') & (host_2 != 'Macrophthalmus_convexus'))

## partial correlation
pcor_node4 <- partial.r(data=node4, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node4



##############
##  node15  ##
##############
## adjust dataset - Varunidae
node15 = subset(dist_df, (host_1_family == 'Varunidae') & (host_2_family == 'Varunidae'))

## partial correlation
pcor_node15 <- partial.r(data=node15, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node15



##############
##  node13  ##
##############
## adjust dataset
node13 = subset(node15, (host_1_genus != 'Metaplax') & (host_2_genus != 'Metaplax'))

## partial correlation
pcor_node13 <- partial.r(data=node13, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node13



##############
##  node55  ##
##############
## adjust dataset - Sesarmidae
node55 = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae'))

## partial correlation
pcor_node55 <- partial.r(data=node55, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node55



##############
##  node53  ##
##############
## adjust dataset
node53 = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
                          & (host_1_genus != 'Lithoselatium') & (host_2_genus != 'Lithoselatium')
						  & (host_1_genus != 'Selatium') & (host_2_genus != 'Selatium')
						  & (host_1_genus != 'Episesarma') & (host_2_genus != 'Episesarma')
						  & (host_1_genus != 'Sarmatium') & (host_2_genus != 'Sarmatium')
						  & (host_1_genus != 'Neosarmatium') & (host_2_genus != 'Neosarmatium')
						  & (host_1_genus != 'Parasesarma') & (host_2_genus != 'Parasesarma')
						  )

## partial correlation
pcor_node53 <- partial.r(data=node53, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node53



##############
##  node51  ##
##############
## adjust dataset
node51 = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus == 'Orisarma') & (host_2_genus == 'Orisarma'))

## partial correlation
pcor_node51 <- partial.r(data=node51, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node51



##############
##  node52  ##
##############
## adjust dataset
node52 = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
                          & (host_1_genus != 'Chiromantes') & (host_2_genus != 'Chiromantes')
						  & (host_1_genus != 'Orisarma') & (host_2_genus != 'Orisarma'))

## partial correlation
pcor_node52 <- partial.r(data=node52, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node52



##############
##  node45  ##
##############
## adjust dataset
node45 = subset(node52, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus != 'Sarmatium') & (host_2_genus != 'Sarmatium')
						  & (host_1_genus != 'Neosarmatium') & (host_2_genus != 'Neosarmatium')
						  & (host_1_genus != 'Parasesarma') & (host_2_genus != 'Parasesarma'))

## partial correlation
pcor_node45 <- partial.r(data=node45, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node45



##############
##  node42  ##
##############
## adjust dataset
node42 = subset(node45, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
                          & (host_1_genus != 'Lithoselatium') & (host_2_genus != 'Lithoselatium')
						  & (host_1_genus != 'Selatium') & (host_2_genus != 'Selatium'))

## partial correlation
pcor_node42 <- partial.r(data=node42, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node42



##############
##  node39  ##
##############
## adjust dataset
node39 = subset(node42, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1 != 'Episesarma_versicolor') & (host_2 != 'Episesarma_versicolor'))

## partial correlation
pcor_node39 <- partial.r(data=node39, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node39



##############
##  node44  ##
##############
## adjust dataset
node44 = subset(node52, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
                          & (host_1_genus != 'Lithoselatium') & (host_2_genus != 'Lithoselatium')
						  & (host_1_genus != 'Selatium') & (host_2_genus != 'Selatium')
						  & (host_1_genus != 'Episesarma') & (host_2_genus != 'Episesarma'))

## partial correlation
pcor_node44 <- partial.r(data=node44, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node44



##############
##  node32  ##
##############
## adjust dataset
node32 = subset(node44, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus != 'Parasesarma') & (host_2_genus != 'Parasesarma'))

## partial correlation
pcor_node32 <- partial.r(data=node32, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node32



##############
##  node33  ##
##############
## adjust dataset - Parasesarma
node33 = subset(node44, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus != 'Sarmatium') & (host_2_genus != 'Sarmatium')
						  & (host_1_genus != 'Neosarmatium') & (host_2_genus != 'Neosarmatium'))

## partial correlation
pcor_node33 <- partial.r(data=node33, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node33

## check node 33 with different subsetting criteria
## node33 = subset(dist_df, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
## 						  & (host_1_genus == 'Parasesarma') & (host_2_genus == 'Parasesarma'))



##############
##  node31  ##
##############
## adjust dataset
node31 = subset(node33, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus == 'Parasesarma') & (host_2_genus == 'Parasesarma')
						  & (host_1 != 'Parasesarma_pictum') & (host_2 != 'Parasesarma_pictum'))

## partial correlation
pcor_node31 <- partial.r(data=node31, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node31



##############
##  node29  ##
##############
## adjust dataset
node29 = subset(node31, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus == 'Parasesarma') & (host_2_genus == 'Parasesarma')
						  & (host_1 != 'Parasesarma_ungulatum') & (host_2 != 'Parasesarma_ungulatum'))

## partial correlation
pcor_node29 <- partial.r(data=node29, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node29



##############
##  node27  ##
##############
## adjust dataset
node27 = subset(node29, (host_1_family == 'Sesarmidae') & (host_2_family == 'Sesarmidae')
						  & (host_1_genus == 'Parasesarma') & (host_2_genus == 'Parasesarma')
						  & (host_1 != 'Parasesarma_affine') & (host_2 != 'Parasesarma_affine'))

## partial correlation
pcor_node27 <- partial.r(data=node27, x=c("Jaccard","divergence"), y="microhabitat_jaccard")
pcor_node27
