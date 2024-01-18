⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
####  Cophylogeny Test  ####
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
## This script was run manually and only the genus names and the crab groups were changed in each set of tests.
## I will use the core bacterial genus Candidatus Bacilloplasma in Sesarmidae as an example


########################################################################
##  Sesarmidae - D_4__Mycoplasmataceae;D_5__Candidatus Bacilloplasma  ##
########################################################################

## Load R
R

## Load library
library(paco)
library(ape)
library(phytools)
library(vegan)

## Prepare the trees of the host and the core bacterial genus, and the association matrix
sesarmidae_tree <- ape::read.tree(file = "../sesarmidae.nwk")
bacteria_tree <- ape::read.tree(file = "Sesarmidae_Candidatus_Bacilloplasma_rooted-tree/tree.nwk")
asso_matrix <- read.table(file = 'Sesarmidae_Candidatus_Bacilloplasma_association-matrix/Sesarmidae_Candidatus_Bacilloplasma_association-matrix.tsv', sep = '\t', header = TRUE, row.names=1, check.names = FALSE)

hosttree <- cophenetic(sesarmidae_tree)
bacteriatree <- cophenetic(bacteria_tree)


#############################
##  Cophylogeny with PACo  ##
#############################
D <- prepare_paco_data(H = hosttree, P = bacteriatree, HP = asso_matrix)
D <- add_pcoord(D, correction = 'cailliez')
D <- PACo(D, nperm = 999, seed = 123456, method = 'quasiswap', symmetric = FALSE, shuffled = FALSE)
D <- paco_links(D)
res <- residuals_paco(D$proc)
print(D$gof)


################################
##  Cophylogeny with ParaFit  ##
################################
parafit(hosttree, bacteriatree, asso_matrix, nperm = 999, test.links = FALSE,
        seed = 123456, correction = "cailliez", silent = FALSE)


####################################
##  Cophylogeny with Mantel test  ##
####################################
jaccard_matrix <- read.table(file = 'Sesarmidae_Candidatus_Bacilloplasma_jaccard/distance-matrix.tsv',
                           sep = '\t', header = TRUE, row.names=1, check.names = FALSE)
						   
mantel(hosttree, jaccard_matrix, method="pearson", permutations=999)

