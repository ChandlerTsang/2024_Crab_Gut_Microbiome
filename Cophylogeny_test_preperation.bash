⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
####  Cophylogeny Test Preparation  ####
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
## This script was run manually and only the genus names and the crab groups were changed in each set of tests.
## I will use the core bacterial genus Candidatus Bacilloplasma in Sesarmidae as an example



##################
##  Sesarmidae  ##
##################

## Trim table to contain Sesarmidae species only
qiime feature-table filter-samples \
  --i-table ../../table_filtered_clean_LC3.qza \
  --m-metadata-file ../../metadata/metadata_qiime_final.tsv \
  --p-where "family='Sesarmidae'" \
  --p-no-exclude-ids \
  --o-filtered-table table_filtered_clean_LC3_Sesarmidae.qza

## note: if working on all crabs, no need to trim the table to having only 1 family,
##       but on the filp side, you can actually choose any samples you want to analysis
##       with the qiime taxa filter-table function below
##           e.g. if you want to look at the cophylogenetic signals in birds that are granivorous,
##                you can substitute the above selection criteria to --p-where "diet='granivore'"



########################################################################
##  Sesarmidae - D_4__Mycoplasmataceae;D_5__Candidatus Bacilloplasma  ##
########################################################################

## Extract table for core microbiome taxa
qiime taxa filter-table \
  --i-table ./table_filtered_clean_LC3_Sesarmidae.qza \
  --i-taxonomy ../../taxonomy.qza \
  --p-include 'D_4__Mycoplasmataceae;D_5__Candidatus Bacilloplasma' \
  --p-mode 'contains' \
  --o-filtered-table table_filtered_clean_Sesarmidae_Candidatus_Bacilloplasma.qza

mkdir Candidatus_Bacilloplasma

## Group count table to crab species level
qiime feature-table group \
  --i-table table_filtered_clean_Sesarmidae_Candidatus_Bacilloplasma.qza \
  --p-axis sample \
  --m-metadata-file ../../metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --p-mode mean-ceiling \
  --o-grouped-table ./Candidatus_Bacilloplasma/Sesarmidae_Candidatus_Bacilloplasma.qza

## Extract sequences identified as Candidatus_Bacilloplasma
qiime feature-table filter-seqs \
  --i-data ../../rep-seqs_filtered_clean_LC3.qza \
  --i-table ./Candidatus_Bacilloplasma/Sesarmidae_Candidatus_Bacilloplasma.qza \
  --p-no-exclude-ids \
  --o-filtered-data ./Candidatus_Bacilloplasma/rep-seqs_Sesarmidae_Candidatus_Bacilloplasma.qza


####################################################
##  Prepare files for doing cophylogeny analysis  ##
####################################################

cd Candidatus_Bacilloplasma

## construct phylogeny of the Candidatus_Bacilloplasma
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_Sesarmidae_Candidatus_Bacilloplasma.qza \
  --p-n-threads 10 \
  --o-alignment rep-seqs_Sesarmidae_Candidatus_Bacilloplasma_aligned.qza \
  --o-masked-alignment rep-seqs_Sesarmidae_Candidatus_Bacilloplasma_masked-aligned.qza \
  --o-tree rep-seqs_Sesarmidae_Candidatus_Bacilloplasma_unrooted-tree.qza \
  --o-rooted-tree rep-seqs_Sesarmidae_Candidatus_Bacilloplasma_rooted-tree.qza

## export Candidatus_Bacilloplasma tree
qiime tools export \
  --input-path rep-seqs_Sesarmidae_Candidatus_Bacilloplasma_rooted-tree.qza \
  --output-path Sesarmidae_Candidatus_Bacilloplasma_rooted-tree

qiime feature-table transpose \
  --i-table Sesarmidae_Candidatus_Bacilloplasma.qza \
  --o-transposed-feature-table Sesarmidae_Candidatus_Bacilloplasma_T.qza 

qiime feature-table presence-absence \
  --i-table Sesarmidae_Candidatus_Bacilloplasma_T.qza \
  --o-presence-absence-table Sesarmidae_Candidatus_Bacilloplasma_T_PA.qza

qiime tools export \
  --input-path Sesarmidae_Candidatus_Bacilloplasma_T_PA.qza \
  --output-path Sesarmidae_Candidatus_Bacilloplasma_association-matrix

biom convert -i Sesarmidae_Candidatus_Bacilloplasma_association-matrix/feature-table.biom \
             -o Sesarmidae_Candidatus_Bacilloplasma_association-matrix/Sesarmidae_Candidatus_Bacilloplasma_association-matrix.tsv --to-tsv

## note: need to remove the header output from biom convert before doing statistical tests in R!!
nano Sesarmidae_Candidatus_Bacilloplasma_association-matrix/Sesarmidae_Candidatus_Bacilloplasma_association-matrix.tsv

## make jaccard distance matrix table from the filtered table
qiime diversity beta \
  --i-table Sesarmidae_Candidatus_Bacilloplasma.qza \
  --p-metric 'jaccard' \
  --p-n-jobs 10 \
  --o-distance-matrix Sesarmidae_Candidatus_Bacilloplasma_jaccard.qza 

## make jaccard distance matrix table from the filtered table
qiime tools export \
  --input-path Sesarmidae_Candidatus_Bacilloplasma_jaccard.qza \
  --output-path Sesarmidae_Candidatus_Bacilloplasma_jaccard

