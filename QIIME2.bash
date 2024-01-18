#######################################
#### ****** QIIME2 PROTOCOL ****** ####
#######################################


⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
**** Part 1 - Import data, cutadapt, denoise ****
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

##############################################
##  Start the conda environment for QIIME2  ##
##############################################
conda activate qiime2-2022.8


######################
##  Import dataset  ##
######################
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./Novogene_Data/rawdata_final \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ./analysis/raw_data.qza


########################
##  cutadapt + DADA2  ##
########################
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences raw_data.qza \
 --p-cores 10 \
 --p-front-f CCTACGGGNGGCWGCAG \
 --p-front-r GACTACHVGGGTATCTAATCC \
 --p-discard-untrimmed True \
 --p-no-indels \
 --o-trimmed-sequences raw_data_cutadapt.qza

qiime dada2 denoise-paired \
--i-demultiplexed-seqs raw_data_cutadapt.qza \
--p-trunc-len-f 226 \
--p-trunc-len-r 222 \
--p-chimera-method 'consensus' \
--p-pooling-method 'pseudo' \
--p-n-reads-learn 10000000 \
--p-n-threads 10 \
--o-table raw_table.qza \
--o-representative-sequences raw_rep-seqs.qza \
--o-denoising-stats raw_denoising-stats.qza




⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
**** Part 2 - Assign taxonomy to the reads and remove unwanted sequences ****
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

##################################################################
##  Assign taxonomy to the reads and remove unwanted sequences  ##
##################################################################
qiime tools import \
 --input-path silva_132_99_16S.fna \
 --output-path silva_132_99_16S.qza \
 --type 'FeatureData[Sequence]'

qiime tools import \
 --input-path consensus_taxonomy_7_levels.txt \
 --output-path consensus_taxonomy_7_levels.qza \
 --type 'FeatureData[Taxonomy]'

qiime feature-classifier extract-reads \
  --i-sequences SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 350 \
  --p-n-jobs 10 \
  --o-reads SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S_V34_ref-seqs_min350.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S_V34_ref-seqs_min350.qza \
  --i-reference-taxonomy SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/consensus_taxonomy_7_levels.qza \
  --o-classifier SILVA_132_classifier_7_levels_min350.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ./SILVA_132_classifier_7_levels_min350.qza \
  --i-reads raw_rep-seqs.qza \
  --p-n-jobs 10 \
  --o-classification taxonomy.qza


##############################################
##  filter table according to the taxonomy  ##
##############################################

## filter features that cannot be classified at least to a phylum level
qiime taxa filter-table \
 --i-table raw_table.qza \
 --i-taxonomy taxonomy.qza \
 --p-include 'D_1__' \
 --p-mode 'contains' \
 --o-filtered-table raw_table_filtered_1.qza

## filter features that are classified as mitochondria
qiime taxa filter-table \
 --i-table raw_table_filtered_1.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Mitochondria' \
 --p-mode 'contains' \
 --o-filtered-table raw_table_filtered_2.qza

## filter features that are unclassified
qiime taxa filter-table \
 --i-table raw_table_filtered_2.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Unassigned' \
 --p-mode 'contains' \
 --o-filtered-table raw_table_filtered_3.qza

## filter features that are classified as chloroplast
qiime taxa filter-table \
 --i-table raw_table_filtered_3.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Chloroplast' \
 --p-mode 'contains' \
 --o-filtered-table raw_table_filtered_4.qza

## filter features that are classified as archaea
qiime taxa filter-table \
 --i-table raw_table_filtered_4.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'D_0__Archaea;' \
 --p-mode 'contains' \
 --o-filtered-table table_filtered_clean.qza

## remove intermediate files
rm raw_table_filtered_1.qza
rm raw_table_filtered_2.qza
rm raw_table_filtered_3.qza
rm raw_table_filtered_4.qza

## remove singleton and doubleton features
qiime feature-table filter-features \
 --i-table table_filtered_clean.qza \
 --p-min-frequency 3 \
 --o-filtered-table table_filtered_clean_LC3.qza


################################################
##  filter rep_seq according to the taxonomy  ##
################################################

## filter sequences that cannot be classified at least to a phylum level
qiime taxa filter-seqs \
 --i-sequences raw_rep-seqs.qza \
 --i-taxonomy taxonomy.qza \
 --p-include 'D_1__' \
 --p-mode 'contains' \
 --o-filtered-sequences raw_rep-seqs_filtered_1.qza

## filter sequences that are classified as mitochondria
qiime taxa filter-seqs \
 --i-sequences raw_rep-seqs_filtered_1.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Mitochondria' \
 --p-mode 'contains' \
 --o-filtered-sequences raw_rep-seqs_filtered_2.qza

## filter sequences that are unclassified
qiime taxa filter-seqs \
 --i-sequences raw_rep-seqs_filtered_2.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Unassigned' \
 --p-mode 'contains' \
 --o-filtered-sequences raw_rep-seqs_filtered_3.qza

## filter sequences that are classified as chloroplast
qiime taxa filter-seqs \
 --i-sequences raw_rep-seqs_filtered_3.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'Chloroplast' \
 --p-mode 'contains' \
 --o-filtered-sequences raw_rep-seqs_filtered_4.qza

## filter sequences that are classified as archaea
qiime taxa filter-seqs \
 --i-sequences raw_rep-seqs_filtered_4.qza \
 --i-taxonomy taxonomy.qza \
 --p-exclude 'D_0__Archaea;' \
 --p-mode 'contains' \
 --o-filtered-sequences rep-seqs_filtered_clean.qza

## remove intermediate files
rm raw_rep-seqs_filtered_1.qza
rm raw_rep-seqs_filtered_2.qza
rm raw_rep-seqs_filtered_3.qza
rm raw_rep-seqs_filtered_4.qza

## filter the rep seq qith the filtered table
qiime feature-table filter-seqs \
 --i-data rep-seqs_filtered_clean.qza \
 --i-table table_filtered_clean_LC3.qza \
 --p-no-exclude-ids \
 --o-filtered-data rep-seqs_filtered_clean_LC3.qza

## Export the rep seq and check number of ASVs
grep -o '>' dna-sequences.fasta | wc -l




⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
**** Part 3 - Calculate alpha-diversity (ASV and goods coverage) ****
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

mkdir alpha-diversity-results_LC3
cd alpha-diversity-results_LC3

## export the number of ASVs
qiime diversity alpha \
  --i-table ../table_filtered_clean_LC3.qza \
  --p-metric 'observed_features' \
  --o-alpha-diversity observed_features_vector.qza

qiime diversity alpha \
  --i-table ../table_filtered_clean_LC3.qza \
  --p-metric 'goods_coverage' \
  --o-alpha-diversity goods_coverage_vector.qza


## Check the differences between the alpha-diversity between families
qiime diversity alpha-group-significance \
  --i-alpha-diversity observed_features_vector.qza \
  --m-metadata-file ../metadata/metadata_qiime_final.tsv \
  --o-visualization observed_features_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity goods_coverage_vector.qza \
  --m-metadata-file ../metadata/metadata_qiime_final.tsv \
  --o-visualization goods_coverage_significance.qzv

cd ../


⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
**** Part 4 - Calculate beta-diversity (Bray-Curtis, Jaccard, Weighted and unweighted UniFrac) ****
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

## SEPP fragment insertion algorithm to insert 16S sequences into a well-curated tree
qiime fragment-insertion sepp \
  --i-representative-sequences rep-seqs_filtered_clean_LC3.qza \
  --i-reference-database sepp-refs-gg-13-8.qza \
  --p-threads 10 \
  --o-tree sepp-insertion-tree_LC3.qza \
  --o-placements insertion-placement_LC3.qza

## 4 basics beta diversity metric (Bray-Curtis, Jaccard, Weighted and unweighted UniFrac)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table table_filtered_clean_LC3.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --output-dir core-metrics-results_LC3






⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛
**** Part 5 - Beta-diversity statistics (PERMANOVA + adonis + PERMDISP) ****
⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛⬛

#################
##  PERMANOVA  ##
#################

## jaccard distance

## PERMANOVA - jaccard - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMANOVA_LC3/jaccard-species-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMANOVA_LC3/jaccard-genus-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMANOVA_LC3/jaccard-family-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMANOVA_LC3/jaccard-microhabitat-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMANOVA_LC3/jaccard-sex-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMANOVA_LC3/jaccard-locality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMANOVA_LC3/jaccard-terrestriality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMANOVA_LC3/jaccard-gutregion-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - jaccard - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMANOVA_LC3/jaccard-diet-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova



## bray_curtis distance

## PERMANOVA - bray_curtis - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMANOVA_LC3/bray_curtis-species-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMANOVA_LC3/bray_curtis-genus-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMANOVA_LC3/bray_curtis-family-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMANOVA_LC3/bray_curtis-microhabitat-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMANOVA_LC3/bray_curtis-sex-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMANOVA_LC3/bray_curtis-locality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMANOVA_LC3/bray_curtis-terrestriality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMANOVA_LC3/bray_curtis-gutregion-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - bray_curtis - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMANOVA_LC3/bray_curtis-diet-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova



## unweighted_unifrac distance

## PERMANOVA - unweighted_unifrac - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-species-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-genus-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-family-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-microhabitat-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-sex-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-locality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-terrestriality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-gutregion-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - unweighted_unifrac - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMANOVA_LC3/unweighted_unifrac-diet-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova



## weighted_unifrac distance

## PERMANOVA - weighted_unifrac - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-species-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-genus-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-family-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-microhabitat-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-sex-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-locality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-terrestriality-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-gutregion-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova

## PERMANOVA - weighted_unifrac - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMANOVA_LC3/weighted_unifrac-diet-significance_permanova.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permanova



################
##  PERMDISP  ##
################

mkdir PERMDISP_LC3

## jaccard

## PERMDISP - jaccard - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMDISP_LC3/jaccard-species-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMDISP_LC3/jaccard-genus-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMDISP_LC3/jaccard-family-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMDISP_LC3/jaccard-microhabitat-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMDISP_LC3/jaccard-sex-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMDISP_LC3/jaccard-locality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMDISP_LC3/jaccard-terrestriality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMDISP_LC3/jaccard-gutregion-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - jaccard - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMDISP_LC3/jaccard-diet-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp



## bray_curtis

## PERMDISP - bray_curtis - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMDISP_LC3/bray_curtis-species-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMDISP_LC3/bray_curtis-genus-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMDISP_LC3/bray_curtis-family-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMDISP_LC3/bray_curtis-microhabitat-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMDISP_LC3/bray_curtis-sex-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMDISP_LC3/bray_curtis-locality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMDISP_LC3/bray_curtis-terrestriality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMDISP_LC3/bray_curtis-gutregion-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - bray_curtis - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMDISP_LC3/bray_curtis-diet-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp



## unweighted_unifrac

## PERMDISP - unweighted_unifrac - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-species-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-genus-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-family-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-microhabitat-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-sex-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-locality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-terrestriality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-gutregion-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - unweighted_unifrac - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMDISP_LC3/unweighted_unifrac-diet-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp



## weighted_unifrac

## PERMDISP - weighted_unifrac - species
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column species \
  --o-visualization PERMDISP_LC3/weighted_unifrac-species-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - genus
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column genus \
  --o-visualization PERMDISP_LC3/weighted_unifrac-genus-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - family
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column family \
  --o-visualization PERMDISP_LC3/weighted_unifrac-family-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - microhabitat
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column microhabitat \
  --o-visualization PERMDISP_LC3/weighted_unifrac-microhabitat-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - sex
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column sex \
  --o-visualization PERMDISP_LC3/weighted_unifrac-sex-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - locality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column locality \
  --o-visualization PERMDISP_LC3/weighted_unifrac-locality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - terrestriality
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column terrestriality \
  --o-visualization PERMDISP_LC3/weighted_unifrac-terrestriality-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - gutregion
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column gutregion \
  --o-visualization PERMDISP_LC3/weighted_unifrac-gutregion-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp

## PERMDISP - weighted_unifrac - diet
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --m-metadata-column diet \
  --o-visualization PERMDISP_LC3/weighted_unifrac-diet-significance_permdisp.qzv \
  --p-pairwise \
  --p-permutations 9999 \
  --p-method permdisp



#############################################
##  PCoA - family with at least 3 species  ##
#############################################

## For Supplementary Figures

mkdir PERMANOVA_LC3-family

## Dotillidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Dotillidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Dotillidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Dotillidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Dotillidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Dotillidae


## Grapsidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Grapsidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Grapsidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Grapsidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Grapsidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Grapsidae


## Macrophthalmidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Macrophthalmidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Macrophthalmidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Macrophthalmidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Macrophthalmidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Macrophthalmidae


## Ocypodidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Ocypodidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Ocypodidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Ocypodidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Ocypodidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Ocypodidae


## Portunidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Portunidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Portunidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Portunidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Portunidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Portunidae


## Sesarmidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Sesarmidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Sesarmidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Sesarmidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Sesarmidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Sesarmidae


## Varunidae
qiime feature-table filter-samples \
  --i-table table_filtered_clean_LC3.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-where "family='Varunidae'" \
  --p-no-exclude-ids \
  --o-filtered-table PERMANOVA_LC3-family/table_filtered_clean_Varunidae.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny sepp-insertion-tree_LC3.qza \
  --i-table PERMANOVA_LC3-family/table_filtered_clean_Varunidae.qza \
  --p-sampling-depth 5000 \
  --p-n-jobs-or-threads 10 \
  --m-metadata-file metadata/metadata_qiime_final_Varunidae.tsv \
  --output-dir PERMANOVA_LC3-family/core-metrics-results_Varunidae




########################################################################################
## Adonis - % variatoin contribution of each factors to the composition of microbiome ##
########################################################################################

mkdir Adonis_LC3

## adonis - jaccard
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-species.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "genus" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-genus.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "family" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-family.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "microhabitat" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-microhabitat.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "sex" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-sex.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "locality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-locality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "terrestriality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-terrestriality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "gutregion" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-gutregion.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-diet.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/jaccard_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species + locality + microhabitat + sex + terrestriality + gutregion + diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-jaccard-all.qzv




## adonis - bray_curtis
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-species.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "genus" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-genus.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "family" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-family.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "microhabitat" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-microhabitat.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "sex" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-sex.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "locality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-locality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "terrestriality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-terrestriality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "gutregion" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-gutregion.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-diet.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species + locality + microhabitat + sex + terrestriality + gutregion + diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-bray_curtis-all.qzv



## adonis - unweighted_unifrac
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-species.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "genus" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-genus.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "family" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-family.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "microhabitat" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-microhabitat.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "sex" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-sex.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "locality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-locality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "terrestriality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-terrestriality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "gutregion" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-gutregion.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-diet.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species + locality + microhabitat + sex + terrestriality + gutregion + diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-unweighted_unifrac-all.qzv

  

## adonis - weighted_unifrac
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-species.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "genus" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-genus.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "family" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-family.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "microhabitat" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-microhabitat.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "sex" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-sex.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "locality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-locality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "terrestriality" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-terrestriality.qzv

qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "gutregion" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-gutregion.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-diet.qzv
  
qiime diversity adonis \
  --i-distance-matrix core-metrics-results_LC3/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata/metadata_qiime_final.tsv \
  --p-formula "species + locality + microhabitat + sex + terrestriality + gutregion + diet" \
  --p-permutations 9999 \
  --p-n-jobs 10 \
  --o-visualization Adonis_LC3/adonis-weighted_unifrac-all.qzv

