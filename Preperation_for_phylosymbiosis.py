######################################
##  Preperation for phylosymbiosis  ##
######################################

## Import library
import numpy as np
import pandas as pd
import skbio as sk
from qiime2 import Artifact
from skbio import TreeNode
from plotnine import *
from os.path import join, abspath
from os import makedirs
from scipy.spatial.distance import squareform, pdist


## Define distance matrices (DMs) file path
dm_dir = abspath('../core-metrics-results_LC3')

dm_fps = {'Unweighted_UniFrac': 'unweighted_unifrac_distance_matrix.qza',
       'Weighted_UniFrac': 'weighted_unifrac_distance_matrix.qza',
       'Jaccard': 'jaccard_distance_matrix.qza',
       'Bray': 'bray_curtis_distance_matrix.qza'}
       
dms = {}


## Load DMs
for metric in dm_fps:
    dm_art = Artifact.load(join(dm_dir, dm_fps[metric]))
    dms[metric] = dm_art.view(sk.DistanceMatrix)


## Load metadata
md_dir = abspath('../metadata')
host_md_fp = join(md_dir, 'metadata_stat_final.tsv')
host_md = pd.read_csv(host_md_fp, sep='\t')


## Subset tree and DMs to same tips
tree_dir = abspath('../host_phylogeny/')
host_tree_fp = join(tree_dir, 'ME2024-tree.newick')
host_tree= sk.read(host_tree_fp, format='newick', into=TreeNode, convert_underscores=False)

host_tips = [x.name for x in host_tree.tips()]

subset = False

if subset:
    host_subset = np.random.choice(host_tips, size=100, replace=False)
else:
    host_subset = host_tips


## Prune full tree to tree with only species with microbiome data 
host_tree_subset = host_tree.shear(host_subset)

len([x.name for x in host_tree_subset.tips()])

host_ids_subset = host_md.loc[host_md['species'].isin(host_subset), 'sample-id']

len(host_ids_subset)


## Using all data
one_per_sp = False

if one_per_sp:
    host_md = host_md.loc[(host_md['sample-id'].isin(host_ids_subset)) &
                          (host_md['sample-id'].isin(unw_dm_subset.ids)),].groupby('species').first()
    host_md =  host_md.loc[(host_md['sample-id'].isin(host_ids_subset)) &
                          (host_md['sample-id'].isin(unw_dm_subset.ids)),].groupby('species').first().reset_index()
    host_ids_subset = list(set(host_ids_subset)) & set(host_md['sample-id'])


## Filter distance matrices:
for dm in dms:
    host_ids_subset_dm = list(set(host_ids_subset) & set(dms[dm].ids))
    dms[dm] = dms[dm].filter(host_ids_subset_dm)


## Parse distance matrices
## Calculate microhabitat Jaccard distances from metadata
microhabitat_cols = ['muddy_mangrove',
             'back_mangrove',
             'open_muddy_mudflat',
             'rocky_mangrove',
             'rocky_shore',
             'sandy_muddy_shore',
             'vegetated_mudflat',
             'terrestrial_environment',
             'freshwater_stream']

host_microhabitat_df = host_md.loc[host_md[microhabitat_cols].sum(axis=1) == 100,
                           ['sample-id'] + microhabitat_cols].dropna()

host_microhabitat_df.set_index('sample-id', inplace=True)

microhabitat_dm = sk.DistanceMatrix(squareform(pdist(host_microhabitat_df.iloc[:, :], metric='jaccard')))
microhabitat_dm.ids = host_microhabitat_df.index




###############################################
##  Build files for phylosymbiosis analysis  ##
###############################################
host_lookup = {x['sample-id']: x['species'] for i, x in host_md.iterrows()}

microhabitat_dm_paired = microhabitat_dm.to_series().reset_index()
microhabitat_dm_paired.columns = ['host_1_id','host_2_id','microhabitat_jaccard']

microhabitat_dm_paired.head()


## Add taxonomic info - Host name
host_lookup = {x['sample-id']: x['species'] for i, x in host_md.iterrows()}

microhabitat_dm_paired['host_1'] = microhabitat_dm_paired['host_1_id'].apply(lambda x: host_lookup[x])
microhabitat_dm_paired['host_2'] = microhabitat_dm_paired['host_2_id'].apply(lambda x: host_lookup[x])

distances = microhabitat_dm_paired.copy()


## Add metdata
host_lookup_genus = {x['sample-id']: x['genus'] for i, x in host_md.iterrows()}
distances['host_1_genus'] = distances['host_1_id'].apply(lambda x: host_lookup_genus[x])
distances['host_2_genus'] = distances['host_2_id'].apply(lambda x: host_lookup_genus[x])


host_lookup_family = {x['sample-id']: x['family'] for i, x in host_md.iterrows()}
distances['host_1_family'] = distances['host_1_id'].apply(lambda x: host_lookup_family[x])
distances['host_2_family'] = distances['host_2_id'].apply(lambda x: host_lookup_family[x])


host_lookup_terrestriality = {x['sample-id']: x['terrestriality'] for i, x in host_md.iterrows()}
distances['host_1_terrestriality'] = distances['host_1_id'].apply(lambda x: host_lookup_terrestriality[x])
distances['host_2_terrestriality'] = distances['host_2_id'].apply(lambda x: host_lookup_terrestriality[x])


host_lookup_sex = {x['sample-id']: x['sex'] for i, x in host_md.iterrows()}
distances['host_1_sex'] = distances['host_1_id'].apply(lambda x: host_lookup_sex[x])
distances['host_2_sex'] = distances['host_2_id'].apply(lambda x: host_lookup_sex[x])


host_lookup_gutregion = {x['SampleID']: x['gutregion'] for i, x in host_md.iterrows()}
distances['host_1_gutregion'] = distances['host_1_id'].apply(lambda x: host_lookup_gutregion[x])
distances['host_2_gutregion'] = distances['host_2_id'].apply(lambda x: host_lookup_gutregion[x])


## check the new table
distances.head()


## Add within/between fields
## For subset graphing
distances['within_family'] = False
distances.loc[distances['host_1_family'] == distances['host_2_family'], 'within_family'] = True

distances['within_species'] = False
distances.loc[distances['host_1'] == distances['host_2'], 'within_species'] = True


## Calcukate patristic distance from tree
patristic_dm = host_tree_subset.tip_tip_distances()


## Define function for data curation
def get_col(df, dm, col1, col2):
    
    values = []

    for i, row in df.iterrows():
        if i % 100000 == 0:
            print(i)
        try:
            values.append(dm[row[col1], row[col2]])
        except:
            values.append(np.nan)

    return(pd.Series(values, index=df.index))


## Add patristic distance to data frame
distances['phylo'] = get_col(distances, patristic_dm, 'host_1', 'host_2')


## Add beta distance to data frame
for metric in dms:
    distances[metric] = get_col(distances, dms[metric], 'host_1_id', 'host_2_id')




##################################################
##  Write files for phylosymbiosis calculation  ##
##################################################

## Jaccard
distances.loc[(distances['host_1'] != distances['host_2']), 'Jaccard'].notna().sum()

out_dir = '../phylosymbiosis_LC3'
makedirs(out_dir, exist_ok=True)
distances.to_csv(join(out_dir, 'phylo.distance_Jaccard_table.csv'))

## Bray
distances.loc[(distances['host_1'] != distances['host_2']), 'Bray'].notna().sum()

out_dir = '../phylosymbiosis_LC3'
makedirs(out_dir, exist_ok=True)
distances.to_csv(join(out_dir, 'phylo.distance_Bray_table.csv'))

## Unweighted_UniFrac
distances.loc[(distances['host_1'] != distances['host_2']), 'Unweighted_UniFrac'].notna().sum()

out_dir = '../phylosymbiosis_LC3'
makedirs(out_dir, exist_ok=True)
distances.to_csv(join(out_dir, 'phylo.distance_Unweighted_UniFrac_table.csv'))

## Weighted_UniFrac
distances.loc[(distances['host_1'] != distances['host_2']), 'Weighted_UniFrac'].notna().sum()

out_dir = '../phylosymbiosis_LC3'
makedirs(out_dir, exist_ok=True)
distances.to_csv(join(out_dir, 'phylo.distance_Weighted_UniFrac_table.csv'))

