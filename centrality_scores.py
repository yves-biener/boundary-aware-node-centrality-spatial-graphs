import squidpy as sq

adata = sq.datasets.merfish()

# connectivity matrix from spatial coordinates
sq.gr.spatial_neighbors(adata)

#print(adata.obs.keys())

# calculate centrality scores
sq.gr.centrality_scores(adata, "Cell_class")
# visualize results
sq.pl.centrality_scores(adata, "Cell_class", save='./merfish/Cell_class')

# --- DATASET: merfish ---
# NOTE: all available keys
# 'Cell_ID', 'Animal_ID', 'Animal_sex', 'Behavior', 'Bregma', 'Centroid_X', 'Centroid_Y', 'Cell_class', 'Neuron_cluster_ID', 'batch'
#
# NOTE: the following keys do not give results:
# 'Bregma', 'Cell_ID'
#
# NOTE: the following keys give results (not sure how ~good~ they are):
# - 'Cell_class'
# - 'batch'
# - 'Neuron_cluster_ID'

# --- DATASET: mibitof ---
# NOTE: all available keys
# 'row_num', 'point', 'cell_id', 'X1', 'center_rowcoord', 'center_colcoord', 'cell_size', 'category', 'donor', 'Cluster', 'batch', 'library_id'
#
# NOTE: the following keys do not give results:
# 'point', 'cell_id', 'cell_size', 'category', 'library_id'
#
# NOTE: the following keys give results (not sure how ~good~ they are):
# - 'Cluster'
# - 'batch'
# - 'donor'
