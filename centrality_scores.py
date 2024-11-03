import squidpy as sq

adata = sq.datasets.mibitof()

# connectivity matrix from spatial coordinates
sq.gr.spatial_neighbors(adata)

# print(adata.obs.keys())
# NOTE: all available keys
# 'row_num', 'point', 'cell_id', 'X1', 'center_rowcoord', 'center_colcoord', 'cell_size', 'category', 'donor', 'Cluster', 'batch', 'library_id'

# NOTE: the following keys do not give results:
# 'point', 'cell_id', 'cell_size', 'category', 'library_id'

# NOTE: the following keys give results (not sure how ~good~ they are):
# - 'batch'
# - 'Cluster'
# - 'donor'

# calculate centrality scores
sq.gr.centrality_scores(adata, "batch")

# visualize results
sq.pl.centrality_scores(adata, "batch", save='./batch')
