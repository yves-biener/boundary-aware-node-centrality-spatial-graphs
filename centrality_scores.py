import numpy as np
import squidpy as sq

# merfish (non-grid dataset)
adata = sq.datasets.mibitof()

# connectivity matrix from spatial coordinates
sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
print(adata)

sq.pl.spatial_scatter(
    adata,
    shape=None,
    connectivity_key="spatial_connectivities",
    save='./mibitof/Spatial_scatter',
)

# seems to be the exact same
sq.pl.spatial_segment(
    adata,
    library_id=["point8"],
    seg_cell_id="cell_id",
    color="Cluster",
    library_key="library_id",
    title=["point8"],
    save='./mibitof/Spacial_segment'
)

# TODO: where is the 'edge'?
# print(adata.obsp["connectivities"])
# print(adata.obsp["distances"])
