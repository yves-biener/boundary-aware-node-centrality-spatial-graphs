import networkx as nx
import squidpy as sq

def merfish():
    adata = sq.datasets.merfish()
    adata = adata[adata.obs.Bregma == -9].copy()
    return adata


def mibitof():
    adata = sq.datasets.mibitof()
    return adata


def visium():
    adata = sq.datasets.visium_hne_adata()
    return adata


def spatial_graph(adata):
    sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
    return nx.from_scipy_sparse_array(adata.obsp['spatial_distances'])
