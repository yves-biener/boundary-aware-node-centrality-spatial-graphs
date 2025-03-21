import networkx as nx
import squidpy as sq

def merfish():
    """
    Merfish dataset from `squidpy`.
    """
    adata = sq.datasets.merfish()
    adata = adata[adata.obs.Bregma == -9].copy()
    return adata


def mibitof():
    """
    Mibitof dataset from `squidpy`.
    """
    adata = sq.datasets.mibitof()
    return adata


def spatial_graph(adata):
    """
    Generate the spatial graph using delaunay for the given `adata`.
    `adata` will contain the calculated spatial graph contents in the keys
    `adata.obps['spatial_distances']` and `adata.obsm['spatial']` afterwards too.
    @return [Graph] generated networkx graph from adata['spatial_distances']
    """
    sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
    return nx.from_scipy_sparse_array(adata.obsp['spatial_distances'])
