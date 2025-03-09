import numpy as np

from src import dataset
from src import metrics
from src import plot
from src import quantification

adata = dataset.merfish()
G = dataset.spatial_graph(adata)

closeness = metrics.closeness(G)
pagerank = metrics.pagerank(G)

pos = adata.obsm['spatial']
edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, closeness, 'Merfish with closeness centrality', './figures/merfish/spatial_graph.svg')

data = plot.Quantification.data(pos, closeness)
plot.Quantification.plot(data, None, None, 'Closeness centrality', './figures/merfish/quantification/closeness.svg')
d = data[:, 0]
C = data[:, 1]
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, plot.Quantification.opt_data(pos, closeness, b_opt), 'Mibitof with closeness centrality', './figures/merfish/spatial_graph_c0.svg')

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Closeness centrality', './figures/merfish/quantification/closeness_fitted.svg')

# pagerank
data = plot.Quantification.data(pos, pagerank)
plot.Quantification.plot(data, None, None, 'Pagerank centrality', './figures/merfish/quantification/pagerank.svg')
d = data[:, 0]
C = data[:, 1]
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Pagerank centrality', './figures/merfish/quantification/pagerank_fitted.svg')
