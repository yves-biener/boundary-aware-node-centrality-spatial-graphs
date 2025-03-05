import numpy as np

from src import dataset
from src import metrics
from src import plot
from src import quantification

adata = dataset.mibitof()
G = dataset.spatial_graph(adata)

closeness = metrics.closeness(G)
pagerank = metrics.pagerank(G)

pos = adata.obsm['spatial']
edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, closeness, 'Mibitof with closeness centrality', './figures/mibitof/spatial_graph.svg')

# closeness
data = plot.Quantification.data(pos, closeness)
d = data[:, 0]
C = data[:, 1]
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Closeness centrality', './figures/mibitof/quantification/closeness.svg')

# pagerank
data = plot.Quantification.data(pos, pagerank)
d = data[:, 0]
C = data[:, 1]
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Pagerank centrality', './figures/mibitof/quantification/pagerank.svg')
