import numpy as np

from src import dataset
from src import metrics
from src import plot
from src import quantification

# NOTE due to the many solvings of the LP (i.e. `quantification.fit_piece_wise_linear`) the runtime for calculating all solutions may be pretty big

adata = dataset.mibitof()
G = dataset.spatial_graph(adata)

pos = adata.obsm['spatial']
edges = plot.edges(adata)

closeness = metrics.closeness(G)

plot.Graph.plot(G, pos, edges, closeness, 'Mibitof with closeness centrality', './figures/mibitof/spatial_graph.svg')

# closeness
data = plot.Quantification.data(pos, closeness)
plot.Quantification.plot(data, None, None, 'Closeness centrality', './figures/mibitof/quantification/closeness.svg')
d = data[:, 0]
C = data[:, 1]
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, plot.Quantification.opt_data(pos, closeness, b_opt), 'Mibitof with closeness centrality', './figures/mibitof/spatial_graph_c0.svg')

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Closeness centrality', './figures/mibitof/quantification/closeness_fitted.svg')

# pagerank
for i in range(0, 10):
    alpha = i / 10.0
    G = dataset.spatial_graph(adata)
    pagerank = metrics.pagerank(G, alpha)
    data = plot.Quantification.data(pos, pagerank)
    plot.Quantification.plot(data, None, None, f"Pagerank centrality (alpha = {alpha})", f"./figures/mibitof/quantification/pagerank_alpha{alpha}.svg")
    d = data[:, 0]
    C = data[:, 1]
    m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

    d_curve = np.linspace(min(d), max(d), 500)
    C_curve = np.piecewise(
        d_curve,
        [d_curve <= b_opt, d_curve > b_opt],
        [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
    )
    plot.Quantification.plot(data, d_curve, C_curve, f"Pagerank centrality (alpha = {alpha})", f"./figures/mibitof/quantification/pagerank_alpha{alpha}_fitted.svg")
