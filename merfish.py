import numpy as np

from src import dataset
from src import metrics
from src import plot
from src import quantification

# NOTE due to the many solvings of the LP (i.e. `quantification.fit_piece_wise_linear`) the runtime for calculating all solutions may be pretty big

adata = dataset.merfish()
G = dataset.spatial_graph(adata)

pos = adata.obsm['spatial']
edges = plot.edges(adata)

# closeness
closeness = metrics.closeness(G)
plot.Graph.plot(G, pos, edges, closeness, 'Merfish with closeness centrality', './figures/merfish/spatial_graph.svg')

# create relationship
data = plot.Quantification.data(pos, closeness)
plot.Quantification.plot(data, None, None, 'Closeness centrality', './figures/merfish/quantification/closeness.svg')
d = data[:, 0]
C = data[:, 1]
# create model
m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

edges = plot.edges(adata)
plot.Graph.plot(G, pos, edges, plot.Quantification.opt_data(pos, closeness, b_opt), 'Merfish with closeness centrality', './figures/merfish/spatial_graph_c0.svg')

# correction of closeness
edges = plot.edges(adata)
corrected_closeness = metrics.Correction.correct(pos, closeness, m_opt, c0_opt, b_opt)
plot.Graph.plot(G, pos, edges, corrected_closeness, 'Merfish with corrected closeness centrality', './figures/merfish/spatial_graph_corrected.svg')
# corrected relationship
corrected_data = plot.Quantification.data(pos, corrected_closeness)
plot.Quantification.plot(corrected_data, None, None, 'Corrected closeness centrality', './figures/merfish/quantification/closeness_corrected.svg')

d_curve = np.linspace(min(d), max(d), 500)
C_curve = np.piecewise(
    d_curve,
    [d_curve <= b_opt, d_curve > b_opt],
    [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
)
plot.Quantification.plot(data, d_curve, C_curve, 'Closeness centrality', './figures/merfish/quantification/closeness_fitted.svg')

# pagerank
# for i in range(0, 10):
#     # for all alpha values 0.0 to 0.9 (increments of 0.1)
#     alpha = i / 10.0
#     G = dataset.spatial_graph(adata)
#     pagerank = metrics.pagerank(G, alpha)
#     data = plot.Quantification.data(pos, pagerank)
#     plot.Quantification.plot(data, None, None, f"Pagerank centrality (alpha = {alpha})", f"./figures/merfish/quantification/pagerank_alpha{alpha}.svg")
#     d = data[:, 0]
#     C = data[:, 1]
#     m_opt, c0_opt, b_opt = quantification.fit_piece_wise_linear(d, C)

#     d_curve = np.linspace(min(d), max(d), 500)
#     C_curve = np.piecewise(
#         d_curve,
#         [d_curve <= b_opt, d_curve > b_opt],
#         [lambda x: m_opt * x + c0_opt, lambda x: m_opt * b_opt + c0_opt]
#     )
#     plot.Quantification.plot(data, d_curve, C_curve, f"Pagerank centrality (alpha = {alpha})", f"./figures/merfish/quantification/pagerank_alpha{alpha}_fitted.svg")
