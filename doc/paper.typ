#import "config.typ": configuration
#show: doc => configuration(
 title: [Boundary-aware node centralities for spatial graphs],
 cols: 2,
 authors: (
  (
   name: "Yves Biener",
   email: "yves.biener@fau.de",
  ),
  (
   name: "Prof. Dr. David Blumenthal",
   email: "david.b.bluementhal@fau.de",
  ),
 ),
 abstract: [
  Centrality measures are used to determine the importance of nodes in a given graph. Usually they are applied to spatial graphs which are by nature a limited section of the actual tissues that are analysed. However the calculation of the importance values is impacted by the boundaries of the spatial graph itself leading to potentially incorrect values, or values which cannot be used to derive meaningful results back to the tissue section. This project attempts to create a model to identify centrality values which are impacted by the boundary and a potential correction of these values using that model.
 ],
 enable_bibliography: true,
 enable_figure_outline: false,
 enable_table_outline : false,
 doc
)

// TODO add missing logo for the Bionet Lehrstuhl
= Motivation
Node centrality measures such as closeness, betweenness, PageRank centrality are standard tools to quantify the importance of individual nodes in a network. Consequently, they are also widely used to analyze spatial graphs derived from tissue sections (e.g., spatial k-nearest neighbor graphs computed based on spatial omics data). However, when using node centrality measures to quantify node importance in such spatial graphs, they tend to prioritize nodes in the center of the graph and de-prioritize nodes that are close to the boundary of the tissue section. This is a problem, because the boundary is very often an arbitrary artifact of the tissue sample collection protocol (e.g, a small skin section was cut out of an arbitrary section of a larger skin area of interest) and should hence not affect node importance quantification.

For this project the proposed approach has been applied for the merfish @merfish and mibitof @mibitof datasets to test if this approach has potential use as a quantification tool. All resulting graphs and diagrams are available in the #link("https://github.com/yves-biener/boundary-aware-node-centrality-spatial-graphs")[project's repository].

= Centrality measurements
A centrality measurement determines a value describing a node's importance in a (di-)graph. This importance may be described by certain properties, like the edge count, shortes paths to other nodes, distance to closest articulation point, etc. Most commonly used is the _closeness centrality_, which describes the average distance to all other nodes of the graph. This can be formalized as followed:

$
C(u) = (n - 1) * (sum_(v in V without {u}) d(v,u))^(-1)\
= ("arithmetic-mean"_(v in V without {v}) d(v,u))^(-1)
$

This property is good because it creates a relationship between every node of the graph making each node comparable to every other node in the graph. However this approach favors the nodes usually located in the center of a spatial graph because these nodes are the most likely ones to have the shortest distances to every other node in the graph. On the other hand nodes close to the boundary are further away from other nodes of the graph earning them a lower score contributing to this imbalance created by the boundary.

#figure(
  image("./figures/merfish/spatial_graph.svg"),
  caption: [Closeness centrality for merfish dataset],
) <closeness-merfish>

@closeness-merfish shows that closeness centrality measure favors the nodes in the center of the spatial graph as nodes in the center of the spatial graph have a higher score than nodes close to the border of the graph.

Another widely used centrality measure is the PageRank @pagerank algorithm. Originally the algorithm ranked web pages by providing each page a score based on the idea that more important pages are likely to receive more links form other pages. This approach can also be applied to any graph, such that the importance of each node is based on the number and quality of connections to that node.

= Approach
Using the relationship between a given centrality measure (in the case of this project: _closeness_ and _pagerank_) and the distance to the bounding box a model can be derived describing the dependency of each node's centrality value to the boundary. The bounding box is determined through a simple gift wrapping algorithm @giftwrapping.

The model itself is based on a _piece-wise_ linear function $h(x) =$ with 
  + $h(b_0) = c$ // add this to the piece-wise linear function description
  + the slope $f(x) = m x + t$
  + the constant $g(x) = c$
to which each point of the relationship is fitted. For the fitting an elipse based approach is used to ensure that the fitted function describes values as close to all nodes corresponding centrality values as possible.

#figure(
  image("./figures/merfish/quantification/closeness_fitted.svg"),
  caption: [Fitted Model for merfish using closeness centrality],
) <fitted-model>

Each point of the scatter plot (blue in @fitted-model) describes a node's centrality value, while the red function plot in @fitted-model describes the piece-wise linear function $h(x)$ for this dataset and centrality measure.

$b_0$ describes the intersection between the linear part ($f(x)$) and the constant part ($g(x)$). The model assumes that every node having a shorter distance to the boundary than the intersection has centrality values associated which are effected by the boundary. On the other hand nodes that have a greater distance to the bounding box than $b_0$ are not effected by the boundary and can be assumed to be trust-worthy.

#figure(
  image("./figures/merfish/spatial_graph_c0.svg"),
  caption: [Every yellow node is affected by the boundary, while nodes highlighted blue are not affected by the boundary.],
)

This would mean that the model can be used to determine if a given centrality algorithm is more or less prone to being affected by the boundary of the spatial graph compared to other algorithms. Additionally the slope of the linear part $f(x)$ can be used to describe the impact of the boundary on the centarity measure itself. A steeper slope would indicate that the boundary has a bigger impact on the centralities than a shallower slope.

Assuming that $b_0$ represents the intersection between affected and unaffected centrality measures, algorithms which place this intersection closer to the boundary describe fewer points affected by the boundary. Depending on the slope $m$ of the linear function $f(x)$, it may be possible to apply corrections to centrality measures or interpret $m$ as a probability indicating how "trustworthy" a given centrality value is for drawing conclusions about the analyses.

// TODO: explain what we did step by step and how we try to improve the corresponding calculated centrality measurements
// - Boundary box (gift wrapping algorithm)
// - Centrality calculation
// - graph creation with the relation for distance to bounding-box <-> centrality of each node
// - function fitting as a (I)LP
//   - slope, breaking point, constant (as variables)
// - resulting correction to be applied to the centrality values of the corresponding nodes from the relationship
// -> new centrality values
// -> detection which nodes are effected by the boundary

= Model-based correction
For the correction each point which is assumed to be affected by the boundary will be corrected. For each point determine the delta of its associated model value (i.e. the point of the slope for the distance to the boundary of that point) to the constant of the model. Effectively for a point with distance $x$ calculate $delta = m * (b_0 - x)$ and add that $delta$ to the associated point's centrality measure. For the shown merfish dataset @merfish this results in the following correction:

#figure(
 grid(
  columns: (1fr, 1fr),
  image("./figures/merfish/quantification/closeness.svg"),
  image("./figures/merfish/quantification/closeness_corrected.svg"),
  image("./figures/merfish/spatial_graph.svg"),
  image("./figures/merfish/spatial_graph_corrected.svg"),
 ),
 caption: [Left: Original relationship and spatial graph; Right: Applied correction based on model for the merfish dataset],
)

Because the correction is applied to every point that is on the left side of the intersection point $b_0$, there are now centrality values which are effectively bigger than any node was before. Due to the fact that the centrality values have been normalized before (hence the max values of the centrality are between $0.0$ and $1.0$) the corrected values go over $1.0$ and are not normalized to leave the centrality values of the points that are on the right hand side of $b_0$ remain unchanged (from the original normalized centrality).

= Future Work
The shown approach may provide a general way to quantify centrality measures for any given (spatial) graph and may even go further by providing a way to correct the centrality values according to the model. However this still requires verification against more datasets and known results. Can the results be quantified using this approach? Where do implications drawn from centrality measures - which are deemed not "trustworthy" by the model - still hold? Can implications from other publications be validated with this model?

As the relationship on which the model is build on top, is relating an algorithm with a given graph (and its boundary), are there any algorithms which are unaffected (or more likely to be unaffected) by the boundary? Are there certain properties that can be derived from such algorithms? Or is this maybe even a property of the graph that applies to specific centrality measures? For example the *Appendix* section contains the calculated models for the PageRank @pagerank centrality measure for both the merfish @merfish and mibitof @mibitof datasets with various dampening factors. All of them have no linear part of the peace-wise function and indicate that they are completely unaffected by the boundary.

#colbreak()
= Appendix <appendix>
This section contains all the generated graphs and images from the mibitof @mibitof and merfish @merfish datasets used for the model generation with their corresponding values.

#figure(
 grid(
  columns: (1fr, 1fr),
  image("./figures/mibitof/quantification/pagerank_alpha0.0_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.1_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.2_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.3_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.4_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.5_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.6_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.7_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.8_fitted.svg"),
  image("./figures/mibitof/quantification/pagerank_alpha0.9_fitted.svg"),
 ),
 caption: [Pagerank algorithm analyses of mibitof dataset],
)

#figure(
 grid(
  columns: (1fr, 1fr),
  image("./figures/merfish/quantification/pagerank_alpha0.0_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.1_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.2_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.3_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.4_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.5_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.6_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.7_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.8_fitted.svg"),
  image("./figures/merfish/quantification/pagerank_alpha0.9_fitted.svg"),
 ),
 caption: [Pagerank algorithm analyses of merfish dataset],
)
