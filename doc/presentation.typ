#let heading_sep_color = rgb("#04316A")

#set page(
 paper: "presentation-16-9",
 margin: (x: 16pt, top: 24pt),
 footer: stack(
  spacing: 0.5em,
  line(length: 100%, stroke: 2pt + heading_sep_color),
  grid(
   columns: (1fr, 1fr),
   align(left)[
    #text(12pt, fill: luma(180))[
     #datetime.today().display() Yves Biener -- Boundary-aware node centralities for spatial graphs
    ]
   ],
   align(right)[
    #text(12pt, fill: luma(180))[
     #context counter(page).display("1")
    ]
   ],
  ),
 ),
)
#set text(size: 22pt)

// FIX: this is in the most recent version bugged
// #show heading.where(outlined: false): heading.with(bookmarked: false)
#show heading: it => [
 #stack(
  spacing: 0.5em,
  grid(
   columns: (4fr, 1fr),
   align(left)[
    #it
   ],
   align(right)[
    #image("/fau.png", height: 28pt)
   ]
  ),
  line(length: 100%, stroke: 2pt + heading_sep_color),
 )
]

// link highlighting
#show link: set text(blue)
#show link: underline

#show raw.where(block: true): box.with(
 fill: luma(240),
 inset: 10pt,
 radius: 4pt,
)

// make text size of figure caption's smaller
#show figure: set text(size: 12pt)

// --------------------------------------------------------------------------------------------------------------------
// NOTE: do not show this heading in the toc
#heading(outlined: false, level: 2)[Boundary-aware node centralities for spatial graphs]
#grid(
  columns: (1.3fr, 1fr),
  image("./figures/merfish/Spatial_scatter.png"),
  image("./figures/merfish/spatial_coordinates.svg")
)
// TODO: - merfish and mibitof images from squidpy

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
// don't show subheadings below the first level
#outline(depth: 1)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Centrality

- Determine a value describing a node's *importance* in a (directed or undirected) graph
- The importance of a node can be described by certain properties, i.e.
  - edge count (in-/out-degree of a node)
  - shortes paths to other nodes (with or without weighted paths)
  - connectivity in a graph (i.e. number of articulation points, k-connectivity graphs, etc.)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Closeness Centrality]

- A node is central if its close (shortest path distance) to all other nodes (apart from itself):
- Normalized *Closeness Centrality*:

$ C(u) = (n - 1) * (sum_(v in V without {u}) d(v,u))^(-1) = ("arithmetic-mean"_(v in V without {v}) d(v,u))^(-1) $

#text(fill: rgb(0, 255, 0))[+] easy to calculate\
#text(fill: rgb(0, 255, 0))[+] easy to understand\
#text(fill: rgb(0, 255, 0))[+] considers all nodes\
#text(fill: rgb(255, 0, 0))[-] is only meaningful in *strongly connected graphs*

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[The Problem of the Boundaries]

// TODO reformulate!
However, when using node centrality measures to quantify node importance in such spatial graphs, they tend to prioritize nodes in the center of the graph and de-prioritize nodes that are close to the boundary of the tissue section. This is a problem, because the boundary is very often an arbitrary artifact of the tissue sample collection protocol (e.g, a small skin section was cut out of an arbitrary section of a larger skin area of interest) and should hence not affect node importance quantification.

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Relation between Boundary and Centrality

// TODO boundary <-> centrality
// - for closeness centrality measurement

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Mibitof]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Mibitof Spatial Graph],
    image("./figures/mibitof/spatial_graph.svg")
  ),
  figure(
    caption: [Mibitof Closeness Function Fitting Result],
    image("./figures/mibitof/quantification/closeness.svg"),
  ),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Merfish]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Mibitof Closeness Function Fitting Result],
    image("./figures/merfish/spatial_graph.svg", height: 80%),
  ),
  figure(
    caption: [Mibitof Closeness Function Fitting Result],
    image("./figures/merfish/quantification/closeness.svg"),
  ),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Function Fitting

- Find a function which describes (an "optimal") relationship between the centrality measure and the distance to the bounding box #footnote(link("https://en.wikipedia.org/wiki/Gift_wrapping_algorithm"))
- *piece-wise* linear function:
  - $f(x) = m * x + t$
  - $g(x) = c$
  - intersection point $c_0$ between the $f(x)$ and $g(x)$

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Mibitof]
// only show the relationship with the result of the ILP
#figure(
  caption: [Mibitof Closeness Function Fitting Result],
  image("./figures/mibitof/quantification/closeness_fitted.svg", height: 75%),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Potentially effected centrality values

- Assume that $c_0$ describes the intersection between effected and non-effected centrality measures
  - centrality measures for nodes that have a shorter distance to the border (and are therefore closer to the boundary) are more likely to be effected by the boundary of the graph

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Mibitof]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Mibitof Spatial Graph with boundary-effected nodes],
    image("./figures/mibitof/spatial_graph_c0.svg")
  ),
  figure(
    caption: [Mibitof Closeness Function Fitting Result],
    image("./figures/mibitof/quantification/closeness_fitted.svg"),
  ),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false)[Merfish]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Merfish Spatial Graph with boundary-effected nodes],
    image("./figures/merfish/spatial_graph_c0.svg"),
  ),
  figure(
    caption: [Merfish Closeness Function Fitting Result],
    image("./figures/merfish/quantification/closeness_fitted.svg"),
  ),
)

// NOTE the corresponding *m* of both fitted functions is pretty similar -> is this a problem of the optimization problem formulation?

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Solution for the fitted piece-wise function

- If the assumption holds that $c_0$ describes the intersection between effected and non-effected centrality measures
  - all points of $g(x)$ are non-effected centrality measures meaning that the intersection point should be as close to the border as possible for any given centrality measure to be independent of the boundary
  - limit the corresponding nodes to a sub-set which we can "trust" for further implications
- Depending on the stepness $m$ of the linear function $f(x)$
  - corresponding corrections could be applied to the centrality measure with probability in relation to the stepness $m$

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false, level: 2)[What about other Centrality Measures? - Mibitof]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Mibitof Closeness Function Fitted],
    image("./figures/mibitof/quantification/closeness_fitted.svg"),
  ),
  figure(
    caption: [Mibitof PageRank Function Fitted],
    image("./figures/mibitof/quantification/pagerank_fitted.svg"),
  ),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
#heading(outlined: false, level: 2)[What about other Centrality Measures? - Merfish]
#grid(
  columns: (1fr, 1fr),
  figure(
    caption: [Merfish Closeness Centrality Fitted],
    image("./figures/merfish/quantification/closeness_fitted.svg"),
  ),
  figure(
    caption: [Merfish PageRank Centrality Fitted],
    image("./figures/merfish/quantification/pagerank_fitted.svg"),
  ),
)

// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Future Work
// - is there a certain property associated with centrality measurement algorithms that cause the results to be either dependent or independent of the graph's boundaries? (i.e. page rank or betweeness)

- Can this approach be used to quantify results extracted from implication made by centrality measures?
  - Is a the linear part ($f(x)$) always a linear function or could it even be logarithmic? Maybe certain centrality measures should use other different function to fit / optimize for depending on their importance indicators.
- Are there certain properties associated with algorithms (i.e. closeness, betweeness, page rank, etc.) that make them more or less dependent or even independent of the boundary?
- And the other way around: Are there certain properties on graphs that indicate better use for particular algorithms?


// --------------------------------------------------------------------------------------------------------------------
#pagebreak()
= Appendix

The following slides are used for specific questions if details are necessary, but are not planned on being used for the actual presentation. However readers may also want to take a look at the appendix slides for further details. See the last slides for a collection of all links.

#pagebreak()
== Centrality Measures

#pagebreak()
=== Degree centrality

#box(
 fill: luma(240),
 inset: 10pt,
 radius: 4pt,
)[
  The *in-degree centrality* of a node in a directed graph is given by the in-degree of each node.

  The *out-degree centrality* (or *degree-centrality*) for undirected graphs are defined analogously.
]

#text(fill: rgb(0, 255, 0))[+] easy to compute and understand\
#text(fill: rgb(255, 0, 0))[-] considers only direct neighbourhood of a node

#pagebreak()
=== PageRank centrality

#box(
 fill: luma(240),
 inset: 10pt,
 radius: 4pt,
)[
  A node is important if it is *highly linked* or if it is *linked from other important nodes* that do not link many other nodes.
]

- Column-based (weighted) adjacency matrix $A$ with $A_(i.j) = "weight of edge from" j "to" i$
- Scale all columns to $A$ to unit sum: $A_(i.j) = A_(i.j) / (sum_(i')A_(i'.j))$
- Add *teleportation edges* ($i,j$) of weight $n^(-1)$ between all nodes $i$ and $j$.
- Random walker #footnote(link("https://en.wikipedia.org/wiki/Random_walk")) walks along $A$ with probability $d in (0,1]$ and along teleportation edges with probability $1 - d$.
- Results in transition probability matrix:
  $ A' = d A + (1 - d) vec(1/n, ..., 1/n, delim: "[") 1^T $

#pagebreak()
=== Betweeness centrality

#box(
 fill: luma(240),
 inset: 10pt,
 radius: 4pt,
)[
  A node is important if it is on the *shortest path* *between* many pairs of *other nodes*.
]

- Non-Normalized *Betweeness Centrality*:
  - $sigma_(s,t)$: Number of shortest paths from $s$ to $t$
  - $sigma_(s,t)(u)$: Number of shortest paths from $s$ to $t$ that pass through $u$
$ B(u) = sum_(s,t in V without {u}) sigma_(s,t)(u) * (sigma_(s,t))^(-1) $

#text(fill: rgb(0, 255, 0))[+] easy to understand\
#text(fill: rgb(0, 255, 0))[+] even though complex possible in $Theta(n m + n^2 + log n)$ run-time complexity #footnote(link("https://doi.org/10.1080/0022250X.2001.9990249"))

#pagebreak()
= Links

- #link("https://doi.org/10.1101/2020.01.17.909796")[MIBI-TOF Dataset]
- #link("https://doi.org/10.1126/science.aau5324")[MERFISH Dataset]
- #link("https://www.bionets.tf.fau.de/files/2024/08/Node-centrality-project.pdf")[Project Description - FAU Bionet]
- #link("https://github.com/yves-biener/boundary-aware-node-centrality-spatial-graphs")[Project Repository - Github]
- #link("https://en.wikipedia.org/wiki/Gift_wrapping_algorithm")[Gift wrapping algorithm - Wikipedia]
- #link("https://en.wikipedia.org/wiki/PageRank")[PageRank algorithm - Wikipedia]
- #link("https://en.wikipedia.org/wiki/Random_walk")[Random Walk - Wikipedia]
- #link("https://doi.org/10.1080/0022250X.2001.9990249")[Fast Betweeness Centrality Algorithm]
