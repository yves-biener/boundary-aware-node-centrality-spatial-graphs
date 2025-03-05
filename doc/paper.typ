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
  #lorem(80)
 ],
 enable_bibliography: true,
 enable_figure_outline: false,
 enable_table_outline : false,
 doc
)

= Motivation
#lorem(80)
// TODO: explain the following aspects:
// - why measurements at the boundaries are usually wrong when compared to the real world
// - improve by adjusting the computated values through the other existing ones
// - working with non-grid datasets (i.e. mibitof and merfish) BIB: both could be cited too

= Centrality measurements

// TODO: what are centrality measurements and how are they calculated (take the we are also plann on improving on)
// - what is the reason then for nodes at the border to become less (or more) potent

== Closeness centrality

== Pagerank centrality
// TODO: Pagerank does not seem to be too effected by the boundary

// TODO:
// BIB: google paper for the pagerank algorithm?

= Correction for boundary

// TODO: explain what we did step by step and how we try to improve the corresponding calculated centrality measurements
// - Boundary box (gift wrapping algorithm)
// - Centrality calculation
// - graph creation with the relation for distance to bounding-box <-> centrality of each node
// - function fitting as a (I)LP
//   - slope, breaking point, constant (as variables)
// - resulting correction to be applied to the centrality values of the corresponding nodes from the relationship
// -> new centrality values
// -> detection which nodes are effected by the boundary

= Usecases

// TODO: What can this approach be used for?
// - Detect how much a metric may be effected by the boundary of a spatial graph
// - Detect the nodes which are effected by the boundary (maybe even with a given threshhold)

= Conclusion

// NOTE: this should serve as a report for the master thesis based on this work
// - what can also be done on top of the existing work?
//   - improvement of algorithms (i.e. runtime improvements - *very unlikely*)
//   - improvement of metric algorithms to use (i.e. which are more or less effected by the boundary)
//   - quantification of metrics (i.e. can a metric be used to make educated guesses about the corresponding cells, etc.)
