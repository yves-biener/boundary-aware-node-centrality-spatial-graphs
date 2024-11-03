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

= Contraction algorithm (Krager's Algorithm)
The following Monte-Carlo algorithm by David Karger describes a simple randomized approach to determine a minimum cut of a connected graph. The main idea for this approach is based on the contraction of edges in the graph to merge both verticies of an edge into a new merged verticy. The corresponding edges to the merged nodes will point to the merged node instead. This results in a multigraph (multiple different edges with the same start and end node). The notation $G\/(u,v)$ is used to denote the contraction of the edge $(u,v)$ in the (multi-)graph $G$.

In the process of the algorithm the number of verticies is reduced until there are only two left, which describe a cut of the original graph. Both potentially merged nodes describe two subsets $S$ and $T$ of the cut-set $C = (S,T)$ with ${(u,v) in | u in S, v in T}$. @contract

```
procedure contract(G = (V,E)):
 while |V| > 2
   choose e = (u,v) ∈ E u.a.r
   /* Contract operation */
  G := G/(u,v)
 C := the only remaining cut in G
 return C
```

There is no guarantee that the determined cut by `contract(G)` is minimal. However it can be proved that `contract(G)` finds a minimal cut with probability of at least $2/n(n-1)$.

$ Pr[C "is minimal cut"] >= Pr[C = K] \
 = product_(i=2)^(n-2) (1- Pr[e_i in K | and_(1<=j<i) e_j in.not K]) \
 >= product_(i=1)^(n-2) (1- 2/(n-i+1)) \
 = product_(i=1)^(n-2) ((n-i-1)/(n-i+1)) \
 = (n-2)/n * (n-3)/(n-1) * (n-4)/(n-2) * ... * 3/5 * 2/4 * 1/3 \
 = 2/(n(n-1)) qed $ <contract-expected-prob>

With @contract-expected-prob the expected number of repetitions until `contract(G)` finds a minimal cut in $G$ is at most $n(n-1)/2$ (_geometric distribution_).

Simple implementations can implement the `contract(G)` algorithm in $O(n^2)$ such that with the repetitions to find a minimal cut the expected runtime in $O(n^4)$.

= Stoer-Wagner algorithm
This deterministic algorithm finds the minimal cut by means of an iterative process. In each iteration phase a cut is determined from which the minimal is returned. The iterative process is started by an arbritary starting node $a$. In each Iteration a subset $A$ containing $a$ in the beginning is created to which the most tightly connected node to $A$ is added until $A = V$. Formally this can be written as:

$z in.not A | w(A,z) = max{w(A,y) | y in.not A}$ where $w(A,y)$ is the sum of the weights of all edges between $A$ and $y$.

From the subset $A$ the cut of the phase is determined by using the last two remaining nodes of the phase. In the end of each iteration these two nodes $s$ and $t$ are merged. When merging the nodes all edges between $s$ and $t$ are removed and all remaining edges are replaced by new edges with the corresponding sum of their weights. This way the algorithm has $n - 1$ iterations.

```
procedure stoer_wagner(G = (V,E), s):
 min_cut = MAX
 while |V| > 1
  /* minimal cut phase */
  A := {s}
  while A != V
   add most tidly connected vertex to A -- this needs to be explained better
  cut_of_phase := cut of the last merge -- this is also not very perceise
  G := G/(A[-1],A[-2])
  if min_cut > cut_of_phase
   min_cut := cut_of_phase
 return min_cut
```

The cut determined by `stoer_wagner(G, s)` is minimal. A minimal cut $C$ separates the graph $G$ into a pair of distinct subsets $A$ and $B$. For the merging of the last two remaining nodes of each phase there are two possible outcomes. Both nodes belong to either $A$ or $B$ or alternatively each node is in either $A$ or $B$. In the first case there is no minimal cut and the replacement of the nodes does not change the minimal cut. On the other hand the second case the minimal cut is found and as the connecting edge, which leads to a minimal cut is removed, a more minimal cut cannot be found.

= Conclusion
By embracing simplicity and harnessing the power of heuristics, randomized algorithms offer a compelling alternative to their deterministic counterparts, paving the way for more efficient and flexible solutions.
