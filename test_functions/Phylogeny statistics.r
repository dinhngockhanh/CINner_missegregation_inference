# phylogeny statistics testing
# install.packages("ape")
# install.packages("phyloTop")

# ==========================================================Package: ape
# ======================================================================
library(ape)
# ------------------------------------------------------input: two trees
# Topology distance between two trees:
# Twice the number of internal branches defining different bipartitions of the tips
# Topology distance: https://rdrr.io/cran/ape/man/dist.topo.html
# Need same number of tips
ta <- rtree(30, rooted = FALSE)
tb <- rtree(30, rooted = FALSE)
dist.topo(ta, ta) # 0
dist.topo(ta, tb) # unlikely to be

# -------------------------------------------------------input: one tree
# Compute depth and Heights of Nodes and Tips
# * It returns vectors, but maybe we could get average of them *
# node.depth: https://rdrr.io/cran/ape/man/node.depth.html
# depth of a node
node.depth(c, method = 1)
# depth of a node by branch lengths
node.depth.edgelength(c)
# heights of nodes and tips as plotted by a phylogram or a cladogram.
node.height(c, clado.style = FALSE)


# =====================================================Package: phyloTop
# ======================================================================
library(phyloTop)
# -------------------------------------------------------input: one tree
# Mean size of ladders in the tree
# avgLadder:
avgLadder(c, normalise = FALSE)

# Number of Pitchforks
# pitchforks: https://rdrr.io/cran/phyloTop/man/pitchforks.html
pitchforks(c, normalise = FALSE)

# Number of Cherries
# cherries: https://rdrr.io/cran/phyloTop/man/cherries.html
cherries(c)

# IL Number
# Computes the number of internal nodes with a single tip child
# https://rdrr.io/cran/phyloTop/man/ILnumber.html
tree <- rtree(10)
plot(tree)
ILnumber(tree)

# Colless Number (Measurement of balance)
# colless.phylo: https://rdrr.io/cran/phyloTop/man/colless.phylo.html
colless.phylo(tree, normalise = TRUE)

# Sackin Index
# measuring the sum of path lengths between leaves and the root
# Index for measuring the balance
# Compute the Sackin index of the tree
# https://rdrr.io/cran/phyloTop/man/sackin.phylo.html
## Sackin index of a random tree with 10 tips:
sackin.phylo(rtree(10))
## normalised Sackin index:
sackin.phylo(rtree(10), normalise = TRUE)

# Stairs:
# return 1: the proportion of subtrees that are imbalanced (i.e. subtrees where the left child has more tip descendants than the right child, or vice versa)
#        2: the average of all the min(l,r)/max(l,r) values of each subtree, where l and r are the number of tips in the left and right children of a subtree.
# https://rdrr.io/cran/phyloTop/man/stairs.html
stairs(rtree(20))

# Width: number of nodes at each depth of the tree
# https://rdrr.io/cran/phyloTop/man/widths.html
## Find the node widths in a random tree with 10 tips:
tree <- rtree(10)
tree$edge.length <- rep(1, 18) # to make it easier to see the width and depths in the plot
plot(tree)
widths(tree)
