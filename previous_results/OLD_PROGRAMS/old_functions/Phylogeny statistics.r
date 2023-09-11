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
plot(ta)
plot(tb)
dist.topo(ta, ta) # 0
dist.topo(ta, tb) # unlikely to be

# -------------------------------------------------------input: one tree
# Compute depth and Heights of Nodes and Tips
# * It returns vectors, but maybe we could get average of them *
# node.depth: https://rdrr.io/cran/ape/man/node.depth.html
# depth of a node
node.depth(ta, method = 1)
# depth of a node by branch lengths
node.depth.edgelength(ta)
mean(node.depth.edgelength(ta))
mean(node.depth.edgelength(tb))
mean(node.depth.edgelength(rtree(30, rooted = FALSE)))
mean(node.depth.edgelength(rtree(30, rooted = FALSE)))
mean(node.depth.edgelength(rtree(33, rooted = FALSE)))
mean(node.depth.edgelength(rtree(30, rooted = FALSE)))
mean(node.depth.edgelength(rtree(30, rooted = FALSE)))
mean(node.depth.edgelength(rtree(30, rooted = FALSE)))


# heights of nodes and tips as plotted by a phylogram or a cladogram.
node.height(ta, clado.style = FALSE)
mean(node.height(ta, clado.style = FALSE))
mean(node.height(tb, clado.style = FALSE))

# =====================================================Package: phyloTop
# ======================================================================
library(phyloTop)
# -------------------------------------------------------input: one tree
# Mean size of ladders in the tree
# avgLadder:
avgLadder(ta, normalise = FALSE)

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
# Inspect the internal nodes, partitioning the tips that descend from them into groups of sizes r and s,
# and computes the sum of absolute values ∣r − s∣ for all nodes.
# colless.phylo: https://rdrr.io/cran/phyloTop/man/colless.phylo.html
colless.phylo(tree, normalise = TRUE)
colless.phylo(tree, normalise = FALSE)
17 / 0.4722222
14 / 0.3888889

tree <- rtree(5)
plot(tree)
sackin.phylo(tree, normalise = TRUE)
sackin.phylo(tree, normalise = FALSE)


tree <- rtree(8)
plot(tree)
sackin.phylo(tree, normalise = TRUE)
sackin.phylo(tree, normalise = FALSE)

tree <- rtree(10)
plot(tree)
sackin.phylo(tree, normalise = TRUE)
sackin.phylo(tree, normalise = FALSE)

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
stairs(rtree(20))[1]
stairs(ta)

# Width: number of nodes at each depth of the tree
# https://rdrr.io/cran/phyloTop/man/widths.html
## Find the node widths in a random tree with 10 tips:
tree <- rtree(10)
tree$edge.length <- rep(1, 18) # to make it easier to see the width and depths in the plot
plot(tree)
widths(tree)
maxWidth(tree)

tree <- ape::read.tree(text = "((((,),),(,)),(((,),),(,)));")
tree$edge.length <- rep(1, 18)
plot(tree)
areaPerPairI(tree)
tree <- ape::read.tree(text = "(((,),(,)),(,),(,));")
colPlaLab(tree, method = "binary")

tree <- rtree(20)
plot(tree)
avgLadder(tree)
ladderSizes(tree)
# note that the ladders can be highlighted in a plot using ladderShow:
ladderShow(tree)
stairs(tree)
19 * 0.579

factorial(2000) * 0.6

tree2 <- rtree(2000)
stairs(tree2)
l <- subtrees(tree2)
l
a <- subtrees(rtree(20))
a[2]

install.packages("apTreeshape")
library(apTreeshape)
tpda <- rtreeshape(1, tip.number = 70)
colless(tpda[[1]]) / 70

merge <- matrix(NA, 3, 2)
merge[, 1] <- c(-3, -1, 2)
merge[, 2] <- c(-4, -2, 1)
tree1 <- treeshape(merge)
merge[, 1] <- c(-1, -3, 1)
merge[, 2] <- c(-2, -4, 2)
tree2 <- treeshape(merge)
plot(tree1, tree2)
all.equal(tree1, tree2)
all.equal(tree1, tree2, height = TRUE)

tree <- rtree(10)
tree$edge.length <- rep(1, 18) # to make it easier to see the width and depths in the plot
plot(tree)
widths(tree)
width.multiPhylo(tree)
summary(ryule(100))

install.packages("treebalance")
library(treebalance)
