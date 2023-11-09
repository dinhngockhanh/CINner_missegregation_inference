# install.packages("ape")



simulate_phylo_with_muts <- function(n_cells, theta_mut, file_name, font_size = 2) {
    library(ape)
    #   Generate a coalescent tree, where each edge ~ Exp(2/(j*(j-1)))
    tree <- rcoal(n_cells)
    #   Generate mutations, where # muts on an edge ~ Poisson(theta_mut * edge_length/2)
    edge_mutcount <- c()
    for (i in 1:length(tree$edge.length)) {
        edge_mutcount[i] <- rpois(1, theta_mut * tree$edge.length[i] / 2)
    }
    edge_labels <- edge_mutcount
    #   Plot the coalescent with mutations
    filename <- paste0(file_name, ".jpeg")
    jpeg(filename, width = 2000, height = 1000)
    plot(tree, show.tip.label = FALSE)
    edgelabels(edge_labels, cex = font_size)
    dev.off()
    #   Return the total mutation count
    return(sum(edge_mutcount))
}

total_mutcount <- simulate_phylo_with_muts(
    n_cells = 20,
    theta_mut = 10,
    file_name = "X",
    font_size = 2
)
print(total_mutcount)
