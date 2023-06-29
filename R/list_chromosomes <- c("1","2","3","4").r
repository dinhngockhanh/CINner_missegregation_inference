list_chromosomes <- c("1", "2", "3", "4")
list_tried_shuffles <- list()
list_tried_shuffles[[1]] <- list_chromosomes
list_chromosomes_new <- list_chromosomes
list_chromosomes_new
list_tried_shuffles
list_chromosomes_new %in% list_tried_shuffles[[length(list_tried_shuffles)]]
list_chromosomes_new <- sample(list_chromosomes, length(list_chromosomes), replace = FALSE)
list_tried_shuffles[[length(list_tried_shuffles) + 1]] <- list_chromosomes_new
list_tried_shuffles

list_chromosomes <- c("1", "2", "3", "4")
list_tried_shuffles <- list()
list_tried_shuffles[[1]] <- list_chromosomes
for (i in 1:10) {
    #   Find a new permutation of chromosomes
    list_chromosomes_new <- list_tried_shuffles[[length(list_tried_shuffles)]]
    while (all(list_chromosomes_new %in% list_tried_shuffles)) {
        list_chromosomes_new <- sample(list_chromosomes, length(list_chromosomes), replace = FALSE)
    }
    list_tried_shuffles[[length(list_tried_shuffles) + 1]] <- list_chromosomes_new
    print(list_chromosomes_new)
}
list_chromosomes_new == list_tried_shuffles[[length(list_tried_shuffles)]]
list_tried_shuffles
list_tried_shuffles[[1]]:list_tried_shuffles[[3]]
1:4 == c(1, 2, 3, 4)
1:4 %in% c(1, 3, 2, 4)
all(1:4 == c(1, 2, 3, 4))
(list_chromosomes_new %in% list_tried_shuffles)[length(list_chromosomes_new)]
list_chromosomes_new <- list()
list_chromosomes_new[[3]] <- c("1", "2", "2", "4")
list_tried_shuffles[[2]] <- c("1", "3", "2", "4")
list_chromosomes_new
list_tried_shuffles

n <- 24
list_chromosomes <- c("1", "2", "3", "4")
list_tried_shuffles <- list()
list_chromosomes_new <- list()
list_tried_shuffles[[1]] <- list_chromosomes
for (i in 1:(n - 1)) {
    #   Find a new permutation of chromosomes
    list_chromosomes_new[[length(list_chromosomes_new) + 1]] <- list_tried_shuffles[[length(list_tried_shuffles)]]
    while ((list_chromosomes_new %in% list_tried_shuffles)[length(list_chromosomes_new)]) {
        list_chromosomes_new[[length(list_chromosomes_new)]] <- sample(list_chromosomes, length(list_chromosomes), replace = FALSE)
    }
    list_tried_shuffles[[length(list_tried_shuffles) + 1]] <- list_chromosomes_new[[length(list_chromosomes_new)]]
    print(list_chromosomes_new[[length(list_chromosomes_new)]])
    print(i)
}
