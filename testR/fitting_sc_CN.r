#---Function to get arm-level CN profiles
#---from bin-level CN profiles
get_arm_CN_profiles <- function(clonal_CN_profiles) {
    clonal_CN_profiles_arm <- vector("list", length(clonal_CN_profiles))
    for (sample in 1:length(clonal_CN_profiles)) {
        clonal_CN_profiles_arm[[sample]] <- vector("list", length(clonal_CN_profiles[[sample]]))
        for (clone in 1:length(clonal_CN_profiles[[sample]])) {
            clone_CN_profile <- clonal_CN_profiles[[sample]][[clone]]
            clone_CN_profile_arm <- data.frame(matrix(ncol = (ncol(clone_CN_profile) + 1), nrow = 0))
            colnames(clone_CN_profile_arm) <- c(colnames(clone_CN_profile), "arm")
            for (chrom in cn_table$Chromosome) {
                for (arm in c("p", "q")) {
                    if (arm == "p") {
                        start <- 1
                        end <- cn_table[cn_table$Chromosome == chrom, ]$Centromere
                    } else if (arm == "q") {
                        start <- 1 + cn_table[cn_table$Chromosome == chrom, ]$Centromere
                        end <- cn_table[cn_table$Chromosome == chrom, ]$Length
                    }
                    locs <- which(
                        clone_CN_profile$chr == chrom & clone_CN_profile$start >= start & clone_CN_profile$end <= end
                    )
                    Min <- as.numeric(names(sort(summary(as.factor(clone_CN_profile[locs, ]$Min)), decreasing = T)[1]))
                    Maj <- as.numeric(names(sort(summary(as.factor(clone_CN_profile[locs, ]$Maj)), decreasing = T)[1]))
                    tmp <- clone_CN_profile[locs, ]
                    state <- Min + Maj
                    copy <- state
                    clone_CN_profile_arm[nrow(clone_CN_profile_arm) + 1, ] <- 0
                    clone_CN_profile_arm$chr[nrow(clone_CN_profile_arm)] <- chrom
                    clone_CN_profile_arm$start[nrow(clone_CN_profile_arm)] <- start
                    clone_CN_profile_arm$end[nrow(clone_CN_profile_arm)] <- end
                    clone_CN_profile_arm$Min[nrow(clone_CN_profile_arm)] <- Min
                    clone_CN_profile_arm$Maj[nrow(clone_CN_profile_arm)] <- Maj
                    clone_CN_profile_arm$state[nrow(clone_CN_profile_arm)] <- state
                    clone_CN_profile_arm$copy[nrow(clone_CN_profile_arm)] <- copy
                    clone_CN_profile_arm$arm[nrow(clone_CN_profile_arm)] <- arm
                }
            }
            clonal_CN_profiles_arm[[sample]][[clone]] <- clone_CN_profile_arm
        }
    }
    return(clonal_CN_profiles_arm)
}

get_each_clonal_CN_profiles <- function(simulations,
                                        arm_level = FALSE,
                                        cn_table = NULL,
                                        bulk = FALSE) {
    #-----------------------------------------Find bin-level CN profiles
    clonal_CN_profiles_all_sims <- list()
    for (i in 1:length(simulations)) {
        simulation <- simulations[[i]]
        clonal_CN_profiles <- simulation$sample$sample_genotype_unique_profile
        sample_genotype_unique <- simulation$sample$sample_genotype_unique
        clonal_CN_populations <- rep(0, length(clonal_CN_profiles))
        for (j in 1:length(clonal_CN_profiles)) {
            clonal_CN_populations[j] <- sum(simulation$sample$sample_clone_ID == sample_genotype_unique[j])
        }
        clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]] <- clonal_CN_profiles
        clonal_CN_profiles_all_sims[["variable=clonal_CN_populations"]][[i]] <- clonal_CN_populations
    }
    #----------------------Convert to arm-level CN profiles if requested
    if (arm_level) {
        clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]] <-
            get_arm_CN_profiles(clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]])
    }
    #--------------------------------Convert to bulk format if necessary
    if (bulk) {
        clonal_CN_profiles_all_sims_new <- list()
        for (i in 1:length(simulations)) {
            clonal_populations <- clonal_CN_profiles_all_sims[["variable=clonal_CN_populations"]][[i]]
            #   Find maximal CN profile (most popular among clones)
            loc_max <- which.max(clonal_populations)
            clonal_CN_profiles_all_sims_new[["variable=maximal_CN_profile"]][[i]] <- clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[loc_max]]
            #   Find average CN profile (among all clones)
            for (clone in 1:length(clonal_populations)) {
                if (clone == 1) {
                    average_CN_profile <- clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[1]]
                    average_CN_profile$Min <- clonal_populations[1] * average_CN_profile$Min
                    average_CN_profile$Maj <- clonal_populations[1] * average_CN_profile$Maj
                } else {
                    average_CN_profile$Min <- average_CN_profile$Min + clonal_populations[clone] * clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[clone]]$Min
                    average_CN_profile$Maj <- average_CN_profile$Maj + clonal_populations[clone] * clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[clone]]$Maj
                }
            }
            average_CN_profile$Min <- round(average_CN_profile$Min / sum(clonal_populations))
            average_CN_profile$Maj <- round(average_CN_profile$Maj / sum(clonal_populations))
            average_CN_profile$state <- average_CN_profile$Min + average_CN_profile$Maj
            average_CN_profile$copy <- average_CN_profile$state
            clonal_CN_profiles_all_sims_new[["variable=average_CN_profile"]][[i]] <- average_CN_profile
        }
        clonal_CN_profiles_all_sims <- clonal_CN_profiles_all_sims_new
    }
    #-----------------------------------Return the clonal CN information
    return(clonal_CN_profiles_all_sims)
}

#---------------------------------Get statistics for each simulation
get_each_statistics <- function(simulations, simulations_clonal_CN, list_targets) {
    library(vegan)
    library(matrixStats)
    library(transport)
    library(ape)
    library(phyloTop)
    library(treebalance)
    simulations_statistics <- list()
    #----------------------Initialize matrices for simulation statistics
    for (stat in list_targets) {
        stat_details <- strsplit(stat, ";")[[1]]
        stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
        stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
        if (stat_target == "genome") {
            simulations_statistics[[stat_ID]] <- matrix(0, nrow = length(simulations), ncol = 1)
        } else if (stat_target == "chromosome") {
            stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
            simulations_statistics[[stat_ID]] <- matrix(0, nrow = length(simulations), ncol = length(stat_chromosome_ID))
        }
    }
    #--------------------------------------Compute simulation statistics
    for (i in 1:length(simulations)) {
        simulation <- simulations[[i]]
        #   Find clonal ancestors for single-cell or bulk average-CN-based event counts
        if (any(grepl("variable=event_count", list_targets) & grepl("data=sc", list_targets)) |
            any(grepl("variable=event_count", list_targets) & grepl("data=bulk", list_targets) & grepl("representative_CN=average_CN", list_targets))) {
            evolution_origin <- simulation$clonal_evolution$evolution_origin
            sample_genotype_unique <- simulation$sample$sample_genotype_unique
            evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
            subclonal_ancestry <- vector("list", length(sample_genotype_unique))
            for (j in 1:length(sample_genotype_unique)) {
                ancestry <- sample_genotype_unique[j]
                while (ancestry[1] != 0) ancestry <- c(evolution_origin[ancestry[1]], ancestry)
                subclonal_ancestry[[j]] <- ancestry
            }
            clonal_ancestry <- find_clonal_ancestry(subclonal_ancestry)
        }
        #   Find clonal ancestors for bulk maximal-CN-based event counts
        if (any(grepl("variable=event_count", list_targets) & grepl("data=bulk", list_targets) & grepl("representative_CN=maximal_CN", list_targets))) {
            #   ....
            #   ....
            #   ....
        }
        #   Get requested statistics for this simulation
        for (stat in list_targets) {
            stat_details <- strsplit(stat, ";")[[1]]
            stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
            stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
            stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
            if (stat_target == "chromosome") {
                stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
                if (stat_variable == "shannon") {
                    #   Get Shannon index
                } else if (stat_variable == "event_count") {
                    #   Get source of data
                    stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
                    #   Get count of clonal/subclonal events of a given type
                    clonal_type <- strsplit(stat_details[grep("type=", stat_details)], "=")[[1]][2]
                    event_type <- strsplit(stat_details[grep("event=", stat_details)], "=")[[1]][2]
                    if (stat_data == "sc") {
                        #   Get event count from single-cell samples
                        if (clonal_type == "clonal") {
                            for (j in 1:length(stat_chromosome_ID)) {
                                simulations_statistics[[stat_ID]][i, j] <- find_event_count(clonal_ancestry, event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                            }
                        } else if (clonal_type == "subclonal") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(stat_chromosome_ID)) {
                                for (k in 1:length(sample_genotype_unique)) {
                                    sample_genotype_event_counts[k] <- find_event_count(subclonal_ancestry[[k]], event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                                }
                                simulations_statistics[[stat_ID]][i, j] <-
                                    sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) /
                                    length(simulation$sample$sample_clone_ID) -
                                    find_event_count(clonal_ancestry, event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                            }
                        } else {
                            stop(paste0("Error: Unknown clonal type: ", stat))
                        }
                    } else if (stat_data == "bulk") {
                        representative_CN_type <- strsplit(stat_details[grep("representative_CN=", stat_details)], "=")[[1]][2]
                        if (representative_CN_type == "average_CN") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(stat_chromosome_ID)) {
                                for (k in 1:length(sample_genotype_unique)) {
                                    sample_genotype_event_counts[k] <- find_event_count(subclonal_ancestry[[k]], event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                                }
                                simulations_statistics[[stat_ID]][i, j] <-
                                    sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) / length(simulation$sample$sample_clone_ID)
                            }
                        } else if (representative_CN_type == "maximal_CN") {
                            simpleError("Error: event count from maximal CN in bulk data not implemented yet")
                        }
                    }
                }
            } else if (stat_target == "genome") {
                if (stat_variable == "shannon") {
                    #   Get Shannon index
                    frequencies <- table(simulation$sample$sample_clone_ID)
                    simulations_statistics[[stat_ID]][i] <- diversity(frequencies, index = "shannon")
                } else if (stat_variable == "event_count") {
                    #   Get source of data
                    stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
                    #   Get count of clonal/subclonal events of a given type
                    clonal_type <- strsplit(stat_details[grep("type=", stat_details)], "=")[[1]][2]
                    event_type <- strsplit(stat_details[grep("event=", stat_details)], "=")[[1]][2]
                    if (stat_data == "sc") {
                        #   Get event count from single-cell samples
                        if (clonal_type == "clonal") {
                            simulations_statistics[[stat_ID]][i] <- find_event_count(clonal_ancestry, event_type, evolution_genotype_changes)
                        } else if (clonal_type == "subclonal") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(sample_genotype_unique)) {
                                sample_genotype_event_counts[j] <- find_event_count(subclonal_ancestry[[j]], event_type, evolution_genotype_changes)
                            }
                            simulations_statistics[[stat_ID]][i] <-
                                sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) /
                                length(simulation$sample$sample_clone_ID) -
                                find_event_count(clonal_ancestry, event_type, evolution_genotype_changes)
                        } else {
                            stop(paste0("Error: Unknown clonal type: ", stat))
                        }
                    } else if (stat_data == "bulk") {
                        #   Get event count from bulk samples
                        representative_CN_type <- strsplit(stat_details[grep("representative_CN=", stat_details)], "=")[[1]][2]
                        if (representative_CN_type == "average_CN") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(sample_genotype_unique)) {
                                sample_genotype_event_counts[j] <- find_event_count(subclonal_ancestry[[j]], event_type, evolution_genotype_changes)
                            }
                            simulations_statistics[[stat_ID]][i] <-
                                sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) / length(simulation$sample$sample_clone_ID)
                        } else if (representative_CN_type == "maximal_CN") {
                            simpleError("Error: event count from maximal CN in bulk data not implemented yet")
                        }
                    }
                } else if (stat_variable == "cherries") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get number of cherries
                    simulations_statistics[[stat_ID]][i] <- cherries(tree, normalise = TRUE)
                } else if (stat_variable == "pitchforks") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get number of pitchforks
                    simulations_statistics[[stat_ID]][i] <- pitchforks(tree, normalise = TRUE)
                } else if (stat_variable == "colless") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get colless index
                    simulations_statistics[[stat_ID]][i] <- colless.phylo(tree, normalise = TRUE)
                } else if (stat_variable == "sackin") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get sackin index
                    simulations_statistics[[stat_ID]][i] <- sackin.phylo(tree, normalise = TRUE)
                } else if (stat_variable == "IL_number") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get IL_number
                    simulations_statistics[[stat_ID]][i] <- ILnumber(tree, normalise = TRUE)
                } else if (stat_variable == "avgLadder") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get avgLadder
                    simulations_statistics[[stat_ID]][i] <- avgLadder(tree, normalise = TRUE)
                } else if (stat_variable == "maxDepth") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get maxDepth
                    simulations_statistics[[stat_ID]][i] <- maxDepth(tree)
                } else if (stat_variable == "stairs") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get stairness
                    simulations_statistics[[stat_ID]][i] <- stairs(tree)[1]
                } else if (stat_variable == "B2") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get B2Index
                    simulations_statistics[[stat_ID]][i] <- B2I(tree, logbase = 2)
                } else if (stat_variable == "clonal_CN") {
                    #   Get clonal CN profiles and their populations
                    simulations_statistics[["variable=clonal_CN_profiles"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_profiles"]][[i]]
                    simulations_statistics[["variable=clonal_CN_populations"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_populations"]][[i]]
                } else if (stat_variable == "maximal_CN") {
                    #   Get maximal CN profile
                    simulations_statistics[["variable=maximal_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=maximal_CN_profile"]][[i]]
                } else if (stat_variable == "average_CN") {
                    #   Get average CN profile
                    simulations_statistics[["variable=average_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=average_CN_profile"]][[i]]
                } else {
                    stop(paste0("Error: Unknown statistic: ", stat))
                }
            }
        }
    }
    return(simulations_statistics)
}

#----------------------------------------------------Local functions
#---Function to find shared ancestral clones between subclones
find_clonal_ancestry <- function(list_subclonal_ancestry) {
    if (length(list_subclonal_ancestry) == 0) {
        clonal_ancestry <- c()
    } else if (length(list_subclonal_ancestry) == 1) {
        clonal_ancestry <- list_subclonal_ancestry[[1]]
    } else {
        clonal_ancestry <- list_subclonal_ancestry[[1]]
        for (i in 2:length(list_subclonal_ancestry)) {
            clonal_ancestry <- intersect(clonal_ancestry, list_subclonal_ancestry[[i]])
        }
    }
    return(clonal_ancestry)
}
#---Function to find the number of events
#---given a list of clones and event type
find_event_count <- function(clone_ancestry, event_type, evolution_genotype_changes, by_chromosome = NULL) {
    event_count <- 0
    for (clone in clone_ancestry) {
        if (clone <= 0) next
        if (length(evolution_genotype_changes[[clone]]) == 0) next
        for (j in 1:length(evolution_genotype_changes[[clone]])) {
            if (evolution_genotype_changes[[clone]][[j]][1] == event_type) {
                if (is.null(by_chromosome)) {
                    event_count <- event_count + 1
                } else {
                    if (evolution_genotype_changes[[clone]][[j]][2] == by_chromosome) {
                        event_count <- event_count + 1
                    }
                }
            }
        }
    }
    return(event_count)
}
#---Function to find distance between two CN profiles
cn_distance <- function(from_CN_profile,
                        to_CN_profile,
                        metric) {
    if (metric == "euclidean") {
        distance <- sqrt(sum((from_CN_profile - to_CN_profile)^2))
    } else if (metric == "wcnd") {
        distance <- 0
    } else if (metric == "correlation") {
        library(energy)
        distance <- dcor(from_CN_profile, to_CN_profile)
    } else if (metric == "mahalanobis") {
        library(MASS)
        distance <- sqrt(t(from_CN_profile - to_CN_profile) %*% ginv(cov(from_CN_profile %*% t(from_CN_profile))) %*% (from_CN_profile - to_CN_profile))
    }
    return(distance)
}

#---Function to find distance between two samples
#---based on clonal population and CN profiles
sample_distance <- function(from_clonal_CN_populations,
                            from_clonal_CN_profiles,
                            to_clonal_CN_populations,
                            to_clonal_CN_profiles,
                            metric) {
    library(transport)
    #   Compute distance between any pair of CN profiles in the samples
    n_sample_from <- length(from_clonal_CN_populations)
    n_sample_to <- length(to_clonal_CN_populations)
    cost_matrix <- matrix(nrow = n_sample_from, ncol = n_sample_to)
    for (i in 1:n_sample_from) {
        for (j in 1:n_sample_to) {
            from_CN_profile <- from_clonal_CN_profiles[[i]]$copy
            to_CN_profile <- to_clonal_CN_profiles[[j]]$copy
            cost_matrix[i, j] <- cn_distance(from_CN_profile, to_CN_profile, metric)
        }
    }
    #   Normalize clonal populations
    from_clonal_CN_populations <- from_clonal_CN_populations / sum(from_clonal_CN_populations)
    to_clonal_CN_populations <- to_clonal_CN_populations / sum(to_clonal_CN_populations)
    #   Compute the Wasserstein distance between the two samples
    wasserstein_dist <- wasserstein(
        from_clonal_CN_populations,
        to_clonal_CN_populations,
        costm = cost_matrix, p = 1
    )
    return(wasserstein_dist)
}
#---Function to find distance between two sample cohorts
#---based on clonal population and CN profiles
cohort_distance <- function(cohort_from, cohort_to, metric, bulk_CN_input = "", bulk = FALSE) {
    library(transport)
    #-----------------------------Define cost matrix between two cohorts
    if (bulk) {
        if (bulk_CN_input == "average") {
            n_cohort_from <- length(cohort_from[["variable=average_CN_profile"]])
            n_cohort_to <- length(cohort_to[["variable=average_CN_profile"]])
        } else if (bulk_CN_input == "maximal") {
            n_cohort_from <- length(cohort_from[["variable=maximal_CN_profile"]])
            n_cohort_to <- length(cohort_to[["variable=maximal_CN_profile"]])
        }
        sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
        for (i in 1:n_cohort_from) {
            for (j in 1:n_cohort_to) {
                if (bulk_CN_input == "average") {
                    from_CN_profile <- cohort_from[["variable=average_CN_profile"]][[i]]$copy
                    to_CN_profile <- cohort_to[["variable=average_CN_profile"]][[j]]$copy
                } else if (bulk_CN_input == "maximal") {
                    from_CN_profile <- cohort_from[["variable=maximal_CN_profile"]][[i]]$copy
                    to_CN_profile <- cohort_to[["variable=maximal_CN_profile"]][[j]]$copy
                }
                sample_distance_mtx[i, j] <- cn_distance(from_CN_profile, to_CN_profile, metric)
            }
        }
    } else {
        #   Compute distance between any pair of samples in the cohorts
        n_cohort_from <- length(cohort_from[["variable=clonal_CN_populations"]])
        n_cohort_to <- length(cohort_to[["variable=clonal_CN_populations"]])
        sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
        for (i in 1:n_cohort_from) {
            for (j in 1:n_cohort_to) {
                from_clonal_CN_populations <- cohort_from[["variable=clonal_CN_populations"]][[i]]
                from_clonal_CN_profiles <- cohort_from[["variable=clonal_CN_profiles"]][[i]]
                to_clonal_CN_populations <- cohort_to[["variable=clonal_CN_populations"]][[j]]
                to_clonal_CN_profiles <- cohort_to[["variable=clonal_CN_profiles"]][[j]]
                sample_distance_mtx[i, j] <- sample_distance(
                    from_clonal_CN_populations = from_clonal_CN_populations,
                    from_clonal_CN_profiles = from_clonal_CN_profiles,
                    to_clonal_CN_populations = to_clonal_CN_populations,
                    to_clonal_CN_profiles = to_clonal_CN_profiles,
                    metric = metric
                )
            }
        }
    }
    #---------------------------------------Normalize clonal populations
    dist_from <- rep(1 / n_cohort_from, n_cohort_from)
    dist_to <- rep(1 / n_cohort_to, n_cohort_to)
    #-----------Compute the Wasserstein distance between the two cohorts
    wasserstein_dist <- wasserstein(
        dist_from,
        dist_to,
        costm = sample_distance_mtx, p = 1
    )
    return(wasserstein_dist)
}

#' @export
get_statistics <- function(list_targets,
                           simulations_sc = NULL,
                           simulations_bulk = NULL,
                           simulations_statistics_sc = NULL,
                           simulations_statistics_bulk = NULL,
                           cn_data_sc = NULL,
                           cn_data_bulk = NULL,
                           arm_level = FALSE,
                           cn_table = NULL,
                           save_sample_statistics = FALSE) {
    #-------------Get clonal CN profiles for all single-cell simulations
    if (is.null(simulations_statistics_sc)) {
        if (any(grepl("variable=clonal_CN", list_targets))) {
            simulations_clonal_CN_sc <- get_each_clonal_CN_profiles(simulations_sc, arm_level, cn_table)
        }
        list_targets_sc <- list_targets[grepl("data=sc", list_targets)]
        simulations_statistics_sc <- get_each_statistics(simulations_sc, simulations_clonal_CN_sc, list_targets_sc)
    }
    #------------Get representative CN profiles for all bulk simulations
    if (is.null(simulations_statistics_bulk)) {
        if ((any(grepl("variable=average_CN", list_targets))) | (any(grepl("variable=maximal_CN", list_targets)))) {
            simulations_clonal_CN_bulk <- get_each_clonal_CN_profiles(simulations_bulk, arm_level, cn_table, bulk = TRUE)
        }
        list_targets_bulk <- list_targets[grepl("data=bulk", list_targets)]
        simulations_statistics_bulk <- get_each_statistics(simulations_bulk, simulations_clonal_CN_bulk, list_targets_bulk)
    }
    #-----------------------------------------Get statistics for fitting
    statistics <- vector("list", length(list_targets))
    for (i in 1:length(list_targets)) {
        stat <- list_targets[i]
        stat_details <- strsplit(stat, ";")[[1]]
        stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
        stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
        stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
        stat_type <- strsplit(stat_details[grepl("statistic=", stat_details)], "=")[[1]][2]
        if (stat_type == "mean") {
            if (stat_data == "sc") {
                statistics[[i]] <- apply(simulations_statistics_sc[[stat_ID]], 2, mean)
            } else if (stat_data == "bulk") {
                statistics[[i]] <- apply(simulations_statistics_bulk[[stat_ID]], 2, mean)
            }
        } else if (stat_type == "var") {
            if (stat_data == "sc") {
                statistics[[i]] <- apply(simulations_statistics_sc[[stat_ID]], 2, var)
            } else if (stat_data == "bulk") {
                statistics[[i]] <- apply(simulations_statistics_bulk[[stat_ID]], 2, var)
            }
        } else if (stat_type == "dist") {
            if (stat_variable == "clonal_CN" & stat_data == "sc") {
                stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                if (is.null(cn_data_sc)) cn_data_sc <- simulations_statistics_sc
                statistics[[i]] <- cohort_distance(
                    cohort_from = simulations_statistics_sc,
                    cohort_to = cn_data_sc,
                    metric = stat_metric
                )
            } else if (stat_variable == "average_CN" & stat_data == "bulk") {
                stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                statistics[[i]] <- cohort_distance(
                    cohort_from = simulations_statistics_bulk,
                    cohort_to = cn_data_bulk,
                    metric = stat_metric,
                    bulk_CN_input = "average",
                    bulk = TRUE
                )
            } else if (stat_variable == "maximal_CN" & stat_data == "bulk") {
                stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                statistics[[i]] <- cohort_distance(
                    cohort_from = simulations_statistics_bulk,
                    cohort_to = cn_data_bulk,
                    metric = stat_metric,
                    bulk_CN_input = "maximal",
                    bulk = TRUE
                )
            }
        } else {
            stop(paste0("Error: Unknown statistic type: ", stat))
        }
    }
    #-----------------------------------------------------Prepare output
    output <- list()
    output$statistics <- statistics
    if (save_sample_statistics) {
        output$simulations_statistics_sc <- simulations_statistics_sc
        output$simulations_statistics_bulk <- simulations_statistics_bulk
    }
    return(output)
}

#---Function to assign parameters to proper positions
assign_paras <- function(model_variables, parameter_IDs, parameters) {
    for (i in 1:length(parameter_IDs)) {
        parameter_ID_input <- parameter_IDs[i]
        parameter_value_input <- parameters[i]
        #   Prepare values for operation on parameter
        if (grepl(":", parameter_ID_input)) {
            parameter_ID <- sub(".*:", "", parameter_ID_input)
            parameter_operator <- sub(":.*", "", parameter_ID_input)
            parameter_operator <- paste0(parameter_operator, "(parameter_value_input)")
        } else {
            parameter_ID <- parameter_ID_input
            parameter_operator <- "parameter_value_input"
        }
        parameter_value <- eval(parse(text = parameter_operator))
        #   Input parameter
        if (parameter_ID %in% model_variables$general_variables$Variable) {
            model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameter_value
        } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
            model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)] <- parameter_value
        }
    }
    return(model_variables)
}
#---Objective function for ABC fitting
func_ABC <- function(parameters,
                     parameter_IDs,
                     model_variables,
                     list_targets,
                     cn_data_sc = NULL,
                     cn_data_bulk = NULL,
                     arm_level = FALSE,
                     save_sample_statistics = FALSE) {
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
    #   Make single-cell simulations
    SIMS_sc <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_simulations_sc,
        stage_final = 3,
        save_simulation = FALSE, report_progress = TRUE,
        lite_memory = TRUE,
        output_variables = c(
            "evolution_origin",
            "evolution_genotype_changes",
            "sample_clone_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile",
            "phylogeny_clustering_truth"
        )
    )
    #   Make bulk simulations
    SIMS_bulk <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_simulations_bulk,
        stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE,
        lite_memory = TRUE,
        output_variables = c(
            "evolution_origin",
            "evolution_genotype_changes",
            "sample_clone_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile"
        )
    )
    #   Get statistics from simulations
    stat <- get_statistics(
        simulations_sc = SIMS_sc,
        simulations_bulk = SIMS_bulk,
        cn_data_sc = cn_data_sc,
        cn_data_bulk = cn_data_bulk,
        list_targets = list_targets,
        arm_level = arm_level,
        cn_table = cn_table,
        save_sample_statistics = save_sample_statistics
    )
    return(stat)
}

#' @export
library_sc_CN <- function(model_name,
                          model_variables,
                          list_parameters,
                          list_targets_library,
                          ABC_simcount = 10000,
                          arm_level = FALSE,
                          cn_table = NULL,
                          ##############################################
                          cn_data_sc = NULL,
                          cn_data_bulk = NULL,
                          n_simulations_sc = 10,
                          n_simulations_bulk = 10,
                          ##############################################
                          n_cores = NULL,
                          library_name,
                          save_sample_statistics = FALSE) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    #---------------------------------List of parameter IDs to be fitted
    parameter_IDs <- list_parameters$Variable
    # =============================================CREATE REFERENCE DATA
    #---------------------------------------Simulate table of parameters
    sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    for (col in 1:ncol(sim_param)) {
        sim_param[, col] <- runif(
            ABC_simcount,
            min = as.numeric(list_parameters$Lower_bound[col]),
            max = as.numeric(list_parameters$Upper_bound[col])
        )
    }
    #-----------------------------------------------Make reference table
    start_time <- Sys.time()
    #   Configure parallel pool
    cl <- makePSOCKcluster(n_cores)
    cat("Creating reference table for ABC...\n")
    sim_param <<- sim_param
    parameter_IDs <<- parameter_IDs
    model_variables <<- model_variables
    func_ABC <<- func_ABC
    assign_paras <<- assign_paras
    list_targets_library <<- list_targets_library
    ####
    ####
    ####
    cn_table <<- cn_table
    save_sample_statistics <<- save_sample_statistics
    arm_level <<- arm_level
    cn_data_sc <<- cn_data_sc
    cn_data_bulk <<- cn_data_bulk
    n_simulations_sc <<- n_simulations_sc
    n_simulations_bulk <<- n_simulations_bulk
    ####
    ####
    ####
    clusterExport(cl, varlist = c(
        "n_simulations_sc", "n_simulations_bulk", "save_sample_statistics", "list_targets_library", "sim_param", "parameter_IDs", "model_variables", "cn_data_sc", "cn_data_sc", "cn_table", "arm_level",
        "func_ABC", "assign_paras", "get_statistics", "get_each_clonal_CN_profiles", "save_sample_statistics", "get_cn_profile", "get_arm_CN_profiles",
        "find_clonal_ancestry", "find_event_count", "cn_distance", "sample_distance", "cohort_distance", "get_each_statistics",
        "vec_CN_block_no", "vec_centromeres",
        "BUILD_driver_library", "simulator_full_program", "one_simulation",
        "SIMULATOR_VARIABLES_for_simulation",
        "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
        "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
        "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
        "SIMULATOR_FULL_PHASE_2_main",
        "SIMULATOR_FULL_PHASE_3_main",
        "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model"
    ))
    e <- new.env()
    e$libs <- .libPaths()
    clusterExport(cl, "libs", envir = e)
    clusterEvalQ(cl, .libPaths(libs))
    #   Create simulated results in parallel
    pbo <- pboptions(type = "txt")
    sim_results_list <- pblapply(cl = cl, X = 1:ABC_simcount, FUN = function(iteration) {
        parameters <- sim_param[iteration, ]
        stat <- func_ABC(
            parameters = parameters, parameter_IDs = parameter_IDs, model_variables = model_variables, list_targets = list_targets_library, arm_level = arm_level,
            cn_data_sc = cn_data_sc, cn_data_bulk = cn_data_bulk,
            save_sample_statistics = save_sample_statistics
        )
        return(stat)
    })
    stopCluster(cl)
    #   Group simulated statistics into one table
    sim_stat <- matrix(0, nrow = ABC_simcount, ncol = length(sim_results_list[[1]]$statistics))
    for (row in 1:ABC_simcount) {
        sim_stat[row, ] <- sim_results_list[[row]]$statistics
    }
    #   Find simulated CN profiles if requested
    if (save_sample_statistics) {
        sim_sample_stat <- list()
        for (row in 1:ABC_simcount) {
            sim_sample_stat[[row]] <- list()
            sim_sample_stat[[row]]$sc <- sim_results_list[[row]]$simulations_statistics_sc
            sim_sample_stat[[row]]$bulk <- sim_results_list[[row]]$simulations_statistics_bulk
        }
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    #---------------------------Save the parameters and their statistics
    print(sim_param)
    print(sim_stat)



    ABC_input <- list()
    ABC_input$model_variables <- model_variables
    ABC_input$parameter_IDs <- parameter_IDs
    ABC_input$sim_param <- sim_param
    ABC_input$sim_stat <- sim_stat
    if (save_sample_statistics) ABC_input$sim_sample_stat <- sim_sample_stat
    filename <- paste0(library_name, "_ABC_input.rda")
    save(ABC_input, file = filename)
}

permutate_chromosome_arms <- function(current_sim_param, current_sim_sample_stat, list_parameters, permutate_chromosome) {
    #--------------------------------Permutate simulation parameters
    locs_old_chrom_arms <- which(list_parameters$Chromosome == permutate_chromosome)
    locs_new_chrom_arms <- locs_old_chrom_arms[length(locs_old_chrom_arms):1]
    new_sim_param <- current_sim_param
    new_sim_param[locs_old_chrom_arms] <- new_sim_param[locs_new_chrom_arms]
    #--------------------------------Permutate simulation statistics
    new_sim_sample_stat <- current_sim_sample_stat
    for (sim in 1:length(current_sim_sample_stat[["variable=clonal_CN_profiles"]])) {
        for (clone in 1:length(current_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]])) {
            current_clonal_CN_profiles <- current_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]][[clone]]
            locs_old_chrom_arms <- which(current_clonal_CN_profiles$chr == permutate_chromosome)
            locs_new_chrom_arms <- locs_old_chrom_arms[length(locs_old_chrom_arms):1]
            new_clonal_CN_profiles <- current_clonal_CN_profiles
            new_clonal_CN_profiles$copy[locs_old_chrom_arms] <- current_clonal_CN_profiles$copy[locs_new_chrom_arms]
            new_clonal_CN_profiles$state[locs_old_chrom_arms] <- current_clonal_CN_profiles$state[locs_new_chrom_arms]
            new_clonal_CN_profiles$Min[locs_old_chrom_arms] <- current_clonal_CN_profiles$Min[locs_new_chrom_arms]
            new_clonal_CN_profiles$Maj[locs_old_chrom_arms] <- current_clonal_CN_profiles$Maj[locs_new_chrom_arms]
            new_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]][[clone]] <- new_clonal_CN_profiles
        }
    }
    #-------------------------------------------------Prepare output
    output <- list()
    output$new_sim_param <- new_sim_param
    output$new_sim_sample_stat <- new_sim_sample_stat
    return(output)
}
permutate_chromosomes <- function(current_sim_param, current_sim_sample_stat, list_parameters, current_chromosomes, new_chromosomes) {
    #------------------------------------Permutate simulation parameters
    new_sim_param <- current_sim_param
    for (i in 1:length(current_chromosomes)) {
        locs_old_chrom <- which(list_parameters$Chromosome == current_chromosomes[i])
        locs_new_chrom <- which(list_parameters$Chromosome == new_chromosomes[i])
        new_sim_param[locs_old_chrom] <- current_sim_param[locs_new_chrom]
    }
    #------------------------------------Permutate simulation statistics
    new_sim_sample_stat <- current_sim_sample_stat
    #   Permutate single-cell CN profiles
    if ("sc" %in% names(current_sim_sample_stat)) {
        current_sim_sample_stat_sc <- current_sim_sample_stat$sc
        for (sim in 1:length(current_sim_sample_stat_sc[["variable=clonal_CN_profiles"]])) {
            for (clone in 1:length(current_sim_sample_stat_sc[["variable=clonal_CN_profiles"]][[sim]])) {
                current_clonal_CN_profiles <- current_sim_sample_stat_sc[["variable=clonal_CN_profiles"]][[sim]][[clone]]
                new_clonal_CN_profiles <- current_clonal_CN_profiles
                for (i in 1:length(current_chromosomes)) {
                    locs_old_chrom <- which(current_clonal_CN_profiles$chr == current_chromosomes[i])
                    locs_new_chrom <- which(current_clonal_CN_profiles$chr == new_chromosomes[i])
                    new_clonal_CN_profiles$copy[locs_old_chrom] <- current_clonal_CN_profiles$copy[locs_new_chrom]
                    new_clonal_CN_profiles$state[locs_old_chrom] <- current_clonal_CN_profiles$state[locs_new_chrom]
                    new_clonal_CN_profiles$Min[locs_old_chrom] <- current_clonal_CN_profiles$Min[locs_new_chrom]
                    new_clonal_CN_profiles$Maj[locs_old_chrom] <- current_clonal_CN_profiles$Maj[locs_new_chrom]
                }
                new_sim_sample_stat$sc[["variable=clonal_CN_profiles"]][[sim]][[clone]] <- new_clonal_CN_profiles
            }
        }
    }
    #   Permutate bulk CN profiles
    if ("bulk" %in% names(current_sim_sample_stat)) {
        current_sim_sample_stat_bulk <- current_sim_sample_stat$bulk
        vec_variable_names <- names(current_sim_sample_stat_bulk)
        if ("variable=average_CN_profile" %in% vec_variable_names) {
            for (sim in 1:length(current_sim_sample_stat_bulk[["variable=average_CN_profile"]])) {
                current_clonal_CN_profiles <- current_sim_sample_stat_bulk[["variable=average_CN_profile"]][[sim]]
                new_clonal_CN_profiles <- current_clonal_CN_profiles
                for (i in 1:length(current_chromosomes)) {
                    locs_old_chrom <- which(current_clonal_CN_profiles$chr == current_chromosomes[i])
                    locs_new_chrom <- which(current_clonal_CN_profiles$chr == new_chromosomes[i])
                    new_clonal_CN_profiles$copy[locs_old_chrom] <- current_clonal_CN_profiles$copy[locs_new_chrom]
                    new_clonal_CN_profiles$state[locs_old_chrom] <- current_clonal_CN_profiles$state[locs_new_chrom]
                    new_clonal_CN_profiles$Min[locs_old_chrom] <- current_clonal_CN_profiles$Min[locs_new_chrom]
                    new_clonal_CN_profiles$Maj[locs_old_chrom] <- current_clonal_CN_profiles$Maj[locs_new_chrom]
                }
                new_sim_sample_stat$bulk[["variable=average_CN_profile"]][[sim]] <- new_clonal_CN_profiles
            }
        }
        if ("variable=maximal_CN_profile" %in% vec_variable_names) {
            for (sim in 1:length(current_sim_sample_stat_bulk[["variable=maximal_CN_profile"]])) {
                current_clonal_CN_profiles <- current_sim_sample_stat_bulk[["variable=maximal_CN_profile"]][[sim]]
                new_clonal_CN_profiles <- current_clonal_CN_profiles
                for (i in 1:length(current_chromosomes)) {
                    locs_old_chrom <- which(current_clonal_CN_profiles$chr == current_chromosomes[i])
                    locs_new_chrom <- which(current_clonal_CN_profiles$chr == new_chromosomes[i])
                    new_clonal_CN_profiles$copy[locs_old_chrom] <- current_clonal_CN_profiles$copy[locs_new_chrom]
                    new_clonal_CN_profiles$state[locs_old_chrom] <- current_clonal_CN_profiles$state[locs_new_chrom]
                    new_clonal_CN_profiles$Min[locs_old_chrom] <- current_clonal_CN_profiles$Min[locs_new_chrom]
                    new_clonal_CN_profiles$Maj[locs_old_chrom] <- current_clonal_CN_profiles$Maj[locs_new_chrom]
                }
                new_sim_sample_stat$bulk[["variable=maximal_CN_profile"]][[sim]] <- new_clonal_CN_profiles
            }
        }
    }
    #-----------------------------------------------------Prepare output
    output <- list()
    output$new_sim_param <- new_sim_param
    output$new_sim_sample_stat <- new_sim_sample_stat
    return(output)
}
########################################################
func_ABC_by_permutation <- function(current_sim_param,
                                    current_sim_sample_stat,
                                    current_chromosomes,
                                    new_chromosomes,
                                    list_parameters,
                                    cn_data_sc = NULL,
                                    cn_data_bulk = NULL,
                                    cn_table = NULL,
                                    arm_level = FALSE) {
    df_permutate_chromosomes <- permutate_chromosomes(
        current_sim_param = current_sim_param,
        current_sim_sample_stat = current_sim_sample_stat,
        list_parameters = list_parameters,
        current_chromosomes = list_chromosomes,
        new_chromosomes = list_chromosomes_new
    )
    sim_stat_new <- get_statistics(
        list_targets = colnames(list_targets)[-1],
        simulations_statistics_sc = df_permutate_chromosomes$new_sim_sample_stat$sc,
        simulations_statistics_bulk = df_permutate_chromosomes$new_sim_sample_stat$bulk,
        cn_data_sc = cn_data_sc,
        cn_data_bulk = cn_data_bulk,
        arm_level = arm_level,
        cn_table = cn_table
    )$statistics
    output <- list()
    output$sim_param <- df_permutate_chromosomes$new_sim_param
    output$sim_stat <- sim_stat_new
    return(output)
}
########################################################


#' @export
fitting_sc_CN <- function(library_name,
                          model_name,
                          copynumber_DATA,
                          parameters_truth = NULL,
                          list_parameters,
                          list_targets,
                          shuffle_num,
                          n_cores = NULL,
                          cn_data_sc = NULL,
                          cn_data_bulk = NULL,
                          cn_table = NULL,
                          arm_level = FALSE,
                          shuffle_chromosome_arms = FALSE,
                          shuffle_chromosomes_by_permutation = FALSE,
                          shuffle_chromosomes_by_moving = FALSE) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    library(data.table)
    library(matrixStats)
    library(combinat)
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    #-----------------------------------------Input simulated CN library
    filename <- paste0(library_name, "_ABC_input.rda")
    load(filename)
    n_samples <- ABC_input$n_samples
    model_variables <- ABC_input$model_variables
    parameter_IDs <- ABC_input$parameter_IDs
    sim_param <- ABC_input$sim_param
    sim_stat <- ABC_input$sim_stat
    sim_sample_stat <- ABC_input$sim_sample_stat
    #--------------------------------Find statistics for each CN heatmap
    DATA_target <- copynumber_DATA$statistics
    # ========================INCREASE SIMULATED LIBRARY VIA PERMUTATION
    if ((shuffle_chromosome_arms | shuffle_chromosomes_by_permutation | shuffle_chromosomes_by_moving) &
        ((any(grepl("sc", names(sim_sample_stat[[1]])))) |
            (any(grepl("bulk", names(sim_sample_stat[[1]]))))
        )) {
        #---Find chromosome and arm for each parameter
        list_parameters$Chromosome <- NA
        list_parameters$Arm <- NA
        for (i in 1:nrow(list_parameters)) {
            if (list_parameters$Type[i] == "Arm_selection_rate") {
                list_parameters$Chromosome[i] <- substr(list_parameters$Variable[i], 1, nchar(list_parameters$Variable[i]) - 1)
                list_parameters$Arm[i] <- substr(list_parameters$Variable[i], nchar(list_parameters$Variable[i]), nchar(list_parameters$Variable[i]))
            }
        }
        list_chromosomes <- unique(list_parameters$Chromosome[!is.na(list_parameters$Chromosome)])
        #---Boost simulated data by permutating chromosome arms
        if (shuffle_chromosome_arms) {
            # sim_param_new <- sim_param
            # sim_stat_new <- sim_stat
            # sim_sample_stat_new <- sim_sample_stat
            # for (chromosome in list_chromosomes) {
            #     if (length(which(list_parameters$Chromosome == chromosome)) <= 1) next
            #     if (length(which(list_parameters$Chromosome == chromosome)) > 2) simpleError("More than two arms for a chromosome")
            #     sim_param_next <- sim_param_new
            #     sim_stat_next <- sim_stat_new
            #     for (sim in 1:nrow(sim_param_new)) {
            #         current_sim_param <- sim_param_new[sim, ]
            #         current_sim_sample_stat <- sim_sample_stat_new[[sim]]
            #         df_permutate_chromosome_arms <- permutate_chromosome_arms(
            #             current_sim_param = current_sim_param,
            #             current_sim_sample_stat = current_sim_sample_stat,
            #             list_parameters = list_parameters,
            #             permutate_chromosome = chromosome
            #         )
            #         new_sim_param <- df_permutate_chromosome_arms$new_sim_param
            #         new_sim_sample_stat <- df_permutate_chromosome_arms$new_sim_sample_stat
            #         sim_param_next <- rbind(sim_param_next, new_sim_param)
            #         sim_stat_next <- rbind(sim_stat_next, get_statistics(
            #             simulations_statistics = new_sim_sample_stat,
            #             list_targets = colnames(list_targets)[-1],
            #             cn_data = cn_data,
            #             arm_level = arm_level,
            #             cn_table = cn_table
            #         )$statistics)
            #     }
            #     sim_param_new <- sim_param_next
            #     sim_stat_new <- sim_stat_next
            # }
            # sim_param <- sim_param_new
            # sim_stat <- sim_stat_new
        }
        #---Boost simulated data by moving chromosomes
        if (shuffle_chromosomes_by_moving) {
            #   Initialize new library
            sim_param_new <- sim_param
            sim_stat_new <- sim_stat
            #   Find maximum number of possible moves
            max_shuffle_count <- length(list_chromosomes)
            shuffle_num <- min(shuffle_num, max_shuffle_count)

            # pb <- txtProgressBar(
            #     min = 2, max = shuffle_num,
            #     style = 3, width = 50, char = "="
            # )
            for (i in 2:shuffle_num) {
                # setTxtProgressBar(pb, i)
                list_chromosomes_new <- c(list_chromosomes[i:length(list_chromosomes)], list_chromosomes[1:(i - 1)])
                #   Permutate chromosomes
                sim_param_next <- sim_param
                sim_stat_next <- sim_stat
                # ==================================================================================================
                start_time <- Sys.time()
                cl <- makePSOCKcluster(n_cores)
                cat("\nShuffling by chromosomes...\n")
                sim_param <<- sim_param
                sim_sample_stat <<- sim_sample_stat
                func_ABC_by_permutation <<- func_ABC_by_permutation
                list_chromosomes <<- list_chromosomes
                list_chromosomes_new <<- list_chromosomes_new
                list_parameters <<- list_parameters
                cn_data_sc <<- cn_data_sc
                cn_data_bulk <<- cn_data_bulk
                arm_level <<- arm_level
                cn_table <<- cn_table
                ###

                clusterExport(cl, varlist = c(
                    "sim_param", "sim_sample_stat", "func_ABC_by_permutation",
                    "list_chromosomes", "list_chromosomes_new", "list_parameters",
                    "cn_data_sc", "cn_data_bulk", "arm_level", "cn_table", "permutate_chromosomes", "get_statistics", "list_targets", "cohort_distance", "cn_distance", "sample_distance"
                ))
                e <- new.env()
                e$libs <- .libPaths()
                clusterExport(cl, "libs", envir = e)
                clusterEvalQ(cl, .libPaths(libs))
                # Create fitted parameters parallel
                pbo <- pboptions(type = "txt")
                sim_results_list <- pblapply(cl = cl, X = 1:nrow(sim_param), FUN = function(iteration) {
                    current_sim_param <- sim_param[iteration, ]
                    current_sim_sample_stat <- sim_sample_stat[[iteration]]
                    output <- func_ABC_by_permutation(
                        current_sim_param = current_sim_param,
                        current_sim_sample_stat = current_sim_sample_stat,
                        current_chromosomes = list_chromosomes,
                        new_chromosomes = list_chromosomes_new,
                        list_parameters = list_parameters,
                        cn_data_sc = cn_data_sc,
                        cn_data_bulk = cn_data_bulk,
                        arm_level = arm_level,
                        cn_table = cn_table
                    )
                    return(output)
                })
                stopCluster(cl)
                #   Group simulated statistics into one table
                for (sim in 1:nrow(sim_param)) {
                    sim_param_next[sim, ] <- sim_results_list[[sim]]$sim_param
                    sim_stat_next[sim, ] <- sim_results_list[[sim]]$sim_stat
                }
                end_time <- Sys.time()
                print(end_time - start_time)
                # ==================================================================================================
                # for (sim in 1:nrow(sim_param)) {
                #     current_sim_param <- sim_param[sim, ]
                #     current_sim_sample_stat <- sim_sample_stat_new[[sim]]
                # df_permutate_chromosomes <- permutate_chromosomes(
                #     current_sim_param = current_sim_param,
                #     current_sim_sample_stat = current_sim_sample_stat,
                #     list_parameters = list_parameters,
                #     current_chromosomes = list_chromosomes,
                #     new_chromosomes = list_chromosomes_new
                # )
                # sim_param_next[sim, ] <- df_permutate_chromosomes$new_sim_param
                # sim_stat_next[sim, ] <- get_statistics(
                #     list_targets = colnames(list_targets)[-1],
                #     simulations_statistics_sc = df_permutate_chromosomes$new_sim_sample_stat$sc,
                #     simulations_statistics_bulk = df_permutate_chromosomes$new_sim_sample_stat$bulk,
                #     cn_data_sc = cn_data_sc,
                #     cn_data_bulk = cn_data_bulk,
                #     arm_level = arm_level,
                #     cn_table = cn_table
                # )$statistics
                #     output <- func_ABC_by_permutation(
                #         current_sim_param = current_sim_param,
                #         current_sim_sample_stat = current_sim_sample_stat,
                #         current_chromosomes = list_chromosomes,
                #         new_chromosomes = list_chromosomes_new,
                #         list_parameters = list_parameters,
                #         cn_data_sc = cn_data_sc,
                #         cn_data_bulk = cn_data_bulk,
                #         arm_level = arm_level,
                #         cn_table = cn_table
                #     )
                #     sim_param_next[sim, ] <- output$sim_param
                #     sim_stat_next[sim, ] <- output$sim_stat
                # }
                sim_param_new <- rbind(sim_param_new, sim_param_next)
                sim_stat_new <- rbind(sim_stat_new, sim_stat_next)
            }
            print(sim_param_new)
            cat("\n")
        }

        #---Boost simulated data by permutating chromosomes
        if (shuffle_chromosomes_by_permutation) {
            #   Initialize new library
            sim_param_new <- sim_param
            sim_stat_new <- sim_stat
            #   Find maximum number of possible permutations
            max_shuffle_count <- factorial(length(list_chromosomes))
            shuffle_num <- min(shuffle_num, max_shuffle_count)
            already_tried_shuffles <- list()
            already_tried_shuffles[[1]] <- list_chromosomes
            for (i in 1:(shuffle_num - 1)) {
                #   Find a new permutation of chromosomes
                vec_comparision <- TRUE
                while (any(vec_comparision)) {
                    list_chromosomes_new <- sample(list_chromosomes, length(list_chromosomes), replace = FALSE)
                    vec_comparision <- rep(FALSE, length(already_tried_shuffles))
                    for (j in 1:length(already_tried_shuffles)) {
                        vec_comparision[j] <- all(list_chromosomes_new == already_tried_shuffles[[j]])
                    }
                }
                already_tried_shuffles[[length(already_tried_shuffles) + 1]] <- list_chromosomes_new
                print("------")
                print(list_chromosomes_new)
                #   Permutate chromosomes
                sim_param_next <- sim_param
                sim_stat_next <- sim_stat
                for (sim in 1:nrow(sim_param)) {
                    current_sim_param <- sim_param[sim, ]
                    current_sim_sample_stat <- sim_sample_stat[[sim]]
                    df_permutate_chromosomes <- permutate_chromosomes(
                        current_sim_param = current_sim_param,
                        current_sim_sample_stat = current_sim_sample_stat,
                        list_parameters = list_parameters,
                        current_chromosomes = list_chromosomes,
                        new_chromosomes = list_chromosomes_new
                    )
                    sim_param_next[sim, ] <- df_permutate_chromosomes$new_sim_param
                    sim_stat_next[sim, ] <- get_statistics(
                        list_targets = colnames(list_targets)[-1],
                        simulations_statistics_sc = df_permutate_chromosomes$new_sim_sample_stat$sc,
                        simulations_statistics_bulk = df_permutate_chromosomes$new_sim_sample_stat$bulk,
                        cn_data_sc = cn_data_sc,
                        cn_data_bulk = cn_data_bulk,
                        arm_level = arm_level,
                        cn_table = cn_table
                    )$statistics
                }
                sim_param_new <- rbind(sim_param_new, sim_param_next)
                sim_stat_new <- rbind(sim_stat_new, sim_stat_next)
            }
            # print("=====")
            # print(sim_param_new)
            sim_param <- sim_param_new
            sim_stat <- sim_stat_new
        }
    }
    # ================================PREPARE SIMULATION LIBRARY FOR ABC
    #   Find ID for each parameter in the prepared library
    sim_param_ID <- list_parameters$Variable
    sim_stat_ID <- colnames(list_targets)[-1]
    #   Prepare the simulation library
    df_sim_param <- data.frame(sim_param)
    colnames(df_sim_param) <- sim_param_ID
    df_sim_stat <- data.frame(sim_stat)
    colnames(df_sim_stat) <- paste0("stat_", 1:ncol(sim_stat))
    # ====================================FITTING WITH ABC RANDOM FOREST
    #---Dataframe for prepared library of parameters
    all_paras <- df_sim_param
    #---Dataframe for prepared library of statistics
    all_data <- df_sim_stat
    #---Dataframe for data observation
    all_obs <- data.frame(matrix(DATA_target, nrow = 1))
    colnames(all_obs) <- paste0("stat_", 1:ncol(sim_stat))
    #---Fit each parameter with ABC-rf
    layout <- matrix(NA, nrow = 7, ncol = ceiling(length(parameter_IDs) / 7))
    gs <- list()
    id <- 0
    for (para in 1:nrow(list_parameters)) {
        start_time <- Sys.time()
        para_ID <- list_parameters$Variable[para]
        para_type <- list_parameters$Type[para]
        cat(paste("\nABC for parameter ", para_ID, " [", para, "/", nrow(list_parameters), "]", "\n", sep = ""))
        #   Prepare observations for this parameter
        mini_obs <- all_obs[, which(list_targets[para, -1] == 1)]
        #   Prepare library of statistics for this parameter
        mini_data <- all_data[, which(list_targets[para, -1] == 1)]
        #   Prepare library of parameters for this parameter
        data_rf <- cbind(all_paras[para_ID], mini_data)
        #   Train the random forest
        colnames(data_rf)[1] <- "para"
        f <- as.formula("para ~.")
        model_rf <- regAbcrf(
            formula = f, data_rf,
            paral = TRUE, ncores = n_cores,
            # ntree = ntree,
            # sampsize = nrow(data_rf),
            # save.memory = TRUE
        )
        #   Predict posterior distribution based on found random forest
        post_rf <- predict(model_rf, mini_obs, data_rf, paral = TRUE, ncores = n_cores)
        #   Choose best values from posterior distribution
        df_dist <- densityPlot_df(model_rf, mini_obs, data_rf)
        post_mean <- weighted.mean(df_dist$x, df_dist$y_posterior)
        post_median <- weightedMedian(df_dist$x, df_dist$y_posterior)
        post_mode <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
        if (!is.null(parameters_truth)) cat("True value: ", parameters_truth$Value[which(parameters_truth$Variable == para_ID)], "\n")
        cat("Posterior mean: ", post_mean, "\n")
        cat("Posterior median: ", post_median, "\n")
        cat("Posterior mode: ", post_mode, "\n")
        #   Save results for fitting this parameter
        ABC_output <- list()
        ABC_output$para_ID <- para_ID
        ABC_output$post_mode <- post_mode
        ABC_output$post_mean <- post_mean
        ABC_output$post_median <- post_median
        filename <- paste0(model_name, "_ABC_output_", para_ID, ".rda")
        save(ABC_output, file = filename)
        #   Plot the prior, posterior and chosen best parameter for all variables
        true_para <- NULL
        if (!is.null(parameters_truth)) {
            true_para <- parameters_truth$Value[which(parameters_truth$Variable == para_ID)]
        }
        id <- id + 1
        row <- id %% 7
        if (row == 0) row <- 7
        col <- ceiling(id / 7)
        layout[row, col] <- id
        gs[[id]] <- plot_parameter_ABC(
            model_rf, all_obs, data_rf,
            protocol = "arm",
            ###
            highlight_values = c(true_para, post_mode, post_mean, post_median),
            highlight_colors = c("black", "#d03333", "#3939ae", "#3da53d"),
            highlight_linetype = c("solid", "dashed", "dashed", "dashed"),
            ###
            fontsize = 20,
            main = para_ID
        )
        end_time <- Sys.time()
        cat(paste0("Best parameter: ", post_mode, "\n"))
        print(end_time - start_time)
        #   Clear memory
        model_rf <- c()
        data_rf <- c()
        post_rf <- c()
        post_mode <- c()
    }
    #   Plot the prior, posterior and chosen best parameter for all variables
    filename <- paste0(model_name, "_ABC_all.jpeg")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
    #---------------------------------Save shuffled simulated CN library
    if (shuffle_chromosome_arms | shuffle_chromosomes_by_permutation | shuffle_chromosomes_by_moving) {
        ABC_input$sim_param <- sim_param
        ABC_input$sim_stat <- sim_stat
        filename <- paste0(library_name, "_ABC_input_shuffled.rda")
        save(ABC_input, file = filename)
    }
}
