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

get_clonal_CN_profiles <- function(simulations,
                                   arm_level = FALSE,
                                   cn_table = NULL) {
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
    return(clonal_CN_profiles_all_sims)
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
find_event_count <- function(clone_group, event_type, evolution_genotype_changes) {
    event_count <- 0
    for (clone in clone_group) {
        if (clone <= 0) next
        if (length(evolution_genotype_changes[[clone]]) == 0) next
        for (j in 1:length(evolution_genotype_changes[[clone]])) {
            if (evolution_genotype_changes[[clone]][[j]][1] == event_type) {
                event_count <- event_count + 1
            }
        }
    }
    return(event_count)
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
            vec1 <- from_clonal_CN_profiles[[i]]$copy
            vec2 <- to_clonal_CN_profiles[[j]]$copy
            if (metric == "euclidean") {
                cost_matrix[i, j] <- sqrt(sum((vec1 - vec2)^2))
            } else if (metric == "wcnd") {
                cost_matrix[i, j] <- 0
            } else if (metric == "correlation") {
                library(energy)
                cost_matrix[i, j] <- dcor(vec1, vec2)
            } else if (metric == "mahalanobis") {
                library(MASS)
                cost_matrix[i, j] <- sqrt(t(vec1 - vec2) %*% ginv(cov(vec1 %*% t(vec1))) %*% (vec1 - vec2))
            }
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
cohort_distance <- function(cohort_from, cohort_to, metric) {
    library(transport)
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
    #   Normalize clonal populations
    dist_from <- rep(1 / n_cohort_from, n_cohort_from)
    dist_to <- rep(1 / n_cohort_to, n_cohort_to)
    #   Compute the Wasserstein distance between the two cohorts
    wasserstein_dist <- wasserstein(
        dist_from,
        dist_to,
        costm = sample_distance_mtx, p = 1
    )
    return(wasserstein_dist)
}

#---------------------------------Get statistics for each simulation
get_statistics_simulations <- function(simulations, simulations_clonal_CN, list_targets) {
    simulations_statistics <- list()
    for (i in 1:length(simulations)) {
        simulation <- simulations[[i]]
        #   Find common clonal ancestors (for event counts)
        if (any(grepl("variable=event_count", list_targets))) {
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
        #   Get requested statistics for this simulation
        for (stat in list_targets) {
            stat_details <- strsplit(stat, ";")[[1]]
            stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
            stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
            if (stat_variable == "shannon") {
                #   Get Shannon index
                frequencies <- table(simulation$sample$sample_clone_ID)
                simulations_statistics[[stat_ID]][i] <- diversity(frequencies, index = "shannon")
            } else if (stat_variable == "event_count") {
                #   Get count of clonal/subclonal events of a given type
                clonal_type <- strsplit(stat_details[grep("type=", stat_details)], "=")[[1]][2]
                event_type <- strsplit(stat_details[grep("event=", stat_details)], "=")[[1]][2]
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
            } else if (stat_variable == "IL_number") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get IL_number
                list_statistics_simulations[[stat_ID]][i] <- ILnumber(tree, normalise = TRUE)
            } else if (stat_variable == "avgLadder") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get avgLadder
                list_statistics_simulations[[stat_ID]][i] <- avgLadder(tree, normalise = TRUE)
            } else if (stat_variable == "maxDepth") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get maxDepth
                list_statistics_simulations[[stat_ID]][i] <- maxDepth(tree)
            } else if (stat_variable == "ColPla") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get ColPla index
                list_statistics_simulations[[stat_ID]][i] <- colPlaLab(tree, method = "binary")
            } else if (stat_variable == "stairs") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get stairness
                list_statistics_simulations[[stat_ID]][i] <- stairs(tree)[1]
            } else if (stat_variable == "B2") {
                #   Get cell phylogeny tree
                tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                #   Get B2Index
                list_statistics_simulations[[stat_ID]][i] <- B2I(tree, logbase = 2)
            } else if (stat_variable == "clonal_CN") {
                #   Get clonal CN profiles and their populations
                simulations_statistics[["variable=clonal_CN_profiles"]][[i]] <-
                    simulations_clonal_CN[["variable=clonal_CN_profiles"]][[i]]
                simulations_statistics[["variable=clonal_CN_populations"]][[i]] <-
                    simulations_clonal_CN[["variable=clonal_CN_populations"]][[i]]
            } else {
                stop(paste0("Error: Unknown statistic: ", stat))
            }
        }
    }
    return(simulations_statistics)
}

#' @export
get_statistics <- function(simulations = NULL,
                           simulations_statistics = NULL,
                           list_targets,
                           cn_data = NULL,
                           arm_level = FALSE,
                           cn_table = NULL,
                           save_sample_statistics = FALSE) {
    library(vegan)
    library(matrixStats)
    library(transport)
    library(ape)
    library(phyloTop)
    #-------------------------Get clonal CN profiles for all simulations
    if (is.null(simulations_statistics)) {
        if (any(grepl("variable=clonal_CN", list_targets))) {
            simulations_clonal_CN <- get_clonal_CN_profiles(simulations, arm_level, cn_table)
        }
        simulations_statistics <- get_statistics_simulations(simulations, simulations_clonal_CN, list_targets)
    }
    #-----------------------------------------Get statistics for fitting
    statistics <- rep(0, length(list_targets))
    for (i in 1:length(list_targets)) {
        stat <- list_targets[i]
        stat_details <- strsplit(stat, ";")[[1]]
        stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
        stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
        stat_type <- strsplit(stat_details[grepl("statistic=", stat_details)], "=")[[1]][2]
        if (stat_type == "mean") {
            statistics[i] <- mean(simulations_statistics[[stat_ID]])
        } else if (stat_type == "var") {
            statistics[i] <- var(simulations_statistics[[stat_ID]])
        } else if (stat_type == "dist") {
            if (stat_variable == "clonal_CN") {
                stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                if (is.null(cn_data)) cn_data <- simulations_statistics
                statistics[i] <- cohort_distance(
                    cohort_from = simulations_statistics,
                    cohort_to = cn_data,
                    metric = stat_metric
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
        output$simulations_statistics <- simulations_statistics
    }
    return(output)
}

#---Function to assign parameters to proper positions
assign_paras <- function(model_variables, parameter_IDs, parameters) {
    for (i in 1:length(parameter_IDs)) {
        parameter_ID <- parameter_IDs[i]
        if (parameter_ID %in% model_variables$general_variables$Variable) {
            model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameters[i]
        } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
            model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)] <- parameters[i]
        }
    }
    return(model_variables)
}
#---Objective function for ABC fitting
func_ABC <- function(parameters,
                     parameter_IDs,
                     model_variables,
                     list_targets,
                     cn_data = NULL,
                     arm_level = FALSE) {
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
    #   Make simulations
    SIMS_chromosome <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_simulations,
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
    #   Get statistics from simulations
    stat <- get_statistics(
        simulations = SIMS_chromosome,
        list_targets = list_targets,
        cn_data = cn_data,
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
                          list_targets,
                          ABC_simcount = 10000,
                          n_simulations = 10,
                          n_cores = NULL,
                          library_name,
                          cn_data = NULL,
                          cn_table = NULL,
                          arm_level = FALSE,
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
    list_targets <<- list_targets
    cn_data <<- cn_data
    ####
    ####
    ####
    cn_table <<- cn_table
    save_sample_statistics <<- save_sample_statistics
    arm_level <<- arm_level
    ####
    ####
    ####
    n_simulations <<- n_simulations
    clusterExport(cl, varlist = c(
        "n_simulations", "list_targets", "sim_param", "parameter_IDs", "model_variables", "cn_data", "cn_table", "arm_level",
        "func_ABC", "assign_paras", "get_statistics", "get_clonal_CN_profiles", "save_sample_statistics", "get_cn_profile", "get_arm_CN_profiles",
        "find_clonal_ancestry", "find_event_count", "cohort_distance", "sample_distance", "get_statistics_simulations",
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
        stat <- func_ABC(parameters, parameter_IDs, model_variables, list_targets, cn_data = cn_data, arm_level = arm_level)
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
            sim_sample_stat[[row]] <- sim_results_list[[row]]$simulations_statistics
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
    #--------------------------------Permutate simulation parameters
    new_sim_param <- current_sim_param
    for (i in 1:length(current_chromosomes)) {
        locs_old_chrom <- which(list_parameters$Chromosome == current_chromosomes[i])
        locs_new_chrom <- which(list_parameters$Chromosome == new_chromosomes[i])
        new_sim_param[locs_old_chrom] <- current_sim_param[locs_new_chrom]
    }
    #--------------------------------Permutate simulation statistics
    new_sim_sample_stat <- current_sim_sample_stat
    for (sim in 1:length(current_sim_sample_stat[["variable=clonal_CN_profiles"]])) {
        for (clone in 1:length(current_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]])) {
            current_clonal_CN_profiles <- current_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]][[clone]]
            new_clonal_CN_profiles <- current_clonal_CN_profiles
            for (i in 1:length(current_chromosomes)) {
                locs_old_chrom <- which(current_clonal_CN_profiles$chr == current_chromosomes[i])
                locs_new_chrom <- which(current_clonal_CN_profiles$chr == new_chromosomes[i])
                new_clonal_CN_profiles$copy[locs_old_chrom] <- current_clonal_CN_profiles$copy[locs_new_chrom]
                new_clonal_CN_profiles$state[locs_old_chrom] <- current_clonal_CN_profiles$state[locs_new_chrom]
                new_clonal_CN_profiles$Min[locs_old_chrom] <- current_clonal_CN_profiles$Min[locs_new_chrom]
                new_clonal_CN_profiles$Maj[locs_old_chrom] <- current_clonal_CN_profiles$Maj[locs_new_chrom]
                new_sim_sample_stat[["variable=clonal_CN_profiles"]][[sim]][[clone]] <- new_clonal_CN_profiles
            }
        }
    }
    #-------------------------------------------------Prepare output
    output <- list()
    output$new_sim_param <- new_sim_param
    output$new_sim_sample_stat <- new_sim_sample_stat
    return(output)
}

#' @export
fitting_sc_CN <- function(library_name,
                          model_name,
                          copynumber_DATA,
                          parameters_truth = NULL,
                          list_parameters,
                          list_targets,
                          n_cores = NULL,
                          cn_data = NULL,
                          cn_table = NULL,
                          arm_level = FALSE,
                          shuffle_chromosome_arms = FALSE,
                          shuffle_chromosomes = FALSE) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    library(data.table)
    library(matrixStats)
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
    if ((shuffle_chromosome_arms | shuffle_chromosomes) & (any(grepl("variable=clonal_CN_profiles", names(sim_sample_stat[[1]]))))) {
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
            sim_param_new <- sim_param
            sim_stat_new <- sim_stat
            sim_sample_stat_new <- sim_sample_stat
            for (chromosome in list_chromosomes) {
                if (length(which(list_parameters$Chromosome == chromosome)) <= 1) next
                if (length(which(list_parameters$Chromosome == chromosome)) > 2) simpleError("More than two arms for a chromosome")
                sim_param_next <- sim_param_new
                sim_stat_next <- sim_stat_new
                for (sim in 1:nrow(sim_param_new)) {
                    current_sim_param <- sim_param_new[sim, ]
                    current_sim_sample_stat <- sim_sample_stat_new[[sim]]
                    df_permutate_chromosome_arms <- permutate_chromosome_arms(
                        current_sim_param = current_sim_param,
                        current_sim_sample_stat = current_sim_sample_stat,
                        list_parameters = list_parameters,
                        permutate_chromosome = chromosome
                    )
                    new_sim_param <- df_permutate_chromosome_arms$new_sim_param
                    new_sim_sample_stat <- df_permutate_chromosome_arms$new_sim_sample_stat
                    sim_param_next <- rbind(sim_param_next, new_sim_param)
                    sim_stat_next <- rbind(sim_stat_next, get_statistics(
                        simulations_statistics = new_sim_sample_stat,
                        list_targets = list_targets,
                        cn_data = cn_data,
                        arm_level = arm_level,
                        cn_table = cn_table
                    )$statistics)
                }
                sim_param_new <- sim_param_next
                sim_stat_new <- sim_stat_next
            }
            sim_param <- sim_param_new
            sim_stat <- sim_stat_new
        }
        #---Boost simulated data by permutating chromosomes
        if (shuffle_chromosomes) {
            sim_param_new <- sim_param
            sim_stat_new <- sim_stat
            sim_sample_stat_new <- sim_sample_stat
            for (i in 2:length(list_chromosomes)) {
                list_chromosomes_new <- list_chromosomes[c(i:length(list_chromosomes), 1:i - 1)]
                sim_param_next <- sim_param
                sim_stat_next <- sim_stat
                for (sim in 1:nrow(sim_param)) {
                    current_sim_param <- sim_param[sim, ]
                    current_sim_sample_stat <- sim_sample_stat_new[[sim]]
                    df_permutate_chromosomes <- permutate_chromosomes(
                        current_sim_param = current_sim_param,
                        current_sim_sample_stat = current_sim_sample_stat,
                        list_parameters = list_parameters,
                        current_chromosomes = list_chromosomes,
                        new_chromosomes = list_chromosomes_new
                    )
                    sim_param_next[sim, ] <- df_permutate_chromosomes$new_sim_param
                    sim_stat_next[sim, ] <- get_statistics(
                        simulations_statistics = df_permutate_chromosomes$new_sim_sample_stat,
                        list_targets = list_targets,
                        cn_data = cn_data,
                        arm_level = arm_level,
                        cn_table = cn_table
                    )$statistics
                }
                sim_param_new <- rbind(sim_param_new, sim_param_next)
                sim_stat_new <- rbind(sim_stat_new, sim_stat_next)
            }
            sim_param <- sim_param_new
            sim_stat <- sim_stat_new
        }
    }
    # ================================PREPARE SIMULATION LIBRARY FOR ABC
    #   Find ID for each parameter in the prepared library
    sim_param_ID <- list_parameters$Variable
    sim_stat_ID <- list_targets
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
        mini_obs <- all_obs
        #   Prepare library of statistics for this parameter
        mini_data <- all_data
        #   Prepare library of parameters for this parameter
        data_rf <- cbind(all_paras[para_ID], all_data)
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
        # filename <- paste0(folder_workplace_tmp, model_name, "_ABC_output_", para_ID, ".rda")
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





    print(df_sim_param)
    print(df_sim_stat)
    print(all_obs)
}
