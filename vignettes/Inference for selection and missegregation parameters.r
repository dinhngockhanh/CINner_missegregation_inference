library(CINner)
library(readxl)
library(parallel)
library(pbapply)
# ==================================================IMPORTANT PARAMETERS
#   Number of single-cell samples in ground-truth data & ABC simulations
N_data_sc <- 50
#   Number of bulk samples in ground-truth data & ABC simulations
N_data_bulk <- 100
#   Bounds for ground-truth selection rates (1/r -> r)
bound_ground_truth_arm_s <- 1.15
#   Bounds for prior distribution of log10(prob_CN_missegregation)
bound_ABC_prob_CN_missegregation_left <- -5
bound_ABC_prob_CN_missegregation_right <- -3
#   Bounds for prior distribution of selection rates
bound_ABC_arm_s <- 1.2
# ===============================================GROUND TRUTH PARAMETERS
cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = c("SA01"), Cell_count = c(1000), Age_sample = c(80))
selection_model <- "chrom-arm-selection"
#---Probabilities of CNA
prob_CN_missegregation <- 2e-4
#---Viability thresholds
bound_maximum_CN <- 8
bound_average_ploidy <- 4.5
bound_homozygosity <- 10
#---Population dynamics
vec_time <- T_0[[1]]:T_end[[1]]
L <- 10000
t_0 <- 20
k <- 0.3
vec_cell_count <- L / (1 + exp(-k * (vec_time - t_0)))
table_population_dynamics <- cbind(vec_time, vec_cell_count)
#---Initialize model variables
model_variables <- BUILD_general_variables(
    cell_lifespan = cell_lifespan,
    T_0 = T_0, T_end = T_end,
    CN_arm_level = TRUE,
    Table_sample = Table_sample,
    prob_CN_missegregation = prob_CN_missegregation,
    selection_model = selection_model,
    bound_average_ploidy = bound_average_ploidy,
    bound_homozygosity = bound_homozygosity,
    table_population_dynamics = table_population_dynamics
)
#---Set up (randomized) chromosome arm selection rates
arm_id <- c(
    paste0(model_variables$cn_info$Chromosome, "p"),
    paste0(model_variables$cn_info$Chromosome, "q")
)
arm_chromosome <- rep(model_variables$cn_info$Chromosome, 2)
arm_start <- c(
    rep(1, length(model_variables$cn_info$Chromosome)),
    model_variables$cn_info$Centromere_location + 1
)
arm_end <- c(
    model_variables$cn_info$Centromere_location,
    model_variables$cn_info$Bin_count
)
arm_s <- rep(1, length(arm_id))
set.seed(1)
for (i in 1:length(arm_s)) {
    if (grepl("q$", arm_id[i])) {
        arm_s[i] <- 1
    } else if (grepl("p$", arm_id[i])) {
        arm_s[i] <- runif(1, 1, bound_ground_truth_arm_s)
        if (runif(1) < 0.5) arm_s[i] <- 1 / arm_s[i]
    }
}
set.seed(NULL)
table_arm_selection_rates <- data.frame(Arm_ID = arm_id, Chromosome = arm_chromosome, Bin_start = arm_start, Bin_end = arm_end, s_rate = arm_s)
model_variables <- BUILD_driver_library(model_variables = model_variables, table_arm_selection_rates = table_arm_selection_rates, )
#---Set up initial cell population
cell_count <- 20
CN_matrix <- BUILD_cn_normal_autosomes(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(model_variables = model_variables, cell_count = cell_count, CN_matrix = CN_matrix, drivers = drivers)
#---Save model variables
model_name <- "Fitting_whole_chroms"
model_variables <- CHECK_model_variables(model_variables)
# ======================================DEFINE LIST OF PARAMETERS TO FIT
list_parameters <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(list_parameters) <- c("Variable", "Title", "Chromosome", "Type", "Lower_bound", "Upper_bound")
list_parameters[nrow(list_parameters) + 1, ] <- c(
    "10^:prob_CN_missegregation", "log10(prob_misseg)", NA, "CNA_probability",
    bound_ABC_prob_CN_missegregation_left, bound_ABC_prob_CN_missegregation_right
)
for (i in 1:nrow(model_variables$chromosome_arm_library)) {
    if (grepl("p$", model_variables$chromosome_arm_library$Arm_ID[i])) {
        list_parameters[nrow(list_parameters) + 1, ] <- c(
            model_variables$chromosome_arm_library$Arm_ID[i], paste0("Selection rate - chromosome ", model_variables$chromosome_arm_library$Chromosome[i]), model_variables$chromosome_arm_library$Chromosome[i], "Selection_rate",
            1 / bound_ABC_arm_s, bound_ABC_arm_s
        )
    }
}
list_chromosomes <- paste(unique(model_variables$chromosome_arm_library$Chromosome), collapse = ",")
# =============DEFINE LIST OF STATISTICS FOR BUILDING SIMULATION LIBRARY
list_targets_library <- c(
    #---Bulk DNA: CN
    "data=bulk;target=genome;statistic=dist;variable=average_CN;metric=euclidean",
    "data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    "data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    paste0("data=bulk;target=chromosome;statistic=dist;variable=average_CN;metric=euclidean;chromosome=", list_chromosomes),
    paste0("data=bulk;target=chromosome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=bulk;target=chromosome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    #---Single-cell DNA: subclonal CN
    "data=sc;target=genome;statistic=dist;variable=clonal_CN;metric=euclidean",
    "data=sc;target=genome;statistic=mean;variable=shannon",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "data=sc;target=genome;statistic=var;variable=shannon",
    "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    paste0("data=sc;target=chromosome;statistic=dist;variable=clonal_CN;metric=euclidean;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=shannon;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=shannon;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    #---Single-cell DNA: phylo stats for tips
    "data=sc;target=genome;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;target=genome;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;target=genome;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;target=genome;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;target=genome;statistic=var;variable=cherries",
    "data=sc;target=genome;statistic=var;variable=pitchforks",
    "data=sc;target=genome;statistic=var;variable=IL_number",
    "data=sc;target=genome;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;target=genome;statistic=var;variable=stairs",
    "data=sc;target=genome;statistic=var;variable=colless",
    "data=sc;target=genome;statistic=var;variable=sackin",
    "data=sc;target=genome;statistic=var;variable=B2",
    "data=sc;target=genome;statistic=var;variable=maxDepth"
)
# ==============GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
cn_table <- model_variables$cn_info
cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
cn_table$Length <- cn_table$Bin_count * cn_bin_length
cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length
vec_CN_block_no <<- model_variables$cn_info$Bin_count
vec_centromeres <<- model_variables$cn_info$Centromere_location
# ================================================MAKE GROUND-TRUTH DATA
#---Make single-cell ground-truth simulations
cat(paste0("\n\n\nMaking ", N_data_sc, " single-cell simulations...\n"))
simulator_full_program(
    model = model_variables, model_prefix = paste0(model_name, "_sc"),
    n_simulations = N_data_sc,
    stage_final = 3,
    compute_parallel = TRUE,
    output_variables = c(
        "evolution_origin",
        "evolution_genotype_changes",
        "sample_clone_ID",
        "sample_genotype_unique",
        "sample_genotype_unique_profile",
        "phylogeny_clustering_truth"
    ),
    R_libPaths = R_libPaths
)
#---Make bulk ground-truth simulations
cat(paste0("\n\n\nMaking ", N_data_bulk, " bulk simulations...\n"))
simulator_full_program(
    model = model_variables, model_prefix = paste0(model_name, "_bulk"),
    n_simulations = N_data_bulk,
    stage_final = 2,
    compute_parallel = TRUE,
    output_variables = c(
        "sample_genotype_unique_profile",
        "sample_genotype_unique",
        "sample_clone_ID"
    ),
    R_libPaths = R_libPaths
)
# ---Print out ground-truth parameters
list_parameters_ground_truth <- list_parameters
list_parameters_ground_truth$Value <- 0
for (row in 1:nrow(list_parameters)) {
    parameter_ID_input <- list_parameters$Variable[row]
    #   Convert parameter operator if necessary
    if (grepl(":", parameter_ID_input)) {
        parameter_ID <- sub(".*:", "", parameter_ID_input)
        parameter_operator <- sub(":.*", "", parameter_ID_input)
    } else {
        parameter_ID <- parameter_ID_input
        parameter_operator <- ""
    }
    if (parameter_ID %in% model_variables$general_variables$Variable) {
        parameter_value_input <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)])
    } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
        parameter_value_input <- as.numeric(model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)])
    }
    if (parameter_operator == "") {
        parameter_value <- parameter_value_input
    } else if (parameter_operator == "10^") {
        parameter_value <- log10(parameter_value_input)
    } else {
        simpleError("Parameter operator not recognized")
    }
    #   Assign parameter value
    list_parameters_ground_truth$Value[row] <- parameter_value
}
write.csv(list_parameters_ground_truth, "parameters_ground_truth.csv")
# ============GET STATISTICS & CN PROFILES FROM GROUND-TRUTH SIMULATIONS
#---Get single-cell statistics & CN profiles
#   Get statistics & clonal CN profiles for each single-cell sample
list_targets_library_sc <- list_targets_library[grepl("data=sc", list_targets_library)]
cat(paste0("Loading ", N_data_sc, " single-cell DNA-seq data sets...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
model_name <<- model_name
clusterExport(cl, varlist = c(
    "model_name", "get_each_clonal_CN_profiles", "get_arm_CN_profiles",
    "cn_table", "get_each_statistics", "list_targets_library_sc", "find_clonal_ancestry", "find_event_count"
))
e <- new.env()
e$libs <- .libPaths()
clusterExport(cl, "libs", envir = e)
clusterEvalQ(cl, .libPaths(libs))
pbo <- pboptions(type = "txt")
ls_cn_sc_ground_truth <- pblapply(cl = cl, X = 1:N_data_sc, FUN = function(i) {
    load(paste0(model_name, "_sc_simulation_", i, ".rda"))
    simulations <- list()
    simulations[[1]] <- simulation
    ls_each_sim <- list()
    ls_each_sim[[1]] <- get_each_clonal_CN_profiles(
        simulations,
        arm_level = TRUE,
        cn_table = cn_table
    )
    ls_each_sim[[2]] <- get_each_statistics(simulations, ls_each_sim[[1]], list_targets_library_sc)
    return(ls_each_sim)
})
stopCluster(cl)
#   Get statistics & clonal CN profiles for entire single-cell cohort
ground_truth_cn_data_sc <- list()
ground_truth_statistics_sc <- list()
for (simulation in 1:N_data_sc) {
    for (statistic in 1:length(ls_cn_sc_ground_truth[[1]][[1]])) {
        if (simulation == 1) {
            ground_truth_cn_data_sc[[statistic]] <- ls_cn_sc_ground_truth[[simulation]][[1]][[statistic]][1]
        } else {
            ground_truth_cn_data_sc[[statistic]] <- c(ground_truth_cn_data_sc[[statistic]], ls_cn_sc_ground_truth[[simulation]][[1]][[statistic]][1])
        }
    }
    names(ground_truth_cn_data_sc) <- names(ls_cn_sc_ground_truth[[1]][[1]])
    for (stat_ID in names(ls_cn_sc_ground_truth[[1]][[2]])) {
        stat_details <- strsplit(stat_ID, ";")[[1]]
        if (simulation == 1) {
            ground_truth_statistics_sc[[stat_ID]] <- ls_cn_sc_ground_truth[[1]][[2]][[stat_ID]]
        } else {
            ground_truth_statistics_sc[[stat_ID]] <- rbind(ground_truth_statistics_sc[[stat_ID]], ls_cn_sc_ground_truth[[simulation]][[2]][[stat_ID]])
        }
    }
    names(ground_truth_statistics_sc) <- names(ls_cn_sc_ground_truth[[1]][[2]])
}
#---Get bulk statistics & CN profiles
#   Get statistics & representative CN profiles for each bulk sample
list_targets_library_bulk <- list_targets_library[grepl("data=bulk", list_targets_library)]
cat(paste0("Loading ", N_data_bulk, " bulk DNA-seq data sets...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
model_name <<- model_name
clusterExport(cl, varlist = c(
    "model_name", "N_data_sc", "get_each_clonal_CN_profiles", "get_arm_CN_profiles",
    "cn_table", "get_each_statistics", "list_targets_library_bulk", "find_clonal_ancestry", "find_event_count"
))
e <- new.env()
e$libs <- .libPaths()
clusterExport(cl, "libs", envir = e)
clusterEvalQ(cl, .libPaths(libs))
pbo <- pboptions(type = "txt")
ls_cn_bulk_ground_truth <- pblapply(cl = cl, X = 1:N_data_bulk, FUN = function(i) {
    load(paste0(model_name, "_bulk_simulation_", i, ".rda"))
    simulations <- list()
    simulations[[1]] <- simulation
    ls_each_sim <- list()
    ls_each_sim[[1]] <- get_each_clonal_CN_profiles(
        simulations,
        arm_level = TRUE,
        cn_table = cn_table,
        bulk = TRUE
    )
    ls_each_sim[[2]] <- get_each_statistics(simulations, ls_each_sim[[1]], list_targets_library_bulk)
    return(ls_each_sim)
})
#   Get statistics & representative CN profiles for entire bulk cohort
ground_truth_cn_data_bulk <- list()
ground_truth_statistics_bulk <- list()
for (simulation in 1:N_data_bulk) {
    for (statistic in 1:length(ls_cn_bulk_ground_truth[[1]][[1]])) {
        if (simulation == 1) {
            ground_truth_cn_data_bulk[[statistic]] <- ls_cn_bulk_ground_truth[[simulation]][[1]][[statistic]][1]
        } else {
            ground_truth_cn_data_bulk[[statistic]] <- c(ground_truth_cn_data_bulk[[statistic]], ls_cn_bulk_ground_truth[[simulation]][[1]][[statistic]][1])
        }
    }
    names(ground_truth_cn_data_bulk) <- names(ls_cn_bulk_ground_truth[[1]][[1]])
    for (stat_ID in names(ls_cn_bulk_ground_truth[[1]][[2]])) {
        stat_details <- strsplit(stat_ID, ";")[[1]]
        if (simulation == 1) {
            ground_truth_statistics_bulk[[stat_ID]] <- ls_cn_bulk_ground_truth[[1]][[2]][[stat_ID]]
        } else {
            ground_truth_statistics_bulk[[stat_ID]] <- rbind(ground_truth_statistics_bulk[[stat_ID]], ls_cn_bulk_ground_truth[[simulation]][[2]][[stat_ID]])
        }
    }
    names(ground_truth_statistics_bulk) <- names(ls_cn_bulk_ground_truth[[1]][[2]])
}
# ===============================================MAKE SIMULATION LIBRARY
ABC_simcount <- 100000
library_simulations(
    library_name = model_name,
    model_variables = model_variables,
    list_parameters = list_parameters,
    list_targets_library = list_targets_library,
    ABC_simcount_start = 0,
    ABC_simcount = ABC_simcount,
    arm_level = TRUE,
    cn_table = cn_table,
    n_simulations_sc = N_data_sc,
    n_simulations_bulk = N_data_bulk
)
# ==================DEFINE LIST OF STATISTICS FOR FITTING EACH PARAMETER
list_targets <- data.frame(matrix(0, ncol = (length(list_targets_library) + 1), nrow = length(list_parameters$Variable)))
colnames(list_targets) <- c("Variable", list_targets_library)
list_targets[, 1] <- list_parameters$Variable
#---Statistics for fitting misseg rate
list_targets_misseg <- c(
    #---Bulk DNA: CN
    "data=bulk;target=genome;statistic=dist;variable=average_CN;metric=euclidean",
    "data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    "data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    #---Single-cell DNA: subclonal CN
    "data=sc;target=genome;statistic=dist;variable=clonal_CN;metric=euclidean",
    "data=sc;target=genome;statistic=mean;variable=shannon",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=shannon",
    "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    #---Single-cell DNA: phylo stats for tips
    "data=sc;target=genome;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;target=genome;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;target=genome;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;target=genome;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;target=genome;statistic=var;variable=cherries",
    "data=sc;target=genome;statistic=var;variable=pitchforks",
    "data=sc;target=genome;statistic=var;variable=IL_number",
    "data=sc;target=genome;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;target=genome;statistic=var;variable=stairs",
    "data=sc;target=genome;statistic=var;variable=colless",
    "data=sc;target=genome;statistic=var;variable=sackin",
    "data=sc;target=genome;statistic=var;variable=B2",
    "data=sc;target=genome;statistic=var;variable=maxDepth"
)
list_targets[1, which(colnames(list_targets) %in% list_targets_misseg)] <- 1
#---Statistics for fitting selection rate
list_targets_selection <- c(
    #---Bulk DNA: CN
    paste0("data=bulk;target=chromosome;statistic=dist;variable=average_CN;metric=euclidean;chromosome=", list_chromosomes),
    paste0("data=bulk;target=chromosome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=bulk;target=chromosome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    #---Single-cell DNA: subclonal CN
    paste0("data=sc;target=chromosome;statistic=dist;variable=clonal_CN;metric=euclidean;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=shannon;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=shannon;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    #---Single-cell DNA: phylo stats for tips
    "data=sc;target=genome;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;target=genome;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;target=genome;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;target=genome;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;target=genome;statistic=var;variable=cherries",
    "data=sc;target=genome;statistic=var;variable=pitchforks",
    "data=sc;target=genome;statistic=var;variable=IL_number",
    "data=sc;target=genome;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;target=genome;statistic=var;variable=stairs",
    "data=sc;target=genome;statistic=var;variable=colless",
    "data=sc;target=genome;statistic=var;variable=sackin",
    "data=sc;target=genome;statistic=var;variable=B2",
    "data=sc;target=genome;statistic=var;variable=maxDepth"
)
for (row in 2:nrow(list_targets)) {
    list_targets[row, which(colnames(list_targets) %in% list_targets_selection)] <- 1
}
# =============================COMPUTE STATISTICS FOR SIMULATION LIBRARY
# library_statistics(
#     library_name = model_name,
#     model_variables = model_variables,
#     list_parameters = list_parameters,
#     list_targets_library = list_targets_library,
#     ABC_simcount_start = 0,
#     ABC_simcount = ABC_simcount,
#     cn_data_sc = ground_truth_cn_data_sc,
#     cn_data_bulk = ground_truth_cn_data_bulk,
#     arm_level = TRUE,
#     cn_table = cn_table
# )
# ==================PLOT CORRELATION OF STATISTICS IN SIMULATION LIBRARY
plot_statistics_correlation(
    filename = "Fitting_whole_chroms_ABC_input.rda"
)
# =========================GET FITTING STATISTICS FROM GROUND-TRUTH DATA
# DLP_stats <- get_statistics(
#     simulations_statistics_sc = ground_truth_statistics_sc,
#     simulations_statistics_bulk = ground_truth_statistics_bulk,
#     list_targets = list_targets_library,
#     cn_data_sc = ground_truth_cn_data_sc,
#     cn_data_bulk = ground_truth_cn_data_bulk,
#     arm_level = TRUE,
#     cn_table = cn_table
# )

# ==============================================FIT PARAMETERS USING ABC
stat_names <- c(
    "Bulk CN distance",
    "Mean(misseg. count in bulk)",
    "Var(misseg. count in bulk)",
    "Bulk CN distance",
    "Mean(misseg. count in bulk)",
    "Var(misseg. count in bulk)",
    #---Single-cell DNA: subclonal CN
    "Single-cell CN distance",
    "Mean(sc Shannon index)",
    "Mean(clonal misseg. count in sc)",
    "Mean(subclonal misseg. count in sc)",
    "Mean(clonal arm-misseg. count in sc)",
    "Mean(subclonal arm-misseg. count in sc)",
    "Var(sc Shannon index)",
    "Var(clonal misseg. count in sc)",
    "Var(subclonal misseg. count in sc)",
    "Var(clonal arm-misseg. count in sc)",
    "Var(subclonal arm-misseg. count in sc)",
    "Single-cell CN distance",
    "Mean(sc Shannon index)",
    "Mean(clonal misseg. count in sc)",
    "Mean(subclonal misseg. count in sc)",
    "Mean(clonal arm-misseg. count in sc)",
    "Mean(subclonal arm-misseg. count in sc)",
    "Var(sc Shannon index)",
    "Var(clonal misseg. count in sc)",
    "Var(subclonal misseg. count in sc)",
    "Var(clonal arm-misseg. count in sc)",
    "Var(subclonal arm-misseg. count in sc)",
    #---Single-cell DNA: phylo stats for tips
    "Mean(cherry count)",
    "Mean(pitchfork count)",
    "Mean(IL number)",
    "Mean(average ladder)",
    "Var(cherry count)",
    "Var(pitchfork count)",
    "Var(IL number)",
    "Var(average ladder)",
    "Mean(stairs)",
    "Mean(Colless index)",
    "Mean(Sackin index)",
    "Mean(B2 index)",
    "Mean(max depth)",
    "Var(stairs)",
    "Var(Colless index)",
    "Var(Sackin index)",
    "Var(B2 index)",
    "Var(max depth)"
)
phylo_tip_stats <- c(
    "Mean(cherry count)",
    "Mean(pitchfork count)",
    "Mean(IL number)",
    "Mean(average ladder)",
    "Var(cherry count)",
    "Var(pitchfork count)",
    "Var(IL number)",
    "Var(average ladder)"
)

phylo_balance_stats <- c(
    "Mean(stairs)",
    "Mean(Colless index)",
    "Mean(Sackin index)",
    "Mean(B2 index)",
    "Mean(max depth)",
    "Var(stairs)",
    "Var(Colless index)",
    "Var(Sackin index)",
    "Var(B2 index)",
    "Var(max depth)"
)

copy_number_stats <- c(
    "Bulk CN distance",
    "Mean(misseg. count in bulk)",
    "Var(misseg. count in bulk)",
    "Bulk CN distance",
    "Mean(misseg. count in bulk)",
    "Var(misseg. count in bulk)",
    #---Single-cell DNA: subclonal CN
    "Single-cell CN distance",
    "Mean(sc Shannon index)",
    "Mean(clonal misseg. count in sc)",
    "Mean(subclonal misseg. count in sc)",
    "Mean(clonal arm-misseg. count in sc)",
    "Mean(subclonal arm-misseg. count in sc)",
    "Var(sc Shannon index)",
    "Var(clonal misseg. count in sc)",
    "Var(subclonal misseg. count in sc)",
    "Var(clonal arm-misseg. count in sc)",
    "Var(subclonal arm-misseg. count in sc)",
    "Single-cell CN distance",
    "Mean(sc Shannon index)",
    "Mean(clonal misseg. count in sc)",
    "Mean(subclonal misseg. count in sc)",
    "Mean(clonal arm-misseg. count in sc)",
    "Mean(subclonal arm-misseg. count in sc)",
    "Var(sc Shannon index)",
    "Var(clonal misseg. count in sc)",
    "Var(subclonal misseg. count in sc)",
    "Var(clonal arm-misseg. count in sc)",
    "Var(subclonal arm-misseg. count in sc)"
)
# ---load ground-truth parameters
parameters_truth <- read.csv("parameters_ground_truth.csv", header = TRUE)
fitting_parameters(
    library_name = model_name,
    copynumber_DATA = DLP_stats,
    parameters_truth = parameters_truth,
    list_parameters = list_parameters,
    list_targets_by_parameter = list_targets,
    plot_ABC_prior_as_uniform = TRUE,
    stat_names = stat_names,
    phylo_tip_stats = phylo_tip_stats,
    phylo_balance_stats = phylo_balance_stats,
    copy_number_stats = copy_number_stats
)
# ===================PLOT CORRELATION BETWEEN INFERENCE AND GROUND TRUTH
# ===================================================FOR SELECTION RATES
parameters_inferred <- read.csv(paste0(model_name, "_para_output.csv"), header = TRUE)
parameters_inferred <- parameters_inferred[which(parameters_inferred$Type == "Selection_rate"), ]
plot_ABC_correlation(
    inference_result = parameters_inferred,
    plot_name = paste0(model_name, "_ABC_correlation.jpeg"),
    title_plot = "",
    value_x = "Ground_truth", title_x = "Ground truth",
    value_y = "Mean", title_y = "Posterior mean +/- std",
    error_y = "Sd",
    color_data = "red",
    plot_Error = TRUE,
    plot_diagonal = TRUE
)
# =============================SENSITIVITY ANALYSIS FOR SIMULATION COUNT
sensitivity_parameter <- "ABC_simcount"
sensitivity_title <- "Simulation count in ABC library"
# sensitivity_values <- c(1000, seq(10000, 100000, by = 10000))
sensitivity_library_statistics(
    library_name = model_name,
    library_sensitivity_name = paste0(model_name, "_simcount"),
    model_variables = model_variables,
    sensitivity_parameter = sensitivity_parameter,
    sensitivity_values = sensitivity_values,
    list_parameters = list_parameters,
    list_targets_library = list_targets_library,
    ABC_simcount_start = 0,
    ABC_simcount = ABC_simcount,
    cn_data_sc = ground_truth_cn_data_sc,
    cn_data_bulk = ground_truth_cn_data_bulk,
    arm_level = TRUE,
    cn_table = cn_table
)
sensitivity_fitting_and_plotting(
    library_name = model_name,
    library_sensitivity_name = paste0(model_name, "_simcount"),
    sensitivity_title = sensitivity_title,
    sensitivity_values = sensitivity_values,
    copynumber_DATA = DLP_stats,
    parameters_truth = parameters_truth,
    list_parameters = list_parameters,
    list_targets_by_parameter = list_targets,
    plot_ABC_prior_as_uniform = TRUE
)
# =============================PLOT FITTING RESULTS FOR DIFFERENT CHOICE
# ======================================================== OF STATISTICS
statistics_title <- "Choice of statistics"
statistics_values <- c("All", "Green", "Red", "Blue", "Red&Green", "Green&Blue", "Red&Blue")
statistics_IDs <- c("All", "CN", "Phylo balance", "Phylo tips", "CN + Phylo balance", "CN + Phylo tips", "Phylo tips + balance")
statistics_fitting_and_plotting(
    library_name = model_name,
    library_statistics_name = paste0(model_name, "_statistics"),
    statistics_title = statistics_title,
    statistics_values = statistics_values,
    statistics_IDs = statistics_IDs,
    Error_targets = c("CNA_probability", "Selection_rate"),
    Error_IDs = c("log10(prob_misseg)", "Selection rate"),
    copynumber_DATA = DLP_stats,
    parameters_truth = parameters_truth,
    list_parameters = list_parameters,
    list_targets_library = list_targets_library,
    plot_ABC_prior_as_uniform = TRUE,
    fontsize = 50,
    fitting = FALSE,
    # fitting = TRUE,
    plot_name = paste0(model_name, "_Choose_Statistics.jpeg")
)
