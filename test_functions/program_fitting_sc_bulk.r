# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
R_workplace <- getwd()
R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/simulation/DLP experiment_ch1&2"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/test/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/DLPfit/test_functions"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/DLPfit/R"



# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(readxl)
library(CancerSimulator)
library(parallel)
library(pbapply)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)
setwd(R_workplace)
# devtools::install_github("dinhngockhanh/CancerSimulator", force = TRUE)
# ======================COUNT OF SAMPLES IN BULK AND SINGLE-CELL COHORTS
N_data_dlp <- 10
N_data_bulk <- 100
# ===============================================GROUND TRUTH PARAMETERS
cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(
    Sample_ID = c("SA01"),
    Cell_count = c(1000),
    Age_sample = c(80)
)
selection_model <- "chrom-arm-selection"
CN_bin_length <- 500000
#---Probabilities of CNA
prob_CN_missegregation <- 2e-4
prob_CN_chrom_arm_missegregation <- 0
#---Viability thresholds
bound_driver <- 3
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
    CN_bin_length = CN_bin_length,
    Table_sample = Table_sample,
    prob_CN_missegregation = prob_CN_missegregation,
    prob_CN_chrom_arm_missegregation = prob_CN_chrom_arm_missegregation,
    selection_model = selection_model,
    bound_driver = bound_driver,
    bound_average_ploidy = bound_average_ploidy,
    bound_homozygosity = bound_homozygosity,
    table_population_dynamics = table_population_dynamics
)
#---Set up (randomized) chromosome arm selection rates
arm_id <- c(
    paste(model_variables$cn_info$Chromosome, "p", sep = ""),
    paste(model_variables$cn_info$Chromosome, "q", sep = "")
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
for (i in 1:length(arm_s)) {
    if (grepl("q$", arm_id[i])) {
        arm_s[i] <- 1
    }
    if (grepl("p$", arm_id[i])) {
        set.seed(i)
        arm_s[i] <- runif(1, 1, 1.2)
        set.seed(i)
        if (runif(1) < 0.5) arm_s[i] <- 1 / arm_s[i]
    }
}
selected_chromosomes <- c(
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
    "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"
)

table_arm_selection_rates <- data.frame(
    Arm_ID = arm_id,
    Chromosome = arm_chromosome,
    Bin_start = arm_start,
    Bin_end = arm_end,
    s_rate = arm_s
)
model_variables <- BUILD_driver_library(
    model_variables = model_variables,
    table_arm_selection_rates = table_arm_selection_rates,
)
#---Set up initial cell population
cell_count <- 20
CN_matrix <- BUILD_cn_normal_XX(model_variables$cn_info)
drivers <- list()
model_variables <- BUILD_initial_population(
    model_variables = model_variables,
    cell_count = cell_count,
    CN_matrix = CN_matrix,
    drivers = drivers
)
#---Save model variables
model_name <- "Simpler_DLP&BULK_DNA"
model_variables <- CHECK_model_variables(model_variables)
SAVE_model_variables(
    model_name = paste0(model_name, "_sc"),
    model_variables = model_variables
)
SAVE_model_variables(
    model_name = paste0(model_name, "_bulk"),
    model_variables = model_variables
)
# ======================================DEFINE LIST OF PARAMETERS TO FIT
list_parameters <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(list_parameters) <- c("Variable", "Type", "Lower_bound", "Upper_bound")
list_parameters[nrow(list_parameters) + 1, ] <- c(
    "10^:prob_CN_missegregation",
    "CNA_probability",
    -5, -3
)
for (i in 1:nrow(model_variables$chromosome_arm_library)) {
    if (grepl("p$", model_variables$chromosome_arm_library$Arm_ID[i])) {
        list_parameters[nrow(list_parameters) + 1, ] <- c(
            model_variables$chromosome_arm_library$Arm_ID[i],
            "Arm_selection_rate",
            1 / 1.2, 1.2
        )
    }
}
# =============DEFINE LIST OF STATISTICS FOR BUILDING SIMULATION LIBRARY
list_targets_library <- c(
    #---Bulk DNA: CN
    "data=bulk;statistic=dist;variable=average_CN;metric=euclidean",
    "data=bulk;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    #---Single-cell DNA: subclonal CN
    "data=sc;statistic=mean;variable=shannon",
    "data=sc;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "data=sc;statistic=var;variable=shannon",
    "data=sc;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "data=sc;statistic=dist;variable=clonal_CN;metric=euclidean",
    #---Single-cell DNA: phylo stats for tips
    "data=sc;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;statistic=var;variable=cherries",
    "data=sc;statistic=var;variable=pitchforks",
    "data=sc;statistic=var;variable=IL_number",
    "data=sc;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;statistic=var;variable=stairs",
    "data=sc;statistic=var;variable=colless",
    "data=sc;statistic=var;variable=sackin",
    "data=sc;statistic=var;variable=B2",
    "data=sc;statistic=var;variable=maxDepth"
)
#---Get statistics for single-cell data and simulations
list_targets_library_sc <- list_targets_library[grepl("data=sc", list_targets_library)]
#---Get statistics for bulk data and simulations
list_targets_library_bulk <- list_targets_library[grepl("data=bulk", list_targets_library)]
# =====================================PRINT OUT GROUND-TRUTH PARAMETERS
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
#---Import ground truth parameters
parameters_truth <- read.csv("parameters_ground_truth.csv", header = TRUE)
# ==============GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
cn_table <- model_variables$cn_info
cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
cn_table$Length <- cn_table$Bin_count * cn_bin_length
cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length
vec_CN_block_no <<- model_variables$cn_info$Bin_count
vec_centromeres <<- model_variables$cn_info$Centromere_location
# ================================================MAKE GROUND-TRUTH DATA
#---Make single-cell ground-truth simulations
# cat(paste0("\n\n\nMaking ", N_data_dlp, " single-cell simulations...\n"))
# tmp <- simulator_full_program(
#     model = paste0(model_name, "_sc"),
#     n_simulations = N_data_dlp,
#     stage_final = 3,
#     compute_parallel = TRUE,
#     output_variables = c(
#         "evolution_origin",
#         "evolution_genotype_changes",
#         "sample_clone_ID",
#         "sample_genotype_unique",
#         "sample_genotype_unique_profile",
#         "phylogeny_clustering_truth"
#     ),
#     R_libPaths = R_libPaths
# )
# tmp <- c()
# #---Make bulk ground-truth simulations
# cat(paste0("\n\n\nMaking ", N_data_bulk, " bulk simulations...\n"))
# tmp <- simulator_full_program(
#     model = paste0(model_name, "_bulk"),
#     n_simulations = N_data_bulk,
#     stage_final = 2,
#     compute_parallel = TRUE,
#     output_variables = c(
#         # "evolution_origin",
#         # "evolution_genotype_changes",
#         "sample_genotype_unique_profile",
#         "sample_genotype_unique",
#         "sample_clone_ID"
#     ),
#     R_libPaths = R_libPaths
# )
# tmp <- c()
# ============GET STATISTICS & CN PROFILES FROM GROUND-TRUTH SIMULATIONS
#---Get single-cell statistics & CN profiles
#   Get statistics & clonal CN profiles for each single-cell sample
cat(paste0("Loading ", N_data_dlp, " single-cell DNA-seq data sets...\n"))
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
ls_cn_sc_ground_truth <- pblapply(cl = cl, X = 1:N_data_dlp, FUN = function(i) {
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
ls_cn_sc_ground_truth_all <- vector("list", length = 2)
for (type in 1:2) {
    #   type = 1 <-> clonal CN profiles; type = 2 <-> statistics
    for (statistic in 1:length(ls_cn_sc_ground_truth[[1]][[type]])) {
        for (simulation in 1:N_data_dlp) {
            if (simulation == 1) {
                ls_cn_sc_ground_truth_all[[type]][[statistic]] <- ls_cn_sc_ground_truth[[simulation]][[type]][[statistic]][1]
            } else {
                ls_cn_sc_ground_truth_all[[type]][[statistic]] <- c(ls_cn_sc_ground_truth_all[[type]][[statistic]], ls_cn_sc_ground_truth[[simulation]][[type]][[statistic]][1])
            }
        }
    }
    names(ls_cn_sc_ground_truth_all[[type]]) <- names(ls_cn_sc_ground_truth[[1]][[type]])
}
#---Get bulk statistics & CN profiles
#   Get statistics & representative CN profiles for each bulk sample
cat(paste0("Loading ", N_data_bulk, " bulk DNA-seq data sets...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
model_name <<- model_name
clusterExport(cl, varlist = c(
    "model_name", "N_data_dlp", "get_each_clonal_CN_profiles", "get_arm_CN_profiles",
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
stopCluster(cl)
#   Get statistics & representative CN profiles for entire bulk cohort
ls_cn_bulk_ground_truth_all <- vector("list", length = 2)
for (type in 1:2) {
    #   type = 1 <-> representative CN profiles; type = 2 <-> statistics
    for (statistic in 1:length(ls_cn_bulk_ground_truth[[1]][[type]])) {
        for (simulation in 1:N_data_bulk) {
            if (simulation == 1) {
                ls_cn_bulk_ground_truth_all[[type]][[statistic]] <- ls_cn_bulk_ground_truth[[simulation]][[type]][[statistic]][1]
            } else {
                ls_cn_bulk_ground_truth_all[[type]][[statistic]] <- c(ls_cn_bulk_ground_truth_all[[type]][[statistic]], ls_cn_bulk_ground_truth[[simulation]][[type]][[statistic]][1])
            }
        }
    }
    names(ls_cn_bulk_ground_truth_all[[type]]) <- names(ls_cn_bulk_ground_truth[[1]][[type]])
}
# ===============================================MAKE SIMULATION LIBRARY
# library_sc_CN(
#     model_name = model_name,
#     model_variables = model_variables,
#     list_parameters = list_parameters,
#     list_targets_library = list_targets_library,
#     ####
#     ####
#     ####
#     ####
#     ####
#     ABC_simcount = 500,
#     arm_level = TRUE,
#     cn_table = cn_table,
#     cn_data_sc = ls_cn_sc_ground_truth_all[[1]],
#     cn_data_bulk = ls_cn_bulk_ground_truth_all[[1]],
#     n_simulations_sc = N_data_dlp,
#     n_simulations_bulk = N_data_bulk,
#     ####
#     ####
#     ####
#     ####
#     ####
#     library_name = model_name,
#     save_sample_statistics = TRUE
# )
# ==================DEFINE LIST OF STATISTICS FOR FITTING EACH PARAMETER
list_targets <- data.frame(matrix(0, ncol = (length(list_targets_library) + 1), nrow = length(list_parameters$Variable)))
colnames(list_targets) <- c("Variable", list_targets_library)
list_targets[, 1] <- list_parameters$Variable
#---Statistics for fitting misseg rate
list_targets_misseg <- c(
    "data=bulk;statistic=dist;variable=average_CN;metric=euclidean",
    "data=bulk;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    # "data=sc;statistic=mean;variable=shannon",
    # "data=sc;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    # "data=sc;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    # "data=sc;statistic=var;variable=shannon",
    # "data=sc;statistic=var;variable=event_count;type=clonal;event=missegregation",
    # "data=sc;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    # "data=sc;statistic=dist;variable=clonal_CN;metric=euclidean",
    # "data=sc;statistic=mean;variable=colless",
    # "data=sc;statistic=mean;variable=sackin",
    # "data=sc;statistic=mean;variable=B2",
    # "data=sc;statistic=mean;variable=maxDepth"
    # "data=sc;statistic=var;variable=B2"
    "data=sc;statistic=mean;variable=shannon",
    "data=sc;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=var;variable=shannon",
    "data=sc;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=dist;variable=clonal_CN;metric=euclidean",
    #---Single-cell DNA: phylo stats for tips
    "data=sc;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;statistic=var;variable=cherries",
    "data=sc;statistic=var;variable=pitchforks",
    "data=sc;statistic=var;variable=IL_number",
    "data=sc;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;statistic=var;variable=stairs",
    "data=sc;statistic=var;variable=colless",
    "data=sc;statistic=var;variable=sackin",
    "data=sc;statistic=var;variable=B2",
    "data=sc;statistic=var;variable=maxDepth"
)
list_targets[1, which(colnames(list_targets) %in% list_targets_misseg)] <- 1
#---Statistics for fitting selection rate
list_targets_selection <- c(
    # "data=bulk;statistic=dist;variable=average_CN;metric=euclidean",
    # "data=sc;statistic=dist;variable=clonal_CN;metric=euclidean",
    # "data=sc;statistic=var;variable=shannon",
    # "data=sc;statistic=mean;variable=shannon",
    # "data=sc;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    # "data=sc;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    # "data=bulk;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    # "data=sc;statistic=mean;variable=colless",
    # "data=sc;statistic=mean;variable=sackin",
    # # "data=sc;statistic=mean;variable=B2",
    # # "data=sc;statistic=var;variable=B2",
    # "data=sc;statistic=mean;variable=maxDepth",
    # "data=sc;statistic=var;variable=event_count;type=clonal;event=missegregation"
    "data=bulk;statistic=dist;variable=average_CN;metric=euclidean",
    "data=bulk;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    #---Single-cell DNA: subclonal CN
    "data=sc;statistic=mean;variable=shannon",
    "data=sc;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=var;variable=shannon",
    "data=sc;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;statistic=dist;variable=clonal_CN;metric=euclidean",
    #---Single-cell DNA: phylo stats for tips
    "data=sc;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "data=sc;statistic=var;variable=cherries",
    "data=sc;statistic=var;variable=pitchforks",
    "data=sc;statistic=var;variable=IL_number",
    "data=sc;statistic=var;variable=avgLadder",
    #---Single-cell DNA: phylo stats for balance
    "data=sc;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;statistic=mean;variable=maxDepth", # height of phylogeny tree
    "data=sc;statistic=var;variable=stairs",
    "data=sc;statistic=var;variable=colless",
    "data=sc;statistic=var;variable=sackin",
    "data=sc;statistic=var;variable=B2",
    "data=sc;statistic=var;variable=maxDepth"
)
for (row in 2:nrow(list_targets)) {
    list_targets[row, which(colnames(list_targets) %in% list_targets_selection)] <- 1
}
# =========================GET FITTING STATISTICS FROM GROUND-TRUTH DATA
DLP_stats <- get_statistics(
    simulations_statistics_sc = ls_cn_sc_ground_truth_all[[2]],
    simulations_statistics_bulk = ls_cn_bulk_ground_truth_all[[2]],
    list_targets = list_targets_library,
    cn_data_sc = ls_cn_sc_ground_truth_all[[1]],
    cn_data_bulk = ls_cn_bulk_ground_truth_all[[1]],
    arm_level = TRUE,
    cn_table = cn_table,
    save_sample_statistics = TRUE
)
# ==============================================FIT PARAMETERS USING ABC
fitting_sc_CN(
    library_name = model_name,
    model_name = model_name,
    copynumber_DATA = DLP_stats,
    parameters_truth = parameters_truth,
    list_parameters = list_parameters,
    list_targets = list_targets,
    shuffle_num = 3,
    cn_data_sc = ls_cn_sc_ground_truth_all[[1]],
    cn_data_bulk = ls_cn_bulk_ground_truth_all[[1]],
    arm_level = TRUE,
    cn_table = cn_table,
    shuffle_chromosome_arms = FALSE,
    shuffle_chromosomes_by_permutation = FALSE,
    shuffle_chromosomes_by_moving = FALSE,
    currentpath
)
