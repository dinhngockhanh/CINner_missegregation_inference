# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
R_workplace <- getwd()
R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/simulation/DLP experiment_ch1&2"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/test/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/DLPfit/test_functions"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/DLPfit/R"



# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
setwd(R_libPaths_extra)
getwd()
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)
setwd(R_workplace)
library(readxl)
library(CancerSimulator)
library(parallel)
library(pbapply)
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
prob_CN_chrom_arm_missegregation <- 2e-4
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



selected_chromosomes <- c("1", "2", "3", "4", "5")
arm_s[which(arm_chromosome %in% selected_chromosomes)] <-
    runif(length(which(arm_chromosome %in% selected_chromosomes)), 1, 1.5)
for (i in 1:length(arm_s)) {
    if (runif(1) < 0.5) arm_s[i] <- 1 / arm_s[i]
}



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
model_name <- "Simpler_DLP_5chr"
model_variables <- CHECK_model_variables(model_variables)
SAVE_model_variables(
    model_name = model_name,
    model_variables = model_variables
)
# =======================================================MAKE "DLP" DATA
####
####
####
####
####
# N_data <- 3
N_data <- 10
####
####
####
####
####
# cat(paste0("\n\n\nMaking ", N_data, " simulations...\n"))
# tmp <- simulator_full_program(
#     model = model_name,
#     n_simulations = N_data,
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
# ======================================DEFINE LIST OF PARAMETERS TO FIT
list_parameters <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(list_parameters) <- c("Variable", "Type", "Lower_bound", "Upper_bound")
# list_parameters[nrow(list_parameters) + 1, ] <- c(
#     "prob_CN_missegregation",
#     "CNA_probability",
#     1e-4, 5e-4
# )
# list_parameters[nrow(list_parameters) + 1, ] <- c(
#     "prob_CN_chrom_arm_missegregation",
#     "CNA_probability",
#     1e-4, 5e-4
# )
list_parameters[nrow(list_parameters) + 1, ] <- c(
    "10^:prob_CN_missegregation",
    "CNA_probability",
    -5, -3
)
list_parameters[nrow(list_parameters) + 1, ] <- c(
    "10^:prob_CN_chrom_arm_missegregation",
    "CNA_probability",
    -5, -3
)
# for (i in 1:nrow(model_variables$chromosome_arm_library)) {
#     if (model_variables$chromosome_arm_library$Chromosome[i] %in% selected_chromosomes) {
#         list_parameters[nrow(list_parameters) + 1, ] <- c(
#             model_variables$chromosome_arm_library$Arm_ID[i],
#             "Arm_selection_rate",
#             1 / 1.5, 1.5
#         )
#     }
# }
# =====================================PRINT OUT GROUND TRUTH PARAMETERS
list_parameters_ground_truth <- list_parameters
list_parameters_ground_truth$Value <- 0
for (row in 1:nrow(list_parameters)) {
    parameter_ID <- list_parameters$Variable[row]
    if (parameter_ID %in% model_variables$general_variables$Variable) {
        list_parameters_ground_truth$Value[row] <- model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)]
    } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
        list_parameters_ground_truth$Value[row] <- model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)]
    }
}
write.csv(list_parameters_ground_truth, "parameters_ground_truth.csv")
# =====================DEFINE STATISTICS FOR BUILDING SIMULATION LIBRARY
list_targets_library <- c(
    "statistic=mean;variable=shannon",
    "statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "statistic=var;variable=shannon",
    "statistic=var;variable=event_count;type=clonal;event=missegregation",
    "statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "statistic=dist;variable=clonal_CN;metric=euclidean",
    #---phylo stats with tips
    "statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    "statistic=var;variable=cherries",
    "statistic=var;variable=pitchforks",
    "statistic=var;variable=IL_number",
    "statistic=var;variable=avgLadder",
    #---phylo stats for balance
    "statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "statistic=mean;variable=colless", # balance index of phylogeny tree
    "statistic=mean;variable=sackin", # balance index of phylogeny tree
    "statistic=mean;variable=B2", # balance index of phylogeny tree
    "statistic=mean;variable=maxDepth", # height of phylogeny tree
    "statistic=var;variable=stairs",
    "statistic=var;variable=colless",
    "statistic=var;variable=sackin",
    "statistic=var;variable=B2",
    "statistic=var;variable=maxDepth"
)
# ==============GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
cn_table <- model_variables$cn_info
cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
cn_table$Length <- cn_table$Bin_count * cn_bin_length
cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length
# ===================================INPUT GROUND TRUTH DATA FOR FITTING
vec_CN_block_no <<- model_variables$cn_info$Bin_count
vec_centromeres <<- model_variables$cn_info$Centromere_location

cat(paste0("Loading ", N_data, " data sets...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
model_name <<- model_name
clusterExport(cl, varlist = c("model_name"))
pbo <- pboptions(type = "txt")
cn_ground_truth <- pblapply(cl = cl, X = 1:N_data, FUN = function(i) {
    load(paste0(model_name, "_simulation_", i, ".rda"))
    return(simulation)
})
stopCluster(cl)
# ======================================GET CLONAL CN PROFILES FROM DATA
data_clonal_CN_profiles <- get_clonal_CN_profiles(
    cn_ground_truth,
    arm_level = TRUE,
    cn_table = cn_table
)
# =======================================FIT PARAMETERS USING "DLP" DATA
#---Produce library of simulations for fitting
n_simulations <- N_data
library_sc_CN(
    model_name = model_name,
    model_variables = model_variables,
    list_parameters = list_parameters,
    list_targets_library = list_targets_library,
    # list_targets = list_targets,
    ####
    ####
    ####
    ####
    ####
    cn_table = cn_table,
    # ABC_simcount = 2,
    arm_level = TRUE,
    ABC_simcount = 10000,
    ####
    ####
    ####
    ####
    ####
    n_simulations = n_simulations,
    library_name = model_name,
    cn_data = data_clonal_CN_profiles,
    save_sample_statistics = TRUE
)
#---Import ground truth parameters
parameters_truth <- read.csv("parameters_ground_truth.csv", header = TRUE)
#---Get statistics from ground truth
# DLP_stats <- get_statistics(
#     simulations = cn_ground_truth,
#     list_targets = list_targets_library,
#     cn_data = data_clonal_CN_profiles,
#     arm_level = TRUE,
#     cn_table = cn_table,
#     save_sample_statistics = TRUE
# )
#---Fit parameters and compare with ground truth
list_targets <- c(
    "statistic=mean;variable=shannon",
    "statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    "statistic=dist;variable=clonal_CN;metric=euclidean",
    #---phylo stats with tips
    "statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    #---phylo stats for balance
    # "statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "statistic=mean;variable=colless" # balance index of phylogeny tree
    # "statistic=mean;variable=sackin" # balance index of phylogeny tree
    # "statistic=mean;variable=B2", # balance index of phylogeny tree
    # "statistic=mean;variable=maxDepth" # height of phylogeny tree
)
# fitting_sc_CN(
#     library_name = model_name,
#     model_name = model_name,
#     copynumber_DATA = DLP_stats,
#     parameters_truth = parameters_truth,
#     list_parameters = list_parameters,
#     list_targets_library = list_targets_library,
#     list_targets = list_targets,
#     shuffle_num = 200,
#     cn_data = data_clonal_CN_profiles,
#     arm_level = TRUE,
#     cn_table = cn_table,
#     shuffle_chromosome_arms = FALSE,
#     shuffle_chromosomes = TRUE
# )
