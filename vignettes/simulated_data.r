# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
R_workplace <- getwd()
R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/simulation/Fitting_experiment_1point1_8500"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh&Zijin - Macmini
# R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/experiment"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/DLPfit/testR"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/test/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/DLPfit/vignettes"
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
# ==================================================IMPORTANT PARAMETERS
#   Number of single-cell samples in ground-truth data & ABC simulations
N_data_sc <- 50
#   Number of bulk samples in ground-truth data & ABC simulations
N_data_bulk <- 100
#   Bounds for ground-truth selection rates (1/r -> r)
bound_ground_truth_arm_s <- 1.1
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
prob_CN_chrom_arm_missegregation <- 0e-4
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
    prob_CN_chrom_arm_missegregation = prob_CN_chrom_arm_missegregation,
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
        arm_s[i] <- runif(1, 1, 1.1)
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
# model_name <- "Chromosome_missegregation"
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
            model_variables$chromosome_arm_library$Arm_ID[i], paste0("Selection rate - chromosome ", model_variables$chromosome_arm_library$Chromosome[i]), model_variables$chromosome_arm_library$Chromosome[i], "Arm_selection_rate",
            1 / bound_ABC_arm_s, bound_ABC_arm_s
        )
    }
}
list_chromosomes <- paste(unique(model_variables$chromosome_arm_library$Chromosome), collapse = ",")
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
sim.sc <- simulator_full_program(
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
sim.bulk <- simulator_full_program(
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
