# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/simulation/Output_experiment"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh&Zijin - Macmini
# R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/experiment"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/DLPfit/testR"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/knd2127/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/knd2127/test/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/DLPfit/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/DLPfit/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(readxl)
# library(CancerSimulator)
library(parallel)
library(pbapply)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)

setwd("/Users/dinhngockhanh/CancerSimulator/R")
files_sources <- list.files(pattern = "*.r$")
sapply(files_sources, source)

setwd(R_workplace)
# devtools::install_github("dinhngockhanh/CancerSimulator", force = TRUE)
# ==================================================IMPORTANT PARAMETERS
#   Number of simulations
n_simulations <- 10
#   Number of cells sampled in each simulation
n_cells <- 500
#   Bounds for ground-truth selection rates (1/r -> r)
bound_ground_truth_arm_s <- 1.15
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
model_name <- "Medicc_testing"
model_variables <- CHECK_model_variables(model_variables)
folder_workplace <- paste0(model_name, "/")
# ======================================================MAKE SIMULATIONS
simulator_full_program(
    model = model_variables, model_prefix = model_name,
    folder_workplace = folder_workplace,
    n_simulations = n_simulations,
    build_cn = TRUE,
    save_cn_profile = TRUE, save_cn_clones = TRUE,
    internal_nodes_cn_info = TRUE,
    save_newick_tree = TRUE,
    compute_parallel = TRUE,
    R_libPaths = R_libPaths
)
# =======================================COMPUTE GROUND-TRUTH STATISTICS
# ....
# ....
# ....
# ....
# ....
# ....
# ....
# =============================================COMPUTE MEDICC STATISTICS
# ....
# ....
# ....
# ....
# ....
# ....
# ....
