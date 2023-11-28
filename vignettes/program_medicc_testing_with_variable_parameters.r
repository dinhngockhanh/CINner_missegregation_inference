# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/simulation/1026_medicc"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/DLPdata/vignettes"
# R_resultPaths <- "/Users/lexie/Documents/DNA/DLPdata/Results/"
# R_outputPaths <- "/Users/lexie/Library/CloudStorage/GoogleDrive-zl3213@columbia.edu/.shortcut-targets-by-id/10rccHeZeICEtbkkvEtKGaMZdk7yOe5PB/2023-10-26.Ground truth for Zhihan/output/"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh&Zijin - Macmini
R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/1128_Medicc_simulations"
R_libPaths <- ""
R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/DLPfit/R"
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
#   Number of simulations
n_simulations <- 1000
#   Number of cells sampled in each simulation
n_cells <- 500
#   Bounds for ground-truth log10(prob_CN_missegregation)
bound_ground_truth_prob_CN_missegregation_left <- -5
bound_ground_truth_prob_CN_missegregation_right <- -3
#   Bounds for ground-truth selection rates (1/r -> r)
bound_ground_truth_arm_s <- 1.15
#   Statistics for the study
list_targets_library <- c(
    #---Single-cell DNA: clonal/subclonal CNAs
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    #---Single-cell DNA: phylo stats for tips
    "data=sc;target=genome;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
    "data=sc;target=genome;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
    "data=sc;target=genome;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
    "data=sc;target=genome;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
    #---Single-cell DNA: phylo stats for balance
    "data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
    "data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
    "data=sc;target=genome;statistic=mean;variable=maxDepth" # height of phylogeny tree
)
# ===============================================GROUND TRUTH PARAMETERS
cell_lifespan <- 30
T_0 <- list(0, "year")
T_end <- list(80, "year")
Table_sample <- data.frame(Sample_ID = c("SA01"), Cell_count = c(1000), Age_sample = c(80))
selection_model <- "chrom-arm-selection"
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
# ==============GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
cn_table <- model_variables$cn_info
cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
cn_table$Length <- cn_table$Bin_count * cn_bin_length
cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length
vec_CN_block_no <<- model_variables$cn_info$Bin_count
vec_centromeres <<- model_variables$cn_info$Centromere_location
# ======================================================MAKE SIMULATIONS
# ===Function to make one simulation, parameters ~ uniform distributions
make_simulation <- function(index) {
    folder_workplace <- paste0(model_name, "_", index, "/")
    #---Sample model parameters from respective uniform distributions
    df_ground_truth <- data.frame(matrix(0, nrow = 0, ncol = 2))
    colnames(df_ground_truth) <- c("Parameter", "Value")
    prob_CN_missegregation <- 10^runif(1, bound_ground_truth_prob_CN_missegregation_left, bound_ground_truth_prob_CN_missegregation_right)
    model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "prob_CN_missegregation")] <- prob_CN_missegregation
    df_ground_truth[nrow(df_ground_truth) + 1, ] <- c("log10_misseg_rate", log10(prob_CN_missegregation))
    for (i in 1:nrow(model_variables$chromosome_arm_library)) {
        if (grepl("q$", model_variables$chromosome_arm_library$Arm_ID[i])) {
            model_variables$chromosome_arm_library$s_rate[i] <- 1
        } else if (grepl("p$", model_variables$chromosome_arm_library$Arm_ID[i])) {
            model_variables$chromosome_arm_library$s_rate[i] <- runif(1, 1, bound_ground_truth_arm_s)
            if (runif(1) < 0.5) model_variables$chromosome_arm_library$s_rate[i] <- 1 / model_variables$chromosome_arm_library$s_rate[i]
            df_ground_truth[nrow(df_ground_truth) + 1, ] <- c(paste0("selection_rate_", model_variables$chromosome_arm_library$Chromosome[i]), model_variables$chromosome_arm_library$s_rate[i])
        }
    }
    #---Produce one simulation
    simulator_full_program(
        model = model_variables, model_prefix = model_name,
        folder_workplace = folder_workplace,
        n_simulations = 1,
        # stage_final = 4,
        build_cn = TRUE,
        save_cn_profile = TRUE,
        save_cn_clones = TRUE,
        internal_nodes_cn_info = TRUE,
        save_newick_tree = TRUE,
        R_libPaths = R_libPaths
    )
    #---Save ground truth parameters
    write.csv(df_ground_truth, paste0(folder_workplace, "/", model_name, "_ground_truth_1.csv"), row.names = FALSE)
}
# ===Make simulations in parallel
simcount_start <- 1
simcount_end <- 1200



numCores <- detectCores()
cl <- makePSOCKcluster(numCores - 1)
if (is.null(R_libPaths) == FALSE) {
    R_libPaths <<- R_libPaths
    clusterExport(cl, varlist = c("R_libPaths"))
    clusterEvalQ(cl = cl, .libPaths(R_libPaths))
}
clusterEvalQ(cl = cl, library(CancerSimulator))
clusterExport(cl, varlist = c(
    "make_simulation",
    "model_name",
    "model_variables",
    "bound_ground_truth_prob_CN_missegregation_left",
    "bound_ground_truth_prob_CN_missegregation_right",
    "bound_ground_truth_arm_s"
))
pbo <- pboptions(type = "txt")
pblapply(cl = cl, X = simcount_start:simcount_end, FUN = function(iteration) {
    make_simulation(index = iteration)
})
stopCluster(cl)
