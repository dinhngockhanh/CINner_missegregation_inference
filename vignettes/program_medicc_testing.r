# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/zx2406/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/simulation/1026_medicc"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/DLPfit/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh&Zijin - Macmini
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

# setwd("/Users/dinhngockhanh/CancerSimulator/R")
# files_sources <- list.files(pattern = "*.r$")
# sapply(files_sources, source)

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
# ==============GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
cn_table <- model_variables$cn_info
cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
cn_table$Length <- cn_table$Bin_count * cn_bin_length
cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length
vec_CN_block_no <<- model_variables$cn_info$Bin_count
vec_centromeres <<- model_variables$cn_info$Centromere_location
# ======================================================MAKE SIMULATIONS
# simulator_full_program(
#     model = model_variables, model_prefix = model_name,
#     folder_workplace = folder_workplace,
#     n_simulations = n_simulations,
#     build_cn = TRUE,
#     save_cn_profile = TRUE, save_cn_clones = TRUE,
#     internal_nodes_cn_info = TRUE,
#     save_newick_tree = TRUE,
#     compute_parallel = TRUE,
#     R_libPaths = R_libPaths
# )
# =========================================================TRANSFER TSVS
# library(reshape2)
# library(tidyr)
# library(readr)
# for (i in 1:n_simulations) {
#     data <- read.csv(paste0(R_workplace, "/", model_name, "/", model_name, "_cn_profiles_long_", i, ".csv"))
#     new_df <- data.frame(matrix(nrow = nrow(data), ncol = 6))
#     colnames(new_df) <- c("sample_id", "chrom", "start", "end", "cn_a", "cn_b")
#     new_df$sample_id <- data$cell_id
#     new_df$chrom <- paste0("chr", data$chr)
#     new_df$start <- as.integer(data$start)
#     new_df$end <- as.integer(data$end)
#     new_df$cn_a <- data$Maj
#     new_df$cn_b <- data$Min
#     new_df <- new_df[grep("Library", new_df$sample_id), ]
#     file_path <- paste0(R_workplace, "/", model_name, "/Sample_", i, ".tsv")
#     write.table(new_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
# }

# =======================================COMPUTE GROUND-TRUTH STATISTICS
#---Get single-cell statistics & CN profiles
#   Get statistics & clonal CN profiles for each single-cell sample
list_targets_library <- c(
    #---Bulk DNA: CN
    # "data=bulk;target=genome;statistic=dist;variable=average_CN;metric=euclidean",
    # "data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    # "data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
    # paste0("data=bulk;target=chromosome;statistic=dist;variable=average_CN;metric=euclidean;chromosome=", list_chromosomes),
    # paste0("data=bulk;target=chromosome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    # paste0("data=bulk;target=chromosome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
    #---Single-cell DNA: subclonal CN
    # "data=sc;target=genome;statistic=dist;variable=clonal_CN;metric=euclidean",
    # "data=sc;target=genome;statistic=mean;variable=shannon",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    # "data=sc;target=genome;statistic=var;variable=shannon",
    "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
    "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
    # paste0("data=sc;target=chromosome;statistic=dist;variable=clonal_CN;metric=euclidean;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=mean;variable=shannon;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=var;variable=shannon;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
    # paste0("data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
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
list_targets_library_sc <- list_targets_library[grepl("data=sc", list_targets_library)]
cat(paste0("Loading ", n_simulations, " single-cell DNA-seq data sets...\n"))
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
ls_cn_sc_ground_truth <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(i) {
    load(paste0(model_name, "_simulation_", i, ".rda"))
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
for (simulation in 1:n_simulations) {
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


# DLP_stats <- get_statistics(
#     simulations_statistics_sc = ground_truth_statistics_sc,
#     # simulations_statistics_bulk = ground_truth_statistics_bulk,
#     list_targets = list_targets_library,
#     cn_data_sc = ground_truth_cn_data_sc,
#     # cn_data_bulk = ground_truth_cn_data_bulk,
#     arm_level = TRUE,
#     cn_table = cn_table
# )
# rbind(ls_cn_sc_ground_truth[[1]][[2]],(ls_cn_sc_ground_truth[[2]][[2]] )


### Add the simulation labels to the data frame
simulation_labels <- c()
for (i in 1:n_simulations) {
    simulation_labels[i] <- paste0("Simulation", i)
}
stats_comb <- data.frame(ls_cn_sc_ground_truth[[1]][[2]])
for (i in 2:n_simulations) {
    stats_comb <- rbind(stats_comb, data.frame(ls_cn_sc_ground_truth[[i]][[2]]))
}
stats_comb <- cbind(simulation = simulation_labels, stats_comb)
write.csv(stats_comb, "Statistics_simulation.csv")
# =============================================COMPUTE MEDICC STATISTICS
# ....
# ....
# ....
# ....
# ....
# ....
# ....
