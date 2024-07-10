library(CINner)
library(CINner_missegregation_inference)
library(readxl)
library(parallel)
library(pbapply)
# ==================================================IMPORTANT PARAMETERS
#   Number of simulations
n_simulations <- 1200
#   Number of cells sampled in each simulation
n_cells <- 500
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
#---Probabilities of CNA
prob_CN_missegregation <- 2e-4
prob_neutral_CN_focal_amplification <- 1e-4
prob_neutral_CN_focal_deletion <- 1e-4

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
    prob_neutral_CN_focal_amplification = prob_neutral_CN_focal_amplification,
    prob_neutral_CN_focal_deletion = prob_neutral_CN_focal_deletion,
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
plot_clonal_phylo(
    model = model_name,
    n_simulations = n_simulations,
    folder_workplace = folder_workplace,
    folder_plots = folder_workplace,
    width = 2000,
    height = 2000,
    compute_parallel = TRUE
)
plot_cn_heatmap(
    model = model_name,
    n_simulations = n_simulations,
    folder_workplace = folder_workplace,
    folder_plots = folder_workplace,
    plotcol = "total-copy",
    CN_data = "TRUTH",
    phylo = "TRUTH",
    width = 1000,
    height = 1000,
    compute_parallel = TRUE
)
# =======================================TRANSFERING CSV TO MEDICC INPUTS
cat(paste0("Transferring ", n_simulations, " csvs to input format of medicc2...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
clusterExport(cl, varlist = c("read.csv", "write.table", "R_workplace", "as.integer", "grep"))
e <- new.env()
e$libs <- .libPaths()
clusterExport(cl, "libs", envir = e)
clusterEvalQ(cl, .libPaths(libs))
pbo <- pboptions(type = "txt")
mediccinput <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(i) {
    R_inputplace <- paste0(R_workplace, "/Medicc_testing_", i)
    data <- read.csv(paste0(R_inputplace, "/Medicc_testing_cn_profiles_long_1", ".csv"))
    new_df <- data.frame(matrix(nrow = nrow(data), ncol = 6))
    colnames(new_df) <- c("sample_id", "chrom", "start", "end", "cn_a", "cn_b")
    new_df$sample_id <- data$cell_id
    new_df$chrom <- paste0("chr", data$chr)
    new_df$start <- as.integer(data$start)
    new_df$end <- as.integer(data$end)
    new_df$cn_a <- data$Maj
    new_df$cn_b <- data$Min
    new_df <- new_df[grep("Library", new_df$sample_id), ]
    file_path <- paste0(R_workplace, "/Simulation", i, ".tsv")
    write.table(new_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
})
stopCluster(cl)
# ==========================TRANSFERING GROUNDTRUTH MISSEG RATE TO ONE DF
cat(paste0("Transferring ", n_simulations, " ground-truth missg rate to missg_rate_df...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
clusterExport(cl, varlist = c("read.csv", "write.table", "R_workplace"))
e <- new.env()
e$libs <- .libPaths()
clusterExport(cl, "libs", envir = e)
clusterEvalQ(cl, .libPaths(libs))
pbo <- pboptions(type = "txt")
misseg_rates <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(i) {
    R_inputplace <- paste0(R_workplace, "/Medicc_testing_", i)
    data <- read.csv(paste0(R_inputplace, "/Medicc_testing_ground_truth_1.csv"))
    misseg_rate <- data$Value[which(data$Parameter == "log10_misseg_rate")]
    return(misseg_rate)
})
stopCluster(cl)
misseg_rate_df <- data.frame(
    Simulation_id = paste0("Simulation_", 1:n_simulations),
    log10_misseg_rate = unlist(misseg_rates)
)
# Write the dataframe to a CSV file
write.csv(misseg_rate_df, file = paste0(R_workplace, "/misseg_rate_df.csv"), row.names = FALSE)
# =======================================COMPUTE GROUND-TRUTH STATISTICS
# ---Get single-cell statistics & CN profiles
# ---Get statistics & clonal CN profiles for each single-cell sample
list_targets_library_sc <- list_targets_library[grepl("data=sc", list_targets_library)]
cat(paste0("Loading ", n_simulations, " single-cell DNA-seq data sets...\n"))
n_cores <- max(detectCores() - 1, 1)
cl <- makePSOCKcluster(n_cores)
model_name <<- model_name
clusterExport(cl, varlist = c(
    "model_name", "R_workplace", "get_each_clonal_CN_profiles", "get_arm_CN_profiles",
    "cn_table", "get_each_statistics", "list_targets_library_sc", "find_clonal_ancestry", "find_event_count"
))
e <- new.env()
e$libs <- .libPaths()
clusterExport(cl, "libs", envir = e)
clusterEvalQ(cl, .libPaths(libs))
pbo <- pboptions(type = "txt")
ls_cn_sc_ground_truth <- pblapply(cl = cl, X = 1:n_simulations, FUN = function(i) {
    R_inputplace <- paste0(R_workplace, "/Medicc_testing_", i)
    load(paste0(R_inputplace, "/", model_name, "_simulation_1", ".rda"))
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
# ---Get statistics & clonal CN profiles for entire single-cell cohort
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
# ---Compute statistics
DLP_stats <- get_statistics(
    simulations_statistics_sc = ground_truth_statistics_sc,
    list_targets = list_targets_library,
    cn_data_sc = ground_truth_cn_data_sc,
    arm_level = TRUE,
    cn_table = cn_table
)
# ---Add the simulation labels to the data frame
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
# ====================================================RUN MEDICC COMMANDS
conda_init_script <- "/Users/xiangzijin/opt/anaconda3/etc/profile.d/conda.sh"
env_name <- "medicc_env"
input_file <- "/Users/xiangzijin/Downloads/medicc_workspace/Simulation1.tsv"
output_folder <- "/Users/xiangzijin/Downloads/medicc_workspace"
# ---Medicc2 command
medicc2_command <- paste(
    "medicc2",
    input_file,
    output_folder,
    "--events --plot none"
)
# ---Create a temporary script
temp_script <- tempfile(fileext = ".sh")
# ---Write the temporary script
writeLines(c(
    paste("source", conda_init_script),
    paste("conda activate", env_name),
    medicc2_command
), temp_script)
# ---Give the script permission
Sys.chmod(temp_script, "0755")
# ---Run the script
system2("bash", c("-c", temp_script))
# ---Clean the script
unlink(temp_script)
# ==============================================COMPUTE MEDICC STATISTICS
library(ape)
library(caper)
library(TreeDist)
# ---Set up list stats
list_targets_library <- c(
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
start_time <- Sys.time()
print(start_time)
##########
setwd(medicc_Paths)
true_stat <- read.csv("Statistics_simulation.csv", row.names = 1)
misseg <- read.csv("misseg_rate_df.csv", row.names = 1)
list_output <- list.dirs(path = medicc_outputPath, full.names = TRUE, recursive = FALSE)
dir_numbers <- as.numeric(gsub("^.*_(\\d+)$", "\\1", basename(list_output)))
list_output <- list_output[order(dir_numbers)]
n_library <- length(list_output)
combined_df <- data.frame()
### ================LOOP FOR SIMULATION
for (i in 1:n_library) {
    if (length(list.files(path = list_output[i])) == 0) {
        next
    }
    newick_file <- list.files(path = list_output[i], pattern = "*.new$", full.names = TRUE)
    medicc_phylogeny <- read.tree(file = newick_file)
    new_tree <- drop.tip(medicc_phylogeny, "diploid")
    groundtruth_newick <- paste0(medicc_Paths, "/Medicc_testing_", i, "/Medicc_testing_cell_phylogeny_1.newick")
    groundtruth_phylogeny <- read.tree(file = groundtruth_newick)
    ### ========== clonal
    events_file <- list.files(path = list_output[i], pattern = "*copynumber_events_df.tsv$", full.names = TRUE)
    events_df <- read.table(events_file, sep = "\t", header = TRUE)
    unique_ids <- unique(events_df$sample_id)
    results <- data.frame(sample_id = character(0), Ncells = numeric(0))
    Ntips <- Ntip(new_tree)
    match_indices <- match(unique_ids, c(new_tree$tip.label, new_tree$node.label))
    for (id in unique_ids) {
        match_number <- match(id, c(new_tree$tip.label, new_tree$node.label))
        descendents <- clade.members(match_number, new_tree, tip.labels = FALSE, include.nodes = FALSE)
        Ncells <- length(descendents)
        results <- rbind(results, data.frame(sample_id = id, Ncells = Ncells))
    }
    merged_data <- merge(events_df, results, by = "sample_id")
    event_count.clonal <- sum(merged_data$Ncells == Ntips)
    event_count.subclonal <- sum(merged_data$Ncells[merged_data$Ncells != 1 & merged_data$Ncells != Ntips]) / Ntips
    ### ========== phylogeny Stat
    medicc_stat <- statistics(sample_phylogeny = new_tree, sample_clone_assignment = None, list_targets = list_targets_library)
    medicc_df <- as.data.frame(medicc_stat)
    medicc_df["data.sc.target.genome.variable.event_count.type.clonal.event.missegregation"] <- event_count.clonal
    medicc_df["data.sc.target.genome.variable.event_count.type.subclonal.event.missegregation"] <- event_count.subclonal
    medicc_df["TreeDistance"] <- TreeDistance(groundtruth_phylogeny, new_tree)
    medicc_df["Simulation"] <- paste0("Simulation", i)
    medicc_df["log10_misseg_rate"] <- misseg[i, 2]
    if (i == 1) {
        combined_df <- medicc_df
    } else {
        combined_df <- rbind(combined_df, medicc_df)
    }
    if (i %% 100 == 0) {
        print(paste0("simulation", i, " finished."))
    }
}
print("Finished.")
# # =====================================================PLOTTING
plot_stats_corr(true_stat, combined_df)
plot_RF(combined_df)
end_time <- Sys.time()
print(end_time)
