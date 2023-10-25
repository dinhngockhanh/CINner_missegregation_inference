#----------------------------------------------------Load the rda
Paths <- "/Users/xiangzijin/Documents/simulation/Fitting_experiment_1point15_100000"
setwd(Paths)
rda_name <- "Fitting_whole_chroms_ABC_input.rda"
load(file = rda_name)
# Paths <- "/Users/khanhngocdinh/Documents/Zijin/experiment"
# setwd(Paths)
# rda_name <- "Simpler_DLP&BULK_DNA_ABC_input_10000.rda"
# load(file = rda_name)
#---------------------------------------------------Get parameters
parameters <- ABC_input$sim_param
# param_names <- list_parameters$Title
param_names <- c(
  "log10(prob_misseg)",
  "Sel.rate(Chrom 1)",
  "Sel.rate(Chrom 2)",
  "Sel.rate(Chrom 3)",
  "Sel.rate(Chrom 4)",
  "Sel.rate(Chrom 5)",
  "Sel.rate(Chrom 6)",
  "Sel.rate(Chrom 7)",
  "Sel.rate(Chrom 8)",
  "Sel.rate(Chrom 9)",
  "Sel.rate(Chrom 10)",
  "Sel.rate(Chrom 11)",
  "Sel.rate(Chrom 12)",
  "Sel.rate(Chrom 13)",
  "Sel.rate(Chrom 14)",
  "Sel.rate(Chrom 15)",
  "Sel.rate(Chrom 16)",
  "Sel.rate(Chrom 17)",
  "Sel.rate(Chrom 18)",
  "Sel.rate(Chrom 19)",
  "Sel.rate(Chrom 20)",
  "Sel.rate(Chrom 21)",
  "Sel.rate(Chrom 22)"
)
# get rid of last column of parameters
colnames(parameters) <- param_names
# prob_armmisseg <- parameters[, "prob_armmisseg"]
prob_misseg <- parameters[, 1]
#---------------------------------------------------Get statistics
correct_misseg_stat_names <- c(
  #-----bulk
  "Wasserstein_dist_bulk_genome",
  "Wasserstein_dist_sc_genome",
  "Mean_Shannon_genome",
  "Mean_Clonal_misseg_count_bulk_genome",
  "Mean_Clonal_misseg_count_sc_genome",
  "Mean_Subclonal_misseg_count_sc_genome",
  "Mean_Clonal_armmisseg_count_sc_genome",
  "Mean_Subclonal_armmisseg_count_sc_genome",
  #-----sc_subclonal CN
  "Var_Shannon_genome",
  "Var_Clonal_misseg_count_bulk_genome",
  "Var_Clonal_misseg_count_sc_genome",
  "Var_Subclonal_misseg_count_sc_genome",
  "Var_Clonal_armmisseg_count_sc_genome",
  "Var_Subclonal_armmisseg_count_sc_genome",
  #-----sc_Phylo Stats
  "Mean_Cherries_genome",
  "Mean_Pitchforks_genome",
  "Mean_IL_number_genome",
  "Mean_AvgLadder_genome",
  "Var_Cherries_genome",
  "Var_Pitchforks_genome",
  "Var_IL_number_genome",
  "Var_AvgLadder_genome",
  # balance
  "Mean_Stairs_genome",
  "Mean_Colless_genome",
  "Mean_Sackin_genome",
  "Mean_B2_genome",
  "Mean_MaxDepth_genome",
  "Var_Stairs_genome",
  "Var_Colless_genome",
  "Var_Sackin_genome",
  "Var_B2_genome",
  "Var_MaxDepth_genome"
)
misseg_stat_names <- c(
  #-----bulk
  "Wasserstein_dist_bulk_genome",
  "Mean_Clonal_misseg_count_bulk_genome",
  "Var_Clonal_misseg_count_bulk_genome",
  "Wasserstein_dist_bulk_chr",
  "Mean_Clonal_misseg_count_bulk_chr",
  "Var_Clonal_misseg_count_bulk_chr",
  #-----sc_subclonal CN
  "Wasserstein_dist_sc_genome",
  "Mean_Shannon_genome",
  "Mean_Clonal_misseg_count_sc_genome",
  "Mean_Subclonal_misseg_count_sc_genome",
  "Mean_Clonal_armmisseg_count_sc_genome",
  "Mean_Subclonal_armmisseg_count_sc_genome",
  "Var_Shannon_genome",
  "Var_Clonal_misseg_count_sc_genome",
  "Var_Subclonal_misseg_count_sc_genome",
  "Var_Clonal_armmisseg_count_sc_genome",
  "Var_Subclonal_armmisseg_count_sc_genome",
  "Wasserstein_dist_sc_chr",
  "Mean_Shannon_chr",
  "Mean_Clonal_misseg_count_sc_chr",
  "Mean_Subclonal_misseg_count_sc_chr",
  "Mean_Clonal_armmisseg_count_sc_chr",
  "Mean_Subclonal_armmisseg_count_sc_chr",
  "Var_Shannon_chr",
  "Var_Clonal_misseg_count_sc_chr",
  "Var_Subclonal_misseg_count_sc_chr",
  "Var_Clonal_armmisseg_count_sc_chr",
  "Var_Subclonal_armmisseg_count_sc_chr",
  #-----sc_Phylo Stats
  "Mean_Cherries_genome",
  "Mean_Pitchforks_genome",
  "Mean_IL_number_genome",
  "Mean_AvgLadder_genome",
  "Var_Cherries_genome",
  "Var_Pitchforks_genome",
  "Var_IL_number_genome",
  "Var_AvgLadder_genome",
  # balance
  "Mean_Stairs_genome",
  "Mean_Colless_genome",
  "Mean_Sackin_genome",
  "Mean_B2_genome",
  "Mean_MaxDepth_genome",
  "Var_Stairs_genome",
  "Var_Colless_genome",
  "Var_Sackin_genome",
  "Var_B2_genome",
  "Var_MaxDepth_genome"
)

sel_stat_names <- c(
  #-----bulk
  "Wasserstein_dist_bulk_chr",
  "Wasserstein_dist_sc_chr",
  "Mean_Shannon_chr",
  "Mean_Clonal_misseg_count_bulk_chr",
  #-----sc_subclonal CN
  "Mean_Clonal_misseg_count_sc_chr",
  "Mean_Subclonal_misseg_count_sc_chr",
  "Mean_Clonal_armmisseg_count_sc_chr",
  "Mean_Subclonal_armmisseg_count_sc_chr",
  "Var_Shannon_chr",
  "Var_Clonal_misseg_count_bulk_chr",
  "Var_Clonal_misseg_count_sc_chr",
  "Var_Subclonal_misseg_count_sc_chr",
  "Var_Clonal_armmisseg_count_sc_chr",
  "Var_Subclonal_armmisseg_count_sc_chr",
  #-----sc_Phylo Stats
  "Mean_Cherries_genome",
  "Mean_Pitchforks_genome",
  "Mean_IL_number_genome",
  "Mean_AvgLadder_genome",
  "Var_Cherries_genome",
  "Var_Pitchforks_genome",
  "Var_IL_number_genome",
  "Var_AvgLadder_genome",
  # balance
  "Mean_Stairs_genome",
  "Mean_Colless_genome",
  "Mean_Sackin_genome",
  "Mean_B2_genome",
  "Mean_MaxDepth_genome",
  "Var_Stairs_genome",
  "Var_Colless_genome",
  "Var_Sackin_genome",
  "Var_B2_genome",
  "Var_MaxDepth_genome"
)
stat_names <- c(
  "Bulk CN distance",
  "Single-cell CN distance",
  "Mean(sc Shannon index)",
  "Mean(misseg. count in bulk)",
  "Mean(clonal misseg. count in sc)",
  "Mean(subclonal misseg. count in sc)",
  "NA",
  "NA",
  "Var(sc Shannon index)",
  "Var(misseg. count in bulk)",
  "Var(clonal misseg. count in sc)",
  "Var(subclonal misseg. count in sc)",
  "NA",
  "NA",
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
# stat_names <- c(
#   #-----bulk
#   "Wasserstein_dist_bulk",
#   "Mean_Clonal_misseg_count_bulk",
#   "Var_Clonal_misseg_count_bulk",
#   #-----sc_subclonal CN
#   "Wasserstein_dist_sc",
#   "Mean_Shannon",
#   "Mean_Clonal_misseg_count_sc",
#   "Mean_Subclonal_misseg_count_sc",
#   "Mean_Clonal_armmisseg_count_sc",
#   "Mean_Subclonal_armmisseg_count_sc",
#   "Var_Shannon",
#   "Var_Clonal_misseg_count_sc",
#   "Var_Subclonal_misseg_count_sc",
#   "Var_Clonal_armmisseg_count_sc",
#   "Var_Subclonal_armmisseg_count_sc",
#   #-----sc_Phylo Stats
#   "Mean_Cherries",
#   "Mean_Pitchforks",
#   "Mean_IL_number",
#   "Mean_AvgLadder",
#   "Var_Cherries",
#   "Var_Pitchforks",
#   "Var_IL_number",
#   "Var_AvgLadder",
#   # balance
#   "Mean_Stairs",
#   "Mean_Colless",
#   "Mean_Sackin",
#   "Mean_B2",
#   "Mean_MaxDepth",
#   "Var_Stairs",
#   "Var_Colless",
#   "Var_Sackin",
#   "Var_B2",
#   "Var_MaxDepth"
# )

#-------------------------------------------Get correlation matrix (misseg)
statistics_misseg <- data.frame(matrix(ncol = 0, nrow = 100000))
for (i in 1:length(correct_misseg_stat_names)) {
  statistics_misseg[, i] <- ABC_input$sim_stat[[which(grepl(correct_misseg_stat_names[i], misseg_stat_names))]]
}
# for (i in 1:length(which(grepl("genome", misseg_stat_names)))) {
#   statistics_misseg[, i] <- ABC_input$sim_stat[which(grepl("genome", misseg_stat_names))[i]]
# }
prob_misseg <- data.frame(ABC_input$sim_param[, 1])
# colnames(statistics_misseg) <- misseg_stat_names[which(grepl("genome", misseg_stat_names))]
corr_mtx_misseg <- cor(y = prob_misseg, x = statistics_misseg)

#-------------------------------------------Get correlation matrix (selection)

corr_mtx_sel_mix <- NULL
for (i in 2:length(param_names)) {
  statistics_sel <- data.frame(matrix(ncol = 0, nrow = 100000))
  for (j in 1:length(sel_stat_names)) {
    if (grepl("genome", sel_stat_names[j])) {
      statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], misseg_stat_names))]]
    } else if (grepl("chr", sel_stat_names[j])) {
      statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], misseg_stat_names))]][, (i - 1)]
    }
  }
  colnames(statistics_sel) <- sel_stat_names
  prob_sel <- parameters[, i]
  corr_mtx_sel <- cor(y = prob_sel, x = statistics_sel)
  corr_mtx_sel_mix <- cbind(corr_mtx_sel_mix, corr_mtx_sel)
}
colnames(corr_mtx_sel_mix) <- param_names[-1]
dim(corr_mtx_sel_mix)
corr_mtx <- cbind(corr_mtx_misseg, corr_mtx_sel_mix)
rownames(corr_mtx) <- stat_names
colnames(corr_mtx) <- param_names

# locs_chosen_stat <- which(total_stat_names %in% stat_names)
# chosen_statistics <- data.frame(statistics[, locs_chosen_stat])
# chosen_statistics <- data.matrix(chosen_statistics)
# Correlation Plot=================================================
# =================================================================
#-------------------------------------------Get correlation matrix
# corr_mtx <- cor(y = parameters[, 1], x = statistics_misseg)
# mean(corr_mtx[, ])
# sort(mean(abs(corr_mtx[, 1])), TRUE)
# sort(abs(corr_mtx[, 1]), TRUE)
print(corr_mtx)
#--------------------------------------------Plot correlation plot

# corr_mtx <- cor(y = prob_misseg, x = statistics_misseg)
corr_plot <- function(corr_mtx) {
  library("ggcorrplot")
  plot <- ggcorrplot(
    corr_mtx,
    method = "circle",
    ggtheme = ggplot2::theme_bw,
    legend.title = "Correlation"
  ) +
    theme(axis.text.x = element_text(colour = c(
      "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
      "#00BA38", "#00BA38", "#00BA38",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
    ), angle = 45), plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "in"))
  # scale_fill_gradient2(limit = c(-0.3, 0.3), low = "blue", high = "red", )
  ggsave(file = "correlation_plot.png", width = 12, height = 10, units = "in", plot = plot, dpi = 300, limitsize = TRUE)
  return(plot)
}

corr_plot(corr_mtx)

# list_targets_library <- c(
#   #---Bulk DNA: CN
#   "statistic_group=Copy number statistics;statistic_ID=Bulk CN distance;data=bulk;target=genome;statistic=dist;variable=average_CN;metric=euclidean",
#   "statistic_group=Copy number statistics;statistic_ID=Mean(misseg. count in bulk);data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=Var(misseg. count in bulk);data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
#   paste0("statistic_group=Copy number statistics;statistic_ID=Bulk CN distance;data=bulk;target=chromosome;statistic=dist;variable=average_CN;metric=euclidean;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Mean(misseg. count in bulk);data=bulk;target=chromosome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Var(misseg. count in bulk);data=bulk;target=chromosome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation;chromosome=", list_chromosomes),
#   #---Single-cell DNA: subclonal CN
#   "statistic_group=Copy number statistics;statistic_ID=Single-cell CN distance;data=sc;target=genome;statistic=dist;variable=clonal_CN;metric=euclidean",
#   "statistic_group=Copy number statistics;statistic_ID=Mean(sc Shannon index);data=sc;target=genome;statistic=mean;variable=shannon",
#   "statistic_group=Copy number statistics;statistic_ID=Mean(clonal misseg. count in sc);data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=Mean(subclonal misseg. count in sc);data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=Var(sc Shannon index);data=sc;target=genome;statistic=var;variable=shannon",
#   "statistic_group=Copy number statistics;statistic_ID=Var(clonal misseg. count in sc);data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=Var(subclonal misseg. count in sc);data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
#   "statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
#   paste0("statistic_group=Copy number statistics;statistic_ID=Single-cell CN distance;data=sc;target=chromosome;statistic=dist;variable=clonal_CN;metric=euclidean;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Mean(sc Shannon index);data=sc;target=chromosome;statistic=mean;variable=shannon;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Mean(clonal misseg. count in sc);data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Mean(subclonal misseg. count in sc);data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=chromosome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=chromosome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Var(sc Shannon index);data=sc;target=chromosome;statistic=var;variable=shannon;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Var(clonal misseg. count in sc);data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=Var(subclonal misseg. count in sc);data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=chromosome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
#   paste0("statistic_group=Copy number statistics;statistic_ID=NA;data=sc;target=chromosome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation;chromosome=", list_chromosomes),
#   #---Single-cell DNA: phylo stats for tips
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Mean(cherry count);data=sc;target=genome;statistic=mean;variable=cherries", # number of internal nodes with 2 tips
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Mean(pitchfork count);data=sc;target=genome;statistic=mean;variable=pitchforks", # number of internal tips with 3 tips
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Mean(IL number);data=sc;target=genome;statistic=mean;variable=IL_number", # number of internal nodes with single tip childs
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Mean(average ladder);data=sc;target=genome;statistic=mean;variable=avgLadder", # mean size of ladder (sequence of internal nodes, each with single tip childs)
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Var(cherry count);data=sc;target=genome;statistic=var;variable=cherries",
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Var(pitchfork count);data=sc;target=genome;statistic=var;variable=pitchforks",
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Var(IL number);data=sc;target=genome;statistic=var;variable=IL_number",
#   "statistic_group=Phylogeny leaf statistics;statistic_ID=Var(average ladder);data=sc;target=genome;statistic=var;variable=avgLadder",
#   #---Single-cell DNA: phylo stats for balance
#   "statistic_group=Phylogeny balance index;statistic_ID=Mean(stairs);data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
#   "statistic_group=Phylogeny balance index;statistic_ID=Mean(Colless index);data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
#   "statistic_group=Phylogeny balance index;statistic_ID=Mean(Sackin index);data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
#   "statistic_group=Phylogeny balance index;statistic_ID=Mean(B2 index);data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
#   "statistic_group=Phylogeny balance index;statistic_ID=Mean(max depth);data=sc;target=genome;statistic=mean;variable=maxDepth", # height of phylogeny tree
#   "statistic_group=Phylogeny balance index;statistic_ID=Var(stairs);data=sc;target=genome;statistic=var;variable=stairs",
#   "statistic_group=Phylogeny balance index;statistic_ID=Var(Colless index);data=sc;target=genome;statistic=var;variable=colless",
#   "statistic_group=Phylogeny balance index;statistic_ID=Var(Sackin index);data=sc;target=genome;statistic=var;variable=sackin",
#   "statistic_group=Phylogeny balance index;statistic_ID=Var(B2 index);data=sc;target=genome;statistic=var;variable=B2",
#   "statistic_group=Phylogeny balance index;statistic_ID=Var(max depth);data=sc;target=genome;statistic=var;variable=maxDepth"
# )
