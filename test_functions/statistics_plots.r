#----------------------------------------------------Load the rda
Paths <- "/Users/xiangzijin/Downloads/"
setwd(Paths)
rda_name <- "Simpler_DLP&BULK_DNA_ABC_input1.rda"
load(file = rda_name)
#---------------------------------------------------Get statistics
statistics <- ABC_input$sim_stat
ncol(statistics)
stat_names <- c(
  "Wasserstein_dist_bulk",
  "Mean_Clonal_misseg_count_bulk",
  "Mean_Shannon",
  "Mean_Clonal_misseg_count_sc",
  "Mean_Subclonal_misseg_count_sc",
  # "Mean_Clonal_armmisseg_count_sc",
  # "Mean_Subclonal_armmisseg_count_sc",
  "Var_Shannon",
  "Var_Clonal_misseg_count_sc",
  "Var_Subclonal_misseg_count_sc",
  # "Var_Clonal_armmisseg_count_sc",
  # "Var_Subclonal_armmisseg_count_sc",
  "Wasserstein_dist_sc",
  "Mean_Cherries",
  "Mean_Pitchforks",
  "Mean_IL_number",
  "Mean_AvgLadder",
  "Var_Cherries",
  "Var_Pitchforks",
  "Var_IL_number",
  "Var_AvgLadder",
  "Mean_Stairs",
  "Mean_Colless",
  "Mean_Sackin",
  "Mean_B2",
  "Mean_MaxDepth",
  "Var_Stairs",
  "Var_Colless",
  "Var_Sackin",
  "Var_B2",
  "Var_MaxDepth"
)
total_stat_names <- c(
  "Wasserstein_dist_bulk",
  "Mean_Clonal_misseg_count_bulk",
  "Mean_Shannon",
  "Mean_Clonal_misseg_count_sc",
  "Mean_Subclonal_misseg_count_sc",
  "Mean_Clonal_armmisseg_count_sc",
  "Mean_Subclonal_armmisseg_count_sc",
  "Var_Shannon",
  "Var_Clonal_misseg_count_sc",
  "Var_Subclonal_misseg_count_sc",
  "Var_Clonal_armmisseg_count_sc",
  "Var_Subclonal_armmisseg_count_sc",
  "Wasserstein_dist_sc",
  "Mean_Cherries",
  "Mean_Pitchforks",
  "Mean_IL_number",
  "Mean_AvgLadder",
  "Var_Cherries",
  "Var_Pitchforks",
  "Var_IL_number",
  "Var_AvgLadder",
  "Mean_Stairs",
  "Mean_Colless",
  "Mean_Sackin",
  "Mean_B2",
  "Mean_MaxDepth",
  "Var_Stairs",
  "Var_Colless",
  "Var_Sackin",
  "Var_B2",
  "Var_MaxDepth"
)
colnames(statistics) <- total_stat_names

locs_chosen_stat <- which(total_stat_names %in% stat_names)
chosen_statistics <- data.frame(statistics[, locs_chosen_stat])
chosen_statistics <- data.matrix(chosen_statistics)
#---------------------------------------------------Get parameters
parameters <- ABC_input$sim_param
param_names <- c(
  "prob_misseg",
  "1p",
  "2p",
  "3p",
  "4p",
  "5p",
  "6p",
  "7p",
  "8p",
  "9p",
  "10p",
  "11p",
  "12p",
  "13p",
  "14p",
  "15p",
  "16p",
  "17p",
  "18p",
  "19p",
  "20p",
  "21p",
  "22p",
  "Xp"
)
# get rid of last column of parameters
colnames(parameters) <- param_names
# prob_armmisseg <- parameters[, "prob_armmisseg"]
prob_misseg <- parameters[, "prob_misseg"]
# Correlation Plot=================================================
# =================================================================
#-------------------------------------------Get correlation matrix
corr_mtx <- cor(y = parameters, x = chosen_statistics)
# mean(corr_mtx[, ])
# sort(mean(abs(corr_mtx[, 1])), TRUE)
# sort(abs(corr_mtx[, 1]), TRUE)
#--------------------------------------------Plot correlation plot


corr_plot <- function(corr_mtx) {
  library("ggcorrplot")
  plot <- ggcorrplot(corr_mtx, ggtheme = ggplot2::theme_gray) +
    theme(axis.text.x = element_text(colour = c(
      "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
      "#00BA38", "#00BA38",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
      # )))
      # theme(axis.text.x = element_text(colour = c(
      # "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
      # "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      # "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
    )))
  # scale_fill_gradient2(limit = c(-0.3, 0.3), low = "blue", high = "red", )
  ggsave(file = "correlation_plot.png", plot = plot, width = 13, height = 10, units = "in", dpi = 300, limitsize = FALSE)
  return(plot)
}
corr_plot(corr_mtx)
min(corr_mtx)
max(corr_mtx)
# P-value Plot=====================================================
# =================================================================
p_value_plot <- function(parameter, statistic, title) {
  #--------------------------------------Set different types of stats
  library(ggplot2)
  CN_stats <- c(
    "Mean_Shannon",
    "Mean_Clonal_misseg_count_sc",
    "Mean_Clonal_misseg_count_bulk",
    "Mean_Subclonal_misseg_count_sc",
    "Mean_Clonal_armmisseg_count_sc",
    "Mean_Subclonal_armmisseg_count_sc",
    "Var_Shannon",
    "Var_Clonal_misseg_count_sc",
    "Var_Subclonal_misseg_count_sc",
    "Var_Clonal_armmisseg_count_sc",
    "Var_Subclonal_armmisseg_count_sc",
    "Wasserstein_dist_sc",
    "Wasserstein_dist_bulk"
  )

  Phylo_tip <- c(
    "Mean_Cherries",
    "Mean_Pitchforks",
    "Mean_IL_number",
    "Mean_AvgLadder",
    "Var_Cherries",
    "Var_Pitchforks",
    "Var_IL_number",
    "Var_AvgLadder"
  )

  Phylo_balance <- c(
    "Mean_Stairs",
    "Mean_Colless",
    "Mean_Sackin",
    "Mean_B2",
    "Mean_MaxDepth",
    "Var_Stairs",
    "Var_Colless",
    "Var_Sackin",
    "Var_B2",
    "Var_MaxDepth"
  )
  #-------------------------------------Get linear regression model
  model <- lm(parameter ~ statistics)
  p_values <- summary(model)$coefficients[, 4][-1]
  orig_names <- names(p_values)
  stats_names <- sub(".*statistics", "", orig_names)
  p_values <- as.numeric(unlist(p_values))
  #------------------------------Get stat names and corres p-values
  p_value_df <- data.frame(stats_names, p_values)
  p_value_df$p_values[which(p_value_df$p_values <= 1e-5)] <- 1e-5
  #------------------------------------------Assign different types
  #--------------------------------------------------Get p-value df
  p_value_df$type <- rep(NA, length(stats_names))
  locs_cn <- which(p_value_df$stats_names %in% CN_stats)
  locs_tip <- which(p_value_df$stats_names %in% Phylo_tip)
  locs_balance <- which(p_value_df$stats_names %in% Phylo_balance)
  p_value_df[locs_cn, "type"] <- "locs_cn"
  p_value_df[locs_tip, "type"] <- "locs_tip"
  p_value_df[locs_balance, "type"] <- "locs_balance"
  #----------------------------------------Plot it and save the plot
  plot <- ggplot(p_value_df, aes(
    x = factor(stats_names, level = stats_names),
    y = p_values, color = factor(type), shape = factor(type),
    group = 1
  )) +
    geom_point(size = 8) +
    geom_line(linewidth = 1, y = log10(1e-2), color = "#a2a0a0") +
    ggtitle(paste0("P_values for ", title)) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    # theme(axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(
      trans = "log10",
      limits = c(1e-5, 1),
      breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1),
      labels = c(expression("" <= 1e-5), "1e-4", "1e-3", "1e-2", "1e-1", "1")
    ) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank()) +
    theme(aspect.ratio = 1) +
    theme(text = element_text(size = 16), plot.title = element_text(size = 20), axis.title.y = element_text(size = 20)) +
    theme(axis.text.x = element_text(colour = c(
      # "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38","#00BA38",
      # "#619CFF", "#619CFF", "#619CFF",
      # "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
      # ))) +
      # theme(axis.text.x = element_text(colour = c(
      "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
      "#00BA38",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
    )))

  ggsave(file = paste0("P_values for ", title, ".png"), plot = plot, width = 10, height = 10, dpi = 300)
  return(plot)
}

p_value_plot(parameters[, 24], statistics, title = "Xp")
prob_sel <- parameters[, "prob_misseg"]
# lm(prob_misseg ~ chosen_statistics)
# summary(lm(prob_misseg ~ statistics))
# DA analysis Plot=====================================================
# =================================================================
library("domir")
