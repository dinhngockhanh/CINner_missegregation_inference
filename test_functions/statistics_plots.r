#----------------------------------------------------Load the rda
# Paths <- "/Users/xiangzijin/Downloads/"
# setwd(Paths)
# rda_name <- "Simpler_DLP_CNA_ABC_input.rda"
# load(file = rda_name)
#---------------------------------------------------Get statistics
statistics <- ABC_input$sim_stat
stat_names <- c(
  "Shannon",
  "Clonal_misseg_count",
  "Subclonal_misseg_count",
  "Clonal_armmisseg_count",
  "Subclonal_armmisseg_count",
  "Wasserstein_dist",
  "Cherries",
  "Pitchforks",
  "IL_number",
  "AvgLadder",
  "Stairs",
  "Colless",
  "Sackin",
  "B2",
  "MaxDepth"
)
colnames(statistics) <- stat_names
#---------------------------------------------------Get parameters
parameters <- ABC_input$sim_param
param_names <- c(
  "prob_misseg",
  "prob_armmisseg"
)
colnames(parameters) <- param_names
prob_armmisseg <- parameters[, "prob_armmisseg"]
prob_misseg <- parameters[, "prob_misseg"]
# Correlation Plot=================================================
# =================================================================
#-------------------------------------------Get correlation matrix
corr_mtx <- cor(y = parameters, x = statistics)
#--------------------------------------------Plot correlation plot
corr_plot <- function(corr_mtx) {
  library("ggcorrplot")
  plot <- ggcorrplot(corr_mtx, ggtheme = ggplot2::theme_gray) +
    theme(axis.text.x = element_text(colour = c(
      "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
      "#619CFF", "#619CFF", "#619CFF", "#619CFF",
      "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
    )))
  ggsave(file = "correlation_plot.png", plot = plot, width = 10, height = 10, dpi = 300)
  return(plot)
}
corr_plot(corr_mtx)
# P-value Plot=====================================================
# =================================================================
p_value_plot <- function(parameter, statistics) {
  #--------------------------------------Set different types of stats
  CN_stats <- c(
    "Shannon",
    "Clonal_misseg_count",
    "Subclonal_misseg_count",
    "Clonal_armmisseg_count",
    "Subclonal_armmisseg_count",
    "Wasserstein_dist"
  )

  Phylo_tip <- c(
    "Cherries",
    "Pitchforks",
    "IL_number",
    "AvgLadder"
  )

  Phylo_balance <- c(
    "Stairs",
    "Colless",
    "Sackin",
    "B2",
    "MaxDepth"
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
    ggtitle("P_values for Prob_misseg") +
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
    theme(text = element_text(size = 16), plot.title = element_text(size = 20), axis.title.y = element_text(size = 20))

  ggsave(file = "p_value.png", plot = plot, width = 10, height = 10, dpi = 300)
  return(plot)
}

p_value_plot(prob_armmisseg, statistics)
