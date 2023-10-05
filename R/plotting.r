densityPlot_df <- function(object,
                           obs,
                           training,
                           add = TRUE,
                           main = "Posterior density",
                           color_prior,
                           chosen_para = NULL,
                           color_posterior,
                           protocol = "",
                           color_vline,
                           log = "",
                           xlim = NULL,
                           ylim = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           paral = FALSE,
                           cutbound = FALSE,
                           lower = NULL,
                           upper = NULL,
                           ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    findweights <- getFromNamespace("findweights", "abcrf")
    ### Checking arguments
    if (!inherits(object, "regAbcrf")) {
        stop("object not of class regAbcrf")
    }

    if (!inherits(training, "data.frame")) {
        stop("training needs to be a data.frame object")
    }

    if (!inherits(obs, "data.frame")) {
        stop("obs needs to be a data.frame object")
    }
    if (nrow(obs) == 0L || is.null(nrow(obs))) {
        stop("no data in obs")
    }
    if (nrow(training) == 0L || is.null(nrow(training))) {
        stop("no simulation in the training reference table (response, sumstat)")
    }

    if ((!is.logical(add)) || (length(add) != 1L)) {
        stop("add should be TRUE or FALSE")
    }
    if ((!is.logical(paral)) || (length(paral) != 1L)) {
        stop("paral should be TRUE or FALSE")
    }
    if (is.na(ncores)) {
        warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
        ncores <- 1
    }

    if (!is.character(log)) {
        stop("log needs to be a character string")
    }
    x <- obs
    if (!is.null(x)) {
        if (is.vector(x)) {
            x <- matrix(x, ncol = 1)
        }
        if (nrow(x) == 0) {
            stop("obs has 0 rows")
        }
        if (any(is.na(x))) {
            stop("missing values in obs")
        }
    }

    # resp and sumsta recover

    mf <- match.call(expand.dots = FALSE)
    mf <- mf[1]
    mf$formula <- object$formula


    mf$data <- training


    training <- mf$data

    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    resp <- model.response(mf)

    obj <- object$model.rf
    inbag <- matrix(unlist(obj$inbag.counts, use.names = FALSE), ncol = obj$num.trees, byrow = FALSE)

    obj[["origNodes"]] <- predict(obj, training, predict.all = TRUE, num.threads = ncores)$predictions
    obj[["origObs"]] <- model.response(mf)

    #####################

    origObs <- obj$origObs
    origNodes <- obj$origNodes

    nodes <- predict(obj, x, predict.all = TRUE, num.threads = ncores)$predictions
    if (is.null(dim(nodes))) nodes <- matrix(nodes, nrow = 1)
    ntree <- obj$num.trees
    nobs <- object$model.rf$num.samples
    nnew <- nrow(x)

    weights <- findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
    weights.std <- weights / ntree

    priorDensity <- density(resp)

    if (add) {
        rangex <- range(priorDensity$x)
        rangey <- range(priorDensity$y)

        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)
            rangex <- range(rangex, postDensity$x)
            rangey <- range(rangey, postDensity$y)
        }

        # plot(priorDensity$x, priorDensity$y,
        #     type = "l", main = main, log = log,
        #     xlim = if (is.null(xlim)) rangex else xlim,
        #     ylim = if (is.null(ylim)) rangey else ylim,
        #     xlab = xlab, ylab = ylab, col = "grey"
        # )
        # for (i in 1:nnew) {
        #     postDensity <- density(resp, weights = weights.std[, i], ...)
        #     points(postDensity$x, postDensity$y, type = "l")
        # }
    } else {
        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)

            # plot(postDensity$x, postDensity$y,
            #     type = "l", main = main, log = log,
            #     xlim = if (is.null(xlim)) range(postDensity$x, priorDensity$x) else xlim,
            #     ylim = if (is.null(ylim)) range(postDensity$y, priorDensity$y) else ylim,
            #     xlab = xlab, ylab = ylab
            # )
            # points(priorDensity$x, priorDensity$y, type = "l", col = "grey")
            # if (nnew > 1 && i < nnew) readline("Press <ENTER> to Continue")
        }
    }



    if (protocol == "TSG") {
        resp <- 1 / resp
    }

    if (cutbound == TRUE) {
        dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)), from = lower, to = upper)
        dist_posterior <- density(resp, weights = weights.std[, i], from = lower, to = upper)
    } else if (cutbound == FALSE) {
        dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)))
        dist_posterior <- density(resp, weights = weights.std[, i])
    }
    df_plot_prior <- data.frame(x = dist_prior$x, y = dist_prior$y)
    # df_plot_posterior <- data.frame(x = dist_posterior$x, y = dist_posterior$y)

    df_plot <- data.frame(x = dist_prior$x, y_prior = dist_prior$y, y_posterior = dist_posterior$y)

    # df_plot <- data.frame(dist_raw = resp)
    # df_plot$weight_prior <- 1 / nrow(df_plot)
    # df_plot$weight_posterior <- weights.std[, i]

    return(df_plot)
}

plot_ABC_inference <- function(object,
                               obs,
                               training,
                               add = TRUE,
                               main = "Posterior density",
                               protocol,
                               color_prior = "lightblue",
                               color_posterior = "darkblue",
                               highlight_values = NULL,
                               highlight_colors = NULL,
                               highlight_linetype = NULL,
                               log = "",
                               xlim = NULL,
                               ylim = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               paral = FALSE,
                               fontsize = 50,
                               plot_ABC_prior_as_uniform = FALSE,
                               cutbound = FALSE,
                               para_lower_bound = NULL,
                               para_upper_bound = NULL,
                               ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    df_plot <- densityPlot_df(
        object = object,
        obs = obs,
        training = training,
        add = add,
        main = main,
        color_prior = color_prior,
        chosen_para = chosen_para,
        color_posterior = color_posterior,
        protocol = protocol,
        color_vline = color_vline,
        log = log,
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        paral = paral,
        ncores = ncores,
        cutbound = cutbound,
        lower = para_lower_bound,
        upper = para_upper_bound
    )

    if (plot_ABC_prior_as_uniform) {
        p_plot <- ggplot(df_plot) +
            annotate("rect",
                xmin = para_lower_bound, xmax = para_upper_bound,
                ymin = 0, ymax = 1 / (para_upper_bound - para_lower_bound),
                color = color_prior, fill = color_prior, alpha = 0.3
            )
    } else {
        p_plot <- ggplot(df_plot) +
            geom_area(aes(x = x, y = y_prior),
                color = color_prior, fill = color_prior, alpha = 0.3
            )
    }
    p_plot <- p_plot +
        geom_area(aes(x = x, y = y_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_prior), color = color_prior, fill = color_prior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        xlab("") +
        ylab("") +
        ggtitle(main) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    if (!is.null(highlight_values)) {
        p_plot <- p_plot +
            geom_vline(
                xintercept = highlight_values,
                colour = highlight_colors,
                linetype = highlight_linetype,
                size = 2
            )
    }
    return(p_plot)
}

#' @export
plot_statistics_correlation <- function(filename,
                                        list_targets_library,
                                        list_parameters) {
    library("ggcorrplot")
    load(file = filename)
    # =========================================Get list of parameter IDs
    parameter_IDs <- list_parameters$Title
    parameter_types <- list_parameters$Type
    # =========================================Get list of statistic IDs
    statistic_IDs <- c()
    statistic_groups <- c()
    for (stat in list_targets_library) {
        stat_details <- strsplit(stat, ";")[[1]]
        stat_ID <- strsplit(stat_details[grep("statistic_ID=", stat_details)], "=")[[1]][2]
        stat_group <- strsplit(stat_details[grep("statistic_group=", stat_details)], "=")[[1]][2]
        statistic_IDs <- c(statistic_IDs, stat_ID)
        statistic_groups <- c(statistic_groups, stat_group)
    }
    unique_statistic_IDs <- unique(statistic_IDs)
    unique_statistic_groups <- unique(statistic_groups)








    # ====================================Compute the correlation matrix
    correlation_matrix <- data.frame(matrix(nrow = 0, ncol = length(parameter_IDs)))
    colnames(correlation_matrix) <- parameter_IDs
    for (i_stat in 1:length(unique_statistic_IDs)) {
        statistic_ID <- unique_statistic_IDs[i_stat]
        statistic_group <- unique_statistic_groups[i_stat]
        for (i_para in 1:length(list_parameters)) {
            parameter_ID <- parameter_IDs[i_para]
            parameter_type <- parameter_types[i_para]
            sim_stat_loc <- which(grepl(statistic_ID, statistic_IDs))
            print(sim_stat_loc)
        }
    }



    # =========================================Initialize the parameters
    parameters <- ABC_input$sim_param
    param_names <- list_parameters$Title
    colnames(parameters) <- param_names
    # ===================================================Get the correlation matrix for misseg rates
    statistics_misseg <- data.frame(matrix(ncol = 0, nrow = ABC_simcount))
    for (i in 1:length(which(grepl("genome", misseg_stat_names)))) {
        statistics_misseg[, i] <- ABC_input$sim_stat[which(grepl("genome", misseg_stat_names))[i]]
    }
    prob_misseg <- parameters[, 1]
    corr_mtx_misseg <- cor(y = prob_misseg, x = statistics_misseg)
    # =====================================================Get the correlation matrix for selections
    # =======================================and combine it with correlation matrix for misseg rates
    corr_mtx_sel_mix <- NULL
    for (i in 2:length(param_names)) {
        statistics_sel <- data.frame(matrix(ncol = 0, nrow = ABC_simcount))
        for (j in 1:length(sel_stat_names)) {
            if (grepl("genome", sel_stat_names[j])) {
                statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], misseg_stat_names))]]
            } else if (grepl("chr", sel_stat_names[j])) {
                statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], misseg_stat_names))]][, (i - 1)]
            }
        }
        prob_sel <- parameters[, i]
        corr_mtx_sel <- cor(y = prob_sel, x = statistics_sel)
        corr_mtx_sel_mix <- cbind(corr_mtx_sel_mix, corr_mtx_sel)
    }
    colnames(corr_mtx_sel_mix) <- param_names[-1]
    corr_mtx <- cbind(corr_mtx_misseg, corr_mtx_sel_mix)
    rownames(corr_mtx) <- stat_names
    colnames(corr_mtx) <- param_names
    # ============================================================Plot the Correlation
    plot <- ggcorrplot(corr_mtx, ggtheme = ggplot2::theme_gray) +
        theme(axis.text.x = element_text(colour = c(
            "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
            "#00BA38", "#00BA38", "#00BA38",
            "#619CFF", "#619CFF", "#619CFF", "#619CFF",
            "#619CFF", "#619CFF", "#619CFF", "#619CFF",
            "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
            "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
        )), plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "in"))
    # scale_fill_gradient2(limit = c(-0.3, 0.3), low = "blue", high = "red", )
    ggsave(file = "correlation_plot.png", width = 10, height = 10, units = "in", plot = plot, dpi = 300, limitsize = TRUE)
    return(plot)
}

#' @export
plot_ABC_correlation <- function(inference_result = parameters_inferred,
                                 plot_name = "ABC_correlation_plot.jpg",
                                 value_x = "Ground_truth",
                                 value_y = NULL,
                                 title_x = "",
                                 title_y = "",
                                 error_x = NULL,
                                 error_y = NULL,
                                 title_plot = "",
                                 color_data = "red",
                                 fontsize = 50,
                                 plot_diagonal = FALSE,
                                 plot_Error = FALSE) {
    library(ggplot2)
    library(ehaGoF)
    if (!is.null(error_y)) {
        parameters_inferred$value_y_min <- parameters_inferred[[value_y]] - parameters_inferred[[error_y]]
        parameters_inferred$value_y_max <- parameters_inferred[[value_y]] + parameters_inferred[[error_y]]
        y_min <- min(parameters_inferred$value_y_min)
        y_max <- max(parameters_inferred$value_y_max)
    } else {
        y_min <- min(parameters_inferred[[value_y]])
        y_max <- max(parameters_inferred[[value_y]])
    }
    if (!is.null(error_x)) {
        parameters_inferred$value_x_min <- parameters_inferred[[value_x]] - parameters_inferred[[error_x]]
        parameters_inferred$value_x_max <- parameters_inferred[[value_x]] + parameters_inferred[[error_x]]
        x_min <- min(parameters_inferred$value_x_min)
        x_max <- max(parameters_inferred$value_x_max)
    } else {
        x_min <- min(parameters_inferred[[value_x]])
        x_max <- max(parameters_inferred[[value_x]])
    }
    corr_plot <- ggplot(parameters_inferred, mapping = aes_string(x = value_x, y = value_y)) +
        geom_point(colour = color_data, size = 10) +
        xlim(min(x_min, y_min), max(x_max, y_max)) +
        ylim(min(x_min, y_min), max(x_max, y_max)) +
        xlab(title_x) +
        ylab(title_y) +
        ggtitle(title_plot) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(aspect.ratio = 1)
    if (!is.null(error_x)) {
        corr_plot <- corr_plot +
            geom_errorbar(aes(ymin = value_x_min, ymax = value_x_max), width = 0, size = 1, colour = color_data)
    }
    if (!is.null(error_y)) {
        corr_plot <- corr_plot +
            geom_errorbar(aes(ymin = value_y_min, ymax = value_y_max), width = 0, size = 1, colour = color_data)
    }
    if (plot_diagonal) {
        corr_plot <- corr_plot + geom_abline(intercept = 0)
    }
    if (plot_Error) {
        Error <- compute_Error(
            results = parameters_inferred,
            ID_actual = value_x, ID_predicted = value_y
        )
        text_x <- min(x_min, y_min) + 0.05 * (max(x_max, y_max) - min(x_min, y_min))
        text_y <- min(x_min, y_min) + 0.95 * (max(x_max, y_max) - min(x_min, y_min))
        corr_plot <- corr_plot +
            annotate(
                "text",
                x = text_x, y = text_y,
                label = paste("RMSE=", round(Error, digits = 3)),
                size = unit(fontsize / 2, "pt"), color = color_data, hjust = 0
            )
    }
    jpeg(plot_name, width = 1500, height = 1500)
    print(corr_plot)
    dev.off()
}
