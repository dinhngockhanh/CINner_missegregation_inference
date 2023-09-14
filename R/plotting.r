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
    dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)))
    dist_posterior <- density(resp, weights = weights.std[, i])

    df_plot_prior <- data.frame(x = dist_prior$x, y = dist_prior$y)
    df_plot_posterior <- data.frame(x = dist_posterior$x, y = dist_posterior$y)

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
                               para_lower_bound = NULL,
                               para_upper_bound = NULL,
                               ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    df_plot <- densityPlot_df(
        object,
        obs,
        training,
        add,
        main,
        color_prior,
        chosen_para,
        color_posterior,
        protocol,
        color_vline,
        log,
        xlim,
        ylim,
        xlab,
        ylab,
        paral,
        ncores
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
plot_statistics_correlation <- function() {

}

#' @export
plot_ABC_correlation <- function(inference_result = parameters_inferred,
                                 library_name = library_name,
                                 value_x = "Ground_truth",
                                 value_y = NULL,
                                 title_x = "",
                                 title_y = "",
                                 error_y = "Sd",
                                 color_data = "red",
                                 color_regression = "blue",
                                 fontsize = 50,
                                 linear_regression = FALSE) {
    df <- parameters_inferred
    xvalue <- df[[value_x]]
    yvalue <- df[[value_y]]
    error <- df[[error_y]]
    corr_plot <- ggplot(df, mapping = aes(x = xvalue, y = yvalue)) +
        geom_errorbar(aes(ymin = yvalue - error, ymax = yvalue + error), width = 0.01, size = 1, colour = color_data) +
        geom_point(colour = color_data, size = 3) +
        xlim(0.89, 1.21) +
        ylim(0.89, 1.21) +
        xlab(title_x) +
        ylab(title_y) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(aspect.ratio = 1)
    if (linear_regression) {
        corr_plot <- corr_plot + stat_smooth(method = "lm", fill = color_regression, colour = color_regression)
    }
    filename <- paste0(library_name, "_ABC_correlation.jpeg")
    jpeg(filename, width = 1500, height = 1500)
    print(corr_plot)
    dev.off()
}
