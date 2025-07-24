Lrnr_screener_name <- R6Class(
  classname = "Lrnr_screener_name",
  inherit = Lrnr_base, portable = TRUE, class = TRUE,
  public = list(
    initialize = function(var_name = NULL) {
      params <- args_to_list()
      private$.var_name <- var_name
      super$initialize(params = params)
    }
  ),
  private = list(
    .properties = c("binomial", "continuous", "categorical", "screener"),
    .var_name = NULL, 
    .train = function(task) {
      outcome_type <- self$get_outcome_type(task)
      X <- task$X
      Y <- outcome_type$format(task$Y)
      covs <- task$nodes$covariates
      
      args <- self$params
      # type <- set_correlation_screener_type(args$type, args$num_screen)
      # method <- args$method
      
      selected <- names(X) == private$.var_name
      
      # list_pvalues <- apply(X, 2, function(x, Y, method) {
      #   ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
      # }, Y = Y, method = method)
      
      # if (type == "rank") {
      #   selected <- (rank(list_pvalues) <= args$num_screen)
      # } else if (type == "threshold") {
      #   selected <- (list_pvalues <= args$pvalue_threshold)
      #   if (sum(selected) < args$min_screen) {
      #     selected[rank(list_pvalues) <= args$min_screen] <- TRUE
      #   }
      # }
      
      selected_names <- names(X)[selected]
      selected_covs <- sapply(covs, function(cov) any(grep(paste0("^", cov, "$"), selected_names)))
      fit_object <- list(selected = covs[selected_covs])
      return(fit_object)
    },
    .predict = function(task) {
      task$data[, private$.fit_object$selected, with = FALSE, drop = FALSE]
    },
    .chain = function(task) {
      return(task$next_in_chain(covariates = private$.fit_object$selected))
    }
  )
)
