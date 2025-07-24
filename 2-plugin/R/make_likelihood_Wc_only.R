make_likelihood_Wc_only <- function (tmle_task, learner_list) {
  Wc_factor <- define_lf(LF_emp, "Wc")
  We_factor <- define_lf(LF_emp, "We")
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1/tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  }  else {
    A_bound <- NULL
  }
  A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], 
                        bound = A_bound)
  Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], 
                        type = "mean")
  factor_list <- list(Wc_factor, We_factor, A_factor, Y_factor)
  if (!is.null(tmle_task$npsem[["Y"]]$censoring_node)) {
    if (is.null(learner_list[["delta_Y"]])) {
      stop("Y is subject to censoring, but no learner was specified for censoring mechanism delta_Y")
    }
    delta_Y_factor <- define_lf(LF_fit, "delta_Y", learner = learner_list[["delta_Y"]], 
                                type = "mean", bound = c(0.025, 1))
    factor_list <- c(factor_list, delta_Y_factor)
  }
  if (!is.null(tmle_task$npsem[["A"]]$censoring_node)) {
    stop("A is subject to censoring, this isn't supported yet")
  }
  if (!is.null(tmle_task$npsem[["Wc"]]$censoring_node)) {
    stop("Wc is subject to censoring, this isn't supported yet")
  }
  if (!is.null(tmle_task$npsem[["We"]]$censoring_node)) {
    stop("We is subject to censoring, this isn't supported yet")
  }
  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}
