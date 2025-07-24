mse <- function(x,y){mean((x-y)^2)}

se_mse = function(x,y){
  loss <- (x-y)^2
  n <- length(loss)
  se <- sd(loss)/sqrt(n)
  return (se)
}

calculate_test_result <- function(fit_object, task_test, learner_name, 
                                  includeBaseAverage = TRUE){
  pred <- fit_object$predict(task_test)
  truth <- task_test$Y
  if ("Stack" %in% class(fit_object)){
    # Base learner prediction
    base_mse <- apply(pred, MARGIN = 2, mse, truth)
    base_se <- apply(pred, MARGIN = 2, se_mse, truth)
    base_names <- names(fit_object$learner_fits)
    result <- 
      data.frame(type = "Base Learner",
                 learner = base_names,
                 mse = base_mse,
                 se = base_se)
    if (includeBaseAverage){
      n_lrnrs <- length(base_names)
      weight <- rep(1/n_lrnrs, n_lrnrs)
      pred_avg <- metalearner_linear(alpha = weight, X = as.matrix(pred))
      result_avg <- 
        data.frame(type = "Super Learner",
                   learner = "Average",
                   mse =  mse(pred_avg, truth),
                   se = se_mse(pred_avg, truth))
      result <- rbind(result_avg, result)
    }
    rownames(result) <- NULL
  } else {
    result <- 
      data.frame(type = "Super Learner",
                 learner = learner_name,
                 mse = mse(pred, truth),
                 se = se_mse(pred, truth))
    
  }
  return (result)
}

runMetaHAL <- function(task, base_lrnrs, max_degree, seed = 2){
  # 0. fast M-HAL-MLE selector
  lrnr_hal <- Lrnr_hal9001$new(max_degree = max_degree,
                               return_lasso = TRUE,
                               fit_control = list(type.measure = "mse",
                                                  cv_select = TRUE))
  meta_hal_fast_cv <- Lrnr_sl$new(learners = base_lrnrs, 
                                  metalearner = lrnr_hal)
  set.seed(seed)
  meta_hal_fast_cv_fit <- meta_hal_fast_cv$train(task)
  
  # 1. Get candidate lambdas
  lambda_star_fast <- meta_hal_fast_cv_fit$metalearner_fit()$lambda_star
  metalearner_fit <- meta_hal_fast_cv_fit$metalearner_fit()$lasso_fit$glmnet.fit
  lambda_list <- metalearner_fit$lambda
  
  # be more conservative with honest lambda selection, and also speed up; around the fast lambda, at most 71 lambdas
  lambda_list <- lambda_list[intersect(which(lambda_star_fast == lambda_list) + (-35):35, seq_along(lambda_list))]
  
  # 2. Pass HAL Params & Generate a meta-HAL for every lambda
  params_honestCV <- meta_hal_fast_cv_fit$params[["metalearner"]]$params
  
  metalrnr_hal_candidates <- Lrnr_hal9001$new(
    max_degree = params_honestCV$max_degree,
    return_lasso = FALSE,
    lambda = lambda_list,
    fit_control = list(type.measure = "mse",
                       cv_select = FALSE))
  
  # 3. Try Honest MetaHAL by discrete SL. It works but slow.
  set.seed(seed) # Control CV in internal folds
  meta_hal_candidates_sl <- Lrnr_sl$new(base_lrnrs, metalearner = metalrnr_hal_candidates)
  meta_hal_candidates_cv <- make_learner(Lrnr_cv, 
                                         meta_hal_candidates_sl,
                                         full_fit=TRUE)
  start_time <- proc.time() # start time
  set.seed(seed)
  meta_hal_candidates_cv_fit <- meta_hal_candidates_cv$train(task)
  
  runtime_sl_fit <- proc.time() - start_time # end time - start time = run time
  print(runtime_sl_fit)
  
  # 4. Find and use lambda_honest to contruct meta_hal_honest_cv
  meta_hal_candidates_cv_pred <- meta_hal_candidates_cv_fit$predict_fold(task, fold_number = "validation")
  
  # keep only the common lambdas
  if (!identical(nrow(meta_hal_candidates_cv_pred), nrow(task$data))) {
    lambdas <- lapply(1:length(task$folds), function(xxx) {
      (meta_hal_candidates_cv_fit$fit_object$fold_fits[[xxx]]$fit_object$cv_meta_fit$fit_object$lambda_star)  
    })
    
    common_lambdas <- lambda_list
    for (xxx in 1:length(task$folds)) {
      common_lambdas <- intersect(common_lambdas, lambdas[[xxx]])
    }
    
    lambda_list <- common_lambdas
    
    params_honestCV <- meta_hal_fast_cv_fit$params[["metalearner"]]$params
    
    metalrnr_hal_candidates <- Lrnr_hal9001$new(
      max_degree = params_honestCV$max_degree,
      return_lasso = FALSE,
      lambda = lambda_list,
      fit_control = list(type.measure = "mse",
                         cv_select = FALSE))
    
    set.seed(seed) # Control CV in internal folds
    meta_hal_candidates_sl <- Lrnr_sl$new(base_lrnrs, metalearner = metalrnr_hal_candidates)
    meta_hal_candidates_cv <- make_learner(Lrnr_cv, 
                                           meta_hal_candidates_sl,
                                           full_fit=TRUE)
    start_time <- proc.time() # start time
    set.seed(seed)
    meta_hal_candidates_cv_fit <- meta_hal_candidates_cv$train(task)
    
    runtime_sl_fit <- proc.time() - start_time # end time - start time = run time
    print(runtime_sl_fit)
    
    meta_hal_candidates_cv_pred <- meta_hal_candidates_cv_fit$predict_fold(task, fold_number = "validation")
  }
  
  
  lambda_cv_risks <- meta_hal_candidates_cv_fit$cv_risk(eval_fun = loss_squared_error)
  lambda_honest <- lambda_list[which.min(lambda_cv_risks$MSE)]
  metalrnr_hal_honest <- Lrnr_hal9001$new(
    max_degree = params_honestCV$max_degree,
    return_lasso = FALSE,
    lambda = lambda_honest, 
    fit_control = list(type.measure = "mse",
                       cv_select = FALSE))
  meta_hal_honest_cv <- Lrnr_sl$new(learners = base_lrnrs, 
                                    metalearner = metalrnr_hal_honest)
  set.seed(seed)
  meta_hal_honest_cv_fit <- meta_hal_honest_cv$train(task)
  
  # 5. fast M-HAL-SL selector
  cv_obj <- meta_hal_fast_cv_fit$fit_object$cv_fit
  outer_folds <- cv_obj$training_task$folds
  outer_fold_fits <- meta_hal_fast_cv_fit$fit_object$cv_fit$fit_object$fold_fits
  # In total V versions have been fitted, and they all can be applied to the data
  fold_preds <- lapply(outer_fold_fits, function(fit_i) {
    fit_i$predict(task)
  })
  
  meta_df     <- as.data.frame(meta_hal_fast_cv_fit$fit_object$cv_meta_task$data)
  outcome     <- "y"
  features    <- names(meta_hal_fast_cv_fit$fit_object$cv_meta_task$X)
  
  list_array <- lapply(1:length(outer_folds), function(v) {
    # v <- 1
    fold_v       <- outer_folds[[v]]
    train_idx_v  <- fold_v$training_set
    valid_idx_v  <- fold_v$validation_set
    
    list_mat <- lapply((1:length(fold_preds))[-v], function(fold_id) {
      fold_preds[[fold_id]][fold_v$validation_set, ]
    })
    
    # training fold data for HAL-MLE
    v_meta_data <- meta_df[train_idx_v, ]
    task_v_meta <- make_sl3_Task(
      data       = v_meta_data,
      outcome    = outcome,
      covariates = features
    )
    hal_candidate <- Lrnr_hal9001$new(
      max_degree   = params_honestCV$max_degree,
      return_lasso = TRUE,        # so we keep the glmnet object
      lambda       = lambda_list, # fix the full sequence
      fit_control  = list(
        type.measure = "mse",
        cv_select    = FALSE      # NO inner CV here
      )
    )
    # training the HAL learner
    hal_fit_v <- hal_candidate$train(task_v_meta)
    
    # predict with the v‐th HAL for averaging -> the v-th version of fast M-HAL-SL
    preds_by_mat <- lapply(list_mat, function(mat_i) {
      X_new_basis <- make_design_matrix(
        X     = as.matrix(mat_i),
        blist = hal_fit_v$fit_object$basis_list
      )
      # include any unpenalized cols if present
      if (!is.null(hal_fit_v$fit_object$X_unpenalized_new)) {
        X_new_basis <- cbind(
          X_new_basis,
          hal_fit_v$fit_object$X_unpenalized_new[valid_idx_v, , drop = FALSE]
        )
      }
      glmnet_obj <- hal_fit_v$fit_object$lasso_fit
      predict(glmnet_obj,
              newx = X_new_basis,
              s    = lambda_list,
              type = "response")
    })   
    avg_mat <- Reduce(`+`, preds_by_mat) / length(preds_by_mat)

    y_v       <- meta_df[[outcome]][valid_idx_v]
    err <- avg_mat - y_v
    
    err  # mse is used for the discrete super learner selecting C
  })
  err_mat <- do.call(rbind, list_array)
  err_by_lambda <- colMeans(err_mat^2)
  
  new_lambda <- lambda_list[which.min(err_by_lambda)]  # selected by fast M-HAL-SL
  
  metalrnr_hal_fastSL <- Lrnr_hal9001$new(
    max_degree = params_honestCV$max_degree,
    return_lasso = FALSE,
    lambda = new_lambda, 
    fit_control = list(type.measure = "mse",
                       cv_select = FALSE))
  meta_hal_fastSL_cv <- Lrnr_sl$new(learners = base_lrnrs, 
                                    metalearner = metalrnr_hal_fastSL)
  set.seed(seed)
  meta_hal_fastSL_cv_fit <- meta_hal_fastSL_cv$train(task)
  
  
  return(list(fast_cv = meta_hal_fast_cv_fit, honest_cv = meta_hal_honest_cv_fit, fastSL_cv = meta_hal_fastSL_cv_fit))
}

trackMetaCompetitor <- function(task, task_test, base_lrnrs, meta_lrnr, name){
  sl <- Lrnr_sl$new(learners = base_lrnrs,
                    metalearner = meta_lrnr)
  sl_fit <- sl$train(task)
  return (calculate_test_result(sl_fit, task_test, name))
}

runCompetition <- function(sim_name, task, task_test, base_lrnrs, seed = 2) {
  # results: SL competitors
  # 09/23/2023: 
  # Correct the configuration! convex = FALSE is nnls, convex = TRUE is convex!
  n <- nrow(task$data)
  meta_lrnr_convex <- Lrnr_nnls$new(convex = TRUE)
  meta_lrnr_nnls <- Lrnr_nnls$new(convex = FALSE)
  meta_lrnr_discrete <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
  
  competitors_results <-
    mapply(trackMetaCompetitor,
           meta_lrnr = c(meta_lrnr_convex, meta_lrnr_nnls, meta_lrnr_discrete),
           name = c("Convex", "NNLS", "Discrete SL"),
           MoreArgs = list(base_lrnrs = base_lrnrs,
                           task = task,
                           task_test = task_test),
           SIMPLIFY = FALSE)
  competitors_results <- data.table::rbindlist(competitors_results)
  
  print(competitors_results)
  
  meta_hal_fit_d1 <- runMetaHAL(task = task, 
                                base_lrnrs = base_lrnrs, 
                                max_degree = 1, seed = seed)
  meta_hal_fit_d2 <- runMetaHAL(task = task,
                                base_lrnrs = base_lrnrs,
                                max_degree = 2, seed = seed)
  
  
  base_lrnrs_fit <- meta_hal_fit_d1[["fast_cv"]]$fit_object$full_fit$learner_fits$Stack
  # base_lrnrs_fit_convex <- meta_lrnr_convex$fit_object$full_fit$learner_fits$Stack
  
  # results: Base Learner
  base_results <- calculate_test_result(base_lrnrs_fit, task_test, "Base Learner")
  
  # results: MetaHAL
  meta_hal_results <- 
    mapply(calculate_test_result, 
           fit_object = c(
             meta_hal_fit_d2[["honest_cv"]],
             meta_hal_fit_d2[["fastSL_cv"]], 
                          meta_hal_fit_d2[["fast_cv"]],
                          meta_hal_fit_d1[["honest_cv"]],
                          meta_hal_fit_d1[["fastSL_cv"]], 
                          meta_hal_fit_d1[["fast_cv"]]),
           learner_name = c(
             "Meta-HAL-d2 (Honest CV)",
             "Meta-HAL-d2 (Fast SL CV)", 
                            "Meta-HAL-d2 (Fast CV)",
                            "Meta-HAL-d1 (Honest CV)",
                            "Meta-HAL-d1 (Fast SL CV)", 
                            "Meta-HAL-d1 (Fast CV)"),
           MoreArgs = list(task_test = task_test),
           SIMPLIFY = FALSE)
  meta_hal_results <- data.table::rbindlist(meta_hal_results)
  
  
  results_all <- cbind(data.frame(sim_name = sim_name,
                                  n = n),
                       rbind(meta_hal_results, 
                             competitors_results,
                             base_results))
  return(results_all)
}

