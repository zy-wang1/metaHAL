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

runMetaHAL <- function(task, base_lrnrs, max_degree){
  # 0. Substantiate metaHAL Super Learner with fast CV
  lrnr_hal <- Lrnr_hal9001$new(max_degree = max_degree,
                               return_lasso = TRUE,
                               fit_control = list(type.measure = "mse",
                                                  cv_select = TRUE))
  meta_hal_fast_cv <- Lrnr_sl$new(learners = base_lrnrs, 
                                  metalearner = lrnr_hal)
  set.seed(2)
  meta_hal_fast_cv_fit <- meta_hal_fast_cv$train(task)
  
  # 1. Get candidate lambdas
  lambda_star_fast <- meta_hal_fast_cv_fit$metalearner_fit()$lambda_star
  metalearner_fit <- meta_hal_fast_cv_fit$metalearner_fit()$lasso_fit$glmnet.fit
  lambda_list <- metalearner_fit$lambda
  
  # 2. Pass HAL Params & Generate a meta-HAL for every lambda
  params_honestCV <- meta_hal_fast_cv_fit$params[["metalearner"]]$params

  metalrnr_hal_candidates <- Lrnr_hal9001$new(
    max_degree = params_honestCV$max_degree,
    return_lasso = FALSE,
    lambda = lambda_list,
    fit_control = list(type.measure = "mse",
                       cv_select = FALSE))
  
  # 3. Try Honest MetaHAL by discrete SL. It works but slow.
  set.seed(2) # Control CV in internal folds
  meta_hal_candidates_sl <- Lrnr_sl$new(base_lrnrs, metalearner = metalrnr_hal_candidates)
  meta_hal_candidates_cv <- make_learner(Lrnr_cv, 
                                        meta_hal_candidates_sl,
                                        full_fit=TRUE)
  start_time <- proc.time() # start time
  set.seed(2)
  meta_hal_candidates_cv_fit <- meta_hal_candidates_cv$train(task)
  
  runtime_sl_fit <- proc.time() - start_time # end time - start time = run time
  print(runtime_sl_fit)
  
  # 4. Find and use lambda_honest to contruct meta_hal_honest_cv
  meta_hal_candidates_cv_pred <- meta_hal_candidates_cv_fit$predict_fold(task, fold_number = "validation")
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
  set.seed(2)
  meta_hal_honest_cv_fit <- meta_hal_honest_cv$train(task)
  return(list(fast_cv = meta_hal_fast_cv_fit, honest_cv = meta_hal_honest_cv_fit))
}

trackMetaCompetitor <- function(task, task_test, base_lrnrs, meta_lrnr, name){
  sl <- Lrnr_sl$new(learners = base_lrnrs,
                    metalearner = meta_lrnr)
  sl_fit <- sl$train(task)
  return (calculate_test_result(sl_fit, task_test, name))
}
  
runCompetition <- function(sim_name, task, task_test, base_lrnrs) {
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
                                max_degree = 1)
  meta_hal_fit_d2 <- runMetaHAL(task = task, 
                                base_lrnrs = base_lrnrs, 
                                max_degree = 2)
  base_lrnrs_fit <- meta_hal_fit_d1[["fast_cv"]]$fit_object$full_fit$learner_fits$Stack
  # base_lrnrs_fit_convex <- meta_lrnr_convex$fit_object$full_fit$learner_fits$Stack
  
  # results: Base Learner
  base_results <- calculate_test_result(base_lrnrs_fit, task_test, "Base Learner")
  
  # results: MetaHAL
  meta_hal_results <- 
    mapply(calculate_test_result, 
           fit_object = c(meta_hal_fit_d2[["honest_cv"]],
                          meta_hal_fit_d2[["fast_cv"]],
                          meta_hal_fit_d1[["honest_cv"]],
                          meta_hal_fit_d1[["fast_cv"]]),
           learner_name = c("Meta-HAL-d2 (Honest CV)",
                            "Meta-HAL-d2 (Fast CV)",
                            "Meta-HAL-d1 (Honest CV)",
                            "Meta-HAL-d1 (Fast CV)"),
           MoreArgs = list(task_test = task_test),
           SIMPLIFY = FALSE)
  meta_hal_results <- data.table::rbindlist(meta_hal_results)
  
  
  results_all <- cbind(data.frame(sim_name = sim_name,
                                  n = n),
                       rbind(meta_hal_results, competitors_results, base_results))
  return(results_all)
}