doHonestCV <- function(meta_hal_dishonest_cv_fit){
  task = meta_hal_dishonest_cv_fit$training_task
  lambda_list = meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$lambda
  lambda_intercept = meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$a0
  lambda_coef = as.matrix(meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$beta)
  L1_list = abs(lambda_intercept) + colSums(abs(lambda_coef))
  
  # Set parameters to do honest CV for metaHAL with lambda candidates
  params_honestCV = meta_hal_dishonest_cv_fit$params[["metalearner"]]$params
  params_honestCV$cv_select = FALSE
  params_honestCV$lambda = lambda_list
  
  meta_hal_candidates = do.call(sl3::Lrnr_hal9001$new, params_honestCV)
  
  base_lrnrs = meta_hal_dishonest_cv_fit$params$learners
  # base_step = make_learner(Lrnr_cv, 
  #                          make_learner(Stack, base_lrnrs), 
  #                          full_fit=TRUE)
  # base_step = make_learner(Lrnr_sl,
  #                          make_learner(Stack, base_lrnrs),
  #                          full_fit=TRUE)
  # 
  # base_hal_pipe = make_learner(Pipeline,
  #                              base_step,
  #                              meta_hal_candidates)
  meta_hal_candidates_sl = make_learner(Lrnr_sl,
                               base_lrnrs,
                               meta_hal_candidates)
  #meta_hal_candidates_external sl = make_learner(Lrnr_sl,meta_hal_candidates_sl,cv_selctor)
  meta_hal_candidates_cv = make_learner(Lrnr_cv, 
                                        meta_hal_candidates_sl,
                                        full_fit=TRUE)
  
  # fix randomness of the internal folds in nested cross validation
  set.seed(2)
  meta_hal_candidates_cv_fit = meta_hal_candidates_cv$train(task)
  
  # But I can not repeat the following code without set.seed()
  # Each time the predictions are different, especially for small lambdas
  meta_hal_candidates_cv_fit_preds = meta_hal_candidates_cv_fit$predict()
  # meta_hal_candidates_cv_fit_preds_2 = meta_hal_candidates_cv_fit$predict()
  # plot(meta_hal_candidates_cv_fit_preds[,100],meta_hal_candidates_cv_fit_preds_2[,100])
  mses <- apply(meta_hal_candidates_cv_fit_preds, 2, 
                function(x,y){mean((x-y)^2)},task$Y)
  
  lambda_star = lambda_list[which(mses==min(mses))]
  #L1_star = colSums(abs(meta_hal_candidates_cv_fit$fit_object$full_fit$learner_fits[[2]]$fit_object$coefs))[which(min(mses) == mses)]
  L1_star = colSums(abs(meta_hal_candidates_cv_fit$fit_object$full_fit$fit_object$full_fit$fit_object$learner_fits[[2]]$fit_object$coefs))[which(min(mses) == mses)]
  
  
  cat("The best meta-HAL has a lambda:", lambda_star, "and L1 norm:", L1_star, "\n")
  max_degree = params_honestCV$max_degree
  try(plot(lambda_list, mses, pch = 16, cex = 0.5, main = paste("metaHAL candidates, max_degree =", max_degree)))
  try(abline(v = lambda_star, col = "red"))
  try(plot(L1_list, mses, pch = 16, cex = 0.5, main = paste("metaHAL candidates, max_degree =", max_degree)))
  try(abline(v = L1_star, col = "red"))
  
  # Pick the best lambda
  params_honestCV$lambda = lambda_star
  lrnr_hal_honest_d1 = do.call(sl3::Lrnr_hal9001$new, params_honestCV)
  meta_hal_honest_cv = sl3::Lrnr_sl$new(learners = base_lrnrs, 
                                        metalearner = lrnr_hal_honest_d1)
  
  return (list(meta_hal_honest_cv = meta_hal_honest_cv, 
               mses = mses, 
               lambda_list = lambda_list,
               L1_list = L1_list))
}

library(zeallot)
fitMetaHAL <- function(task, task_test, base_lrnrs) {
  lrnr_hal_d1 = sl3::make_learner(sl3::Lrnr_hal9001, 
                                  max_degree = 1, 
                                  keep = TRUE, 
                                  type.measure = "mse")
  lrnr_hal_d2 = sl3::make_learner(sl3::Lrnr_hal9001, 
                                  max_degree = 2, 
                                  keep = TRUE, 
                                  type.measure = "mse")
  
  # metaHAL, max_degree = 1
  meta_hal_dishonest_cv_d1 = sl3::Lrnr_sl$new(learners = base_lrnrs, 
                                           metalearner = lrnr_hal_d1)
  
  ## Fix randomness when cv.glmnet is making folds.
  set.seed(2)
  meta_hal_dishonest_cv_d1_fit = meta_hal_dishonest_cv_d1$train(task)
  lambda_l_d1 = meta_hal_dishonest_cv_d1_fit$metalearner_fit()$lambda_star
  L1_u_d1 = sum(abs(meta_hal_dishonest_cv_d1_fit$metalearner_fit()$coefs))
  cat('metaHAL_d1:\n Upper bound of L1 norm:', L1_u_d1,'\n Lower bound of lambda:', lambda_l_d1,'\n')
  zeallot::`%<-%`(c(meta_hal_honest_cv_d1, 
                    meta_hal_honest_cv_d1_mses,
                    lambda_list_d1,
                    L1_list_d1),
                  doHonestCV(meta_hal_dishonest_cv_d1_fit))
  meta_hal_honest_cv_d1_fit = meta_hal_honest_cv_d1$train(task)
  
  # metaHAL, max_degree = 2
  meta_hal_dishonest_cv_d2 = sl3::Lrnr_sl$new(learners = base_lrnrs, 
                                              metalearner = lrnr_hal_d2)
  set.seed(2) 
  meta_hal_dishonest_cv_d2_fit = meta_hal_dishonest_cv_d2$train(task)
  lambda_l_d2 = meta_hal_dishonest_cv_d2_fit$metalearner_fit()$lambda_star
  L1_u_d2 = sum(abs(meta_hal_dishonest_cv_d2_fit$metalearner_fit()$coefs))
  cat('metaHAL_d2:\n Upper bound of L1 norm:', L1_u_d2,'. Lower bound of lambda:', lambda_l_d2,'\n')
  zeallot::`%<-%`(c(meta_hal_honest_cv_d2, 
                    meta_hal_honest_cv_d2_mses,
                    lambda_list_d2,
                    L1_list_d2),
                  doHonestCV(meta_hal_dishonest_cv_d2_fit))
  meta_hal_honest_cv_d2_fit = meta_hal_honest_cv_d2$train(task)
  
  # convex super learner
  sl_convex = sl3::Lrnr_sl$new(learners = base_lrnrs)
  set.seed(2)
  sl_convex_fit = sl_convex$train(task)
  # set.seed(2) # Fix the internal folds of nested Lrnr_cv
  # cv_sl_convex = make_learner(Lrnr_cv,
  #                             sl_convex_fit,
  #                             full_fit=TRUE)
  # cv_sl_convex_fit = cv_sl_convex$train(task)
  # cv_sl_convex_mse = cv_sl_convex_fit$cv_risk(loss_squared_error)$risk
  
  # nnls super learner
  set.seed(2)
  meta_nnls = sl3::make_learner(Lrnr_nnls, convex = TRUE)
  sl_nnls = sl3::Lrnr_sl$new(learners = base_lrnrs,
                             metalearner = meta_nnls)
  sl_nnls_fit = sl_nnls$train(task)

  # set.seed(2)
  # cv_sl_nnls = make_learner(Lrnr_cv,
  #                             sl_nnls_fit,
  #                             full_fit=TRUE)
  # cv_sl_nnls_fit = cv_sl_nnls$train(task)
  # cv_sl_nnls_mse = cv_sl_nnls_fit$cv_risk(loss_squared_error)$risk
  
  
  # average super learner
  set.seed(2)
  meta_average = sl3::make_learner(Lrnr_average)
  sl_average = sl3::Lrnr_sl$new(learners = base_lrnrs,
                             metalearner = meta_average)
  sl_average_fit = sl_average$train(task)
  # set.seed(2)
  # cv_sl_average = make_learner(Lrnr_cv,
  #                              sl_average_fit,
  #                           full_fit=TRUE)
  # cv_sl_average_fit = cv_sl_average$train(task)
  # cv_sl_average_mse = cv_sl_average_fit$cv_risk(loss_squared_error)$risk
  

  #  discrete super learner
  # meta_discrete = make_learner(Lrnr_cv_selector) 
  # sl_discrete = make_learner(Lrnr_cv, 
  #                            c(#base_lrnrs, 
  #                              meta_hal_honest_cv_d2, 
  #                              meta_hal_dishonest_cv_d2, 
  #                              meta_hal_honest_cv_d1, 
  #                              meta_hal_dishonest_cv_d1, 
  #                              sl_convex,
  #                              sl_nnls,
  #                              sl_average),
  #                              metalearner = meta_discrete)
  # sl_discrete_fit = sl_discrete$train(task)
  
  meta_hal_honest_cv_d1_test_pred = meta_hal_honest_cv_d1_fit$predict(task_test)
  meta_hal_dishonest_cv_d1_test_pred = meta_hal_dishonest_cv_d1_fit$predict(task_test)
  meta_hal_honest_cv_d2_test_pred = meta_hal_honest_cv_d2_fit$predict(task_test)
  meta_hal_dishonest_cv_d2_test_pred = meta_hal_dishonest_cv_d2_fit$predict(task_test)
  
  # basic_HAL_test_pred = base_HAL_fit$predict(task_test)
  convex_test_pred = sl_convex_fit$predict(task_test)
  nnls_test_pred = sl_nnls_fit$predict(task_test)
  average_test_pred = sl_average_fit$predict(task_test)
  # discrete_test_pred = sl_discrete_fit$predict(task_test)
  # which_discrete = which.min(c(min(meta_hal_honest_cv_d2_mses),
  #                              min(meta_hal_honest_cv_d1_mses),
  #                              cv_sl_convex_mse,
  #                              cv_sl_nnls_mse,
  #                              cv_sl_average_mse
  #                              ))
  # discrete_test_pred = list(meta_hal_honest_cv_d2_test_pred,
  #                        meta_hal_honest_cv_d1_test_pred,
  #                        convex_test_pred,
  #                        nnls_test_pred,
  #                        average_test_pred
  #                        )[[which_discrete]]

  base_lrnrs_test_pred = sl_convex_fit$fit_object$full_fit$fit_object$learner_fits$Stack$predict(task_test)
  # average_test_pred = metalearner_linear(X = as.matrix(base_lrnrs_test_pred), alpha = rep(1/length(base_lrnrs),length(base_lrnrs)))
  
  # calculate mse on test dataset
  calculate_test_mse = function(fit){mean((fit-task_test$Y)^2)}
  calculate_se = function(pred){
    loss = (pred - task_test$Y)^2
    n = length(loss)
    return (sqrt(1/n)*sd(loss))
  }
  
  test_preds = list(meta_lrnrs_pred = 
                      list("metaHAL_Honest_CV_d2" = meta_hal_honest_cv_d2_test_pred, 
                           "metaHAL_Dishonest_CV_d2" = meta_hal_dishonest_cv_d2_test_pred,
                           "metaHAL_Honest_CV_d1" = meta_hal_honest_cv_d1_test_pred, 
                           "metaHAL_Dishonest_CV_d1" = meta_hal_dishonest_cv_d1_test_pred,
                           "Convex" = convex_test_pred,
                           "NNLS" = nnls_test_pred,
                           "Average" = average_test_pred
                           # ,
                           # "Discrete" = discrete_test_pred
                           ),
                    base_lrnrs_pred = base_lrnrs_test_pred)
  
  summaryResult = data.frame(
    metalearner = c("metaHAL_Honest_CV_d2", 
                    "metaHAL_Dishonest_CV_d2", 
                    "metaHAL_Honest_CV_d1", 
                    "metaHAL_Dishonest_CV_d1", 
                    "Convex", 
                    "NNLS", 
                    "Average", 
                    # "Discrete",
                    names(base_lrnrs_test_pred)),
    test_mse = c(unlist(sapply(list(meta_hal_honest_cv_d2_test_pred, 
                                    meta_hal_dishonest_cv_d2_test_pred,
                                    meta_hal_honest_cv_d1_test_pred, 
                                    meta_hal_dishonest_cv_d1_test_pred,
                                    convex_test_pred,
                                    nnls_test_pred,
                                    average_test_pred
                                    # ,
                                    # discrete_test_pred
                                    ),
                               FUN = calculate_test_mse)),
                 sapply(base_lrnrs_test_pred, calculate_test_mse)),
    se = c(unlist(sapply(list(meta_hal_honest_cv_d2_test_pred, 
                              meta_hal_dishonest_cv_d2_test_pred,
                              meta_hal_honest_cv_d1_test_pred, 
                              meta_hal_dishonest_cv_d1_test_pred,
                              convex_test_pred,
                              nnls_test_pred,
                              average_test_pred
                              # ,
                              # discrete_test_pred
                              ),
                         FUN = calculate_se)),
           sapply(base_lrnrs_test_pred, calculate_se)))
  result = list(lambda_l_d1 = lambda_l_d1,
                L1_u_d1 = L1_u_d1,
                lambda_list_d1 = lambda_list_d1,
                L1_list_d1 = L1_list_d1,
                lambda_star_d1 = meta_hal_honest_cv_d1_fit$metalearner_fit()$lambda, 
                meta_hal_honest_cv_d1_mses = meta_hal_honest_cv_d1_mses,
                lambda_l_d2 = lambda_l_d2,
                L1_u_d2 = L1_u_d2,
                lambda_list_d2 = lambda_list_d2,
                L1_list_d2 = L1_list_d2,
                lambda_star_d2 = meta_hal_honest_cv_d2_fit$metalearner_fit()$lambda, 
                meta_hal_honest_cv_d2_mses = meta_hal_honest_cv_d2_mses,
                base_lrnrs = base_lrnrs,
                task_test_Y = task_test$Y,
                test_preds = test_preds, 
                summaryResult = summaryResult)
  return(result)
}


Lrnr_average <- R6::R6Class(
  classname = "Lrnr_cv_selector",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(loss_function = loss_squared_error,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c(
      "continuous", "binomial", "categorical", "weights"
    ),
    
    .train = function(task) {
      verbose <- getOption("sl3.verbose")
      params <- self$params
      loss_function <- params$loss_function
      outcome_type <- self$get_outcome_type(task)
      
      # specify data
      X <- as.matrix(task$X)
      Y <- outcome_type$format(task$Y)
      weights <- task$weights
      
      risk <- function(preds) {
        loss <- loss_function(preds, Y)
        risk <- weighted.mean(loss, weights)
        return(risk)
      }
      risks <- apply(X, 2, risk)
      
      fit_object <- list()
      fit_object$name <- colnames(task$X)
      # Set coefficients to be equal
      coef <- rep(1/ncol(X),ncol(X))
      fit_object$coefficients <- as.numeric(coef)
      return(fit_object)
    },
    
    .predict = function(task = NULL) {
      verbose <- getOption("sl3.verbose")
      X <- as.matrix(task$X)
      predictions <- X[, self$fit_object$name] %*% private$.fit_object$coefficients
      return(predictions)
    }
  )
)
