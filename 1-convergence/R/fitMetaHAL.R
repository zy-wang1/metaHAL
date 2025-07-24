fitMetaHAL <- function(task, task_test, base_lrnrs) {
  # Get an interval of lambda candidates for metaHAL
  # base_HAL = sl3::make_learner(sl3::Lrnr_hal9001, max_degree = 1, return_lasso = TRUE)
  
  # base_HAL_fit = base_HAL$train(task)
  # lambda_hal = base_HAL_fit$fit_object$lambda_star
  # L1_hal = sum(abs(base_HAL_fit$fit_object$coefs))
  # cat('L1 norm of HAL estimator on all covariates:', 
  #     L1_hal, '. Lower bound of lambda:', lambda_l,'\n')
  
  lrnr_hal_d1 = sl3::make_learner(sl3::Lrnr_hal9001, 
                                  max_degree = 1, 
                                  keep = TRUE, 
                                  type.measure = "mse")
  meta_hal_dishonest_cv = sl3::Lrnr_sl$new(learners = base_lrnrs, 
                                           metalearner = lrnr_hal_d1)
  
  # Lower bound of L1, Upper bound of lambda
  meta_hal_dishonest_cv_fit = meta_hal_dishonest_cv$train(task)
  lambda_l = meta_hal_dishonest_cv_fit$metalearner_fit()$lambda_star
  # cat("Upper bound of L1 without intercept:",sum(abs(meta_hal_dishonest_cv_fit$metalearner_fit()$coefs)[-1]),"\n")
  L1_u = sum(abs(meta_hal_dishonest_cv_fit$metalearner_fit()$coefs))
  cat('Upper bound of L1 norm:',
      L1_u,'. Lower bound of lambda:',
      lambda_l,'\n')
  
  ################ Comment out meta HAL cv part here. Only try overfitted metaHAL #####
  ######  if (FALSE) {
  #n_lambda = 20
  lambda_list = meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$lambda
  lambda_intercept = meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$a0
  lambda_coef = as.matrix(meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$beta)
  L1_list = abs(lambda_intercept) + colSums(abs(lambda_coef))
  # lambda_dev_ratio = meta_hal_dishonest_cv_fit$metalearner_fit()$glmnet_lasso$dev.ratio
  # plot(lambda_list, L1_list)
  # plot(lambda_list, lambda_dev_ratio)
  # plot(L1_list, lambda_dev_ratio, 'l')
  ## Generate metaHAL candidates
  # lambda_list_metahal = lambda_list[lambda_list >= lambda_l]
  # L1_list_metahal = L1_list[lambda_list >= lambda_l]
  meta_hal_candidates = sl3::make_learner(sl3::Lrnr_hal9001,
                                          max_degree = 1,
                                          lambda = lambda_list,
                                          cv_select = FALSE)
  # meta_hal_candidates_all = sl3::make_learner(sl3::Lrnr_hal9001,
  #                             max_degree = 1,
  #                             lambda = lambda_list,
  #                             cv_select = FALSE)
  
  base_step = make_learner(Lrnr_cv, 
                           make_learner(Stack, base_lrnrs), 
                           full_fit=TRUE)
  base_hal_pipe = make_learner(Pipeline,
                               base_step,
                               meta_hal_candidates)
  meta_hal_candidates_cv = make_learner(Lrnr_cv, 
                                        base_hal_pipe,
                                        full_fit=TRUE)
  meta_hal_candidates_cv_fit = meta_hal_candidates_cv$train(task)
  meta_hal_candidates_cv_fit_preds = meta_hal_candidates_cv_fit$predict()
  mses <- apply(meta_hal_candidates_cv_fit_preds, 2, 
                function(x,y){mean((x-y)^2)},task$Y)
  
  lambda_star = lambda_list[which(mses==min(mses))]
  L1_star = colSums(abs(meta_hal_candidates_cv_fit$fit_object$full_fit$learner_fits[[2]]$fit_object$coefs))[which(min(mses) == mses)]
  cat("The best meta-HAL has a lambda =", lambda_star, "and L1 norm =", L1_star, "\n")
  try(plot(lambda_list, mses, pch = 16, cex = 0.5))
  try(abline(v = lambda_star, col = "red"))
  try(plot(L1_list, mses, pch = 16, cex = 0.5))
  try(abline(v = L1_star, col = "red"))
  
  
  if (lambda_star == lambda_l){
    meta_hal_honest_cv = meta_hal_dishonest_cv
    meta_hal_honest_cv_fit = meta_hal_dishonest_cv_fit
  } else {
    lrnr_hal_d1_star = sl3::make_learner(sl3::Lrnr_hal9001, 
                                         max_degree = 1, 
                                         lambda = c(lambda_star),
                                         cv_select = FALSE,
                                         type.measure = "mse")
    meta_hal_honest_cv = sl3::Lrnr_sl$new(learners = base_lrnrs, 
                                          metalearner = lrnr_hal_d1_star)
    meta_hal_honest_cv_fit = meta_hal_honest_cv$train(task)
  }
  
  # convex super learner
  sl_convex = sl3::Lrnr_sl$new(learners = base_lrnrs)
  sl_convex_fit = sl_convex$train(task)
  
  # nnls super learner
  meta_nnls = sl3::make_learner(Lrnr_nnls, convex = TRUE)
  sl_nnls = sl3::Lrnr_sl$new(learners = base_lrnrs,
                             metalearner = meta_nnls)
  sl_nnls_fit = sl_nnls$train(task)
  
  # discrete super learner
  meta_discrete = make_learner(Lrnr_cv_selector) 
  sl_discrete = make_learner(Lrnr_sl, base_lrnrs, meta_discrete)
  sl_discrete_fit = sl_discrete$train(task)
  
  meta_hal_honest_cv_test_pred = meta_hal_honest_cv_fit$predict(task_test)
  meta_hal_dishonest_cv_test_pred = meta_hal_dishonest_cv_fit$predict(task_test)
  # basic_HAL_test_pred = base_HAL_fit$predict(task_test)
  convex_test_pred = sl_convex_fit$predict(task_test)
  nnls_test_pred = sl_nnls_fit$predict(task_test)
  discrete_test_pred = sl_discrete_fit$predict(task_test)
  
  base_lrnrs_test_pred = sl_convex_fit$fit_object$full_fit$fit_object$learner_fits$Stack$predict(task_test)
  average_test_pred = metalearner_linear(X = as.matrix(base_lrnrs_test_pred), alpha = rep(1/length(base_lrnrs),length(base_lrnrs)))
  
  # calculate mse on test dataset
  calculate_test_mse = function(fit){mean((fit-task_test$Y)^2)}
  calculate_se = function(pred){
    loss = (pred - task_test$Y)^2
    n = length(loss)
    return (sqrt(1/n)*sd(loss))
  }
  
  test_preds = list(meta_lrnrs_pred = 
                      list("metaHAL_Honest_CV" = meta_hal_honest_cv_test_pred, 
                           "metaHAL_Dishonest_CV" = meta_hal_dishonest_cv_test_pred,
                           "Convex" = convex_test_pred,
                           "NNLS" = nnls_test_pred,
                           "Average" = average_test_pred,
                           "Discrete" = discrete_test_pred),
                    base_lrnrs_pred = base_lrnrs_test_pred)
  
  summaryResult = data.frame(
    metalearner = c("metaHAL_Honest_CV", 
                    "metaHAL_Dishonest_CV", 
                    "Convex", 
                    "NNLS", 
                    "Average", 
                    "Discrete",
                    names(base_lrnrs_test_pred)),
    test_mse = c(unlist(sapply(list(meta_hal_honest_cv_test_pred, 
                                    meta_hal_dishonest_cv_test_pred,
                                    convex_test_pred,
                                    nnls_test_pred,
                                    average_test_pred,
                                    discrete_test_pred),
                               FUN = calculate_test_mse)),
                 sapply(base_lrnrs_test_pred, calculate_test_mse)),
    se = c(unlist(sapply(list(meta_hal_honest_cv_test_pred, 
                              meta_hal_dishonest_cv_test_pred,
                              convex_test_pred,
                              nnls_test_pred,
                              average_test_pred,
                              discrete_test_pred),
                         FUN = calculate_se)),
           sapply(base_lrnrs_test_pred, calculate_se)))
  result = list(lambda_l = lambda_l,
                L1_u = L1_u,
                lambda_list = lambda_list,
                L1_list = L1_list,
                lambda_star = lambda_star, 
                L1_star = L1_star, 
                mses = mses,
                base_lrnrs = base_lrnrs,
                task_test_Y = task_test$Y,
                test_preds = test_preds, 
                summaryResult = summaryResult)
  return(result)
}