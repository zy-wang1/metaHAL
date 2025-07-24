makeTask <- function(data, Vfold = 10){
  covariates = paste0("x",seq(1:20))
  outcome = "y"
  outcome_type = "continuous"
  task = sl3::make_sl3_Task(data = data$data,
                            covariates = data$covariates,
                            folds = origami::make_folds(data$data, 
                                                        fold_fun = origami::folds_vfold,
                                                        V = Vfold),
                                outcome = data$outcome,
                                outcome_type = data$outcome_type)
  return (task)
}

generateBaselrnrs <- function(covariates, 
                              type = c("CART", "OLS", "Multiple"),
                              subset = TRUE,
                              n_subsets = 4,
                              order = TRUE){
  if (type == "Multiple"){
    base_lrnrs <- list(Lrnr_mean = sl3::make_learner(sl3::Lrnr_mean),
                               Lrnr_glm = sl3::make_learner(sl3::Lrnr_glm),
                               Lrnr_xgboost = make_learner(Lrnr_xgboost),
                               Lrnr_svm = make_learner(Lrnr_svm),
                               Lrnr_rf = make_learner(Lrnr_ranger))
    
  } else{
    base_lrnr_for_subset <- switch(type, 
                                   CART = make_learner(Lrnr_rpart),
                                   OLS = make_learner(Lrnr_glm))
    if (subset){
      chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
      # four groups
      if (order){
        sub_covariates <- chunk(seq(1:length(covariates)), n_subsets) 
      } else{
        sub_covariates <- split(seq(1:length(covariates)), 1:n_subsets)
      }
      base_lrnrs <- list()
      for (index in sub_covariates){
        subset <- covariates[index]
        lrnr_subset <- Lrnr_subset_covariates$new(covariates = subset)
        learner <- make_learner(Pipeline, lrnr_subset, base_lrnr_for_subset)
        base_lrnrs[[paste0(subset,collapse ="+")]] <- learner
      }
    } else{
      base_lrnrs <- base_lrnr_for_subset
    }
  }
  return (base_lrnrs)
}
