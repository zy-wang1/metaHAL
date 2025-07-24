library(dplyr)
library(stringr)
library(here)
library(tmle3)
library(sl3)
library(data.table)
library(parallel)
setwd(here::here())
library(R6)
library(xlsx)
library(glmnet)
library(sl3)
library(Matrix)

n_iter <- 500  # number of iterations
set.seed(123)
vec_seed <- sample(10^9, n_iter)  # generate seed

doc_path <- rstudioapi::getActiveDocumentContext()$path 
file_name <- basename(doc_path)
n_cores <- 1  # set to >1 to parallel
if (!dir.exists(paste0("2-plugin/Results/", file_name))) paste0("2-plugin/Results/", file_name) %>% dir.create

source("2-plugin/R/Lrnr_glmnet_overfit_sparse.R")
source("2-plugin/R/helpers.R")
source("2-plugin/R/Lrnr_screener_name.R")
source("2-plugin/R/Lrnr_identity.R")
source("2-plugin/R/Lrnr_glm_insert.R")
source("2-plugin/R/make_likelihood_Wc_only.R")
source("2-plugin/R/Lrnr_nnet_value.R")

p_cov <- 4  # clinical cov dimension
sample_size <- 200  # n
p_image <- 4  # additional cov dimension; 1 * 4 here for low-d
ratio_overfit = 1  # no overfit in glm
ratio_nonzero <- 1  # lowD setting
truth <- trt_coef <- 1

# To regenerate outputs: delete ./plugin/Results 
list_result <- mclapply(1:n_iter, mc.cores = n_cores, FUN = function(i_iter) {
  
  coef_A <- 0.2
  int_A <- - round(p_image * ratio_nonzero) * coef_A * 0.5
  coef_Y <- 0.6
  
  seed <- vec_seed[i_iter]
  data <- generate_data(n = sample_size, p_cov = p_cov, p_image = p_image, ratio_nonzero = ratio_nonzero, seed = seed, trt_coef = trt_coef, if_simpleG = F
                                 , coef_A = coef_A, coef_Y = coef_Y, if_binary = T, int_A = int_A
  ) %>% as.data.frame

  # get true functions
  p_nonzero <- round(ratio_nonzero * p_image)
  G0 <- expit(rowSums(data[, 1:(p_cov + p_nonzero)] * coef_A) + int_A)
  Q1 <- rowSums(data[, 1:(p_cov + p_nonzero)] * coef_Y) + trt_coef
  Q0 <- rowSums(data[, 1:(p_cov + p_nonzero)] * coef_Y)
  IC0 <- data$trt / G0 * (data$Y - Q1) - 
    (1 - data$trt) / (1 - G0) * (data$Y - Q0)
  n <- nrow(data)
  var_est_0 <- var(IC0) / n
  
  output_target <- paste0("2-plugin/Results/", file_name, "/", "n_", sample_size, "_pimage_", p_image, "_ratio_nonzero_", ratio_nonzero, "_ratio_overfit_", ratio_overfit, "_raw_results_", i_iter, ".rds")
  if (file.exists(output_target)) {
    record <- readRDS(output_target)
  } else {
    
    record <- list()
    
    # completely no adjustment
    est_noadj <- mean(data$Y[data$trt == 1]) - mean(data$Y[data$trt == 0])
    record$noadj <- c(est_noadj, NA, NA, NA)  # no update, no CI
    
    # full-data tmle
    {
      data_use <- data
      # setup
      {
        node_list <- list(
          W = names(data_use) %>% head(-2),
          A = "trt",
          Y = "Y"
        )
        # choose base learners
        lrnr_glm <- make_learner(Lrnr_glm)
        lrnr_G <- lrnr_glm
        lrnr_Q <- lrnr_glm
      }
      {
          learner_list <- list(A = lrnr_G, Y = lrnr_Q)
        
        ate_spec <- tmle_ATE(
          treatment_level = 1,
          control_level = 0
        )
        {
          tmle_spec <- ate_spec
          start_time <- proc.time()
          set.seed(seed)  # a random draw of folds will follow
          tmle_task <- tmle_spec$make_tmle_task(data_use, node_list
                                                # , folds = folds_backup
          )
          task_time <- proc.time()
          initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, 
                                                                  learner_list)
          likelihood_time <- proc.time()
          updater <- tmle3_Update$new(cvtmle = F)
          targeted_likelihood <- tmle_spec$make_targeted_likelihood(initial_likelihood, 
                                                                    updater)
          tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
          updater$tmle_params <- tmle_params
          params_time <- proc.time()
          fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, 
                           updater)
          fit_time <- proc.time()
          fit$set_timings(start_time, task_time, likelihood_time, params_time, 
                          fit_time)
        }
        fit_summary_backup <- fit$summary
        folds_backup <- tmle_task$folds
      }
    }
    fit_summary_backup
    fit$timings
    record$tmle <- c(fit_summary_backup$init_est, fit_summary_backup$tmle_est, fit_summary_backup$lower, fit_summary_backup$upper, 
                     fit_summary_backup$se^2, fit$ED)
    

    # metaHAL initial plugin
    {
      data_use <- data
      
      # setup
      {
        node_list <- list(
          W = names(data_use) %>% head(-2),
          A = "trt",
          Y = "Y"
        )
        # choose base learners
        lrnr_glm <- make_learner(Lrnr_glm)
        
        lrnr_1 <- make_learner(Pipeline, Lrnr_screener_name$new(var_name = paste0("X", 1)), 
                               Lrnr_identity$new())  
        lrnr_2 <- make_learner(Pipeline, Lrnr_screener_name$new(var_name = paste0("X", 2)), 
                               Lrnr_identity$new())  
        lrnr_3 <- make_learner(Pipeline, Lrnr_screener_name$new(var_name = paste0("X", 3)), 
                               Lrnr_identity$new())  
        lrnr_4 <- make_learner(Pipeline, Lrnr_screener_name$new(var_name = paste0("X", 4)), 
                               Lrnr_identity$new())  
        lrnr_trt <- make_learner(Pipeline, Lrnr_screener_name$new(var_name = paste0("trt")), 
                                 Lrnr_identity$new())
        
        lrnr_Q_1 <- Lrnr_glm_insert$new(trt_name = "trt", trt_value = 1)
        lrnr_Q_0 <- Lrnr_glm_insert$new(trt_name = "trt", trt_value = 0)
        
        lrnr_G <- Lrnr_glm_insert_G$new(force_DV = "trt")  # family binomial is forced in the customized learner
        
          base_lrnrs <- list(lrnr_1, lrnr_2, lrnr_3, lrnr_4, lrnr_trt, lrnr_Q_1, lrnr_Q_0
                             , lrnr_G
          )
        
        lrnr_hal_d1 = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0)
        meta_hal_dishonest_cv_d1 = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1)
        
        task = sl3::make_sl3_Task(data = data_use,
                                  covariates = c(node_list$W, node_list$A
                                  ),
                                  folds = folds_backup,
                                  outcome = "Y",
                                  outcome_type = "continuous")
        meta_hal_dishonest_cv_d1_fit = meta_hal_dishonest_cv_d1$train(task)
        lambda_l_d1 = meta_hal_dishonest_cv_d1_fit$metalearner_fit()$lambda_star

        data1 <- data_use
        data1$trt <- 1
        task1 = sl3::make_sl3_Task(data = data1,
                                   covariates = c(node_list$W, node_list$A
                                   ),
                                   folds = folds_backup,
                                   outcome = "Y",
                                   outcome_type = "continuous")
        
        data0 <- data_use
        data0$trt <- 0
        task0 = sl3::make_sl3_Task(data = data0,
                                   covariates = c(node_list$W, node_list$A
                                   ),
                                   folds = folds_backup,
                                   outcome = "Y",
                                   outcome_type = "continuous")
        
      }
    }
    record$meta_init <- c(meta_hal_dishonest_cv_d1_fit$predict(task1) %>% mean - meta_hal_dishonest_cv_d1_fit$predict(task0) %>% mean, 
                          var_est_0
    )
    
    # undersmooth with metaG
    {
      # fit metaG
      {
        lrnr_treatment <- lrnr_glm
        
          base_lrnrs_metaG <- list(lrnr_1, lrnr_2, lrnr_3, lrnr_4 
                                   , lrnr_treatment
          )
        
        lrnr_hal_d1_metaG = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0)
        meta_hal_dishonest_cv_d1_metaG = sl3::Lrnr_sl$new(learners = base_lrnrs_metaG, metalearner = lrnr_hal_d1_metaG)
        
        task_metaG = sl3::make_sl3_Task(data = data_use,
                                        covariates = c(node_list$W
                                        ),
                                        folds = folds_backup,
                                        outcome = "trt",
                                        outcome_type = "binomial")
        meta_hal_dishonest_cv_d1_fit_metaG = meta_hal_dishonest_cv_d1_metaG$train(task_metaG)
      }
      
      {
        current_step_size <- 1/2
        current_lambda <- lambda_l_d1
        current_step <- 1
        record_lambda <- c()
        record_lambda[current_step] <- current_lambda
        vec_if_good_lambda <- c()
        vec_if_good_lambda[current_step] <- T
        vec_eq <- c()
        vec_threshold <- c()

        meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_fit

        IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
          (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
        
        n <- nrow(data_use)
        se_Dstar <- sqrt(var(IC) / n)
        var_est_init <- se_Dstar^2  # track
        ED_threshold <- se_Dstar / min(log(n), 10)
        ED_init <- mean(IC)
        
        record_est_ci <- list()
        temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
        record_est_ci[[current_step]] <- c(
          temp_est, 
          temp_est - 1.96 * se_Dstar, 
          temp_est + 1.96 * se_Dstar, 
          current_lambda
        )
        
        vec_eq[current_step] <- abs(mean(IC))
        vec_threshold[current_step] <- ED_threshold
        
        max_solve_step <- 15
        while(current_step <= max_solve_step) {
          if (abs(mean(IC)) <= ED_threshold) break()
          current_step <- current_step + 1
          current_lambda <- current_lambda * (1 - current_step_size)
          record_lambda[current_step] <- current_lambda
          vec_if_good_lambda[current_step] <- T  # consider it as good for now
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
          
          vec_eq[current_step] <- abs(mean(IC))
          vec_threshold[current_step] <- ED_threshold
          
          temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
          record_est_ci[[current_step]] <- c(
            temp_est, 
            temp_est - 1.96 * se_Dstar, 
            temp_est + 1.96 * se_Dstar, 
            current_lambda
          )
          
          if (abs(mean(IC)) <= ED_threshold) {
            # consider record this successful lambda
            break()
          }  
        }
      }
      
      IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold
      

      record$meta_undersmoothed_1_metaG <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                             var_est_init, 
                                             var_est_update, ED_init, ED_update, if_update_solved
      )
      
      {
        if (length(which(vec_if_good_lambda)) > 0) {
          solver_choice_lambda <- current_reduce_lambda <- record_lambda[which.min(vec_eq[unique(c(1, which(vec_if_good_lambda)))])]
        } else {
          solver_choice_lambda <- current_reduce_lambda <- current_lambda
        }
        
        # try smaller steps
        if_solved <- (vec_eq < vec_threshold)[last(which(vec_if_good_lambda))]  # see if it has been solved
        if (length(if_solved) == 0) if_solved <- F  # if haven't been solved yet in the end, then it is not solved
        if_last_two <- length(which(vec_if_good_lambda)) >= 2
        if (if_solved & if_last_two) {
          current_reduce_step <- 1
          loc_last_two <- which(vec_if_good_lambda) %>% tail(2)
          previous_reduce_lambda <- record_lambda[loc_last_two][1]  # the last good lambda not solving it
          list_reduce_lambda <- c()
          
          max_reduce_step <- 15
          while (current_reduce_step <= max_reduce_step) {
            current_reduce_lambda <- mean(c(previous_reduce_lambda, current_reduce_lambda))  # try if somewhere in between can still solve it
            list_reduce_lambda[current_reduce_step] <- current_reduce_lambda
            
            # see if the equation is still being solved here
            {
              lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
              meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
              meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
              
              IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
                (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
              n <- nrow(data_use)
              se_Dstar <- sqrt(var(IC) / n)
              ED_threshold <- se_Dstar / min(log(n), 10)
            }
            
            # if still solved, continue
            if (abs(mean(IC)) <= ED_threshold) {
              current_reduce_step <- current_reduce_step + 1
              current_reduce_lambda <- last(list_reduce_lambda)
            } else {  # if not, good sign, the last one is good enough
              current_reduce_lambda <- ifelse(length(list_reduce_lambda) >= 2, 
                                              list_reduce_lambda[current_reduce_step - 1], 
                                              solver_choice_lambda  # revert to the previous choice, if the first step immediately stops solving it
              )
              break  # o.w. stop here
            } 
          }
        }
        
        {
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
        }
      }
      
      IC <- data_use$trt / meta_hal_dishonest_cv_d1_fit_metaG$predict(task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - meta_hal_dishonest_cv_d1_fit_metaG$predict(task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))                
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold
      
      # try smaller steps
      record$meta_undersmoothed_2_metaG <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                             var_est_init, 
                                             var_est_update, ED_init, ED_update, if_update_solved
      )
      
    }
    
    # undersmooth with known G
    {
      {
        current_step_size <- 1/2
        current_lambda <- lambda_l_d1
        current_step <- 1
        record_lambda <- c()
        record_lambda[current_step] <- current_lambda
        vec_if_good_lambda <- c()
        vec_if_good_lambda[current_step] <- T
        vec_eq <- c()
        vec_threshold <- c()
        
        
        meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_fit
        p_nonzero <- round(ratio_nonzero * p_image)
        # clinical columns, and then nonzero img columns
        G0 <- expit(rowSums(task1$data[, 1:(p_cov + p_nonzero)] * coef_A) + int_A)
        
        IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
          (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
        n <- nrow(data_use)
        se_Dstar <- sqrt(var(IC) / n)
        var_est_init <- se_Dstar^2  # track
        ED_threshold <- se_Dstar / min(log(n), 10)
        ED_init <- mean(IC)
        
        record_est_ci <- list()
        temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
        record_est_ci[[current_step]] <- c(
          temp_est, 
          temp_est - 1.96 * se_Dstar, 
          temp_est + 1.96 * se_Dstar, 
          current_lambda
        )
        
        vec_eq[current_step] <- abs(mean(IC))
        vec_threshold[current_step] <- ED_threshold
        
        max_solve_step <- 15
        while(current_step <= max_solve_step) {
          if (abs(mean(IC)) <= ED_threshold) break()
          current_step <- current_step + 1
          current_lambda <- current_lambda * (1 - current_step_size)
          record_lambda[current_step] <- current_lambda
          vec_if_good_lambda[current_step] <- T  # consider it as good for now
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
          
          vec_eq[current_step] <- abs(mean(IC))
          vec_threshold[current_step] <- ED_threshold
          
          temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
          record_est_ci[[current_step]] <- c(
            temp_est, 
            temp_est - 1.96 * se_Dstar, 
            temp_est + 1.96 * se_Dstar, 
            current_lambda
          )
          
          if (abs(mean(IC)) <= ED_threshold) {
            # consider record this successful lambda
            break()
          }  
        }
      }
      
      IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold

      record$meta_undersmoothed_1 <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                       var_est_init, 
                                       var_est_update, ED_init, ED_update, if_update_solved
      )
      
      {
        if (length(which(vec_if_good_lambda)) > 0) {
          solver_choice_lambda <- current_reduce_lambda <- record_lambda[which.min(vec_eq[unique(c(1, which(vec_if_good_lambda)))])]
        } else {
          solver_choice_lambda <- current_reduce_lambda <- current_lambda
        }
        
        # try smaller steps
        if_solved <- (vec_eq < vec_threshold)[last(which(vec_if_good_lambda))]  # see if it has been solved
        if (length(if_solved) == 0) if_solved <- F  # if haven't been solved yet in the end, then it is not solved
        if_last_two <- length(which(vec_if_good_lambda)) >= 2
        if (if_solved & if_last_two) {
          current_reduce_step <- 1
          loc_last_two <- which(vec_if_good_lambda) %>% tail(2)
          previous_reduce_lambda <- record_lambda[loc_last_two][1]  # the last good lambda not solving it
          list_reduce_lambda <- c()
          
          max_reduce_step <- 15
          while (current_reduce_step <= max_reduce_step) {
            current_reduce_lambda <- mean(c(previous_reduce_lambda, current_reduce_lambda))  # try if somewhere in between can still solve it
            list_reduce_lambda[current_reduce_step] <- current_reduce_lambda
            
            # see if the equation is still being solved here
            {
              lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
              meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
              meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
              
              IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
                (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
              n <- nrow(data_use)
              se_Dstar <- sqrt(var(IC) / n)
              ED_threshold <- se_Dstar / min(log(n), 10)
            }
            
            # if still solved, continue
            if (abs(mean(IC)) <= ED_threshold) {
              current_reduce_step <- current_reduce_step + 1
              current_reduce_lambda <- last(list_reduce_lambda)
            } else {  # if not, good sign, the last one is good enough
              current_reduce_lambda <- ifelse(length(list_reduce_lambda) >= 2, 
                                              list_reduce_lambda[current_reduce_step - 1], 
                                              solver_choice_lambda  # revert to the previous choice, if the first step immediately stops solving it
              )
              break  # o.w. stop here
            } 
            
            
          }
        }
        
        
        {
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
        }
      }
      
      IC <- data_use$trt / G0 * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - G0) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold
      
      # try smaller steps
      record$meta_undersmoothed_2 <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                       var_est_init, 
                                       var_est_update, ED_init, ED_update, if_update_solved
      )
    }
    
    # undersmooth with overfitted initial G
    {
      {
        current_step_size <- 1/2
        current_lambda <- lambda_l_d1
        current_step <- 1
        record_lambda <- c()
        record_lambda[current_step] <- current_lambda
        vec_if_good_lambda <- c()
        vec_if_good_lambda[current_step] <- T
        vec_eq <- c()
        vec_threshold <- c()
        
        
        meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_fit
        
        IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
          (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
        
        n <- nrow(data_use)
        se_Dstar <- sqrt(var(IC) / n)
        var_est_init <- se_Dstar^2  # track
        ED_threshold <- se_Dstar / min(log(n), 10)
        ED_init <- mean(IC)
        
        record_est_ci <- list()
        temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
        record_est_ci[[current_step]] <- c(
          temp_est, 
          temp_est - 1.96 * se_Dstar, 
          temp_est + 1.96 * se_Dstar, 
          current_lambda
        )
        
        vec_eq[current_step] <- abs(mean(IC))
        vec_threshold[current_step] <- ED_threshold
        
        max_solve_step <- 15
        while(current_step <= max_solve_step) {
          if (abs(mean(IC)) <= ED_threshold) break()
          current_step <- current_step + 1
          current_lambda <- current_lambda * (1 - current_step_size)
          record_lambda[current_step] <- current_lambda
          vec_if_good_lambda[current_step] <- T  # consider it as good for now
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
          
          vec_eq[current_step] <- abs(mean(IC))
          vec_threshold[current_step] <- ED_threshold
          
          temp_est <- mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))              
          record_est_ci[[current_step]] <- c(
            temp_est, 
            temp_est - 1.96 * se_Dstar, 
            temp_est + 1.96 * se_Dstar, 
            current_lambda
          )
          
          if (abs(mean(IC)) <= ED_threshold) {
            # consider record this successful lambda
            break()
          }  
        }
      }
      
      IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold
      
      record$meta_undersmoothed_1_overfitG <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                                var_est_init, 
                                                var_est_update, ED_init, ED_update, if_update_solved
      )
      
      
      {
        if (length(which(vec_if_good_lambda)) > 0) {
          solver_choice_lambda <- current_reduce_lambda <- record_lambda[which.min(vec_eq[unique(c(1, which(vec_if_good_lambda)))])]
        } else {
          solver_choice_lambda <- current_reduce_lambda <- current_lambda
        }
        
        # try smaller steps
        if_solved <- (vec_eq < vec_threshold)[last(which(vec_if_good_lambda))]  # see if it has been solved
        if (length(if_solved) == 0) if_solved <- F  # if haven't been solved yet in the end, then it is not solved
        if_last_two <- length(which(vec_if_good_lambda)) >= 2
        if (if_solved & if_last_two) {
          current_reduce_step <- 1
          loc_last_two <- which(vec_if_good_lambda) %>% tail(2)
          previous_reduce_lambda <- record_lambda[loc_last_two][1]  # the last good lambda not solving it
          list_reduce_lambda <- c()
          
          max_reduce_step <- 15
          while (current_reduce_step <= max_reduce_step) {
            current_reduce_lambda <- mean(c(previous_reduce_lambda, current_reduce_lambda))  # try if somewhere in between can still solve it
            list_reduce_lambda[current_reduce_step] <- current_reduce_lambda
            
            # see if the equation is still being solved here
            {
              lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
              meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
              meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
              
              IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
                (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
              n <- nrow(data_use)
              se_Dstar <- sqrt(var(IC) / n)
              ED_threshold <- se_Dstar / min(log(n), 10)
            }
            
            # if still solved, continue
            if (abs(mean(IC)) <= ED_threshold) {
              current_reduce_step <- current_reduce_step + 1
              current_reduce_lambda <- last(list_reduce_lambda)
            } else {  # if not, good sign, the last one is good enough
              current_reduce_lambda <- ifelse(length(list_reduce_lambda) >= 2, 
                                              list_reduce_lambda[current_reduce_step - 1], 
                                              solver_choice_lambda  # revert to the previous choice, if the first step immediately stops solving it
              )
              break  # o.w. stop here
            } 
         }
        }
        
        {
          lrnr_hal_d1_undersmoothed = Lrnr_hal9001$new(max_degree = 2, smoothness_orders = 0, lambda = current_reduce_lambda, fit_control = list(cv_select = F))
          meta_hal_dishonest_cv_d1_undersmoothed = sl3::Lrnr_sl$new(learners = base_lrnrs, metalearner = lrnr_hal_d1_undersmoothed)
          meta_hal_dishonest_cv_d1_undersmoothed_fit <- meta_hal_dishonest_cv_d1_undersmoothed$train(task)
          
          IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
            (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
          n <- nrow(data_use)
          se_Dstar <- sqrt(var(IC) / n)
          ED_threshold <- se_Dstar / min(log(n), 10)
        }
      }
      
      IC <- data_use$trt / initial_likelihood$factor_list$A$learner$predict(task = task1) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - 
        (1 - data_use$trt) / (1 - initial_likelihood$factor_list$A$learner$predict(task = task1)) * (data_use$Y - meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0))
      n <- nrow(data_use)
      se_Dstar <- sqrt(var(IC) / n)
      var_est_update <- se_Dstar^2  # track
      ED_update <- mean(IC)
      ED_threshold <- se_Dstar / min(log(n), 10)
      if_update_solved <- abs(ED_update) < ED_threshold
      
      # try smaller steps
      record$meta_undersmoothed_2_overfitG <- c(mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task1)) - mean(meta_hal_dishonest_cv_d1_undersmoothed_fit$predict(task0)), 
                                                var_est_init, 
                                                var_est_update, ED_init, ED_update, if_update_solved
      )
    }
    
    record %>% saveRDS(output_target)
    print("This iteration is now finished")
  }
  
  return(record)
})


