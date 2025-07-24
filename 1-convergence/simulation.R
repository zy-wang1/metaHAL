library(parallel)
library(sl3)
library(future)
library(xgboost)
library(nnls)
source("1-convergence/R/DGD_updated.R")
source("1-convergence/R/helpers_updated.R")
source("1-convergence/R/fitMetaHAL_updated.R")
source("1-convergence/R/sim_config_updated.R")
library(hal9001)
source("1-convergence/R/fitMetaHAL_addCV.R")

# i_scenario = 1
for (i_scenario in c(1, 3, 2, 6)) {
sim_name <- names(sim_data_lrnr_pairs)[i_scenario]
# the following i_scenario and sim_name are needed
# 1: scenario_1_1
# 3: scenario_2_1
# 2: scenario_1_2
# 6: scenario_2_4

simData <- sim_data_lrnr_pairs[[sim_name]]$simData
base_lrnrs <- sim_data_lrnr_pairs[[sim_name]]$base_lrnrs
results_dir <- file.path("convergence", "Results",sim_name)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = T)

set.seed(123)
seed_vec_10k <- sample(10^7, 10^4, replace = F)

nCores <- 1  # optional: nCores>1 to run multiple threads in parallel
B <- 100
seed_vec <- seed_vec_10k[1:B]

list_results <- mclapply(X = 1:B, mc.cores = nCores, FUN = function(rep_id) {
  seed <- seed_vec[rep_id]
  
  # Generate data
  set.seed(seed)
  data_train_200 = simData(200)
  data_train_500 = simData(500)
  data_train_1000 = simData(1000)
  data_test = simData(5000)
  data_train_2000 = simData(2000)
  
  # Make Tasks
  task_200 <- makeTask(data_train_200)
  task_500 <- makeTask(data_train_500)
  task_1000 <- makeTask(data_train_1000)
  task_test <- makeTask(data_test)
  task_2000 <- makeTask(data_train_2000)
  
  task_test <- sl3::make_sl3_Task(data = data_test$data,
                                  covariates = data_test$covariates,
                                  folds = origami::make_folds(data_test$data,
                                                              fold_fun = origami::folds_vfold,
                                                              V = 10),
                                  outcome = data_test$outcome,
                                  outcome_type = data_test$outcome_type)
  
  results_200 <- runCompetition(sim_name, task_200, task_test, base_lrnrs, seed = seed)  
  results_500 <- runCompetition(sim_name, task_500, task_test, base_lrnrs, seed = seed)
  results_1000 <- runCompetition(sim_name, task_1000, task_test, base_lrnrs, seed = seed)
  results_2000 <- runCompetition(sim_name, task_2000, task_test, base_lrnrs, seed = seed)
  
  rep_output <- list(results_200, results_500, results_1000, results_2000)
  saveRDS(rep_output, file.path(paste0(results_dir, "/z_", rep_id, ".rds" )))
  
  return(rep_output)  
})
}