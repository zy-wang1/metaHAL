library(dplyr)
library(oro.nifti)  # readNIfTI
library(Matrix)
library(glmnet)
library(biglasso)
library(bigmemory)
library(bigmemoryExt)  # cbindBM
library(origami)  # folds_vfold
library(hal9001)  # fit_hal
library(nnet)
library(parallel)
source("3-imaging_application/R/mediation_helper.R")

# scenarios
# if_true_effect <- F  # T will use real features; F will use noise

for (if_true_effect in c(T, F)) {

# computation setup
n_cores = 1  # set to >1 to parallel
n_iter <- 100

feature_sets <- list(10, 18, 34 , 50, 101, 152, 200,
  c(10, 18, 34, 50, 101, 152, 200)
)
length(feature_sets)

# imaging features
data_dir <- "data/MediationPain/"
mask <- readNIfTI(file.path(data_dir, "brainmask.nii")) %>% as.vector
label_exclude <- c("IE", "EXP", "SCEBL")  # excluded study labels; psychological interventions and potentially different pathways
data_trial_total <- read.csv(file.path(data_dir, "data/nifti/trialdata.csv"))
data_trial_total <- data_trial_total[!(data_trial_total$studyLabel %in% label_exclude), ]
cut_values <- rep(median(data_trial_total$temperature), 2)
data_trial_total$A <- (data_trial_total$temperature > cut_values[2]) * 1  # use binary temperature intervention: if temperature high
data_trial_total$Y <- data_trial_total$rate101/100  # normalize subjective pain ratings
n_trial_total <- nrow(data_trial_total)

# make seed sequence
set.seed(123)
vec_seeds <- sample(10E5, 10000)

doc_path <- rstudioapi::getActiveDocumentContext()$path 
file_name <- basename(doc_path)
task_label <- paste0(file_name, "_", if_true_effect)
if (!dir.exists("3-imaging_application/Results")) dir.create("3-imaging_application/Results")
if (!dir.exists(paste0("3-imaging_application/Results/", task_label))) dir.create(paste0("3-imaging_application/Results/", task_label))

# set parameters here
N <- 10000  # sample size of each iteration
V <- 5  # V fold cross-validation

for (loc_feature in seq_along(feature_sets)) {
  feature_set <- feature_sets[[loc_feature]]
  feature_label <- ifelse(length(feature_set) == 1, feature_set, "ensemble")
  vec_exist <- sapply(1:n_iter, function(i_iter) {
    output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", if_true_effect, "_feature_", feature_label, "_", i_iter, ".rds")  
    file.exists(output_target)
  })
  
  list_result <- mclapply((1:n_iter)[!vec_exist], mc.cores = n_cores, FUN = function(i_iter) {
    output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", if_true_effect, "_feature_", feature_label, "_", i_iter, ".rds")  
    
    SEED <- vec_seeds[i_iter]  # random seed for generating subsamples
    set.seed(SEED)
    
    n_trial <- N
    loc_trial <- sample(n_trial_total, n_trial, replace = T)  # bootstrap sample mimicing a true distribution as the empirical distribution of training samples
    data_trial <- data_trial_total[loc_trial, ]
    
    # load meta covariates
    list_raw_meta_data <- list()
    for (algorithm in paste0("feature_nnet_MN_", feature_set
    )) {
      data_nifti <- list()
      
      for (i in 1:n_trial) {
        data_nifti[[i]] <-       (paste0(
          "data/extracted_features_", 
          strsplit(algorithm, "_")[[1]][-(1:2)] %>% paste0(collapse = "_"), 
          "/study0", data_trial$studyNum[i], "/", 
          "subj", data_trial$subjectNum[i] %>% sprintf(fmt = "%02d"), "_", 
          "trial", data_trial$trial[i] %>% sprintf(fmt = "%02d"), ".nii", ".csv"
        ) %>% read.csv()) [, 1]
      }
      mat_nifti <- do.call(rbind, data_nifti)
      rm(data_nifti)
      
      # force completely random values for fake signals
      if (!if_true_effect) {
        mat_nifti <- rnorm(nrow(mat_nifti) * ncol(mat_nifti)) %>% matrix(nrow = n_trial)
      }
      
      rownames(mat_nifti) <- NULL
      colnames(mat_nifti) <- paste0("Z_", 1:ncol(mat_nifti))
      mat_W <- model.matrix(~ studyLabel - 1, data = data_trial)[, -1] %>% as.matrix
      ncol_W <- ncol(mat_W)
      data_WZ <- cbind(
        mat_W, 
        mat_nifti
      )
      data_WAZ <- cbind(
        mat_W, 
        data.frame(A = data_trial$A) %>% as.matrix, 
        mat_nifti
      )
      data_W1Z <- data_WAZ
      data_W1Z[, ncol_W + 1] <- 1  # double check column loc of A
      
      data_WA <- data.frame(
        mat_W, 
        A = data_trial$A)
      data_WAY <- data.frame(
        mat_W, 
        A = data_trial$A, Y = data_trial$Y)
      data_W1Y <- data_W0Y <- data_WAY
      data_W1Y[, "A"] <- 1
      data_W0Y[, "A"] <- 0
      if (grep("_enet", algorithm) %>% length == 1) {
        data_WZ <- data_WZ %>% funMakeItSparse
        data_WAZ <- data_WAZ %>% funMakeItSparse
        data_W1Z <- data_W1Z %>% funMakeItSparse
      }
      
      # create meta-level dataset
      set.seed(SEED)
      folds <- folds_vfold(nrow(data_trial), V = V)
      
      # Wr, A, Zr, Y; 6 dimensions in total
      # Wr: A given W, QZ of W, A = 0
      # Zr: A given Z and W, QY of Z, W, A = 1
      list_meta_data <- lapply(1:V, function(i) {
        # v-th training set
        i_train <- folds[[i]]$training_set
        i_validation <- folds[[i]]$validation_set
        train_model_A_W <- glm(A ~ ., data_WA[i_train, ], family = "binomial")
        train_model_A_WZ <- nnet(A ~ ., data = data.frame(data_WZ[i_train, ] %>% as.matrix, A = as.factor(data_trial$A)[i_train]), size = 3, maxit = 100, MaxNWts = 10000)
        train_model_Y_WAZ <- nnet(data_WAZ[i_train, ], data_trial$Y[i_train], size = 3, maxit = 100, MaxNWts = 10000
        )
        # get training sample QY
        train_QY <- predict(train_model_Y_WAZ, newx = data_W1Z[i_train, ])
        train_model_QY_WA <- lm(QY ~ ., data = data.frame(QY = as.vector(train_QY), data_WA[i_train, ]))
        
        pred_validation_EA <- predict(object = train_model_A_W, newdata = data_WAY[i_validation, ], type = "response")
        pred_validation_EA_WZ <- predict(train_model_A_WZ, newdata = data_WZ[i_validation, ], type = "raw")  # this gives the conditional mean estimates; conditional density depends on A
        # pred_validation_EA_WZ <- (1 - pred_validation_EA_WZ) / pred_validation_EA_WZ
        pred_validation_QY <- predict(train_model_Y_WAZ, newdata = data_W1Z[i_validation, ])
        pred_validation_QZ <- predict(train_model_QY_WA, newdata = data_W0Y[i_validation, ])
        
        temp <- data.frame(EA = pred_validation_EA %>% as.vector, QZ = pred_validation_QZ %>% as.vector, 
                           A = data_trial$A[i_validation], 
                           EA_WZ = pred_validation_EA_WZ %>% as.vector, QY = pred_validation_QY %>% as.vector, 
                           Y = data_trial$Y[i_validation]
        )
        return(temp)  # folds[[i]]$validation_set are the corresponding rows
      })
      meta_data <- list_meta_data %>% do.call(what = rbind)
      meta_data <- meta_data[lapply(folds, function(fold) fold$validation_set) %>% unlist %>% order, ]  # based on the permutation, find the way to move it back to the original order
      rownames(meta_data) <- NULL
      
      list_raw_meta_data[[algorithm]] <- meta_data
    }
    
    # TMLE on ATE using full-data models (same for all methods)
    {
      model_Y_WA <- lm(Y ~ ., data = data_WAY)
      model_A_W <- glm(A ~ ., data_WA, family = "binomial")
      
      EA <- predict(object = model_A_W, newdata = data_WAY, type = "response")  # prob(A = 1 | W)
      TMLE_Y1 <- funTMLE(observed = data_trial$Y, 
                         initial = predict(model_Y_WA, newdata = data_WAY) %>% bound(0.005), 
                         H = data_trial$A / EA
      )
      TMLE_Y0 <- funTMLE(observed = data_trial$Y, 
                         initial = predict(model_Y_WA, newdata = data_WAY) %>% bound(0.005), 
                         H = (1 - data_trial$A) / (1 - EA)
      )
      
      
      TMLE_Y1_updated <- submodel(TMLE_Y1$epsilon, predict(model_Y_WA, newdata = data_W1Y) %>% bound(0.005), 
                                  matrix(1 / EA)
      )
      TMLE_Y0_updated <- submodel(TMLE_Y0$epsilon, predict(model_Y_WA, newdata = data_W0Y) %>% bound(0.005), 
                                  matrix((1 - 0) / (1 - EA))
      )
    }
    
    
    
    
    if (length(feature_set) == 1) {
      
      # simple glm and regular lasso as initial models    
      model_A_WZ <- nnet(A ~ ., data = data.frame(data_WZ %>% as.matrix, A = as.factor(data_trial$A)), size = 3, maxit = 100, MaxNWts = 10000)
      model_Y_WAZ <- nnet(data_WAZ, data_trial$Y, size = 3, maxit = 100, MaxNWts = 10000, linout = TRUE)
      
      # seq reg
      {
        if (class(model_Y_WAZ) == "nnet") {
          QY <- predict(model_Y_WAZ, data_W1Z) 
        } else {
          QY <- predict(model_Y_WAZ, newx = data_W1Z, s = "lambda.min") 
        }
        
        model_QY_WA <- lm(QY ~ ., data = data.frame(QY = as.vector(QY), data_WA))
        
        est_vec_seq_reg <- c(
          predict(object = model_Y_WA, newdata = data_W0Y) %>% mean, 
          predict(model_QY_WA, newdata = data_W0Y) %>% mean, 
          predict(object = model_Y_WA, newdata = data_W1Y) %>% mean
        )
        est_vec_seq_reg %>% getPercent
      }
      
      # seq reg TMLE
      {
        # TMLE on ATE
        EA <- predict(object = model_A_W, newdata = data_WAY, type = "response")  # prob(A = 1 | W)
        TMLE_Y1 <- funTMLE(observed = data_trial$Y, 
                           initial = predict(model_Y_WA, newdata = data_WAY) %>% bound(0.005), 
                           H = data_trial$A / EA
        )
        TMLE_Y0 <- funTMLE(observed = data_trial$Y, 
                           initial = predict(model_Y_WA, newdata = data_WAY) %>% bound(0.005), 
                           H = (1 - data_trial$A) / (1 - EA)
        )
        
        
        TMLE_Y1_updated <- submodel(TMLE_Y1$epsilon, predict(model_Y_WA, newdata = data_W1Y) %>% bound(0.005), 
                                    matrix(1 / EA)
        )
        TMLE_Y0_updated <- submodel(TMLE_Y0$epsilon, predict(model_Y_WA, newdata = data_W0Y) %>% bound(0.005), 
                                    matrix((1 - 0) / (1 - EA))
        )
        
        # TMLE on NIE
        temp <- predict(model_A_WZ, newx = data_WZ, type = "raw") %>% bound(0.001)  # this gives the conditional mean estimates; conditional density depends on A
        EA <- predict(object = model_A_W, newdata = data_WAY, type = "response")
        TMLE_QY <- funTMLE(observed = data_trial$Y, 
                           initial =  predict(model_Y_WAZ, newx = data_WAZ) %>% as.vector %>% bound(0.005),
                           H = data_trial$A/(1 - EA)*(1-temp)/temp
        )
        QY_updated <- submodel(TMLE_QY$epsilon, predict(model_Y_WAZ, newx = data_W1Z) %>% as.vector %>% bound(0.005), 
                               matrix(1/(1 - EA)*(1-temp)/temp)
        )
        data_WAQY_updated <- data.frame(data_WA, QY_updated = as.vector(QY_updated))
        model_QY_WA_TMLE_initial <- lm(QY_updated ~ ., data_WAQY_updated)
        TMLE_QZ <- funTMLE(observed = QY_updated %>% as.vector, 
                           initial =  predict(model_QY_WA_TMLE_initial, newdata = data_WAY) %>% as.vector %>% bound(0.005),
                           H = (1 - data_trial$A)/(1 - EA)
        )
        QZ_updated <- submodel(TMLE_QZ$epsilon, predict(model_QY_WA_TMLE_initial, newdata = data_W0Y) %>% as.vector %>% bound(0.005), 
                               matrix((1 - 0)/(1 - EA))
        )
        
        est_vec_seq_reg_TMLE <- c(
          TMLE_Y0_updated %>% mean,    
          QZ_updated %>% mean, 
          TMLE_Y1_updated %>% mean 
        )
        est_vec_seq_reg_TMLE %>% getPercent
      }
      
      # seq reg IPW
      {
        est_vec_seq_reg_IPW <- c(
          mean(data_trial$Y * (1 - data_trial$A) / (1 - EA)), 
          mean(data_trial$Y * data_trial$A / EA * (1 - temp)/temp), 
          mean(data_trial$Y * data_trial$A / EA)
        )
        est_vec_seq_reg_IPW %>% getPercent
      }
    }
    
    if (length(feature_set) > 1) {
      est_vec_seq_reg <- est_vec_seq_reg_TMLE <- est_vec_seq_reg_IPW <- c(NA, NA, NA)
    }
    
    # ensemble metaHAL
    {
      ensemble_data <- data.frame(
        mat_W, 
        lapply(list_raw_meta_data, function(mat) {
          mat[, 1:2]
        }) %>% do.call(what = cbind), 
        A = data_WAY$A, 
        lapply(list_raw_meta_data, function(mat) {
          mat[, 4:5]
        }) %>% do.call(what = cbind), 
        Y = data_WAY$Y
      )
      
      meta_data0 <- meta_data1 <- meta_data <- ensemble_data
      meta_data0$A <- 0
      meta_data1$A <- 1

      loc_var_A <- which(colnames(meta_data) == "A")
      loc_var_Y <- which(colnames(meta_data) == "Y")
      
      # undersmoothing by equation solving
      {
        meta_model_Y_WAZ_candidate <- fit_hal(X = meta_data[, 1:(loc_var_Y - 1)] %>% as.matrix, Y = meta_data$Y, 
                                              # max_degree = 1,  # use default; 2 way for more than 20 dimensions
                                              smoothness_orders = 1, fit_control = list(cv_select = FALSE))
        meta_model_A_W <- glm(A ~ ., meta_data[, 1:loc_var_A], family = "binomial")  # W is simple study label info; use glm here
        meta_model_A_WZ <- fit_hal(X = meta_data[, c(1:(loc_var_A-1), (loc_var_A+1):(loc_var_Y-1))], Y = meta_data$A, 
                                   # max_degree = 3,  # use default for now
                                   smoothness_orders = 1, family = "binomial")      
        # QY update
        {
          temp <- predict(meta_model_A_WZ, new_data = meta_data[, c(1:(loc_var_A-1), (loc_var_A+1):(loc_var_Y-1))], type = "response")
          meta_EA <- predict(object = meta_model_A_W, meta_data[, 1:(loc_var_A-1)], type = "response")
          
          # initial is a matrix now
          candidates <- predict(meta_model_Y_WAZ_candidate, new_data = meta_data[, 1:(loc_var_Y - 1)] %>% as.matrix) %>% bound(0.001)
          loc_select <- funUS(observed = meta_data$Y,
                              initial =  candidates,
                              H = (meta_data$A/(1 - meta_EA)*(1-temp)/temp))
          updates <- predict(meta_model_Y_WAZ_candidate, new_data = meta_data1[, 1:(loc_var_Y - 1)] %>% as.matrix) %>% bound(0.001)
          meta_QY_updated <- updates[, loc_select]
        }
        
        # QZ update
        {
          meta_model_QY_WA_TMLE_candidate <- fit_hal(meta_data[, 1:loc_var_A], meta_QY_updated,
                                                     smoothness_orders = 1, fit_control = list(cv_select = FALSE))
          candidates <- predict(meta_model_QY_WA_TMLE_candidate, new_data = meta_data[, 1:loc_var_A])
          loc_select <- funUS(observed = meta_QY_updated,
                              initial =  candidates,
                              H = (1 - meta_data$A)/(1 - meta_EA)
          )
          updates <- predict(meta_model_QY_WA_TMLE_candidate, new_data = meta_data0[, 1:loc_var_A])
          meta_QZ_updated <- updates[, loc_select]
        }
        meta_est_vec_seq_reg_US <- c(
          TMLE_Y0_updated %>% mean,    
          meta_QZ_updated %>% mean(na.rm = TRUE), 
          TMLE_Y1_updated %>% mean 
        )
        meta_est_vec_seq_reg_US %>% getPercent
      }
    }
    
    # metaHAL plug-in
    {
      # meta_data0 <- meta_data1 <- meta_data <- ensemble_data
      # meta_data0$A <- 0
      # meta_data1$A <- 1
      # 
      # loc_var_A <- which(colnames(meta_data) == "A")
      # loc_var_Y <- which(colnames(meta_data) == "Y")
      meta_model_Y_WAZ <- fit_hal(X = meta_data[, 1:(loc_var_Y - 1)] %>% as.matrix, Y = meta_data$Y,
                                  # max_degree = 3,  # use default for now
                                  smoothness_orders = 1)
      meta_QY <- predict(meta_model_Y_WAZ, meta_data1[, 1:(loc_var_Y-1)])
      meta_model_QY_WA <- fit_hal(meta_data[, 1:(loc_var_A)], meta_QY, 
                                  smoothness_orders = 1)
      meta_QZ <- predict(meta_model_QY_WA, new_data = meta_data0[, 1:(loc_var_A)])
      
      meta_est_vec_seq_reg <- c(
        predict(object = model_Y_WA, newdata = data_W0Y) %>% mean, 
        meta_QZ %>% mean(na.rm = TRUE), 
        predict(object = model_Y_WA, newdata = data_W1Y) %>% mean
      )
      meta_est_vec_seq_reg %>% getPercent
    }
    
    result <- 
      c(
        est_vec_seq_reg %>% getPercent,
        est_vec_seq_reg_TMLE %>% getPercent,
        est_vec_seq_reg_IPW %>% getPercent,
        meta_est_vec_seq_reg %>% getPercent, 
        meta_est_vec_seq_reg_US %>% getPercent
      )
    result
    
    saveRDS(result, 
            output_target
    )
  })
}

}