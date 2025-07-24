expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x/(1-x))
use_bounds <- function(x, bounds) {
  x[x > bounds[2]] <- bounds[2]
  x[x < bounds[1]] <- bounds[1]
  return(x)
}

# generate data partition as a list
get_partition <- function(n, V = 10, seed = NULL) {
  # n <- 101
  # V <- 11
  # seed <- 123
  
  if (!is.null(seed)) set.seed(seed)
  
  f <- rep(1:V, each = n/V)
  if (n%%V != 0) f <- c(f, 1:(n%%V))  # assign the leftover
  cv_list <- base::split(sample(1:n), f)
  names(cv_list) <- 1:V
  
  return(cv_list)
}
# examples
# get_partition(100, 10)
# sapply(get_partition(100, 13), length)

get_v_th_training_ids <- function(partition, v) {
  which_fold_val <- which(names(partition) %in% v)  # character or numerical v
  sort(as.vector(unlist(partition[-which_fold_val])))
}

# library(data.table)
# get Q_{n, v_{i}} for each subject according to a cv partition
# currently for WAY, with three dimensional Qnv
get_meta_level_data <- function(data, partition, node_list, learners = "glm", bias = NULL, g_bounds = c(0.001, 0.999)) {
  # start with glm/hal/mean learners; later to consider SuperLearner or sl3 package or a list
  # start with W, A, Y problems; should be the name of the node_list
  n_temp <- nrow(data)
  V_temp <- length(partition)
  V_vec <- names(partition)  # training sample in fast selection might have it not equal to 1:V_temp
  n_output <- sum(sapply(partition, length))
  if (n_temp != n_output) warning("partial partition is provided")
  
  # create id and v columns
  meta_data <- data.table(id = 1:n_temp, v = rep("0", n_temp))  # left 0 rows will be discarded in the end
  for (v_temp in V_vec) {
    meta_data[partition[[v_temp]], v := v_temp]
  }
  
  if (learners == "glm") {
    # train V models; can parallel
    V_lists_coefs <- lapply(V_vec, function(v_temp) {
      v_th_training_set <- get_v_th_training_ids(partition, v_temp)
      # glm type    
      Q_model <- glm(formula = paste0(node_list$Y, " ~ ", node_list$A, " + ", paste(node_list$W, collapse = " + ")), 
                     family = gaussian, 
                     data = data[v_th_training_set, ]
      )
      g_model <- glm(formula = paste0(node_list$A, " ~ ", paste(node_list$W, collapse = " + ")),
                     family = binomial, 
                     data = data[v_th_training_set, ]
      )
      list_coefs <- list()
      list_coefs$Q <- Q_model$coefficients
      list_coefs$g <- g_model$coefficients
      return(list_coefs)
    })
    names(V_lists_coefs) <- names(partition)
    
    # for each subjects, use the list of V list_coefs, calculate the needed Qnv/Wr
    meta_data[, Q1 := 0]
    meta_data[, Q0 := 0]
    meta_data[, G := 0]
    data_impute <- copy(data)
    for (v_temp in V_vec) {
      rows_id <- partition[[v_temp]]
      meta_data[rows_id, 
                Q1 := as.vector(tcrossprod(as.matrix(data_impute[, (node_list$A) := 1][rows_id, c(node_list$W, node_list$A), with = F]), 
                                           t(V_lists_coefs[[v_temp]]$Q[-1])) + V_lists_coefs[[v_temp]]$Q[1])
      ]
      meta_data[rows_id, 
                Q0 := as.vector(tcrossprod(as.matrix(data_impute[, (node_list$A) := 0][rows_id, c(node_list$W, node_list$A), with = F]), 
                                           t(V_lists_coefs[[v_temp]]$Q[-1])) + V_lists_coefs[[v_temp]]$Q[1])
      ]
      temp_G <- expit(as.vector(tcrossprod(as.matrix(data_impute[rows_id, c(node_list$W), with = F]),
                                           t(V_lists_coefs[[v_temp]]$g[-1])) + V_lists_coefs[[v_temp]]$g[1]))
      temp_G[is.nan(temp_G)] <- 1
      temp_G <- use_bounds(temp_G, g_bounds)
      meta_data[rows_id,
                G := temp_G
      ]
    }
  } else if (learners == "hal") {
    
  } else if (learners == "mean") {
    # train V models; can parallel
    V_lists_coefs <- lapply(V_vec, function(v_temp) {
      v_th_training_set <- get_v_th_training_ids(partition, v_temp)
      # mean
      list_coefs <- list()
      list_coefs$Q <- mean(unlist(data[v_th_training_set, node_list$Y, with = F]))
      list_coefs$g <- mean(unlist(data[v_th_training_set, node_list$A, with = F]))
      return(list_coefs)
    })
    names(V_lists_coefs) <- names(partition)
    
    # for each subjects, use the list of V list_coefs, calculate the needed Qnv/Wr
    meta_data[, Q1 := 0]
    meta_data[, Q0 := 0]
    meta_data[, G := 0]
    data_impute <- copy(data)
    for (v_temp in V_vec) {
      rows_id <- partition[[v_temp]]
      meta_data[rows_id, 
                Q1 := V_lists_coefs[[v_temp]]$Q
      ]
      meta_data[rows_id, 
                Q0 := V_lists_coefs[[v_temp]]$Q
      ]
      temp_G <- V_lists_coefs[[v_temp]]$g
      temp_G[is.nan(temp_G)] <- 1
      temp_G <- use_bounds(temp_G, g_bounds)
      meta_data[rows_id,
                G := temp_G
      ]
    }
    
  } else stop("learners other than glm/hal/mean; under development")
  
  
  meta_data[, (node_list$A) := data[[node_list$A]]]
  meta_data[, (node_list$Y) := data[[node_list$Y]]]
  
  # can add bias to learners
  
  return(list(meta_data = meta_data, 
              partition = partition))
}
# examples
# partition <- get_partition(1000, 10)
# node_list <- list(
#   # Wc = names(data_use) %>% head(-2) %>% head(2),
#   # We = names(data_use) %>% head(-2) %>% tail(-2),
#   W = names(data) %>% head(-2),
#   A = "treat",
#   Y = "y"
# )
# setDT(data)
# meta_level_data <- get_meta_level_data(data, partition, node_list, learners = "glm")
# data.table(meta_level_data$meta_data, A = data[[node_list$A]], Y = data[[node_list$Y]])

v_to_partition <- function(v_) {
  # v_ <- meta_data$v
  v_names <- unique(v_) %>% sort
  temp_partition <- lapply(v_names, function(each_v) {
    which(v_ == each_v)
  })
  names(temp_partition) <- v_names
  return(temp_partition)
}


estimator_m_hal_sl <- function(data, meta_level_data, option = "fast", lambda) {
  if (!(option %in% c("fast", "honest"))) stop("only support fast or honest selector")
  
  partition <- meta_level_data$partition
  
  if (option == "fast") {
    for (fold in names(partition)) {
      meta_level_data$meta_data[v != fold, ]
      # train a HAL with this data
      
      
    }
  } else if (option == "honest") {
    for (fold in names(partition)) {
      # create meta data of this fold
      temp_data <- data[meta_level_data$meta_data$v != fold, ]
      temp_partition <- meta_level_data$meta_data[v != fold, ]$v %>% v_to_partition
      temp_meta_level_data <- get_meta_level_data(temp_data, partition = temp_partition, node_list = node_list, learners = "glm")
      temp_meta_level_data$meta_data
      # train a HAL with this data
    }
  }
  
  
}


generate_data <- function(n = 250, p_cov = 10, p_image = 1000, ratio_nonzero = 0.5 , if_image = T, seed = NULL, image_noise_function = NULL, alpha = NULL, value_treat = NULL, trt_coef = -2, if_simpleG = F, 
                                   coef_A = 0.1, coef_Y = 0.1, if_binary = T, int_A = 0
) {
  if (p_cov < 2) stop("need more than two clinical covariates.")
  if (p_cov %% 2 == 1) stop("need even number p_cov; half normal half Bernoulli.")
  if (p_image %% 4 != 0) stop("need p_image divided by 4; nonzero ratio decides true confounding")
  
  if (if_image) p_nonzero <- round(ratio_nonzero * p_image)
  p_nonzero_A <- p_nonzero_Y <- 0
  if (p_nonzero_A + p_nonzero_Y + p_nonzero > p_image)  stop("too many nonzero coefficients")
  
  if (length(coef_A) == 1) if (if_image) {
    coef_A <- list(cov = rep(coef_A, p_cov), 
                   img = rep(coef_A, p_nonzero_A + p_nonzero))  # add more trt related covariates
  } 
  if (length(coef_Y) == 1) if (if_image) {
    coef_Y <- list(cov = rep(coef_Y, p_cov), 
                   img = rep(coef_Y, p_nonzero_Y + p_nonzero)) 
  } 
  
  if (length(coef_A$cov) != p_cov) stop("wrong cov dimension")
  # if (if_image) if (length(coef_A$img) != p_nonzero) stop("wrong img dimension")
  
  if (!is.null(seed)) set.seed(seed)
  
  # if (is.null(alpha)) {
  #   alpha <- rep(0.2, p_cov)
  # }  
  # alpha_inter <- 1
  
  # clinical covariates  
  # x1 = rnorm(n*p_cov/2) %>% matrix(nrow = n)
  if (if_binary) {
    x1 = rbinom(n*p_cov,1,0.5) %>% matrix(nrow = n)    
  } else {
    x1 = rnorm(n*p_cov) %>% matrix(nrow = n)    
  }
  X <- cbind(x1)
  trt_prob <- int_A + X %*% coef_A$cov
  outcome_prob <- X %*% coef_Y$cov
  
  if (if_image) { 
    
    if (if_binary) {
      img_both <- rbinom(n*p_nonzero, 1, 0.5) %>% matrix(nrow = n)
      img_trt <- rbinom(n*p_nonzero_A, 1, 0.5) %>% matrix(nrow = n)
      img_out <- rbinom(n*p_nonzero_Y, 1, 0.5) %>% matrix(nrow = n)
      img_noise <- rbinom(n*(p_image - p_nonzero - p_nonzero_A - p_nonzero_Y), 1, 0.5) %>% matrix(nrow = n)
    } else {
      img_both <- rnorm(n*p_nonzero) %>% matrix(nrow = n)
      img_trt <- rnorm(n*p_nonzero_A) %>% matrix(nrow = n)
      img_out <- rnorm(n*p_nonzero_Y) %>% matrix(nrow = n)
      img_noise <- rnorm(n*(p_image - p_nonzero - p_nonzero_A - p_nonzero_Y)) %>% matrix(nrow = n)
    }
    
    img_mat <- cbind(img_both, img_trt, img_out, img_noise)
    # img_mat <- img_mat[, sample(ncol(img_mat))]
    
    X <- cbind(X, img_mat)
    if (!is.null(image_noise_function)) {
      X[, -(1:p_cov)] <- apply(X[, -(1:p_cov)], 2, image_noise_function)
    }
    if (!if_simpleG) {
      trt_prob <- trt_prob + cbind(img_both, img_trt) %*% coef_A$img
    }  
    # if if_simpleG == T, then G does not depend on images; default is F
    outcome_prob <- outcome_prob + cbind(img_both, img_out) %*% coef_Y$img
  }
  
  trt_prob <- trt_prob %>% expit
  
  if (is.null(value_treat)) value_treat <- rbinom(n = n, size = 1, trt_prob)  # otherwise keep the imputed value
  X_output <- cbind(X, value_treat)
  
  colnames(X_output)[1:p_cov] <- paste0("X", 1:p_cov)
  colnames(X_output)[(1:p_image) + p_cov] <- paste0("img", 1:p_image)
  colnames(X_output)[p_cov + p_image + 1] <- "trt"
  
  outcome_prob <- outcome_prob + trt_coef * value_treat
  # outcome_prob <- outcome_prob %>% expit
  # value_outcome <- rbinom(n = n, size = 1, outcome_prob)
  value_outcome <- outcome_prob + rnorm(length(outcome_prob))
  X_output <- data.frame(X_output, Y = value_outcome)
  
  X_output %>% return
}




call_with_args <- function(fun, args, other_valid = list(), keep_all = FALSE,
                           silent = FALSE, ignore = c()) {
  
  # drop ignore args
  args <- args[!(names(args) %in% ignore)]
  if (!keep_all) {
    # catch arguments to be kept
    formal_args <- names(formals(fun))
    all_valid <- c(formal_args, other_valid)
    
    # find invalid arguments based on combination of formals and other_valid
    invalid <- names(args)[which(!(names(args) %in% all_valid))]
    
    # subset arguments to pass
    args <- args[which(names(args) %in% all_valid)]
    
    # return warnings when dropping arguments
    if (!silent & length(invalid) > 0) {
      message(sprintf(
        "Learner called function %s with unknown args: %s. These will be dropped.\nCheck the params supported by this learner.",
        as.character(substitute(fun)), paste(invalid, collapse = ", ")
      ))
    }
  }
  do.call(fun, args)
}
