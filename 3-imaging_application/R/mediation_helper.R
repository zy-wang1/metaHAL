funMakeItSparse <- function(x) {
  Matrix(x, sparse = T)
}
funMakeItBig <- function(x) {
  as.big.matrix(as.matrix(x))
}
getPercent <- function(est_vec) {
  (est_vec[3] - est_vec[2])/(est_vec[3] - est_vec[1])
}
bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}
submodel = function(epsilon, initial, H) {
  plogis(qlogis(initial) + H %*% epsilon)
}
apply_submodel = function(submodel_data, epsilon) {
  submodel(epsilon, submodel_data$initial, submodel_data$H)
}
funTMLE <- function(observed, initial, H) {
  submodel_data <- list(
    observed = observed, 
    initial = initial, 
    H = H 
  )
  if (is.null(ncol(submodel_data$H))) submodel_data$H <- matrix(submodel_data$H, ncol = 1)
  submodel_fit <- glm(observed ~ H - 1, submodel_data,
                      offset = qlogis(submodel_data$initial),
                      family = binomial(),
                      start = rep(0, ncol(submodel_data$H))
  )
  epsilon <- coef(submodel_fit)
  update <- apply_submodel(submodel_data, epsilon)
  return(list(update = update, epsilon = epsilon))
}
funUS <- function(observed,
                  initial,
                  H) {
  n_update <- length(H)
  # check |P_n D*| / SE(D*) =< max(1/log(n), 1/10)
  level_solved <- apply(initial, 2, function(each_col) {
    mean(H * (observed - each_col))  
  }) %>% unlist
  vec_threshold <- apply(initial, 2, function(each_col) {
    sqrt(var(H * (observed - each_col))/n_update) / min(log(n_update), 10)
  }) %>% unlist
  loc_select <- which(abs(level_solved) < vec_threshold)[1]
  if (is.na(loc_select)) loc_select <- ncol(initial)
  return(loc_select)
}
