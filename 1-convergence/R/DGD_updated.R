makeData5 <- function(n){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  
  y <- -1*as.numeric(x1< -3)*x3 + 0.5*as.numeric(x1 > -2) - 
    1*as.numeric(x1>0) + 2*as.numeric(x1>2)*x3 - 3*as.numeric(x1>3) + 
    1.5*as.numeric(x2 > -1) - 
    5*as.numeric(x2>1)*x3 + 2*as.numeric(x2>3) + as.numeric(x4 < 0)*2 - 
    as.numeric(x5 > 5)*1 - as.numeric(x4<0)*as.numeric(x1<0) + 2*x3 + 
    rnorm(n)
  data <- data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5, y=y)
  covariates <- paste0("x",seq(1:5))
  outcome <- "y"
  outcome_type <- "continuous"
  return(list(data = data, 
              covariates = covariates,
              outcome = outcome,
              outcome_type = outcome_type))
}


makeData20 <- function(n, interaction){
  x1 <- runif(n,-4,4)
  x2 <- runif(n,-4,4)
  x3 <- rbinom(n, 1, 0.5)
  x4 <- rnorm(n)
  x5 <- rgamma(n, 2, 1)
  x6 <- rpois(n, 2)
  x7 <- rexp(n, 3)
  x8 <- rbeta(n, 1, 1)
  x9 <- rchisq(n, 2)
  x10 <- rgeom(n, 0.6)
  ## same
  x11 <- runif(n,-4,4)
  x12 <- runif(n,-4,4)
  x13 <- rbinom(n, 1, 0.5)
  x14 <- rnorm(n)
  x15 <- rgamma(n, 2, 1)
  x16 <- rpois(n, 1)
  x17 <- rexp(n, 1)
  x18 <- rbeta(n, 2, 1)
  x19 <- rchisq(n, 1)
  x20 <- rgeom(n, 0.8)  
  
  # x1 ~ x5
  g1 <- -1*as.numeric(x1< -3)*x3 + 
    0.5*as.numeric(x1 > -2) - 
    1*as.numeric(x1>0) + 
    2*as.numeric(x1>2)*x3 - 
    3*as.numeric(x1>3) + 
    1.5*as.numeric(x2 > -1) - 
    5*as.numeric(x2>1)*x3 + 
    2*as.numeric(x2>3) + 
    as.numeric(x4 < 0)*x2 - 
    as.numeric(x5 > 5)*x1 - 
    as.numeric(x4<0)*as.numeric(x1<0) + 
    2*x3 
  
  # x6 ~ 10
  g2 <- -1*as.numeric(x10 > 3)*x8 + 
    0.5*as.numeric(x10 > 2) - 
    1*as.numeric(x10>0) + 
    2*as.numeric(x10>2)*x8 - 
    3*as.numeric(x10>3) + 
    1.5*as.numeric(x9 > 5) - 
    5*as.numeric(x9>1)*x8 + 
    2*as.numeric(x9>3) + 
    as.numeric(x7 > 4)*x9 - 
    as.numeric(x6 > 5)*x10 - 
    as.numeric(x7 > 1)*as.numeric(x10 < 2) + 
    2*x8 
  
  # x11 ~ x15
  g3 <- -1*as.numeric(x11< -3)*x13 + 
    0.5*as.numeric(x11 > -2) - 
    1*as.numeric(x11>0) + 
    2*as.numeric(x11>2)*x13 - 
    3*as.numeric(x11>3) + 
    1.5*as.numeric(x12 > -1) - 
    5*as.numeric(x12>1)*x13 + 
    2*as.numeric(x12>3) + 
    as.numeric(x14 < 0)*x12 - 
    as.numeric(x15 > 5)*x11 - 
    as.numeric(x14 < 0)*as.numeric(x11<0) + 
    2*x13 
  # x16 ~ x20
  g4 <- -1*as.numeric(x19 > 3)*x17 + 
    0.5*as.numeric(x19 > 2) - 
    1*as.numeric(x19>0) + 
    2*as.numeric(x19>2)*x18 - 
    3*as.numeric(x19>3) + 
    1.5*as.numeric(x16 > 5) - 
    5*as.numeric(x16>1)*x18 + 
    2*as.numeric(x16>3) + 
    as.numeric(x16 > 4)*x19 - 
    as.numeric(x16 > 5)*x20 - 
    as.numeric(x17 > 1)*as.numeric(x20 < 2) + 
    2 * x18
  
  if (interaction){
    y <- g1 + g2 + g3 + g4 + g1*g3 - g1*g2 + g2*g4 - g3*g4 + rnorm(n)
  } else{
    y <- g1 + g2 + g3 + g4 + rnorm(n)
  }
  
  data <- data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10,
                     x11=x11,x12=x12,x13=x13,x14=x14,x15=x15,x16=x16,x17=x17,x18=x18,
                     x19=x19, x20=x20, g1=g1, g2=g2, g3=g3, g4=g4, y=y)
  covariates <- paste0("x",seq(1:20))
  outcome <- "y"
  outcome_type <- "continuous"
  return(list(data = data, 
              covariates = covariates,
              outcome = outcome,
              outcome_type = outcome_type))
}
