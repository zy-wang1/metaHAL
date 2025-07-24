sim_data_lrnr_pairs <-
  list(
    "scenario_1_1" = list(simData = makeData5,
                          base_lrnrs = generateBaselrnrs(type = "Multiple")),
    "scenario_1_2" = list(
      simData = makeData5,
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:5),
        type = "OLS",
        subset = TRUE,
        n_subsets = 5,
        order = TRUE
      )
    ),
    "scenario_2_1" = list(
      simData = function(n)
        return(makeData20(n, interaction = FALSE)),
      base_lrnrs = generateBaselrnrs(type = "Multiple")
    ),
    "scenario_2_2" = list(
      simData = function(n)
        return(makeData20(n, interaction = FALSE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "CART",
        subset = TRUE,
        n_subsets = 4,
        order = TRUE
      )
    ),
    "scenario_2_3" = list(
      simData = function(n)
        return(makeData20(n, interaction = FALSE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "CART",
        subset = TRUE,
        n_subsets = 4,
        order = FALSE
      )
    ),
    "scenario_2_4" = list(
      simData = function(n)
        return(makeData20(n, interaction = FALSE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "OLS",
        subset = TRUE,
        n_subsets = 20,
        order = TRUE
      )
    ),
    "scenario_3_1" = list(
      simData = function(n)
        return(makeData20(n, interaction = TRUE)),
      base_lrnrs = generateBaselrnrs(type = "Multiple")
    ),
    "scenario_3_2" = list(
      simData = function(n)
        return(makeData20(n, interaction = TRUE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "CART",
        subset = TRUE,
        n_subsets = 4,
        order = TRUE
      )
    ),
    "scenario_3_3" = list(
      simData = function(n)
        return(makeData20(n, interaction = TRUE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "CART",
        subset = TRUE,
        n_subsets = 4,
        order = FALSE
      )
    ),
    "scenario_3_4" = list(
      simData = function(n)
        return(makeData20(n, interaction = TRUE)),
      base_lrnrs = generateBaselrnrs(
        covariates = paste0("x", 1:20),
        type = "OLS",
        subset = TRUE,
        n_subsets = 20,
        order = TRUE
      )
    )
  )
