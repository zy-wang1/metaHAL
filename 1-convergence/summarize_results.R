library(dplyr)
library(purrr)
# create latex outputs
library(knitr)
library(kableExtra)

# run each scenario to regenerate the summary tables
# i_scenario <- 1
# i_scenario <- 3
# i_scenario <- 2
# i_scenario <- 6

for (i_scenario in c(1, 3, 2, 6)) {

vec_sim_names <- c("scenario_1_1", "scenario_1_2", "scenario_2_1", "scenario_2_2", 
                   "scenario_2_3", "scenario_2_4", "scenario_3_1", "scenario_3_2", 
                   "scenario_3_3", "scenario_3_4")
sim_name <- vec_sim_names[i_scenario]
print(sim_name)

results_dir <- file.path("1-convergence/Results/", sim_name)
vec_files <- list.files(results_dir)
vec_match <- grep("^z_.*\\.rds$", vec_files, value = TRUE)
length(vec_match)  # number of completed runs

all_files <- paste0("z_", 1:100, ".rds")
missing_files <- setdiff(all_files, vec_match)
missing_numbers <- as.numeric(gsub("z_(\\d+)\\.rds", "\\1", missing_files))
length(missing_numbers)  # number of missing runs

list_results <- lapply(vec_match, function(x) {
  readRDS(file.path(results_dir, x))
})

key_cols <- c("sim_name", "n", "type", "learner")

# function: averaging mse and relative mse
average_df_position <- function(position) {
  dfs <- map(list_results, function(x) {
    df <- x[[position]]
    mse_ref <- df$mse[7]
    df <- df %>%
      mutate(relative_mse = mse / mse_ref)
    return(df)
  })
  combined_df <- bind_rows(dfs)

  # get original row order
  original_order <- list_results[[1]][[position]] %>%
    select(all_of(key_cols)) %>%
    mutate(.row_id = row_number())
  
  # summarize metrics
  summarized <- combined_df %>%
    group_by(across(all_of(key_cols))) %>%
    summarize(
      mean_mse = mean(mse),
      mean_se = mean(se),
      var_of_mse = var(mse),
      relative_mse = mean(relative_mse),
      .groups = "drop"
    )
  
  # create outputs
  final <- original_order %>%
    left_join(summarized, by = key_cols) %>%
    arrange(.row_id) %>%
    select(-.row_id) %>%
    mutate(across(c(mean_mse, mean_se, 
                    # var_of_mse, 
                    relative_mse
                    ), ~ round(.x, 3)))
  
  return(final)
}

# apply to different sample sizes
average_results <- list(
  average_df_position(1),
  average_df_position(2), 
  average_df_position(3),
  average_df_position(4)
)
names(average_results) <- paste0("result_", seq_along(average_results))

average_results

trimmed_results <- lapply(seq_along(average_results), function(i) {
  df <- average_results[[i]]
  suffix <- paste0("_", i)
  df_subset <- df[, c("sim_name", "type", "learner", "n", "mean_mse", 
                      # "var_of_mse", 
                      "relative_mse"
                      )]
  colnames(df_subset)[4:6] <- paste0(colnames(df_subset)[4:6], suffix)
  return(df_subset)
})

combined_result <- Reduce(function(x, y) full_join(x, y, by = c("sim_name", "type", "learner")), trimmed_results)
combined_result <- combined_result %>% select(-sim_name)

combined_result



table_data <- combined_result[, c("learner", 
                                  "mean_mse_1", "relative_mse_1", 
                                  "mean_mse_2", "relative_mse_2",
                                  "mean_mse_3", "relative_mse_3", 
                                  "mean_mse_4", "relative_mse_4")]
colnames(table_data) <- c("metalearner", 
                          "mse", "relative_mse", 
                          "mse", "relative_mse",
                          "mse", "relative_mse",
                          "mse", "relative_mse")
table_data[, -1] <- round(table_data[, -1], 2)

latex_table <- kable(table_data, format = "latex", booktabs = TRUE, align = "lcccccccc", linesep = "") %>%
  add_header_above(c(" " = 1, 
                     "n = 200" = 2, 
                     "n = 500" = 2, 
                     "n = 1000" = 2, 
                     "n = 2000" = 2)) %>%
  kable_styling(latex_options = NULL, full_width = FALSE, position = "left") %>%
  row_spec(0, bold = FALSE)

cat(capture.output(latex_table), file = file.path("table", paste0(sim_name, ".tex")), sep = "\n")
}