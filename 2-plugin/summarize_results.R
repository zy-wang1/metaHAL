library(dplyr)
library(purrr)

# summarize the two folders generate tables and figures for lowD and highD settings
# file_name <- "plugin_glm.R"
# file_name <- "plugin_lasso.R"

for (file_name in c("plugin_glm.R", "plugin_lasso.R")) {

n_iter <- 500
sample_size <- 200
truth <- 1
if (file_name == "plugin_glm.R") {
  p_image <- 4
  ratio_nonzero <- vec_ratio_overfit <- 1
} else {
  p_image <- 200 * 4  # dim of additional covariates
  ratio_nonzero <- 0.005  # ratio of true non-zero coefficients from the additional covariates 
  vec_ratio_overfit <- c(1, 0.001)
}

list_report <- list()
for (loc_overfit in seq_along(vec_ratio_overfit)) {
  ratio_overfit <- vec_ratio_overfit[loc_overfit]
  vec_exist <- sapply(1:n_iter, function(i_iter) {
    output_target <- paste0("2-plugin/Results/", file_name, "/", "n_", sample_size, "_pimage_", p_image, "_ratio_nonzero_", ratio_nonzero, "_ratio_overfit_", ratio_overfit, "_raw_results_", i_iter, ".rds")  
    file.exists(output_target)
  })
  table(vec_exist)
  
  vec_files <- paste0("2-plugin/Results/", file_name, "/", "n_", sample_size, "_pimage_", p_image, "_ratio_nonzero_", ratio_nonzero, "_ratio_overfit_", ratio_overfit, "_raw_results_", 1:n_iter, ".rds")
  list_result <- lapply(vec_files, function(x) if(file.exists(x)) readRDS(x))
  
  list_est <- lapply(list_result, function(raw_iter) {
    vec_iter <- c(noadj = raw_iter$noadj[1],
                  # tmle_init = raw_iter$tmle[1],
                  tmle_full = raw_iter$tmle[2],
                  meta_init = raw_iter$meta_init[1],
                  meta_us = raw_iter$meta_undersmoothed_2_overfitG[1]
                  # ,meta_us_metaG = raw_iter$meta_undersmoothed_2_metaG[1]
    )
    return(vec_iter)
  })
  
  tab_est <- list_est %>% (purrr::compact) %>% do.call(what = rbind)
  colnames(tab_est) <- c("noadj",
                         # "tmle_init",
                         "tmle", "meta_init",
                         "meta_undersmoothed")
  col_order <- colnames(tab_est)
  
  
  library(ggplot2)
  
  df_long <- tab_est %>% as.data.frame %>%
    tidyr::pivot_longer(everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Column = factor(Column, levels = col_order))
  
  mse <- function(x) {
    mean((x - truth)^2)
  }
  
  pp <- ggplot(df_long, aes(x = Column, y = Value)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, color = "red") +  # Add horizontal line at y = 1 with red color
    labs(x = "", y = "Estimate") +
    theme_bw()
  pp %>% ggsave(filename = paste0("figure/summary_p_image_",  p_image, "_overfit_", ratio_overfit, ".png"), 
                units = "in", width = 8, height = 6, dpi = 300, device = "png"
                )
  
  
  tab_measures <- cbind(
    MSE = tab_est %>% apply(2, function(ests) {
      mean((ests - truth)^2)
    }),
    Bias = tab_est %>% apply(2, function(ests) {
      mean(ests) - truth
    }),
    SD = tab_est %>% apply(2, function(ests) {
      sd(ests)
    })
  )
  tab_measures <- as.data.frame(tab_measures)
  tab_measures <- tab_measures %>% mutate(Ratio = Bias/SD) %>% round(3)
  
  list_report[[loc_overfit]] <- tab_measures
}

list_report



library(knitr)
library(kableExtra)

if (length(list_report) == 1) {
  table_data <- list_report[[1]]
  latex_table <- kable(table_data, format = "latex", row.names = T, 
                       col.names = c(" ", colnames(table_data)),
                       align = "rrrrr", 
                       booktabs = TRUE,
                       linesep = "")
} else {
  noadj_row <- list_report[[1]][1, ]
  blocks <- lapply(list_report, function(x) x[-1, ])
  all_blocks <- c(list(noadj_row), blocks)
  block_names <- c("", "$\\lambda_{\\text{cv}}$", "$0.1\\% \\lambda_{\\text{cv}}$")
  
  latex_lines <- c("\\begin{tabular}{rrrrr}",
             "  \\hline")
  for (i in seq_along(all_blocks)) {
    latex_lines <- c(latex_lines, paste0(block_names[i], " & MSE & Bias & SD & Ratio \\\\"), 
                     "\\hline")  
    block_latex <- kable(all_blocks[[i]], format = "latex", row.names = TRUE, col.names = NULL)
    block_body <- strsplit(block_latex, "\n")[[1]]
    content_lines <- block_body[grepl("&", block_body)]
    latex_lines <- c(latex_lines, content_lines, "  \\hline")
  }
  latex_lines <- c(latex_lines, "\\end{tabular}")
  latex_table <- latex_lines
}
cat(latex_table, file = paste0("table/summary_p_image_",  p_image, ".tex"), sep = "\n")

}