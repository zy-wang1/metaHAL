library(dplyr)
library(tidyr)  # gather
library(ggplot2)

# T will use real features; F will use noise
# if_true_effect <- F
if_true_effect <- T

n_iter <- 100
feature_sets <- list(10, 18, 34 , 50, 101, 152, 200,
                     c(10, 18, 34, 50, 101, 152, 200)
)
file_name <- "plugin_mediation.R"
task_label <- paste0(file_name, "_", if_true_effect)
if (!dir.exists("3-imaging_application/Results")) dir.create("3-imaging_application/Results")
if (!dir.exists(paste0("3-imaging_application/Results/", task_label))) dir.create(paste0("3-imaging_application/Results/", task_label))
N <- 10000  # sample size of each iteration
V <- 5  # V fold cross-validation

feature_set <- feature_sets[[1]]
feature_label <- ifelse(length(feature_set) == 1, feature_set, "ensemble")
vec_exist <- sapply(1:n_iter, function(i_iter) {
  output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", if_true_effect, "_feature_", feature_label, "_", i_iter, ".rds")  
  file.exists(output_target)
})

list_vec <- list()
for (loc_feature in seq_along(feature_sets)) {
  # loc_feature <- 1
  feature_set <- feature_sets[[loc_feature]]
  feature_label <- ifelse(length(feature_set) == 1, feature_set, "ensemble")
  task_label <- paste0(file_name, "_", if_true_effect)
  list_results <- lapply(1:n_iter, function(i_iter) {
    output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", if_true_effect, "_feature_", feature_label, "_", i_iter, ".rds")
    if(file.exists(output_target))
      readRDS(output_target)
  })
  vec1 <- list_results %>% lapply(function(x) if(is.character(x)) return(NULL) else return(x)) %>% do.call(what = rbind) %>% as.data.frame
  vec <- vec1[, 5] %>% matrix(nrow = nrow(vec1))
  vec %>% colMeans
  vec %>% apply(2, function(x) quantile(x, c(0.025, 0.975)))
  vec <- data.frame(vec)
  
  colnames(vec) <- paste0("ResNet_", feature_label)
  
  list_vec <- c(list_vec, vec)
}
tab_vec <- do.call(cbind, list_vec)
df_vec <- data.frame(tab_vec)
names(df_vec)[ncol(df_vec)] <- "Ensemble"

long_df <- gather(df_vec, factor_key = TRUE)

# Calculate the quantiles
quantiles_df <- long_df %>%
  group_by(key) %>%
  summarize(lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

# Plot
temp <- ggplot(quantiles_df, aes(x = key, ymin = lower, ymax = upper, group = key, color = key)) +
  geom_errorbar(width = 0.2) +
  geom_point(aes(y = lower)) +
  geom_point(aes(y = upper)) +
  # scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0.05)) + # Display y-axis labels as percentages
  labs(y = "Percentage Mediated", title = "", x = "", color = "")+
  theme_bw() + 
  geom_hline(yintercept = 0, linetype = "dotted")
temp
temp %>% ggsave(filename = paste0("figure/imaging_ensemble.png"), device = "png", width = 8, height = 6, units = "in", dpi = 300)




# comparison plot
loc_feature <- 1  # can set to 1-7 for 10, 18, 34, 50, 101, 152, 200 layers
feature_set <- feature_sets[[loc_feature]]
feature_label <- ifelse(length(feature_set) == 1, feature_set, "ensemble")
task_label <- paste0(file_name, "_", "TRUE")
list_results_T <- lapply(1:n_iter, function(i_iter) {
  output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", "TRUE", "_feature_", feature_label, "_", i_iter, ".rds")
  if(file.exists(output_target))
    readRDS(output_target)
})
task_label <- paste0(file_name, "_", "FALSE")
list_results_F <- lapply(1:n_iter, function(i_iter) {
  output_target <- paste0("3-imaging_application/Results/", task_label, "/", "n_", N, "_V_", V, "_real_feature_", "FALSE", "_feature_", feature_label, "_", i_iter, ".rds")
  if(file.exists(output_target))
    readRDS(output_target)
})

vec1 <- list_results_T %>% lapply(function(x) if(is.character(x)) return(NULL) else return(x)) %>% do.call(what = rbind) %>% as.data.frame
vec2 <- list_results_F %>% lapply(function(x) if(is.character(x)) return(NULL) else return(x)) %>% do.call(what = rbind) %>% as.data.frame

list_quantiles <- lapply(c(3, 2, 5), function(iii) {
  min_nrow <- min(nrow(vec1), nrow(vec2))
  vec <- cbind(vec1[1:min_nrow, iii], vec2[1:min_nrow, iii])
  vec %>% apply(2, function(x) quantile(x, c(0.025, 0.975)))
  vec <- data.frame(vec)
  
  colnames(vec) <- c("True Samples", "Null Samples")
  long_df <- gather(vec, factor_key = TRUE)
  
  quantiles_df <- long_df %>%
    group_by(key) %>%
    summarize(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975))
  
  quantiles_df
})

combined_df <- do.call(rbind, list_quantiles)
combined_df$type <- rep(c("IPW", "TMLE", "metaHAL"), c(2, 2, 2))
combined_df$type <- factor(combined_df$type, levels = c("IPW", "TMLE", "metaHAL"))


# ignore ipw here too biased
temp <- ggplot(combined_df[-(1:2), ]
               , aes(x = key, ymin = lower, ymax = upper, group = key, color = key)) +
  facet_grid(. ~ type, scales = "free"
             # , labeller = label_both
  ) + 
  geom_errorbar(width = 0.2, show.legend = FALSE) +
  geom_point(aes(y = lower), show.legend = FALSE) +
  geom_point(aes(y = upper), show.legend = FALSE) +
  labs(y = "Percentage Mediated", title = "", x = "", color = "")+
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")
temp
temp %>% ggsave(filename = paste0("figure/imaging_compare_ResNet_", feature_label, ".png"), device = "png", width = 7, height = 4, units = "in", dpi = 300)
