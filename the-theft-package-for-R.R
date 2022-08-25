#----------------------------------------------
# This script produces the results and graphics 
# presented in https://arxiv.org/abs/2208.06146
#----------------------------------------------

#----------------------------------------
# Author: Trent Henderson, 25 August 2022
#----------------------------------------

library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(R.matlab)
library(theft)
library(scales)
library(tidyr)
library(reticulate)

# Set up Python environment
# NOTE: This will be different for others depending where they installed tsfresh, Kats, and TSFEL

init_theft("~/opt/anaconda3/bin/python")

# Load in data

tmp <- process_hctsa_file("/Users/trenthenderson/Downloads/INP_Bonn_EEG.mat") %>% # Users will need to set this to wherever they download it
  mutate(id = unlist(id),
         group = unlist(group))

#------------------- Time series plot -----------------------

set.seed(123)

#' Function to sample IDs by class
#' @param data the dataframe containing time-series data
#' @param group_name string specifying the class to filter by
#' @param n number of samples to generate
#' @return object of class vector
#' @author Trent Henderson
#' 

draw_samples <- function(data, group_name, n){
  
  samps <- data %>%
    filter(group == group_name) %>%
    dplyr::select(id) %>%
    distinct() %>%
    pull(id) %>%
    sample(size = n) %>%
    unlist()
  
  return(samps)
}

ids <- unique(tmp$group) %>%
  purrr::map(~ draw_samples(data = tmp, group_name = .x, n = 1)) %>%
  unlist()

tsplot <- tmp %>%
  filter(id %in% ids) %>%
  mutate(id = gsub(".dat", "\\1", id)) %>%
  ggplot(aes(x = timepoint, y = values, colour = group)) +
  geom_line(size = 0.3) +
  labs(title = "Raw time series samples from all five classes",
       x = "Time",
       y = "Value",
       colour = NULL) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  facet_wrap(~id, ncol = 2, scales = "free_y")

#------------------- theft functionality --------------------

#--------------------------------
# Calculate features for all sets
#--------------------------------

# Raw data version

all_features <- calculate_features(data = tmp, 
                                   id_var = "id", 
                                   time_var = "timepoint", 
                                   values_var = "values", 
                                   group_var = "group",
                                   feature_set = c("catch22", "feasts", "tsfeatures", "tsfresh", "TSFEL", "Kats"),
                                   seed = 123)

# z-scored version

z_tmp <- tmp %>%
  group_by(id) %>%
  mutate(values = (values - mean(values, na.rm = TRUE)) / sd(values, na.rm = TRUE)) %>%
  ungroup()

all_features_z <- calculate_features(data = z_tmp, 
                                   id_var = "id", 
                                   time_var = "timepoint", 
                                   values_var = "values", 
                                   group_var = "group",
                                   feature_set = c("catch22", "feasts", "tsfeatures", "tsfresh", "TSFEL", "Kats"),
                                   seed = 123)

#--------------------
# Plot quality matrix
#--------------------

p <- plot_quality_matrix(all_features)

#-------------------
# Normalise features
#-------------------

x <- all_features %>%
  filter(names == "DN_HistogramMode_5") %>%
  pull(values)

xnormed <- normalise_feature_vector(x, method = "RobustSigmoid")

normed <- normalise_feature_frame(all_features, 
                                  names_var = "names", 
                                  values_var = "values", 
                                  method = "RobustSigmoid")

#----------------------------
# Plot feature matrix heatmap
#----------------------------

p2 <- plot_all_features(all_features, 
                        is_normalised = FALSE,
                        id_var = "id", 
                        method = "RobustSigmoid",
                        clust_method = "average",
                        interactive = FALSE) +
  theme(text = element_text(size = 14),
        legend.key.width = unit(1.5, "cm"))

#--------------------------------
# Plot low dimensional projection
#--------------------------------

p3 <- plot_low_dimension(all_features_z, 
                         is_normalised = FALSE, 
                         id_var = "id",
                         group_var = "group", 
                         method = "MinMax", 
                         low_dim_method = "t-SNE", 
                         perplexity = 15,
                         plot = TRUE,
                         seed = 123) +
  theme(text = element_text(size = 12))

#---------------------------
# Top feature classification
#---------------------------

# Run the function

top_feature_results <- compute_top_features(all_features, 
                                            id_var = "id", 
                                            group_var = "group",
                                            num_features = 40, 
                                            normalise_violin_plots = FALSE,
                                            method = "RobustSigmoid",
                                            cor_method = "spearman",
                                            test_method = "svmLinear",
                                            clust_method = "average",
                                            use_balanced_accuracy = FALSE,
                                            use_k_fold = TRUE,
                                            num_folds = 10,
                                            use_empirical_null =  TRUE,
                                            null_testing_method = "ModelFreeShuffles",
                                            num_permutations = 1000,
                                            p_value_method = "gaussian",
                                            pool_empirical_null = FALSE,
                                            seed = 123)

# z-scored

top_feature_results_z <- compute_top_features(all_features_z, 
                                              id_var = "id", 
                                              group_var = "group",
                                              num_features = 40, 
                                              normalise_violin_plots = FALSE,
                                              method = "RobustSigmoid",
                                              cor_method = "spearman",
                                              test_method = "svmLinear",
                                              clust_method = "average",
                                              use_balanced_accuracy = FALSE,
                                              use_k_fold = TRUE,
                                              num_folds = 10,
                                              use_empirical_null =  TRUE,
                                              null_testing_method = "ModelFreeShuffles",
                                              num_permutations = 1000,
                                              p_value_method = "gaussian",
                                              pool_empirical_null = FALSE,
                                              seed = 123)

# Generate table summary for article

top_table_summary <- top_feature_results$ResultsTable %>%
  dplyr::select(c(feature, accuracy, p_value_accuracy)) %>%
  arrange(p_value_accuracy) %>%
  mutate(feature_set = case_when(
          grepl("catch22", feature)    ~ "catch22",
          grepl("feasts", feature)     ~ "feasts",
          grepl("tsfeatures", feature) ~ "tsfeatures",
          grepl("kats", feature)       ~ "Kats",
          grepl("tsfel", feature)      ~ "TSFEL",
          grepl("tsfresh", feature)    ~ "tsfresh")) %>%
  mutate(feature = gsub("catch22_", "\\1", feature),
         feature = gsub("feasts_", "\\1", feature),
         feature = gsub("tsfeatures_", "\\1", feature),
         feature = gsub("kats_", "\\1", feature),
         feature = gsub("tsfel_", "\\1", feature),
         feature = gsub("tsfresh_", "\\1", feature)) %>%
  dplyr::select(c(feature, feature_set, accuracy, p_value_accuracy)) %>%
  mutate(accuracy = accuracy * 100,
         p_value_accuracy = ifelse(p_value_accuracy < .001, "p < .001", p_value_accuracy))

#------------------------------------------------------------------
# Replicate theft functionality to style correlation plot for paper
#------------------------------------------------------------------

# Copy theft function

draw_top_feature_plot <- function(data, cor_method, clust_method, num_features){
  
  # Wrangle dataframe
  
  cor_dat <- data %>%
    dplyr::select(c(.data$id, .data$names, .data$values)) %>%
    tidyr::drop_na() %>%
    tidyr::pivot_wider(id_cols = .data$id, names_from = .data$names, values_from = .data$values) %>%
    dplyr::select(-c(.data$id))
  
  # Calculate correlations
  
  result <- abs(stats::cor(cor_dat, method = cor_method))
  
  # Perform clustering
  
  row.order <- stats::hclust(stats::dist(result, method = "euclidean"), method = clust_method)$order # Hierarchical cluster on rows
  col.order <- stats::hclust(stats::dist(t(result), method = "euclidean"), method = clust_method)$order # Hierarchical cluster on columns
  dat_new <- result[row.order, col.order] # Re-order matrix by cluster outputs
  cluster_out <- reshape2::melt(as.matrix(dat_new)) # Turn into dataframe
  
  # Define a nice colour palette consistent with RColorBrewer in other functions
  
  mypalette <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
  
  # Draw plot
  
  FeatureFeatureCorrelationPlot <- cluster_out %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$Var1, y = .data$Var2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
    ggplot2::labs(title = paste0("Pairwise correlation matrix of top ", num_features, " features"),
                  x = NULL,
                  y = NULL,
                  fill = "Absolute correlation coefficient") +
    ggplot2::scale_fill_stepsn(n.breaks = 6, colours = rev(mypalette),
                               show.limits = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(FeatureFeatureCorrelationPlot)
}

# Filter the main dataset and mimic theft name cleaning

all_features_top_40 <- all_features %>%
  mutate(names = gsub("__", "_", names),
         names = gsub(" ", "_", names),
         names = gsub('"', "", names),
         names = gsub("\\\\", "_", names),
         names = paste0(method, "_", names),
         names = stringr::str_to_lower(names),
         names = ifelse(names == "catch22_sc_fluctanal_2_rsrangefit_50_1_logi_prop_r1", 
                        "catch22_sc_fluct_anal_2_rsrangefit_50_1_logi_prop_r1", names)) %>% # Cleaning up to replicate {theft} processing
  filter(names %in% unique(top_feature_results$ResultsTable$feature)) %>%
  mutate(names = gsub("catch22", "CATCH22", names),
         names = gsub("feasts", "FEASTS", names),
         names = gsub("tsfeatures", "TSFEATURES", names),
         names = gsub("kats", "KATS", names),
         names = gsub("tsfresh", "TSFRESH", names),
         names = gsub("tsfel", "TSFEL", names)) # For visual clarity

# Draw the plot

FeatureFeatureCorrelationPlot <- draw_top_feature_plot(data = all_features_top_40,
                                                       cor_method = "spearman",
                                                       clust_method = "average",
                                                       num_features = 40) +
  theme(title = element_text(size = 14))

#--------------------------------------------------------------
# Replicate theft functionality to style violin plots for paper
#--------------------------------------------------------------

# Copy and slightly modify theft function for paper

plot_feature_discrimination <- function(data, id_var = "id", group_var = "group",
                                        normalise = FALSE,
                                        method = c("z-score", "Sigmoid", "RobustSigmoid", "MinMax")){
  
  #------------- Normalise data -------------------
  
  if(normalise){
    
    normed <- data %>%
      dplyr::select(c(.data$id, .data$names, .data$values, .data$group)) %>%
      tidyr::drop_na() %>%
      dplyr::group_by(.data$names) %>%
      dplyr::mutate(values = normalise_feature_vector(.data$values, method = method)) %>%
      dplyr::ungroup() %>%
      tidyr::drop_na()
    
    if(nrow(normed) != nrow(data)){
      message("Filtered out rows containing NaNs.")
    }
  } else{
    normed <- data
  }
  
  #------------- Produce plots --------------------
  
  # Draw plot
  
  p <- normed %>%
    dplyr::mutate(group = as.factor(.data$group)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$group, y = .data$values, colour = .data$group)) +
    ggplot2::geom_violin()
  
  if(length(unique(normed$names)) > 8){
    p <- p +
      ggplot2::geom_point(size = 0.7, alpha = 0.9, position = ggplot2::position_jitter(w = 0.05))
  } else{
    p <- p +
      ggplot2::geom_point(size = 1, alpha = 0.9, position = ggplot2::position_jitter(w = 0.05))
  }
  
  p <- p +  
    ggplot2::labs(title = "Class discrimination for sample of top performing features",
                  x = "Class",
                  y = "Value") +
    ggplot2::scale_colour_brewer(palette = "Dark2") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90))
  
  if(normalise){
    p <- p +
      ggplot2::facet_wrap(~names, ncol = 2)
  } else{
    p <- p +
      ggplot2::facet_wrap(~names, ncol = 2, scales = "free_y")
  }
  
  return(p)
}

# Filter and prepare core dataframe to mimic theft name cleaning

all_features_filt <- all_features %>%
  filter(names %in% c("values__autocorrelation__lag_5", "0_Standard deviation")) %>%
  mutate(names = gsub("__", "_", names),
         names = gsub(" ", "_", names),
         names = stringr::str_to_lower(names),
         names = paste0(stringr::str_to_upper(method), "_", names), # Cleaning up to replicate {theft} processing
         names = factor(names, levels = c("TSFRESH_values_autocorrelation_lag_5", "TSFEL_0_standard_deviation"), ordered = TRUE)) # Order according to results table for matrix

# Draw plot

ViolinPlots <- plot_feature_discrimination(all_features_filt, id_var = "id", group_var = "group",
                                           normalise = FALSE, method = "RobustSigmoid") +
  theme(text = element_text(size = 14))

#-------------------------
# Multi-feature classifier
#-------------------------

# Run function

mf_results <- fit_multi_feature_classifier(all_features_z, 
                                           id_var = "id", 
                                           group_var = "group",
                                           by_set = TRUE, 
                                           test_method = "svmLinear", 
                                           use_balanced_accuracy = FALSE,
                                           use_k_fold = TRUE, 
                                           num_folds = 10, 
                                           use_empirical_null = TRUE, 
                                           null_testing_method = "ModelFreeShuffles",
                                           p_value_method = "gaussian", 
                                           num_permutations = 1000, 
                                           seed = 123)

FeatureSetResultsPlot <- mf_results$FeatureSetResultsPlot  +
  theme(text = element_text(size = 14)) +
  scale_y_continuous(limits = c(50, 80),
                     breaks = c(50, 60, 70, 80),
                     labels = function(x)paste0(x, "%")) +
  labs(subtitle = "Number of features is indicated in parentheses. Error bars are +/- 1 times SD")

#--------------------------------
# Multi-feature classifier sample
#--------------------------------

# Generate subset of IDs

set.seed(123)

ids_mf <- c("hippocampus", "epileptogenic") %>%
  purrr::map(~ draw_samples(data = tmp, group_name = .x, n = 14)) %>%
  unlist()

all_features_z_sample <- all_features_z %>%
  filter(id %in% ids_mf)

# Fit classifiers

mf_results_sample <- fit_multi_feature_classifier(all_features_z_sample, 
                                                  id_var = "id", 
                                                  group_var = "group",
                                                  by_set = TRUE, 
                                                  test_method = "svmLinear", 
                                                  use_balanced_accuracy = FALSE,
                                                  use_k_fold = TRUE, 
                                                  num_folds = 10, 
                                                  use_empirical_null = TRUE, 
                                                  null_testing_method = "ModelFreeShuffles",
                                                  p_value_method = "gaussian", 
                                                  num_permutations = 10000, 
                                                  seed = 123)

TestStatistics_sample <- mf_results_sample$TestStatistics %>%
  dplyr::select(c(method, accuracy, p_value_accuracy)) %>%
  mutate(accuracy = accuracy * 100,
         p_value_accuracy_adj = p.adjust(p_value_accuracy, method = "holm")) %>%
  arrange(-accuracy)
