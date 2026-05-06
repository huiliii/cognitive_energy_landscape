rm(list = ls())
library(gamlss)
library(dplyr)
library(ggplot2)
library(tidyr)
library(paletteer)
library(splines)

# Functions ---------------------------------------------------------------
get_prediction = function(model, data, xvar, yvar){
  # For models with random effects, manually calculate z-scores
  fitted_mu = fitted(model, what = "mu")
  fitted_sigma = fitted(model, what = "sigma")
  zscores = (data[[yvar]] - fitted_mu) / fitted_sigma
  
  ex_ind = is.na(zscores) | zscores > 1000 | zscores < -1000
  zscores = zscores[!ex_ind]
  data_clean = data[!ex_ind,]
  cents = pnorm(zscores)
  centpreds = quantile(data_clean[[yvar]], probs = cents)
  
  # Simplified prediction - just get mu
  termpreds = data.frame(mu = fitted_mu[!ex_ind])
  
  # Skip slope calculation - not needed for AIC/BIC comparison
  slopes = data.frame(mu = numeric(0))
  
  return(list("zscores" = zscores,
              "centpreds" = centpreds,
              "termpreds" = termpreds,
              "excluded" = ex_ind,
              "slopes" = slopes,
              "cscores" = cents))
}

assess_gamlss_fit = function(model, prediction, actual_data, yvar) {
  # Calculate residuals and predictions
  residuals = residuals(model)
  predicted = prediction$termpreds$mu
  
  # Get actual values, excluding those that were excluded in prediction
  actual = actual_data[[yvar]][!prediction$excluded]
  
  # Ensure lengths match
  min_len = min(length(actual), length(predicted), length(residuals))
  actual = actual[1:min_len]
  predicted = predicted[1:min_len]
  residuals = residuals[1:min_len]
  
  # 1. Shapiro-Wilk test for normality of residuals (sample if too large)
  if (length(residuals) > 5000) {
    residuals_sample = sample(residuals, 5000)
    sw_test = shapiro.test(residuals_sample)
  } else if (length(residuals) >= 3) {
    sw_test = shapiro.test(residuals)
  } else {
    sw_test = list(p.value = NA)
  }
  sw_pvalue = sw_test$p.value
  
  # 2. Kolmogorov-Smirnov test
  ks_test = ks.test(residuals, "pnorm")
  ks_pvalue = ks_test$p.value
  
  # 3. Basic error metrics
  mae = mean(abs(actual - predicted), na.rm = TRUE)
  rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  mape = mean(abs((actual - predicted)/actual), na.rm = TRUE) * 100
  r_squared = 1 - (sum((actual - predicted)^2, na.rm = TRUE)/sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE))
  
  # 4. Z-score coverage
  z_scores = prediction$zscores
  within_1sd = mean(abs(z_scores) <= 1, na.rm = TRUE) * 100
  within_2sd = mean(abs(z_scores) <= 2, na.rm = TRUE) * 100
  within_3sd = mean(abs(z_scores) <= 3, na.rm = TRUE) * 100
  
  # 5. Correlation test between actual and predicted values
  if (length(actual) > 2 && length(predicted) > 2 && 
      sum(!is.na(actual)) > 2 && sum(!is.na(predicted)) > 2) {
    cor_test = cor.test(actual, predicted)
    cor_pvalue = cor_test$p.value
  } else {
    cor_pvalue = NA
  }
  
  # 6. Get model fit statistics (AIC and BIC)
  # Ensure all are scalar values
  gof = ifelse(is.null(model$aic), NA, model$aic[1])
  bic = ifelse(is.null(model$sbc), NA, model$sbc[1])
  dev = ifelse(is.null(model$deviance), NA, sum(model$deviance))  # Sum deviance if it's a vector
  
  # Create results list with p-values
  results = list(
    error_metrics = data.frame(
      Metric = c("MAE", "RMSE", "MAPE", "R-squared"),
      Value = c(mae, rmse, mape, r_squared),
      stringsAsFactors = FALSE
    ),
    coverage = data.frame(
      SD_Range = c("Within 1 SD", "Within 2 SD", "Within 3 SD"),
      Percentage = c(within_1sd, within_2sd, within_3sd),
      stringsAsFactors = FALSE
    ),
    statistical_tests = data.frame(
      Test = c("Shapiro-Wilk", "Kolmogorov-Smirnov", "Correlation"),
      P_value = c(sw_pvalue, ks_pvalue, cor_pvalue),
      stringsAsFactors = FALSE
    ),
    model_fit = data.frame(
      Metric = c("AIC", "BIC", "Deviance"),
      Value = c(gof, bic, dev),
      stringsAsFactors = FALSE
    ),
    model_summary = summary(model)
  )
  
  return(results)
}

# Function to fit model based on family type
fit_gamlss_model = function(yvar, family_dist, family_name, data) {
  # Define formula
  base_formula = as.formula(paste0("`", yvar, "` ~ bs(pseudo_log_age, df=3) + random(study)"))
  
  # Different families require different parameter specifications
  # 2-parameter families: NO, TF, PE, GA, IG, LOGNO
  # 3-parameter families: TF (can use nu)
  # 4-parameter families: BCT, BCPE, JSU, ST3
  
  if (family_name %in% c("NO", "PE", "GA", "IG", "LOGNO")) {
    # 2-parameter families: only mu and sigma
    model <- gamlss(base_formula,
                    sigma.formula = ~bs(pseudo_log_age, df=3),
                    family = family_dist,
                    data = data,
                    control = gamlss.control(n.cyc = 100, trace = FALSE))
  } else if (family_name == "TF") {
    # TF can have 2 or 3 parameters
    model <- gamlss(base_formula,
                    sigma.formula = ~bs(pseudo_log_age, df=3),
                    nu.formula = ~1,
                    family = family_dist,
                    data = data,
                    control = gamlss.control(n.cyc = 100, trace = FALSE))
  } else if (family_name %in% c("BCT", "BCPE", "JSU", "ST3")) {
    # 4-parameter families: mu, sigma, nu, tau
    model <- gamlss(base_formula,
                    sigma.formula = ~bs(pseudo_log_age, df=3),
                    nu.formula = ~1,
                    tau.formula = ~1,
                    family = family_dist,
                    data = data,
                    control = gamlss.control(n.cyc = 100, trace = FALSE))
  } else {
    # Default: try with all parameters
    model <- gamlss(base_formula,
                    sigma.formula = ~bs(pseudo_log_age, df=3),
                    nu.formula = ~1,
                    tau.formula = ~1,
                    family = family_dist,
                    data = data,
                    control = gamlss.control(n.cyc = 100, trace = FALSE))
  }
  
  return(model)
}

# Script ------------------------------------------------------------------
dataf <- read.csv("/Users/huilisun/Library/CloudStorage/OneDrive-PennO365/ControlCog/ce_baseline_keywords_site.csv")

# Define pseudo-log transformation function
pseudo_log_scale <- function(x, breakpoint = 0) {
  result <- ifelse(x <= breakpoint,
                   x,
                   breakpoint + log1p(x - breakpoint))
  return(result)
}

# Create pseudo-log age variable
dataf$pseudo_log_age <- pseudo_log_scale(dataf$age, breakpoint = 0)

# Ensure study is a factor
dataf$study <- as.factor(dataf$study)

# Create output directory if it doesn't exist
output_dir <- "/Users/huilisun/Library/CloudStorage/OneDrive-PennO365/ControlCog/Output/figs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define list of family distributions to test
family_list <- list(
  NO = NO,          # Normal (2 parameters: mu, sigma)
  TF = TF,          # t Family (3 parameters: mu, sigma, nu)
  PE = PE,          # Power Exponential (2 parameters: mu, sigma)
  LOGNO = LOGNO,    # Log Normal (2 parameters: mu, sigma)
  GA = GA,          # Gamma (2 parameters: mu, sigma)
  IG = IG,          # Inverse Gaussian (2 parameters: mu, sigma)
  BCT = BCT,        # Box-Cox t (4 parameters: mu, sigma, nu, tau)
  BCPE = BCPE,      # Box-Cox Power Exponential (4 parameters)
  JSU = JSU,        # Johnson's Su (4 parameters)
  ST3 = ST3         # Skew t type 3 (4 parameters)
)

# Get all node_energy_mean columns (columns 4 to 104)
node_energy_cols <- names(dataf)[4:103]
node_energy_cols <- node_energy_cols[grepl("^node_energy_mean_", node_energy_cols)]

# Create a summary dataframe to store results
summary_results <- data.frame(
  Variable = character(),
  Family = character(),
  AIC = numeric(),
  BIC = numeric(),
  Deviance = numeric(),
  MAE = numeric(),
  RMSE = numeric(),
  R_squared = numeric(),
  Converged = logical(),
  Error_Message = character(),
  stringsAsFactors = FALSE
)

# Create a dataframe to store best model for each variable
best_models <- data.frame(
  Variable = character(),
  Best_Family_AIC = character(),
  Best_AIC = numeric(),
  Best_Family_BIC = character(),
  Best_BIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each variable
cat("Processing", length(node_energy_cols), "variables with", length(family_list), "family distributions...\n\n")

for (i in seq_along(node_energy_cols)) {
  yvar <- node_energy_cols[i]
  
  cat(sprintf("\n========== Processing %d/%d: %s ==========\n", i, length(node_energy_cols), yvar))
  
  # Store results for this variable across all families
  variable_results <- data.frame(
    Family = character(),
    AIC = numeric(),
    BIC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each family distribution
  for (family_name in names(family_list)) {
    family_dist <- family_list[[family_name]]
    
    cat(sprintf("  Testing family: %s ... ", family_name))
    
    # Initialize tracking variables
    converged <- TRUE
    error_msg <- ""
    
    tryCatch({
      # Fit the model using appropriate specification for this family
      model <- fit_gamlss_model(yvar, family_dist, family_name, dataf)
      
      # Get predictions
      prediction <- get_prediction(model, data = dataf,
                                   xvar = "pseudo_log_age", yvar = yvar)
      
      # Assess fit
      fit_assessment <- assess_gamlss_fit(model, prediction, dataf, yvar = yvar)
      
      # Extract AIC and BIC
      aic_val <- ifelse(is.null(model$aic), NA, model$aic)
      bic_val <- ifelse(is.null(model$sbc), NA, model$sbc)
      
      # Store in variable results
      variable_results <- rbind(variable_results, data.frame(
        Family = family_name,
        AIC = aic_val,
        BIC = bic_val,
        stringsAsFactors = FALSE
      ))
      
      # Store summary results
      new_row <- data.frame(
        Variable = yvar,
        Family = family_name,
        AIC = aic_val,
        BIC = bic_val,
        Deviance = ifelse(is.null(model$deviance), NA, model$deviance),
        MAE = ifelse(length(fit_assessment$error_metrics$Value) >= 1,
                     fit_assessment$error_metrics$Value[1], NA),
        RMSE = ifelse(length(fit_assessment$error_metrics$Value) >= 2,
                      fit_assessment$error_metrics$Value[2], NA),
        R_squared = ifelse(length(fit_assessment$error_metrics$Value) >= 4,
                           fit_assessment$error_metrics$Value[4], NA),
        Converged = TRUE,
        Error_Message = "",
        stringsAsFactors = FALSE
      )
      summary_results <- rbind(summary_results, new_row)
      
      cat("✓ AIC:", round(aic_val, 2), "BIC:", round(bic_val, 2), "\n")
      
    }, error = function(e) {
      error_msg <- paste(as.character(e$message), collapse = " ")
      error_trace <- paste(deparse(e$call), collapse = " ")
      
      cat("✗ Error:", substr(error_msg, 1, 80), "\n")
      cat("   Trace:", substr(error_trace, 1, 80), "\n")
      
      # Store error in summary
      new_row <- data.frame(
        Variable = yvar,
        Family = family_name,
        AIC = NA_real_,
        BIC = NA_real_,
        Deviance = NA_real_,
        MAE = NA_real_,
        RMSE = NA_real_,
        R_squared = NA_real_,
        Converged = FALSE,
        Error_Message = substr(error_msg, 1, 200),
        stringsAsFactors = FALSE
      )
      summary_results <<- rbind(summary_results, new_row)
    })
  }
  
  # Find best model for this variable based on AIC and BIC
  if (nrow(variable_results) > 0) {
    valid_results <- variable_results[!is.na(variable_results$AIC) & !is.na(variable_results$BIC), ]
    
    if (nrow(valid_results) > 0) {
      best_aic_idx <- which.min(valid_results$AIC)
      best_bic_idx <- which.min(valid_results$BIC)
      
      best_models <- rbind(best_models, data.frame(
        Variable = yvar,
        Best_Family_AIC = valid_results$Family[best_aic_idx],
        Best_AIC = valid_results$AIC[best_aic_idx],
        Best_Family_BIC = valid_results$Family[best_bic_idx],
        Best_BIC = valid_results$BIC[best_bic_idx],
        stringsAsFactors = FALSE
      ))
      
      cat(sprintf("\n  → Best by AIC: %s (%.2f)\n", 
                  valid_results$Family[best_aic_idx], 
                  valid_results$AIC[best_aic_idx]))
      cat(sprintf("  → Best by BIC: %s (%.2f)\n", 
                  valid_results$Family[best_bic_idx], 
                  valid_results$BIC[best_bic_idx]))
    }
  }
}

# Save summary results
summary_filename <- file.path(output_dir, "gamlss_family_comparison_summary.csv")
write.csv(summary_results, summary_filename, row.names = FALSE)

# Save best models summary
best_models_filename <- file.path(output_dir, "gamlss_best_models_summary.csv")
write.csv(best_models, best_models_filename, row.names = FALSE)

# Generate summary statistics
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Analysis Complete!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Total variable-family combinations processed:", nrow(summary_results), "\n")
cat("Successfully converged:", sum(summary_results$Converged), "\n")
cat("Failed to converge:", sum(!summary_results$Converged), "\n")

output_dir <- "/Users/huilisun/Library/CloudStorage/OneDrive-PennO365/ControlCog/Output/gamlss_family_selection"
summary_results <- read.csv(file.path(output_dir, "gamlss_family_comparison_summary.csv"))
best_models <- read.csv(file.path(output_dir, "gamlss_best_models_summary.csv"))
# Summary by family
cat("\n--- Convergence by Family ---\n")
family_summary <- summary_results %>%
  group_by(Family) %>%
  summarise(
    Total = n(),
    Converged = sum(Converged),
    Failed = sum(!Converged),
    Mean_AIC = mean(AIC, na.rm = TRUE),
    Mean_BIC = mean(BIC, na.rm = TRUE)
  )
print(family_summary)

# Most frequent best families
if (nrow(best_models) > 0) {
  cat("\n--- Most Frequently Selected Families ---\n")
  cat("By AIC:\n")
  print(table(best_models$Best_Family_AIC))
  cat("\nBy BIC:\n")
  print(table(best_models$Best_Family_BIC))
}

cat("\nResults saved to:", summary_filename, "\n")
cat("Best models saved to:", best_models_filename, "\n")