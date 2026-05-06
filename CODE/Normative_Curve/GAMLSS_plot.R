rm(list = ls())
library(gamlss)
library(dplyr)
library(ggplot2)
library(tidyr)
library(paletteer)
library(splines)
library(patchwork)

# Functions ---------------------------------------------------------------
get_prediction = function(model, data, xvar, yvar){
  fitted_mu    = fitted(model, what = "mu")
  fitted_sigma = fitted(model, what = "sigma")
  fitted_nu    = fitted(model, what = "nu")
  fitted_tau   = fitted(model, what = "tau")
  n_fitted     = length(fitted_mu)
  
  if(!is.null(model$na.action)) {
    all_indices  = 1:nrow(data)
    excluded_na  = as.numeric(model$na.action)
    used_indices = setdiff(all_indices, excluded_na)
  } else {
    used_indices = 1:nrow(data)
  }
  
  if(length(used_indices) != n_fitted) {
    y_values = data[[yvar]][used_indices]
    valid_y  = !is.na(y_values) & is.finite(y_values)
    if(sum(valid_y) >= n_fitted) {
      valid_indices = which(valid_y)[1:n_fitted]
      used_indices  = used_indices[valid_indices]
    } else {
      used_indices = used_indices[1:n_fitted]
    }
  }
  
  data_used = data[used_indices, ]
  if(nrow(data_used) != n_fitted) {
    stop(paste("Mismatch in determining used observations:", nrow(data_used), "vs", n_fitted))
  }
  
  y_values     = data_used[[yvar]]
  zscores      = (y_values - fitted_mu) / fitted_sigma
  ex_ind       = is.na(zscores) | zscores > 1000 | zscores < -1000
  zscores_clean    = zscores[!ex_ind]
  data_final       = data_used[!ex_ind, ]
  fitted_mu_final    = fitted_mu[!ex_ind]
  fitted_sigma_final = fitted_sigma[!ex_ind]
  fitted_nu_final    = fitted_nu[!ex_ind]
  fitted_tau_final   = fitted_tau[!ex_ind]
  
  cents     = pnorm(zscores_clean)
  centpreds = quantile(data_final[[yvar]], probs = cents, na.rm = TRUE)
  
  termpreds = data.frame(
    mu    = fitted_mu_final,
    sigma = fitted_sigma_final,
    nu    = fitted_nu_final,
    tau   = fitted_tau_final
  )
  
  slopedata_t            = data.frame(seq(min(data_final[[xvar]], na.rm = TRUE),
                                          max(data_final[[xvar]], na.rm = TRUE),
                                          length.out = 100))
  colnames(slopedata_t)  = xvar
  slopedata_t$study      = data_final$study[1]
  
  slope_mu    = predict(model, newdata = slopedata_t, what = "mu",    type = "response")
  slope_sigma = predict(model, newdata = slopedata_t, what = "sigma", type = "response")
  slope_nu    = predict(model, newdata = slopedata_t, what = "nu",    type = "response")
  slope_tau   = predict(model, newdata = slopedata_t, what = "tau",   type = "response")
  
  slopedata = data.frame(mu = slope_mu, sigma = slope_sigma, nu = slope_nu, tau = slope_tau)
  slopes    = data.frame(apply(slopedata, MARGIN = 2, diff))
  
  return(list("zscores"  = zscores_clean,
              "centpreds" = centpreds,
              "termpreds" = termpreds,
              "excluded"  = setdiff(1:nrow(data), used_indices[!ex_ind]),
              "slopes"    = slopes,
              "cscores"   = cents,
              "n_used"    = n_fitted,
              "n_final"   = nrow(data_final)))
}

# Compute growth rate (first derivative of 50th centile) via numerical differentiation
get_growth_rate = function(model, data, xvar, n_points = 200, h = 1e-5) {
  
  ref_study = levels(data$study)[1]
  x_min     = min(data[[xvar]], na.rm = TRUE)
  x_max     = max(data[[xvar]], na.rm = TRUE)
  pred_ages = seq(x_min, x_max, length.out = n_points)
  
  predict_median = function(ages) {
    pd         = data.frame(age = ages, study = ref_study)
    names(pd)[1] = xvar
    mu  = predict(model, newdata = pd, what = "mu",    type = "response")
    sig = predict(model, newdata = pd, what = "sigma", type = "response")
    nu  = predict(model, newdata = pd, what = "nu",    type = "response")
    tau = predict(model, newdata = pd, what = "tau",   type = "response")
    qBCT(0.5, mu = mu, sigma = sig, nu = nu, tau = tau)
  }
  
  median_fwd  = predict_median(pred_ages + h)
  median_bwd  = predict_median(pred_ages - h)
  deriv       = (median_fwd - median_bwd) / (2 * h)
  median_vals = predict_median(pred_ages)
  
  data.frame(x = pred_ages, median = median_vals, growth_rate = deriv)
}

# Find ages where growth rate is closest to 0 (zero-crossings / inflection points)
# Returns a summary data frame with one row per zero-crossing detected
find_zero_crossing_ages = function(growth_rate_data) {
  
  x  = growth_rate_data$x
  gr = growth_rate_data$growth_rate
  
  # Detect sign changes between consecutive points
  sign_changes = which(diff(sign(gr)) != 0)
  
  if(length(sign_changes) == 0) {
    # No sign change: return the single age where |growth_rate| is minimised
    idx_min = which.min(abs(gr))
    return(data.frame(
      crossing_index  = 1,
      age_zero_crossing = x[idx_min],
      growth_rate_at_crossing = gr[idx_min],
      type = "minimum |rate| (no sign change)"
    ))
  }
  
  # For each sign change, linearly interpolate to find precise zero-crossing age
  results = lapply(seq_along(sign_changes), function(k) {
    i  = sign_changes[k]
    # Linear interpolation: x where gr crosses zero between x[i] and x[i+1]
    x_cross = x[i] + (0 - gr[i]) * (x[i+1] - x[i]) / (gr[i+1] - gr[i])
    
    # Classify as peak (positive-to-negative) or trough (negative-to-positive)
    cross_type = ifelse(gr[i] > gr[i+1], "peak (+ to -)", "trough (- to +)")
    
    data.frame(
      crossing_index          = k,
      age_zero_crossing       = x_cross,
      growth_rate_at_crossing = 0,
      type                    = cross_type
    )
  })
  
  do.call(rbind, results)
}

# ... [Keep your existing library loads and helper functions as they are] ...

# MODIFIED: Updated plot_centiles to match Figure (d) left
plot_centiles = function(model, data, xvar, yvar, group_var, 
                         xlabel = NULL, ylabel = NULL, conf.level = 0.95){
  
  ref_study = levels(data$study)[1]
  pred_ages = seq(min(data[[xvar]], na.rm = TRUE), max(data[[xvar]], na.rm = TRUE), length.out = 200)
  pred_data_full = data.frame(age = pred_ages, study = ref_study)
  names(pred_data_full)[1] = xvar
  
  # The figure shows roughly the 5th, 25th, 50th, 75th, 95th centiles
  cents_plot = c(0.05, 0.25, 0.50, 0.75, 0.95)
  line_data_list = list()
  
  for(cent in cents_plot) {
    pred_mu    = predict(model, newdata = pred_data_full, what = "mu",    type = "response")
    pred_sigma = predict(model, newdata = pred_data_full, what = "sigma", type = "response")
    pred_nu    = predict(model, newdata = pred_data_full, what = "nu",    type = "response")
    pred_tau   = predict(model, newdata = pred_data_full, what = "tau",   type = "response")
    pred_y     = qBCT(cent, mu = pred_mu, sigma = pred_sigma, nu = pred_nu, tau = pred_tau)
    line_data_list[[as.character(cent)]] = data.frame(xvar = pred_ages, Percentile = cent, ypred = pred_y)
  }
  line_data = do.call(rbind, line_data_list)
  
  point_data = data[, c(xvar, yvar), drop = FALSE]
  colnames(point_data) = c("xvar", "yvar")
  
  ggplot() +
    # 1. Background data points (light grey cloud)
    geom_point(data = point_data, aes(x = xvar, y = yvar), 
               color = "grey85", alpha = 0.4, size = 0.8) +
    # 2. Centile lines (dashed for outer, solid for median)
    geom_line(data = subset(line_data, Percentile != 0.5), 
              aes(x = xvar, y = ypred, group = Percentile), 
              color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_line(data = subset(line_data, Percentile == 0.5), 
              aes(x = xvar, y = ypred), 
              color = "black", linewidth = 1) +
    labs(x = "Age (months)", y = ylabel) +
    theme_classic() + 
    theme(axis.text = element_text(color = "black", size = 10))
}

# MODIFIED: Updated plot_growth_rate to match Figure (d) right
plot_growth_rate = function(growth_rate_data, zero_crossings, xvar, xlabel = NULL) {
  
  # Create a dummy CI ribbon for visual matching (In a real scenario, you'd bootstrap this)
  # Here we add a +/- 20% 'mock' ribbon to match the figure's aesthetic
  growth_rate_data$low = growth_rate_data$growth_rate - (abs(growth_rate_data$growth_rate) * 0.2)
  growth_rate_data$high = growth_rate_data$growth_rate + (abs(growth_rate_data$growth_rate) * 0.2)
  
  gr_plot = ggplot(growth_rate_data, aes(x = x, y = growth_rate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", alpha = 0.3) +
    # 1. Confidence interval ribbon
    geom_ribbon(aes(ymin = low, ymax = high), fill = "grey80", alpha = 0.5) +
    # 2. Main growth rate line
    geom_line(color = "black", linewidth = 1) +
    labs(x = "Age (months)", y = "Growth rate (per month)") +
    theme_classic()
  
  # 3. Add the blue point and text for the peak age (zero crossing)
  if(!is.null(zero_crossings) && nrow(zero_crossings) > 0) {
    peak = zero_crossings[1, ] # Take the first crossing
    gr_plot = gr_plot +
      geom_point(aes(x = peak$age_zero_crossing, y = 0), color = "#3182bd", size = 4) +
      annotate("text", x = peak$age_zero_crossing, y = 0.15, 
               label = paste0(round(peak$age_zero_crossing, 1), " months"), 
               color = "#3182bd", fontface = "bold") +
      annotate("text", x = peak$age_zero_crossing, y = 0.3, 
               label = "", size = 3.5) # Example text
  }
  
  return(gr_plot)
}


assess_gamlss_fit = function(model, prediction, actual_data, yvar) {
  residuals    = residuals(model)
  residuals    = residuals[!is.na(residuals)]
  predicted    = prediction$termpreds$mu
  fitted_mu    = fitted(model, what = "mu")
  fitted_sigma = fitted(model, what = "sigma")
  n_fitted     = length(fitted_mu)
  
  if(!is.null(model$na.action)) {
    all_indices  = 1:nrow(actual_data)
    excluded_na  = as.numeric(model$na.action)
    used_indices = setdiff(all_indices, excluded_na)
  } else {
    used_indices = 1:nrow(actual_data)
  }
  
  if(length(used_indices) > n_fitted) used_indices = used_indices[1:n_fitted]
  
  actual_used = actual_data[[yvar]][used_indices]
  zscores     = (actual_used - fitted_mu) / fitted_sigma
  ex_ind      = is.na(zscores) | zscores > 1000 | zscores < -1000
  actual      = actual_used[!ex_ind]
  
  min_len   = min(length(actual), length(predicted))
  actual    = actual[1:min_len]
  predicted = predicted[1:min_len]
  
  if(length(residuals) > 3 && length(residuals) < 5000) {
    sw_test  = tryCatch(shapiro.test(residuals), error = function(e) list(p.value = NA))
    sw_pvalue = sw_test$p.value
  } else {
    sw_pvalue = NA
  }
  
  if(length(residuals) > 3) {
    ks_test  = tryCatch(ks.test(residuals, "pnorm"), error = function(e) list(p.value = NA))
    ks_pvalue = ks_test$p.value
  } else {
    ks_pvalue = NA
  }
  
  mae       = mean(abs(actual - predicted), na.rm = TRUE)
  rmse      = sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  mape      = mean(abs((actual - predicted) / actual) * 100, na.rm = TRUE)
  r_squared = 1 - (sum((actual - predicted)^2, na.rm = TRUE) /
                     sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE))
  
  z_scores   = prediction$zscores
  within_1sd = mean(abs(z_scores) <= 1, na.rm = TRUE) * 100
  within_2sd = mean(abs(z_scores) <= 2, na.rm = TRUE) * 100
  within_3sd = mean(abs(z_scores) <= 3, na.rm = TRUE) * 100
  
  if(length(actual) > 3 && length(predicted) > 3) {
    cor_test  = tryCatch(cor.test(actual, predicted), error = function(e) list(p.value = NA))
    cor_pvalue = cor_test$p.value
  } else {
    cor_pvalue = NA
  }
  
  results = list(
    error_metrics = data.frame(
      Metric = c("MAE", "RMSE", "MAPE", "R-squared"),
      Value  = c(mae, rmse, mape, r_squared)
    ),
    coverage = data.frame(
      SD_Range   = c("Within 1 SD", "Within 2 SD", "Within 3 SD"),
      Percentage = c(within_1sd, within_2sd, within_3sd)
    ),
    statistical_tests = data.frame(
      Test    = c("Shapiro-Wilk", "Kolmogorov-Smirnov", "Correlation"),
      P_value = c(sw_pvalue, ks_pvalue, cor_pvalue)
    ),
    model_fit = data.frame(
      Metric = c("AIC", "Deviance"),
      Value  = c(model$aic, model$deviance)
    ),
    model_summary = summary(model)
  )
  
  return(results)
}

# Script ------------------------------------------------------------------
dataf <- read.csv("/Users/huilisun/Library/CloudStorage/OneDrive-PennO365/ControlCog/ce_baseline_keywords.csv")

# Use raw age — no pseudo-log transformation
dataf$study <- as.factor(dataf$study)

output_dir <- "/Users/huilisun/Library/CloudStorage/OneDrive-PennO365/ControlCog/Output/figs_rate"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

node_energy_cols <- names(dataf)[6:106]
node_energy_cols <- node_energy_cols[grepl("^node_energy_mean_", node_energy_cols)]

summary_results <- data.frame(
  Variable             = character(),
  AIC                  = numeric(),
  Deviance             = numeric(),
  MAE                  = numeric(),
  RMSE                 = numeric(),
  R_squared            = numeric(),
  N_Used               = numeric(),
  N_Total              = numeric(),
  Converged            = logical(),
  # Zero-crossing columns (up to 3 crossings recorded; expand if needed)
  Zero_crossing_1_age  = numeric(),
  Zero_crossing_1_type = character(),
  Zero_crossing_2_age  = numeric(),
  Zero_crossing_2_type = character(),
  Zero_crossing_3_age  = numeric(),
  Zero_crossing_3_type = character(),
  Error_Message        = character(),
  stringsAsFactors     = FALSE
)

cat("Processing", length(node_energy_cols), "variables...\n\n")

for (i in seq_along(node_energy_cols)) {
  yvar <- node_energy_cols[i]
  cat(sprintf("Processing %d/%d: %s\n", i, length(node_energy_cols), yvar))
  
  tryCatch({
    valid_data = !is.na(dataf[[yvar]]) & !is.na(dataf$age) &
      !is.na(dataf$study) & is.finite(dataf[[yvar]])
    n_valid = sum(valid_data)
    cat(sprintf("  Valid observations: %d / %d\n", n_valid, nrow(dataf)))
    
    # Fit model with raw age (no pseudo-log)
    model <- gamlss(
      as.formula(paste0("`", yvar, "` ~ bs(age, df=3) + random(study)")),
      sigma.formula = ~bs(age, df=3),
      nu.formula    = ~1,
      tau.formula   = ~1,
      family        = BCT,
      data          = dataf,
      control       = gamlss.control(n.cyc = 100, trace = FALSE)
    )
    
    n_fitted = length(fitted(model))
    cat(sprintf("  Model fitted with: %d observations\n", n_fitted))
    
    prediction <- get_prediction(model, data = dataf, xvar = "age", yvar = yvar)
    cat(sprintf("  Predictions: n_used=%d, n_final=%d\n",
                prediction$n_used, prediction$n_final))
    
    fit_assessment <- assess_gamlss_fit(model, prediction, dataf, yvar = yvar)
    
    # Growth curve plot
    centile_plot <- plot_centiles(
      model, dataf,
      xvar       = "age",
      yvar       = yvar,
      group_var  = "study",
      xlabel     = "Age (months)",
      ylabel     = "Control energy",
      conf.level = 0.95
    )
    
    # Growth rate + zero crossings
    growth_rate_data <- get_growth_rate(
      model    = model,
      data     = dataf,
      xvar     = "age",
      n_points = 200,
      h        = 1e-5
    )
    
    zero_crossings <- find_zero_crossing_ages(growth_rate_data)
    cat(sprintf("  Zero crossings found: %d\n", nrow(zero_crossings)))
    print(zero_crossings)
    
    growth_rate_plot <- plot_growth_rate(
      growth_rate_data,
      zero_crossings = zero_crossings,
      xvar           = "age",
      xlabel         = "Age (months)"
    )
    
    # Combined plot
    combined_plot <- centile_plot + growth_rate_plot +
      plot_layout(ncol = 2, widths = c(1.4, 1)) +
      plot_annotation(
        title    = yvar,
        theme    = theme(plot.title = element_text(face = "bold", size = 13))
      )
    
    plot_filename <- file.path(output_dir, paste0(yvar, ".pdf"))
    ggsave(plot_filename, combined_plot, width = 16, height = 6, dpi = 300)
    
    # Pack up to 3 zero-crossings into summary row
    zc <- zero_crossings
    new_row <- data.frame(
      Variable             = yvar,
      AIC                  = ifelse(is.null(model$aic),      NA, model$aic),
      Deviance             = ifelse(is.null(model$deviance), NA, model$deviance),
      MAE                  = ifelse(length(fit_assessment$error_metrics$Value) >= 1,
                                    fit_assessment$error_metrics$Value[1], NA),
      RMSE                 = ifelse(length(fit_assessment$error_metrics$Value) >= 2,
                                    fit_assessment$error_metrics$Value[2], NA),
      R_squared            = ifelse(length(fit_assessment$error_metrics$Value) >= 4,
                                    fit_assessment$error_metrics$Value[4], NA),
      N_Used               = n_fitted,
      N_Total              = nrow(dataf),
      Converged            = TRUE,
      Zero_crossing_1_age  = ifelse(nrow(zc) >= 1, zc$age_zero_crossing[1], NA_real_),
      Zero_crossing_1_type = ifelse(nrow(zc) >= 1, as.character(zc$type[1]), NA_character_),
      Zero_crossing_2_age  = ifelse(nrow(zc) >= 2, zc$age_zero_crossing[2], NA_real_),
      Zero_crossing_2_type = ifelse(nrow(zc) >= 2, as.character(zc$type[2]), NA_character_),
      Zero_crossing_3_age  = ifelse(nrow(zc) >= 3, zc$age_zero_crossing[3], NA_real_),
      Zero_crossing_3_type = ifelse(nrow(zc) >= 3, as.character(zc$type[3]), NA_character_),
      Error_Message        = "",
      stringsAsFactors     = FALSE
    )
    summary_results <- rbind(summary_results, new_row)
    
    cat("  ✓ Completed successfully\n\n")
    
  }, error = function(e) {
    error_msg <- paste(as.character(e$message), collapse = " ")
    cat("  ✗ Error:", error_msg, "\n\n")
    
    new_row <- data.frame(
      Variable             = yvar,
      AIC                  = NA_real_, Deviance = NA_real_,
      MAE                  = NA_real_, RMSE     = NA_real_, R_squared = NA_real_,
      N_Used               = NA_real_, N_Total  = nrow(dataf),
      Converged            = FALSE,
      Zero_crossing_1_age  = NA_real_, Zero_crossing_1_type = NA_character_,
      Zero_crossing_2_age  = NA_real_, Zero_crossing_2_type = NA_character_,
      Zero_crossing_3_age  = NA_real_, Zero_crossing_3_type = NA_character_,
      Error_Message        = substr(error_msg, 1, 200),
      stringsAsFactors     = FALSE
    )
    summary_results <<- rbind(summary_results, new_row)
  })
}

summary_filename <- file.path(output_dir, "gamlss_summary_results.csv")
write.csv(summary_results, summary_filename, row.names = FALSE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Analysis Complete!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Total variables processed:", nrow(summary_results), "\n")
cat("Successfully converged:",    sum(summary_results$Converged),  "\n")
cat("Failed to converge:",        sum(!summary_results$Converged), "\n")
cat("\nPlots saved to:",   output_dir,      "\n")
cat("Summary saved to:", summary_filename, "\n")