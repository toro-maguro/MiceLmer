#' @title Extract random effects from lmer with multiple imputation datasets using mice
#' @description \code{ExtractRandomEffect} Return dataframe including random effects at lmer results with mice
#'
#' @param imputed_datasets mice output.
#' @param lmer_formula Formula for lmer. The class of this parameter has to be character.
#' @param group_v1 Clustering variable. If you applying three-level multilevel model, group_v1 has to be a level2 variable.
#' @param group_v2 Clustering variable at three-level multilevel models. group_v2 has to be a level3 variable. If you applying two-level multilevel model, group_v2 has to be NA.
#' @return data.frame including random effects at lmer results with mice.
#' @export
#' @name MiceLmer
#' @import tidyverse dplyr stats lme4 lmerTest broom.mixed

ExtractRandomEffect <- function(imputed_datasets, lmer_formula, group_v1, group_v2 = NA){
  formula <- as.formula(lmer_formula)
  lmer_results <- lapply(1:imputed_datasets$m, function(i){
    dfm <- complete(imputed_datasets, action=i)
    lmer(formula = formula,
         data = dfm,
         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
         REML = TRUE)
  })

  n_imputation <- imputed_datasets$m
  m <- c(1:n_imputation)
  res_group_v1 <- numeric(n_imputation)
  res_group_v2 <- numeric(n_imputation)
  res_residual <- numeric(n_imputation)
  res_total_variance <- numeric(n_imputation)

  if (is.na(group_v2) == TRUE){
    for (i in 1:n_imputation){
      lmer_summary <- summary(lmer_results[[i]])
      var_summary <- as.data.frame(lmer_summary$varcor)

      d_temp_group_v1 <- dplyr::filter(var_summary, grp == group_v1 & var1 == "(Intercept)" & is.na(var2) == TRUE)
      d_temp_residual <- dplyr::filter(var_summary, grp == "Residual" & is.na(var1) == TRUE & is.na(var2) == TRUE)

      res_group_v1[i] <- d_temp_group_v1[,5]
      res_residual[i] <- d_temp_residual[,5]
      res_total_variance[i] <- sqrt(sum(var_summary$vcov)) # total varianceを取ってSD計算
    }
    df_group_v1 <- data.frame(m, group = group_v1, sd_random_effect = res_group_v1)
    df_residual <- data.frame(m, group = "Residual", sd_random_effect = res_residual)
    df_total_variance <- data.frame(m, group = "SD_TotalVariance", sd_random_effect = res_total_variance)

    df <- rbind(df_group_v1, df_residual, df_total_variance)
  } else {
    for (i in 1:n_imputation){
      lmer_summary <- summary(lmer_results[[i]])
      var_summary <- as.data.frame(lmer_summary$varcor)

      d_temp_group_v1 <- dplyr::filter(var_summary, grp == group_v1 & var1 == "(Intercept)" & is.na(var2) == TRUE) # レベルの低い方から
      d_temp_group_v2 <- dplyr::filter(var_summary, grp == group_v2 & var1 == "(Intercept)" & is.na(var2) == TRUE) # レベルの高い (より広範囲の) 変数へ
      d_temp_residual <- dplyr::filter(var_summary, grp == "Residual" & is.na(var1) == TRUE & is.na(var2) == TRUE)

      res_group_v1[i] <- d_temp_group_v1[,5]
      res_group_v2[i] <- d_temp_group_v2[,5]
      res_residual[i] <- d_temp_residual[,5]
      res_total_variance[i] <- sqrt(sum(var_summary$vcov)) # total varianceを取ってSD計算
    }
    df_group_v1 <- data.frame(m, group = group_v1, sd_random_effect = res_group_v1)
    df_group_v2 <- data.frame(m, group = group_v2, sd_random_effect = res_group_v2)
    df_residual <- data.frame(m, group = "Residual", sd_random_effect = res_residual)
    df_total_variance <- data.frame(m, group = "SD_TotalVariance", sd_random_effect = res_total_variance)

    df <- rbind(df_group_v1, df_group_v2, df_residual, df_total_variance)
  }

  return(df)
}

#' @title Plot distributions of random effects
#' @description \code{PlotRandomEffectDistribution} To confirm the parameters are following normal distributions (Rubin's rule), this function shows distributions of random effects parameters
#'
#' @param data output of ExtractRandomEffect function.
#' @param group_v1 Clustering variable. If you applying three-level multilevel model, group_v1 has to be a level2 variable.
#' @param group_v2 Clustering variable at three-level multilevel models. group_v2 has to be a level3 variable. If you applying two-level multilevel model, group_v2 has to be NA.
#' @return distributions of random effects
#' @export
#' @name MiceLmer
#' @import tidyverse patchwork dplyr stats

PlotRandomEffectDistribution <- function(data, group_v1, group_v2 = NA){
  d_group_v1 <- data %>%
    dplyr::filter(group == group_v1)
  p1 <- ggplot(data = d_group_v1, aes(x = sd_random_effect)) +
    geom_histogram() +
    theme_bw() +
    xlab(group_v1)

  d_residual <- data %>%
    dplyr::filter(group == "Residual")
  p_residual <- ggplot(data = d_residual, aes(x = sd_random_effect)) +
    geom_histogram() +
    theme_bw() +
    xlab("Residual")

  d_total <- data %>%
    dplyr::filter(group == "SD_TotalVariance")
  p_total <- ggplot(data = d_total, aes(x = sd_random_effect)) +
    geom_histogram() +
    theme_bw() +
    xlab("SD_TotalVariance")

  if (is.na(group_v2) == TRUE){
    plot <- p1 / (p_residual | p_total) +
      plot_layout(ncol = 1) +
      plot_annotation(title = "distributions of standard deviation in random effects (intercept)")

  } else {
    d2 <- data %>% dplyr::filter(group == group_v2)
    p2 <- ggplot(data = d2, aes(x = sd_random_effect)) +
      geom_histogram() +
      theme_bw() +
      xlab(group_v2)

    plot <- (p1 | p2) / (p_residual | p_total) +
      plot_layout(ncol = 1) +
      plot_annotation(title = "distributions of standard deviation in random effects (intercept)")

  }
  return(plot)
}

#' @title Combine random effect parameters calculated from multiply imputed datasets based on Rubin's rule
#' @description \code{GetProportionOfRandomEffect} Return random effect parameters (variance and standard deviation) and proportions of the total variance
#'
#' @param data output of ExtractRandomEffect function.
#' @return dataframe including random effect parameters (variance and standard deviation) and proportions of the total variance
#' @export
#' @name MiceLmer
#' @import tidyverse dplyr

GetProportionOfRandomEffect <- function(data){
  d <- data %>%
    dplyr::filter(sd_random_effect > 0.00001) %>%
    group_by(group) %>%
    summarise(sd_mean_value = mean(sd_random_effect)) %>%
    mutate(variance_pooled = sd_mean_value^2) %>%
    ungroup() %>%
    mutate(sum_variance = max(variance_pooled)) %>%
    mutate(proportion_variance = variance_pooled / sum_variance) %>%
    arrange(proportion_variance)
  return(d)
}
