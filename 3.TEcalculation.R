############################################################
## 02_calculate_TE_components.R
## Calculate TE decomposition components using the simulated
## city data and coefficient draws.

FOLDER_NAME<-"YOUR PATH"
############################################################
## 1. Read simulated data and coefficient draws
############################################################

sim_city <- read.csv(
  paste0(FOLDER_NAME,"/Data/sim_city_daily.csv")
)

BTcoef <- read.csv(
  paste0(FOLDER_NAME,"/Output/BTcoef_city.csv")
)

sim_city$date <- as.Date(sim_city$date)

cat("Number of rows in simulated data:", nrow(sim_city), "\n")
cat("Number of coefficient draws:", nrow(BTcoef), "\n")
cat("Number of HW days:", sum(sim_city$HW), "\n")


############################################################
## 2. Define variables used in the Supplement equations
############################################################

## Binary heat-wave variable
HWval <- sim_city$HW

## Reference heat-wave value
hw.ref <- 0

## Temperature indicator.
## Following the empirical projection code, the ozone-mediated
## part is only allowed when temperature is >= 25 C.
Tempval <- sim_city$temp >= 25

## Predicted logO3 under actual HW condition
logO3_pred_HW <- sim_city$logO3

## Predicted logO3 under no-HW / reference condition
logO3_pred_NHW <- sim_city$logO3_nohw_above25 

## Set logO3 variables to zero when temperature is below 25 C.
## This makes PIE and mediated interaction components inactive
## below 25 C.

### YOU SHOULD PICK EITHER OF THESE TWO SET BELOW!

#If you wish to estimate TE for T>=25 vs T=25....
logO3_pred_HW <- ifelse(Tempval == 1, sim_city$logO3, 0)  ## Natural O3
logO3_pred_NHW <- ifelse(Tempval == 1, sim_city$logO3_base, 0) ## Natural O3 - HW's effect on O3

#If you wish to estimate TE for HW=1 vs HW=0....
#logO3_pred_HW <- ifelse(sim_city$HW==1, sim_city$logO3, 0)  ## Natural O3
#logO3_pred_NHW <- ifelse(sim_city$HW==1, sim_city$logO3_nohw_above25, 0) ## Natural O3 - HW's effect on O3

##In this example, you will get the same PIE estimates in the end,
##because (logO3_pred_HW - logO3_pred_NHW) was set to be equal (because HW additionally increases O3, "independent of the temperature value")



## If the heat-mediated ozone contrast is negative, set both
logO3_diff <- logO3_pred_HW - logO3_pred_NHW
logO3_pred_HW <- ifelse(logO3_diff < 0, 0, logO3_pred_HW)
logO3_pred_NHW <- ifelse(logO3_diff < 0, 0, logO3_pred_NHW)

## Reference logO3 value   
### YOU SHOULD PICK EITHER OF THESE TWO BELOW!
o3.ref <- sim_city$log_o3_ref ## if you use this, you estimate TE for T>25 vs T=25 (regarding the direct effect of T=25, this effect is equal to the effect of HW in this simulation setting)
#o3.ref <- mean(log(sim_city[sim_city$HW==0,"O3"])) # or median(log(sim_city[sim_city$HW==0,"O3"])) ## if you use this, you estimate TE for HW=0 vs HW=1



## Prediction uncertainty of logO3
O3pred.sigma <- sim_city$log_o3_sigma

## Observed simulated daily deaths.
## For a real-data application, this would be the observed mortality count.
death <- sim_city$death

## For a smoother purely simulated example, one could use true_mu instead:
## death <- sim_city$true_mu
############################################################
## 3. Create blank output objects
############################################################

n_sim <- nrow(BTcoef)

results_total <- data.frame(
  sim = 1:n_sim,
  beta_HW = NA,
  beta_O3 = NA,
  beta_INT = NA,
  death_total = NA,
  AN_CDE = NA,
  AN_INTref = NA,
  AN_INTmed = NA,
  AN_PIE = NA,
  AN_TE = NA,
  mean_RR_CDE = NA,
  mean_RR_INT = NA,
  mean_RR_PIE = NA,
  mean_RR_TE = NA
)

results_year <- data.frame(
  sim = rep(1:n_sim, each = 10),
  year = rep(2011:2020, times = n_sim),
  death_total = NA,
  AN_CDE = NA,
  AN_INT = NA,
  AN_PIE = NA,
  AN_TE = NA
)


############################################################
## 4. Loop over coefficient draws and calculate components
############################################################

for (nn in 1:n_sim) {
  
  ##########################################################
  ## 4.1. Pick one coefficient draw
  ##########################################################
  
  HW.e  <- BTcoef$beta_HW[nn]
  O3.e  <- BTcoef$beta_O3[nn]
  Int.e <- BTcoef$beta_INT[nn]
  
  
  ##########################################################
  ## 4.2. Controlled direct effect
  ##########################################################
  
  RR_CDE <- exp((HW.e + Int.e * o3.ref) * HWval)
  
  
  ##########################################################
  ## 4.3. Kappa term
  ##########################################################
  
  kappa <- exp(
    O3.e * o3.ref +
      Int.e * hw.ref * o3.ref -
      (O3.e + Int.e * hw.ref) * logO3_pred_NHW -
      0.5 * (O3.e + Int.e * hw.ref)^2 * O3pred.sigma^2
  )
  
  
  ##########################################################
  ## 4.4. Reference interaction
  ##########################################################
  
  RR_INTref <-
    exp(
      HW.e * HWval -
        O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * HWval) * logO3_pred_NHW +
        0.5 * (O3.e + Int.e * HWval)^2 * O3pred.sigma^2
    ) -
    exp(
      -O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * hw.ref) * logO3_pred_NHW +
        0.5 * (O3.e + Int.e * hw.ref)^2 * O3pred.sigma^2
    ) -
    exp(
      (HW.e + Int.e * o3.ref) * HWval
    ) +
    1
  
  RR_INTref<-ifelse(HWval==1,RR_INTref,0) ## Conceptually, no interaction exists when HW==0
  
  ##########################################################
  ## 4.5. Mediated interaction
  ##########################################################
  
  RR_INTmed <-
    exp(
      HW.e * HWval -
        O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * HWval) * logO3_pred_HW +
        0.5 * (O3.e + Int.e * HWval)^2 * O3pred.sigma^2
    ) -
    exp(
      -O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * hw.ref) * logO3_pred_HW +
        0.5 * (O3.e + Int.e * hw.ref)^2 * O3pred.sigma^2
    ) -
    exp(
      HW.e * HWval -
        O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * HWval) * logO3_pred_NHW +
        0.5 * (O3.e + Int.e * HWval)^2 * O3pred.sigma^2
    ) +
    exp(
      -O3.e * o3.ref -
        Int.e * hw.ref * o3.ref +
        (O3.e + Int.e * hw.ref) * logO3_pred_NHW +
        0.5 * (O3.e + Int.e * hw.ref)^2 * O3pred.sigma^2
    )
  RR_INTmed<-ifelse(HWval==1,RR_INTmed,0) ## Conceptually, no interaction exists when HW==0
  
  
  ##########################################################
  ## 4.6. Pure indirect effect
  ##########################################################
  
  ## Supplement equation:
  ## RR_PIE does not multiply by HWval.
  ## The mediator contrast is:
  ##   logO3_pred_HW - logO3_pred_NHW
  
  RR_PIE <- exp(
    (O3.e + Int.e * hw.ref) *
      (logO3_pred_HW - logO3_pred_NHW)
  )
  
  
  ##########################################################
  ## 4.7. Total effect
  ##########################################################
  
  RR_TE <- kappa * (
    (RR_CDE - 1) +
      RR_INTref +
      RR_INTmed
  ) + RR_PIE
  
  
  ##########################################################
  ## 4.8. Attributable number decomposition
  ##########################################################
  
  ## Total attributable number:
  ##   AN_TE = death * (RR_TE - 1) / RR_TE
  ##
  ## Component-specific attributable numbers:
  ##   These add up to AN_TE because:
  ##
  ##   RR_TE - 1 =
  ##     kappa * (RR_CDE - 1)
  ##     + kappa * RR_INTref
  ##     + kappa * RR_INTmed
  ##     + (RR_PIE - 1)
  
  AF_CDE_daily <- (RR_CDE - 1) / RR_CDE
  AF_PIE_daily <- (RR_PIE - 1) / RR_PIE
  AF_TE_daily <- (RR_TE - 1) / RR_TE
  AF_INF_daily <- AF_TE_daily-AF_PIE_daily-AF_CDE_daily

  AN_CDE_daily <- death * AF_CDE_daily
  AN_PIE_daily <- death * AF_PIE_daily
  AN_TE_daily <- death * AF_TE_daily
  AN_INT_daily <- death * AF_INF_daily
  
  
  ##########################################################
  ## 4.9. Save total 2011-2020 results
  ##########################################################
  
  results_total$beta_HW[nn]  <- HW.e
  results_total$beta_O3[nn]  <- O3.e
  results_total$beta_INT[nn] <- Int.e
  
  results_total$death_total[nn] <- sum(death, na.rm = TRUE)
  
  results_total$AN_CDE[nn]    <- sum(AN_CDE_daily, na.rm = TRUE)
  results_total$AN_INT[nn] <- sum(AN_INT_daily, na.rm = TRUE)
  results_total$AN_PIE[nn]    <- sum(AN_PIE_daily, na.rm = TRUE)
  results_total$AN_TE[nn]     <- sum(AN_TE_daily, na.rm = TRUE)
  
  ##########################################################
  ## 4.10. Save year-specific results
  ##########################################################
  
  row_start <- (nn - 1) * 10 + 1
  row_end   <- row_start + 9
  
  temp_year_result <- data.frame(
    sim = nn,
    year = 2011:2020,
    death_total = NA,
    AN_CDE = NA,
    AN_INT = NA,
    AN_PIE = NA,
    AN_TE = NA
  )
  
  for (yy in 2011:2020) {
    
    pick <- which(sim_city$year == yy)
    row_pick <- which(temp_year_result$year == yy)
    
    temp_year_result$death_total[row_pick] <- sum(death[pick], na.rm = TRUE)
    
    temp_year_result$AN_CDE[row_pick] <-
      sum(AN_CDE_daily[pick], na.rm = TRUE)
    
    temp_year_result$AN_INT[row_pick] <-
      sum(AN_INT_daily[pick], na.rm = TRUE)
    
    temp_year_result$AN_PIE[row_pick] <-
      sum(AN_PIE_daily[pick], na.rm = TRUE)
    
    temp_year_result$AN_TE[row_pick] <-
      sum(AN_TE_daily[pick], na.rm = TRUE)
  }
  
  results_year[row_start:row_end, ] <- temp_year_result
  
  
  ##########################################################
  ## 4.11. Print progress occasionally
  ##########################################################
  
  if (nn %% 100 == 0) {
    cat("Completed simulation draw:", nn, "of", n_sim, "\n")
  }
}


############################################################
## 5. Summarize uncertainty intervals
############################################################

summary_total <- data.frame(
  component = c("CDE", "INT", "PIE", "TE"),
  mean = NA,
  lower_2.5 = NA,
  upper_97.5 = NA
)

summary_total$mean[summary_total$component == "CDE"] <-
  mean(results_total$AN_CDE, na.rm = TRUE)
summary_total$lower_2.5[summary_total$component == "CDE"] <-
  quantile(results_total$AN_CDE, probs = 0.025, na.rm = TRUE)
summary_total$upper_97.5[summary_total$component == "CDE"] <-
  quantile(results_total$AN_CDE, probs = 0.975, na.rm = TRUE)

summary_total$mean[summary_total$component == "INT"] <-
  mean(results_total$AN_INT, na.rm = TRUE)
summary_total$lower_2.5[summary_total$component == "INT"] <-
  quantile(results_total$AN_INT, probs = 0.025, na.rm = TRUE)
summary_total$upper_97.5[summary_total$component == "INT"] <-
  quantile(results_total$AN_INT, probs = 0.975, na.rm = TRUE)

summary_total$mean[summary_total$component == "PIE"] <-
  mean(results_total$AN_PIE, na.rm = TRUE)
summary_total$lower_2.5[summary_total$component == "PIE"] <-
  quantile(results_total$AN_PIE, probs = 0.025, na.rm = TRUE)
summary_total$upper_97.5[summary_total$component == "PIE"] <-
  quantile(results_total$AN_PIE, probs = 0.975, na.rm = TRUE)

summary_total$mean[summary_total$component == "TE"] <-
  mean(results_total$AN_TE, na.rm = TRUE)
summary_total$lower_2.5[summary_total$component == "TE"] <-
  quantile(results_total$AN_TE, probs = 0.025, na.rm = TRUE)
summary_total$upper_97.5[summary_total$component == "TE"] <-
  quantile(results_total$AN_TE, probs = 0.975, na.rm = TRUE)

############################################################
## 6. Save outputs
############################################################

dir.create(
  paste0(FOLDER_NAME,"/Output"),
  showWarnings = FALSE
)

write.csv(
  results_total,
  file = paste0(FOLDER_NAME,"/Output/TE_components_total_draws_city.csv"),
  row.names = FALSE
)

write.csv(
  results_year,
  file = paste0(FOLDER_NAME,"/Output/TE_components_year_draws_city.csv"),
  row.names = FALSE
)

write.csv(
  summary_total,
  file = paste0(FOLDER_NAME,"/Output/TE_components_total_summary_city.csv"),
  row.names = FALSE
)


############################################################
## 7. Print checks
############################################################

cat("\nSaved files:\n")
cat("Output/TE_components_total_draws_city.csv\n")
cat("Output/TE_components_year_draws_city.csv\n")
cat("Output/TE_components_total_summary_city.csv\n\n")

cat("Total attributable deaths summary:\n")
print(summary_total)

cat("\nCheck additivity using posterior/draw means:\n")

mean_CDE    <- mean(results_total$AN_CDE, na.rm = TRUE)
mean_INT <- mean(results_total$AN_INT, na.rm = TRUE)
mean_PIE    <- mean(results_total$AN_PIE, na.rm = TRUE)
mean_TE     <- mean(results_total$AN_TE, na.rm = TRUE)

cat("CDE + INT  + PIE =",
    mean_CDE + mean_INT + mean_PIE, "\n")
cat("TE =", mean_TE, "\n")
cat("Difference =",
    mean_TE - (mean_CDE + mean_INT + mean_PIE), "\n")
