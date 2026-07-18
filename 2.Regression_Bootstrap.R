############################################################
## 01_fit_model_make_BTcoef.R
## Fit the outcome model and create coefficient draws
## for the TE decomposition example.
##
## This script uses the simulated city-level daily data
## created by 00_simulate_city_data.R.
##
## No personal exposure variables are used.
############################################################
FOLDER_NAME<-"YOUR PATH"
set.seed(20300110)


############################################################
## 1. Read simulated data
############################################################

sim_city <- read.csv(
  paste0(FOLDER_NAME,"/Data/sim_city_daily.csv")
)

sim_city$date <- as.Date(sim_city$date)

sim_city$month <- as.factor(sim_city$month)
sim_city$dow   <- as.factor(sim_city$dow)

cat("Number of rows:", nrow(sim_city), "\n")
cat("Date range:", as.character(min(sim_city$date)), "to", as.character(max(sim_city$date)), "\n")
cat("Number of HW days:", sum(sim_city$HW), "\n")
cat("Mean daily deaths:", round(mean(sim_city$death), 2), "\n")


############################################################
## 2. Fit outcome model
############################################################

## Outcome model:
##
##   log E(Y) =
##     log(population)
##     + beta_HW * HW
##     + beta_O3 * logO3
##     + beta_INT * HW * logO3
##     + month indicators
##     + day-of-week indicators
##
## This is the simplified binary-HW + logO3 model used for
## the hands-on TE decomposition example.

fit_city <- glm(
  death ~ HW + logO3 + HW:logO3 + month + dow,
  family = quasipoisson(link = "log"),
  offset = log(population),
  data = sim_city
)

print(summary(fit_city))


############################################################
## 3. Extract coefficients used in the TE decomposition
############################################################

coef_all <- coef(fit_city)
vcov_all <- vcov(fit_city)

beta_HW_hat  <- coef_all["HW"]
beta_O3_hat  <- coef_all["logO3"]
beta_INT_hat <- coef_all["HW:logO3"]

cat("\nEstimated coefficients used for decomposition:\n")
cat("beta_HW :", beta_HW_hat, "\n")
cat("beta_O3 :", beta_O3_hat, "\n")
cat("beta_INT:", beta_INT_hat, "\n")

beta_hat <- c(beta_HW_hat, beta_O3_hat, beta_INT_hat)
names(beta_hat) <- c("beta_HW", "beta_O3", "beta_INT")

vcov_beta <- vcov_all[
  c("HW", "logO3", "HW:logO3"),
  c("HW", "logO3", "HW:logO3")
]

print(vcov_beta)


############################################################
## 4. Create coefficient draws
############################################################

## The original empirical analysis used non-parametric bootstrap.
## For this simulated example, we generate
## multivariate normal coefficient draws based on the fitted
## model estimate and covariance matrix for simplicity. This is parametric. 
##
## BTcoef[,1] = beta_HW
## BTcoef[,2] = beta_O3
## BTcoef[,3] = beta_INT

n_sim <- 1000

z_mat <- matrix(rnorm(n_sim * 3), nrow = n_sim, ncol = 3)

## Small ridge term for numerical stability
vcov_beta_stable <- vcov_beta + diag(0.000000000001, 3)

chol_vcov <- chol(vcov_beta_stable)

BTcoef <- z_mat %*% chol_vcov

BTcoef[, 1] <- BTcoef[, 1] + beta_hat[1]
BTcoef[, 2] <- BTcoef[, 2] + beta_hat[2]
BTcoef[, 3] <- BTcoef[, 3] + beta_hat[3]

BTcoef <- as.data.frame(BTcoef)

colnames(BTcoef) <- c("beta_HW", "beta_O3", "beta_INT")

cat("\nSummary of coefficient draws:\n")
print(summary(BTcoef))


############################################################
## 5. Save outputs
############################################################

dir.create(
  "FOLDER_NAME/Output",
  showWarnings = FALSE
)

write.csv(
  BTcoef,
  file = paste0(FOLDER_NAME,"/Output/BTcoef_city.csv"),
  row.names = FALSE
)

coef_summary <- data.frame(
  term = c("beta_HW", "beta_O3", "beta_INT"),
  estimate = c(beta_HW_hat, beta_O3_hat, beta_INT_hat),
  se = c(
    sqrt(vcov_beta[1, 1]),
    sqrt(vcov_beta[2, 2]),
    sqrt(vcov_beta[3, 3])
  )
)

write.csv(
  coef_summary,
  file = paste0(FOLDER_NAME,"/Output/model_coef_summary_city.csv"),
  row.names = FALSE
)

save(
  fit_city,
  beta_hat,
  vcov_beta,
  BTcoef,
  file = paste0(FOLDER_NAME,"/Output/fitted_model_city.RData")
)


############################################################
## 6. Simple checks
############################################################

cat("\nSaved files:\n")
cat("Output/BTcoef_city.csv\n")
cat("Output/model_coef_summary_city.csv\n")
cat("Output/fitted_model_city.RData\n")

cat("\nMean of BTcoef draws:\n")
print(colMeans(BTcoef))

cat("\nOriginal coefficient estimates:\n")
print(beta_hat)