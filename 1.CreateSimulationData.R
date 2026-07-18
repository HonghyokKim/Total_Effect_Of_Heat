############################################################
## 00_simulate_city_data.R
## Simulated daily city summer data for TE decomposition example
## No personal exposure variables are generated.
## Data structure:
##   temperature -> heat wave -> heat-mediated ozone -> mortality
############################################################
FOLDER_NAME<-"YOUR PATH"
set.seed(12423)

############################################################
## 1. Create date sequence: city, 10 summers, 2011-2020
############################################################

dates <- seq(as.Date("2011-06-01"), as.Date("2020-08-31"), by = "day")

dat <- data.frame(date = dates)

dat$year  <- as.integer(format(dat$date, "%Y"))
dat$month <- as.integer(format(dat$date, "%m"))
dat$day   <- as.integer(format(dat$date, "%d"))

## Keep summer months only: June-August
dat <- dat[dat$month %in% c(6, 7, 8), ]

dat$city <- "city"

## Day of summer: June 1 = 1
dat$day_summer <- as.integer(dat$date - as.Date(paste0(dat$year, "-06-01"))) + 1

## Day of week: 0 = Sunday, 1 = Monday, ..., 6 = Saturday
dat$dow <- as.POSIXlt(dat$date)$wday
dat$dow_name <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")[dat$dow + 1]


############################################################
## 2. Generate synthetic city summer temperature
############################################################

## Hard-coded annual anomalies.
## These are NOT real city values.
## They just create year-to-year variability.
dat$year_anomaly <- NA

dat$year_anomaly[dat$year == 2011] <- -0.6
dat$year_anomaly[dat$year == 2012] <-  0.2
dat$year_anomaly[dat$year == 2013] <-  0.4
dat$year_anomaly[dat$year == 2014] <- -0.3
dat$year_anomaly[dat$year == 2015] <-  0.1
dat$year_anomaly[dat$year == 2016] <-  0.6
dat$year_anomaly[dat$year == 2017] <-  0.3
dat$year_anomaly[dat$year == 2018] <-  1.0
dat$year_anomaly[dat$year == 2019] <-  0.5
dat$year_anomaly[dat$year == 2020] <-  0.2

## Slow warming trend across the 10-year synthetic period
dat$warming_trend <- 0.04 * (dat$year - 2011)

## Seasonal mean temperature.
## This peaks around late July / early August.
dat$temp_seasonal_mean <- 25 +
  3.5 * sin(2 * pi * (dat$day_summer - 37) / 122)

## Autocorrelated daily weather noise
dat$temp_weather_noise <- NA

for (i in 1:nrow(dat)) {
  
  if (i == 1) {
    dat$temp_weather_noise[i] <- rnorm(1, mean = 0, sd = 2.0)
  }
  
  if (i > 1) {
    
    ## If the previous row is the previous calendar day, keep autocorrelation.
    ## If not, this is the beginning of a new summer.
    if (as.numeric(dat$date[i] - dat$date[i - 1]) == 1) {
      dat$temp_weather_noise[i] <- 0.75 * dat$temp_weather_noise[i - 1] +
        rnorm(1, mean = 0, sd = 1.6)
    }
    
    if (as.numeric(dat$date[i] - dat$date[i - 1]) != 1) {
      dat$temp_weather_noise[i] <- rnorm(1, mean = 0, sd = 2.0)
    }
  }
}

dat$temp <- dat$temp_seasonal_mean +
  dat$year_anomaly +
  dat$warming_trend +
  dat$temp_weather_noise


############################################################
## 3. Define city-specific heat wave threshold and HW variable
############################################################

## city-specific 95th percentile threshold from the simulated data
hw_threshold <- as.numeric(quantile(dat$temp, probs = 0.95, na.rm = TRUE))

dat$hw_threshold <- hw_threshold

## Daily exceedance
dat$temp_above_hw_threshold <- 0
dat$temp_above_hw_threshold[dat$temp > hw_threshold] <- 1

## Heat wave: at least two consecutive days above threshold
dat$HW <- 0

for (i in 2:nrow(dat)) {
  
  if (as.numeric(dat$date[i] - dat$date[i - 1]) == 1) {
    
    if (dat$temp_above_hw_threshold[i] == 1 &
        dat$temp_above_hw_threshold[i - 1] == 1) {
      
      dat$HW[i] <- 1
      dat$HW[i - 1] <- 1
    }
  }
}


############################################################
## 4. Generate synthetic ozone
############################################################

## Autocorrelated ozone noise
dat$o3_noise <- NA

for (i in 1:nrow(dat)) {
  
  if (i == 1) {
    dat$o3_noise[i] <- rnorm(1, mean = 0, sd = 3.0)
  }
  
  if (i > 1) {
    
    if (as.numeric(dat$date[i] - dat$date[i - 1]) == 1) {
      dat$o3_noise[i] <- 0.60 * dat$o3_noise[i - 1] +
        rnorm(1, mean = 0, sd = 2.8)
    }
    
    if (as.numeric(dat$date[i] - dat$date[i - 1]) != 1) {
      dat$o3_noise[i] <- rnorm(1, mean = 0, sd = 3.0)
    }
  }
}

## Month effects for ozone.
## These are synthetic and only used to create plausible seasonality.
dat$o3_month_effect <- 0
dat$o3_month_effect[dat$month == 6] <- 2
dat$o3_month_effect[dat$month == 7] <- 4
dat$o3_month_effect[dat$month == 8] <- 1


############################################################
## 4.1. Baseline ozone
############################################################

## Baseline ozone before temperature/HW-related increments.
## Unit here is ppm, but this is simulated.
dat$O3_base <- 31 +
  dat$o3_month_effect +
  dat$o3_noise

dat$O3_base <- pmax(dat$O3_base, 10)
dat$O3_base <-dat$O3_base *0.001 #ppm

############################################################
## 4.2. General high-temperature ozone increment
############################################################

## This component applies when temperature is >= 25 C,
## regardless of heat-wave status.
##
## Therefore, this component exists for:
##   HW == 0 and temp >= 25
##   HW == 1 and temp >= 25
##
## This is NOT the heat-wave-mediated ozone contrast by itself.
dat$o3_temp25_increment_ppb <- 0

dat$o3_temp25_increment_ppb[dat$HW == 0 & dat$temp >= 25] <-
 1 * (dat$temp[dat$HW == 0 & dat$temp >= 25] - 25)

dat$o3_temp25_increment_ppb <-dat$o3_temp25_increment_ppb *0.001 #ppm
############################################################
## 4.3. HW-only ozone increment
############################################################

## This component applies whenever HW == 1.
## It represents an additional ozone increase associated with
## the binary heat-wave condition.
dat$o3_HW_increment_ppb <- 0

dat$o3_HW_increment_ppb[dat$HW == 1] <- 10
dat$o3_HW_increment_ppb <-dat$o3_HW_increment_ppb *0.001 #ppm

############################################################
## 4.4. Counterfactual and observed ozone
############################################################

## Counterfactual ozone under no heat-wave-mediated increment.
##
## This still includes the general temp >= 25 C ozone increase.
## Therefore, HW == 0 and temp >= 25 days can have elevated O3_nohw.
dat$O3_nohw_above25 <- dat$O3_base +
  dat$o3_temp25_increment_ppb

## Observed ozone under actual heat-wave conditions.
##
## This includes:
##   1) baseline ozone
##   2) general temp >= 25 C increment
##   3) HW-only increment
dat$O3 <- dat$O3_base +
  dat$o3_temp25_increment_ppb +
  dat$o3_HW_increment_ppb


############################################################
## 4.5. Heat-wave-mediated ozone contrast
############################################################
dat$HW_o3_increment_ppb <- ifelse(dat$HW==1,dat$O3-dat$O3_base,0)
dat$HW_o3_increment_ppb <-dat$HW_o3_increment_ppb *0.001 #ppm

############################################################
## 4.6. Log ozone variables
############################################################

## Log ozone variables used for the TE component calculations
dat$logO3_nohw_above25<- log(dat$O3_nohw_above25)
dat$logO3 <- log(dat$O3)
dat$logO3_base <- log(dat$O3_base)

############################################################
## 4.7. Reference ozone level
############################################################

## Reference ozone level.
## No heat wave days + temperature <25¡ÆC
log_o3_ref_value <- median(log(dat$O3_base)[dat$HW == 0], na.rm = TRUE)
dat$log_o3_ref <- log_o3_ref_value


############################################################
## 4.8. Prediction uncertainty for log ozone
############################################################
## Constant value for the example.
dat$log_o3_sigma <- 1


############################################################
## 5. Generate synthetic city population
############################################################

## Synthetic city population.
## Roughly declining over time, but not intended to be real.
dat$population <- NA

dat$population[dat$year == 2011] <- 10400000
dat$population[dat$year == 2012] <- 10350000
dat$population[dat$year == 2013] <- 10300000
dat$population[dat$year == 2014] <- 10250000
dat$population[dat$year == 2015] <- 10200000
dat$population[dat$year == 2016] <- 10150000
dat$population[dat$year == 2017] <- 10100000
dat$population[dat$year == 2018] <- 10050000
dat$population[dat$year == 2019] <- 10000000
dat$population[dat$year == 2020] <-  9950000


############################################################
## 6. Generate synthetic daily mortality counts
############################################################

## Synthetic daily baseline mortality rate:
## approximately 13 deaths per 1,000,000 persons per day.
baseline_daily_mortality_rate <- 13 / 1000000

## True coefficients used only to generate the simulated outcome.
## These are not intended to be real epidemiologic estimates.
beta_HW_true  <- 0.2
beta_O3_true  <- 0.005
beta_INT_true <- 0.01

## Month effects for mortality
dat$mort_month_effect <- 0
dat$mort_month_effect[dat$month == 6] <- 0.00
dat$mort_month_effect[dat$month == 7] <- 0.01
dat$mort_month_effect[dat$month == 8] <- 0.02

## Day-of-week effects for mortality
dat$mort_dow_effect <- 0
dat$mort_dow_effect[dat$dow == 1] <-  0.01   ## Monday
dat$mort_dow_effect[dat$dow == 2] <-  0.00   ## Tuesday
dat$mort_dow_effect[dat$dow == 3] <-  0.00   ## Wednesday
dat$mort_dow_effect[dat$dow == 4] <-  0.00   ## Thursday
dat$mort_dow_effect[dat$dow == 5] <- -0.01   ## Friday
dat$mort_dow_effect[dat$dow == 6] <- -0.02   ## Saturday
dat$mort_dow_effect[dat$dow == 0] <- -0.02   ## Sunday

## Intercept adjusted so that the baseline rate is interpretable
alpha0 <- log(baseline_daily_mortality_rate) -
  beta_O3_true * log_o3_ref_value

## Log expected daily deaths
dat$log_mu <- log(dat$population) +
  alpha0 +
  dat$mort_month_effect +
  dat$mort_dow_effect +
  beta_HW_true * dat$HW +
  beta_O3_true * dat$logO3 +
  beta_INT_true * dat$HW * dat$logO3

dat$true_mu <- exp(dat$log_mu)

## Simulated mortality count
dat$death <- rpois(nrow(dat), lambda = dat$true_mu)


############################################################
## 7. Keep public-use columns and save
############################################################

sim_city <- dat[, c(
  "city",
  "date",
  "year",
  "month",
  "day",
  "dow",
  "dow_name",
  "day_summer",
  "population",
  "temp",
  "hw_threshold",
  "temp_above_hw_threshold",
  "HW",
  "O3_nohw_above25",
  "O3",
  "logO3_nohw_above25",
  "logO3",
  "logO3_base",
  "HW_o3_increment_ppb",
  "log_o3_ref",
  "log_o3_sigma",
  "true_mu",
  "death"
)]

dir.create("data", showWarnings = FALSE)

write.csv(
  sim_city,
  file = paste0(FOLDER_NAME,"/Data/sim_city_daily.csv"),
  row.names = FALSE
)


############################################################
## 8. Simple checks
############################################################

cat("Number of rows:", nrow(sim_city), "\n")
cat("Date range:", as.character(min(sim_city$date)), "to", as.character(max(sim_city$date)), "\n")
cat("HW threshold:", round(hw_threshold, 2), "C\n")
cat("Number of HW days:", sum(sim_city$HW), "\n")
cat("Mean temperature:", round(mean(sim_city$temp), 2), "C\n")
cat("Mean O3:", round(mean(sim_city$O3), 2), "ppb\n")
cat("Mean daily deaths:", round(mean(sim_city$death), 2), "\n")

print(summary(sim_city[, c("temp", "O3_nohw_above25", "O3", "logO3_nohw_above25", "logO3", "death")]))
print(table(sim_city$year, sim_city$HW))