# ============================================================
#Project 2 
#Nisa Ekmekcioglu 12225999
# ============================================================
rm(list = ls())
setwd("/Users/nisaekmekcioglu/Desktop/Zeitreihenanalyse")

# ============================
#Figure 2.1: (a) Periodogram of quarterly log INDPRO
#(b) Cumulative periodogram of quarterly log INDPRO
#(c) Periodogram of first differences of quarterly log INDPRO
#(d) Cumulative periodogram of first differences
# ============================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Converting date column to Date class
df$observation_date <- as.Date(df$observation_date)

#Converting series to numeric (important: spectrum needs numeric)
df$INDPRO <- as.numeric(df$INDPRO)

#Dropping missing / invalid values
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample to avoid extreme historical shocks (WWII + Covid)
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]


#Building the series we need

#(Level) log series
y <- log(df$INDPRO)

#(Growth) first differences of logs: Δlog(INDPRO)
#diff() shortens the series by 1 observation
dy <- diff(y)


#Computing periodograms
spec_y  <- spectrum(y,  plot = FALSE, demean = TRUE)
spec_dy <- spectrum(dy, plot = FALSE, demean = TRUE)

#Converting frequencies to radians (0..pi) like in many textbooks:
#omega = 2*pi*f
omega_y  <- spec_y$freq  * 2*pi
omega_dy <- spec_dy$freq * 2*pi

#Cumulative periodogram = cumulative sum of spectral mass / total mass
cum_y  <- cumsum(spec_y$spec)  / sum(spec_y$spec)
cum_dy <- cumsum(spec_dy$spec) / sum(spec_dy$spec)


#Ploting
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))     

#(a) Periodogram of log(INDPRO)
plot(omega_y, spec_y$spec,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Periodogram",
     main = "(a)")

#(b) Cumulative periodogram of log(INDPRO)
plot(omega_y, cum_y,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Cumulative periodogram",
     main = "(b)")

#(c) Periodogram of first differences Δlog(INDPRO)
plot(omega_dy, spec_dy$spec,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Periodogram",
     main = "(c)")

#(d) Cumulative periodogram of first differences Δlog(INDPRO)
plot(omega_dy, cum_dy,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Cumulative periodogram",
     main = "(d)")


# ============================
#Table 2.1: Testing the differenced log INDPRO for white noise
# ============================

#Preparing the data (same as before)
df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series and first differences
y  <- log(df$INDPRO)
dy <- diff(y)

#Computing periodogram
spec_dy <- spectrum(dy, plot = FALSE, demean = TRUE)

I_w   <- spec_dy$spec          #periodogram values
freq  <- spec_dy$freq * 2*pi   #frequencies in radians

#Sample variance of dy
s2 <- var(dy)


j_values <- c(10, 20, 30, 40, 50)

#Creating storage vectors
freq_j  <- numeric(length(j_values))
T_j     <- numeric(length(j_values))
crit95  <- numeric(length(j_values))
p_val   <- numeric(length(j_values))


#Computing test statistic
for (i in seq_along(j_values)) {
  
  j <- j_values[i]
  
  #j-th frequency
  freq_j[i] <- freq[j]
  
  #Test statistic:
  #T_j = (4*pi / s^2) * sum of first j periodogram ordinates
  T_j[i] <- (4*pi / s2) * sum(I_w[1:j])
  
  #0.95 quantile of chi^2(2j)
  crit95[i] <- qchisq(0.95, df = 2*j)
  
  #p-value
  p_val[i] <- 1 - pchisq(T_j[i], df = 2*j)
}


#Creating the table
result_table <- data.frame(
  j = j_values,
  frequency = round(freq_j, 6),
  T_statistic = round(T_j, 5),
  chi2_0.95 = round(crit95, 5),
  p_value = round(p_val, 6)
)

print(result_table)


# ============================
#Figure 2.2: (a) Periodogram of daily log returns of AAPL
#(b) Cumulative periodogram of daily log returns of AAPL
# ============================

#Reading APPL data
library(quantmod)

getSymbols("AAPL", src = "yahoo", from = "1980-01-01")

#Creating a data frame with Adjusted Close
aapl <- data.frame(
  Date = index(AAPL),
  Close = as.numeric(Ad(AAPL))   #adjusted close
)

write.csv(aapl, "AAPL_full.csv", row.names = FALSE)

aapl <- read.csv("AAPL_full.csv")

names(aapl)
head(aapl)

nrow(aapl)

aapl$Date  <- as.Date(aapl$Date)
aapl$Close <- as.numeric(aapl$Close)

aapl <- aapl[!is.na(aapl$Date) & is.finite(aapl$Close), ]

#Sorting by date 
aapl <- aapl[order(aapl$Date), ]

#Computing daily log returns
log_price <- log(aapl$Close)

#r_t = log(P_t) - log(P_{t-1})
returns <- diff(log_price)

#Computing periodogram
spec_ret <- spectrum(returns, plot = FALSE, demean = TRUE)

omega <- spec_ret$freq * 2*pi
I_w   <- spec_ret$spec

#Cumulative periodogram
cum_ret <- cumsum(I_w) / sum(I_w)


#Ploting
par(mfrow = c(2,1), mar = c(4,4,2,1))

#(a) Periodogram
plot(omega, I_w,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Periodogram",
     main = "(a)")

#(b) Cumulative periodogram
plot(omega, cum_ret,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Cumulative periodogram",
     main = "(b)")


# ============================
#Table 2.2: Testing the daily log returns of AAPL for white noise
# ============================

#Preparing the data (same as before)
aapl <- read.csv("AAPL_full.csv", stringsAsFactors = FALSE)

aapl$Date  <- as.Date(aapl$Date)
aapl$Close <- as.numeric(aapl$Close)

aapl <- aapl[!is.na(aapl$Date) & is.finite(aapl$Close), ]

aapl <- aapl[order(aapl$Date), ]


#Computing daily log returns
log_price <- log(aapl$Close)
returns   <- diff(log_price)


#Computing the periodogram
spec_ret <- spectrum(returns, plot = FALSE, demean = TRUE)

I_w  <- spec_ret$spec
freq <- spec_ret$freq * 2*pi

#Sample variance
s2 <- var(returns)

j_values <- c(1000, 2000, 3000, 4000, 5000)

#Storage
freq_j  <- numeric(length(j_values))
T_j     <- numeric(length(j_values))
crit95  <- numeric(length(j_values))
p_val   <- numeric(length(j_values))


#Computing test statistics
for (i in seq_along(j_values)) {
  
  j <- j_values[i]
  
  freq_j[i] <- freq[j]
  
  T_j[i] <- (4*pi / s2) * sum(I_w[1:j])
  
  crit95[i] <- qchisq(0.95, df = 2*j)
  
  p_val[i] <- 1 - pchisq(T_j[i], df = 2*j)
}


#Creating the results table
result_table <- data.frame(
  j = j_values,
  frequency = round(freq_j, 6),
  T_statistic = round(T_j, 3),
  chi2_0.95 = round(crit95, 3),
  p_value = round(p_val, 5)
)

print(result_table)


# ============================
#Figure 2.3: Periodogram of monthly log returns of AAPl
# ============================

#Preparing the data (same as before)
aapl <- read.csv("AAPL_full.csv", stringsAsFactors = FALSE)

aapl$Date  <- as.Date(aapl$Date)
aapl$Close <- as.numeric(aapl$Close)

aapl <- aapl[!is.na(aapl$Date) & is.finite(aapl$Close), ]
aapl <- aapl[order(aapl$Date), ]


#Converting to monthly prices

#Taking the last trading day of each month
aapl$YearMonth <- format(aapl$Date, "%Y-%m")

monthly_data <- aapl[!duplicated(aapl$YearMonth, fromLast = TRUE), ]

monthly_data <- monthly_data[order(monthly_data$Date), ]


#Computing monthly log returns
log_price_m <- log(monthly_data$Close)
returns_m   <- diff(log_price_m)


#Computing the periodogram
spec_m <- spectrum(returns_m, plot = FALSE, demean = TRUE)

omega <- spec_m$freq * 2*pi
I_w   <- spec_m$spec


#Seasonal frequencies (monthly)
seasonal_freq <- 2*pi*(1:6)/12


#Ploting
par(mfrow = c(1,1))

plot(omega, I_w,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Periodogram",
     main = "")

#red circles: seasonal frequencies
points(seasonal_freq,
       approx(omega, I_w, seasonal_freq)$y,
       col = "red",
       pch = 16)


# ============================
#Table 2.3: Testing the monthly log returns of AAPL for seasonal patterns
# ============================

#Preparing the data (same as before)
aapl <- read.csv("AAPL_full.csv", stringsAsFactors = FALSE)

aapl$Date  <- as.Date(aapl$Date)
aapl$Close <- as.numeric(aapl$Close)

aapl <- aapl[!is.na(aapl$Date) & is.finite(aapl$Close), ]
aapl <- aapl[order(aapl$Date), ]

#last trading day each month
aapl$YearMonth <- format(aapl$Date, "%Y-%m")
monthly_data <- aapl[!duplicated(aapl$YearMonth, fromLast = TRUE), ]
monthly_data <- monthly_data[order(monthly_data$Date), ]

#Monthly log returns
log_price_m <- log(monthly_data$Close)
returns_m   <- diff(log_price_m)


#Periodogram
spec_m <- spectrum(returns_m, plot = FALSE, demean = TRUE)

omega <- spec_m$freq * 2*pi
I_w   <- spec_m$spec

s2 <- var(returns_m)


#Seasonal frequencies
seasonal_freq <- 2*pi*(1:6)/12

#Finding the closest frequencies in spectrum
idx <- sapply(seasonal_freq,
              function(w) which.min(abs(omega - w)))


j_values <- c(1, 3, 6)

T_j    <- numeric(length(j_values))
df_j   <- numeric(length(j_values))
p_val  <- numeric(length(j_values))

for (i in seq_along(j_values)) {
  
  j <- j_values[i]
  
  #Sum of first j seasonal frequencies
  selected_idx <- idx[1:j]
  
  T_j[i] <- (4*pi / s2) * sum(I_w[selected_idx])
  
  #Degrees of freedom:
  #each seasonal frequency gives 2 df
  #but (k=6) contributes only 1 df
  if (j < 6) {
    df_j[i] <- 2*j
  } else {
    df_j[i] <- 2*j - 1
  }
  
  p_val[i] <- 1 - pchisq(T_j[i], df = df_j[i])
}

#Result table
result_table <- data.frame(
  j = j_values,
  T_statistic = round(T_j, 6),
  degrees_of_freedom = df_j,
  p_value = round(p_val, 5)
)

print(result_table)


# ============================
#Figure 2.3: Seasonal periodogram of monthly log returns of AAPl
# ============================

#Preparing the data (same as before)
aapl <- read.csv("AAPL_full.csv", stringsAsFactors = FALSE)

aapl$Date  <- as.Date(aapl$Date)
aapl$Close <- as.numeric(aapl$Close)

aapl <- aapl[!is.na(aapl$Date) & is.finite(aapl$Close), ]
aapl <- aapl[order(aapl$Date), ]


#Monthly prices
aapl$YearMonth <- format(aapl$Date, "%Y-%m")

monthly_data <- aapl[!duplicated(aapl$YearMonth, fromLast = TRUE), ]
monthly_data <- monthly_data[order(monthly_data$Date), ]


#Monthly log returns
log_price_m <- log(monthly_data$Close)
r <- diff(log_price_m)

n <- length(r)

#Making sure the number of returns divisible by 12
n12 <- n - (n %% 12)      #removing the remainder
r12 <- r[(n - n12 + 1):n] #keeping the last n12 obs


#Manual periodogram via FFT
m12 <- floor(n12 / 2)

#Frequencies in radians
f12 <- (2*pi/n12) * (1:m12)

#Fourier transform
ft12 <- fft(r12)

#Keeping only positive frequencies
ft12 <- ft12[2:(m12+1)]

#Periodogram
pg12 <- (1/(2*pi*n12)) * (Mod(ft12))^2


#Ploting
par(mar=c(4,4,2,1))

plot(f12, pg12,
     type="l",
     xlab="Frequency (radians)",
     ylab="Periodogram",
     main="Monthly log returns of AAPL")

#Seasonal frequencies
n.years <- n12 / 12               #number of full years
s <- n.years * (1:6)              #indices of seasonal frequencies

points(f12[s], pg12[s],
       pch=20,
       col="red")


# ============================
#Table 1: j=3: First three seasonal frequencies (excl. pi)
# ============================

j <- 3

#Test statistic
T3 <- (4*pi/s2) * sum(pg12[s[1:3]])

#0.95-quantile of chi^2(6)
crit3 <- qchisq(0.95, df = 2*j)

#p-value
p3 <- 1 - pchisq(T3, df = 2*j)

#Displaying the results
cat("Test statistic:", T3, "\n")
cat("0.95-quantile (chi^2(6)):", crit3, "\n")
cat("p-value:", p3, "\n")

#Test statistic: 4.492324 < 0.95-quantile (chi^2(6)): 12.59159 => H0 is not rejected.

# ============================
#Table 2: j=6: All seasonal frequencies (incl. pi)
# ============================

#Test statistic
T6 <- (4*pi/s2) * sum(pg12[s[1:5]]) +
  (2*pi/s2) * pg12[s[6]]

#0.95-quantile of chi^2(11)
crit6 <- qchisq(0.95, df = 11)

#p-value
p6 <- 1 - pchisq(T6, df = 11)

#Displaying the results
cat("Test statistic:", T6, "\n")
cat("0.95-quantile (chi^2(11)):", crit6, "\n")
cat("p-value:", p6, "\n")

#Test statistic: 25.82234 > 0.95-quantile (chi^2(6)): 19.67514 => H0 is rejected.
