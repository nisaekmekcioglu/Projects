# ============================================================
#Project 3
#Nisa Ekmekcioglu 12225999
# ============================================================
rm(list = ls())
setwd("/Users/nisaekmekcioglu/Desktop/Zeitreihenanalyse")

# ============================
#Tables 3.1 and 3.2 (AIC & BIC model selection)
# ============================

#Setting up

#Preparing the data (same as before)
df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Monthly growth rates (first differences of log series)
y <- log(df$INDPRO)
r <- diff(y)                         

#Demeaned growth rates
r_dm <- r - mean(r, na.rm = TRUE)    

n <- length(r_dm)

#Grid search ARMA(p,q), p<=3, q<=3
max_p <- 3
max_q <- 3


AIC <- matrix(NA_real_, nrow = max_p + 1, ncol = max_q + 1)
BIC <- matrix(NA_real_, nrow = max_p + 1, ncol = max_q + 1)



for (p in 0:max_p) {
  for (q in 0:max_q) {
    
    fit <- tryCatch(
      arima(r_dm,
            order = c(p, 0, q),
            include.mean = FALSE,
            method = "ML"),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      ok_mat[p + 1, q + 1] <- FALSE
      next
    }
    
    #log-likelihood is stored in fit$loglik
    L <- as.numeric(fit$loglik)
    
    #k = number of estimated parameters
    #Here + 1 (innovation variance)
    k <- p + q + 1
    
   
    #AIC(p,q) = -2L + 2k
    #BIC(p,q) = -2L + k*log(n)
    AIC[p + 1, q + 1] <- -2 * L + 2 * k
    BIC[p + 1, q + 1] <- -2 * L + k * log(n)
  }
}


#Adding row/col names (p on rows, q on columns)
dimnames(AIC) <- list(p = 0:max_p, q = 0:max_q)
dimnames(BIC) <- list(p = 0:max_p, q = 0:max_q)


print(round(AIC, 3))

print(round(BIC, 3))


#Identifying minimum AIC / minimum BIC
which_min_AIC <- which(AIC == min(AIC, na.rm = TRUE), arr.ind = TRUE)[1, ]
which_min_BIC <- which(BIC == min(BIC, na.rm = TRUE), arr.ind = TRUE)[1, ]

p_AIC <- as.integer(rownames(AIC)[which_min_AIC[1]])
q_AIC <- as.integer(colnames(AIC)[which_min_AIC[2]])

p_BIC <- as.integer(rownames(BIC)[which_min_BIC[1]])
q_BIC <- as.integer(colnames(BIC)[which_min_BIC[2]])


cat("\nSelected by AIC:  ARMA(", p_AIC, ",", q_AIC, ")\n", sep = "")
cat("Selected by BIC:  ARMA(", p_BIC, ",", q_BIC, ")\n", sep = "")


# ============================================================
#Figure 3.1: the (demeaned) growth rates of the
#quarterly INDPRO with the residuals from the ARMA(3,0) BIC model (from Table 3.2)
# ============================================================

#Data prep (same as before)
df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Monthly growth rates (first differences of log series)
y <- log(df$INDPRO)
r <- diff(y)                              

#Demeaned growth rates
r_dm <- r - mean(r, na.rm = TRUE)         

#Dates for r_dm (diff removes first date)
d <- df$observation_date
d_r <- d[-1]


#Fitting the BIC-selected ARMA model (from Table 3.2)
h <- arima(r_dm,
           order = c(p_BIC, 0, q_BIC),
           include.mean = FALSE,
           method = "ML",
           transform.pars = TRUE)

#Residuals from the fitted ARMA model
res <- as.numeric(h$residuals)

#Fitted values = data - residuals
fit_vals <- r_dm - res


#Plotting
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

#Using the same y-limits across panels 
ylim_all <- range(c(r_dm, res, fit_vals), na.rm = TRUE)

#Demeaned growth rates
plot(d_r, r_dm,
     type = "l",
     xlab = "Time",
     ylab = "",
     main = "",
     ylim = ylim_all)

#Residuals
plot(d_r, res,
     type = "l",
     xlab = "Time",
     ylab = "",
     main = "",
     ylim = ylim_all)

#Fitted values
plot(d_r, fit_vals,
     type = "l",
     xlab = "Time",
     ylab = "",
     main = "",
     ylim = ylim_all)

#This ARMA model explains very little. 
#The residuals are almost of the same size as the data.


# ============================================================
#Figure 3.2: The periodogram of the (demeaned) growth
#rates of the quarterly INDPRO together with the spectral
#densities implied by fitted ARMA models of order (3,0)
#and (3,0), which are the best ARMA models according to
#BIC and AIC, respectively (from Table 3.2). 
# ============================================================

#Data prep (same as before)
df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Monthly growth rates (first differences of log series)
y <- log(df$INDPRO)
r <- diff(y)

#Demeaned growth rates
r_dm <- r - mean(r, na.rm = TRUE)

n <- length(r_dm)


#Getting the raw periodogram of r_dm 
h_pg <- spec.pgram(r_dm, taper = 0, detrend = FALSE, fast = FALSE, plot = FALSE)

#Frequencies in radians 
f  <- 2*pi*h_pg$freq
pg <- h_pg$spec / (2*pi)

#minimum AIC / minimum BIC
p_AIC ; q_AIC
p_BIC ; q_BIC

# Fit the two selected ARMA models
fit_AIC <- arima(r_dm, order = c(p_AIC, 0, q_AIC), include.mean = FALSE, method = "ML")
fit_BIC <- arima(r_dm, order = c(q_BIC, 0, q_BIC), include.mean = FALSE, method = "ML")


#ARMA spectral density at frequencies f (radians)
arma.sp <- function(f, phi, theta, sigma2) {
  nf <- length(f)
  sp <- rep(sigma2 / (2*pi), nf)
  
  q <- length(theta)
  p <- length(phi)
  
  if (q > 0) {
    for (i in 1:nf) {
      num <- 1 + sum(theta * exp(-1i * f[i] * (1:q)))
      sp[i] <- sp[i] * (Mod(num)^2)
    }
  }
  
  if (p > 0) {
    for (i in 1:nf) {
      den <- 1 - sum(phi * exp(-1i * f[i] * (1:p)))
      sp[i] <- sp[i] / (Mod(den)^2)
    }
  }
  
  return(sp)
}

#Extracting AR and MA coefficients 
phi_AIC <- if (p_AIC > 0) as.numeric(fit_AIC$coef[grep("^ar", names(fit_AIC$coef))]) else numeric(0)
theta_AIC  <- if (q_AIC > 0) as.numeric(fit_AIC$coef[grep("^ma", names(fit_AIC$coef))]) else numeric(0)

phi_BIC <- if (p_BIC > 0) as.numeric(fit_BIC$coef[grep("^ar", names(fit_BIC$coef))]) else numeric(0)
theta_BIC  <- if (q_BIC > 0) as.numeric(fit_BIC$coef[grep("^ma", names(fit_BIC$coef))]) else numeric(0)

#sigma^2 estimate (innovation variance)
s2_AIC <- as.numeric(fit_AIC$sigma2)
s2_BIC <- as.numeric(fit_BIC$sigma2)

#Spectral densities implied by the fitted ARMA models
sp_AIC <- arma.sp(f, phi_AIC, theta_AIC, s2_AIC)
sp_BIC <- arma.sp(f, phi_BIC, theta_BIC, s2_BIC)


#Plotting 
par(mfrow = c(1,1))

#Periodogram 
plot(f, pg,
     type = "l",
     xlab = "Frequency (radians)",
     ylab = "Periodogram / Spectral density",
     main = "")

#Adding ARMA spectra (BIC/AIC)
lines(f, sp_BIC, col = "green")
lines(f, sp_AIC, col = "red")


legend("topright",
       legend = c("Periodogram", paste0("BIC: ARMA(", p_BIC, ",", q_BIC, ")"),
                  paste0("AIC: ARMA(", p_AIC, ",", q_AIC, ")")),
       col = c("black", "green", "red"),
       lty = 1, bty = "n")


# ============================================================
#Figure 3.3: Stationarity check for log(INDPRO) via slope d of 
#log periodogram vs  -2*log(sin(omega/2)) near frequency 0
# ============================================================

#Data prep (same as before)
df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log level series 
y <- log(df$INDPRO)


#Manual periodogram for y using FFT 
N <- length(y)
M <- floor(N/2)

#Frequencies in radians: omega_k = (2*pi/N)*k, k=1,...,M
F <- (2*pi/N) * (1:M)

#Fourier transform: keeping only frequencies 1,...,M (excluding k=0 and negatives)
FT <- fft(y)[2:(M+1)]

#Periodogram ordinates
PG <- (1/(2*pi*N)) * (Mod(FT))^2


#Plotting

#Only the 15 lowest frequencies
k_use <- 15

x <- -2 * log( sin(F[1:k_use] / 2) )
y <- log(PG[1:k_use])


par(mfrow = c(1,1), mar = c(4,4,2,1))

plot(x, y,
     type = "o",
     pch = 20,
     xlab = expression(-2*log(sin(omega/2))),
     ylab = expression(log(I(omega))),
     main = "")

#Several reference lines with slope 0.5 
for (a in seq(-10, 10, 0.5)) abline(a = a, b = 0.50, col = "gray")



#Estimating the slope from a simple regression
#log(PG) = alpha + d * [-2 log(sin(omega/2))] + error
fit <- lm(y ~ x)
cat("\nEstimated slope d (using first 15 frequencies):", coef(fit)[2], "\n")

#The log INDPRO is nonstationary because the slope is much greater than 0.5.



