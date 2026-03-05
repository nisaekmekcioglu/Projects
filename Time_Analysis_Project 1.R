# ============================================================
#Project 1 
#Nisa Ekmekcioglu 12225999
# ============================================================
rm(list = ls())
setwd("/Users/nisaekmekcioglu/Desktop/Zeitreihenanalyse")


library(quantmod)

getSymbols("AAPL", src = "yahoo", from = "1980-01-01")

#Creating a data frame with Adjusted Close
df <- data.frame(
  Date = index(AAPL),
  Close = as.numeric(Ad(AAPL))   #adjusted close
)

write.csv(df, "AAPL_full.csv", row.names = FALSE)

# ============================
#Figure 1.1: Historical daily quotes of AAPL
# ============================

df <- read.csv("AAPL_full.csv")

#Converting Date to Date format
df$Date <- as.Date(df$Date)

#Sorting by date 
df <- df[order(df$Date), ]

#Using Close price
P <- df$Close

#Removing invalid observations
ok <- is.finite(P) & P > 0 & !is.na(df$Date)

df <- df[ok, ]
P  <- df$Close

#Computing log prices
logP <- log(P)

#Computing log returns
# r_t = log(P_t) - log(P_{t-1})
r <- diff(logP)

#Corresponding dates for returns
r_dates <- df$Date[-1]

#Squared and 4th power returns
r2 <- r^2
r4 <- r^4


#Ploting
par(mfrow = c(2, 2),mar = c(4, 4, 3, 1))

#Log prices
plot(df$Date, logP,type = "l",main = "Log Prices",xlab = "Time",ylab = "log(P)")

#Log returns
plot(r_dates, r,type = "l",main = "Log Returns",xlab = "Time",ylab = "r[t]")

#Squared log returns
plot(r_dates, r2,type = "l",main = "Squared Log Returns",xlab = "Time",ylab = expression(r[t]^2))

#4th power log returns
plot(r_dates, r4,type = "l",main = "4th Power of Log Returns",xlab = "Time",ylab = expression(r[t]^4))


# ==============================
#Figure 1.2: Polynomial trends of degree 1 (a), 2 (b), 3 (c), 4 (d)
#fitted to the quarterly (seasonally adjusted real) log INDPRO
# ==============================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Parsing and cleaning
df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series and time index
y <- log(df$INDPRO)
t <- seq_along(y)

#Polynomial fits degree 1-4
fit1 <- lm(y ~ poly(t, 1, raw = TRUE))
fit2 <- lm(y ~ poly(t, 2, raw = TRUE))
fit3 <- lm(y ~ poly(t, 3, raw = TRUE))
fit4 <- lm(y ~ poly(t, 4, raw = TRUE))

#y^-
yhat <- list(
  fitted(fit1),
  fitted(fit2),
  fitted(fit3),
  fitted(fit4))


#Ploting
par(mfrow = c(2, 2),mar = c(4, 4, 2, 1))

panel_labels <- c("1 (a)", "2 (b)", "3 (c)", "4 (d)")

for (i in 1:4) {
  
  plot(df$observation_date, y,
       type = "l",
       xlab = "Time",
       ylab = "log(INDPRO)",
       main = panel_labels[i])
  
  lines(df$observation_date, yhat[[i]],
        col = "red",
        lwd = 2)
}


# ==========================================
#Figure 1.3: Deviations of quarterly log INDPRO from
#polynomial trends of degree 1 (a), 2 (b), 3 (c), 4 (d) 
# ==========================================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Parsing and cleaning
df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series and time index
y <- log(df$INDPRO)
t <- seq_along(y)

#Polynomial fits degree 1-4
fit1 <- lm(y ~ poly(t, 1, raw = TRUE))
fit2 <- lm(y ~ poly(t, 2, raw = TRUE))
fit3 <- lm(y ~ poly(t, 3, raw = TRUE))
fit4 <- lm(y ~ poly(t, 4, raw = TRUE))

#Deviations (residuals): u_t = y_t - yhat_t
u1 <- resid(fit1)
u2 <- resid(fit2)
u3 <- resid(fit3)
u4 <- resid(fit4)

u_list <- list(u1, u2, u3, u4)


#Ploting
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

panel_labels <- c("1 (a)", "2 (b)", "3 (c)", "4 (d)")

for (i in 1:4) {
  plot(df$observation_date, u_list[[i]],
       type = "p",
       pch  = 20,         
       cex  = 0.6,
       xlab = "Time",
       ylab = "Deviation",
       main = panel_labels[i])
  abline(h = 0, lty = 2)  # zero line 
}


# ==========================================
#Figure 1.4: Broken linear trends fitted to the quarterly log INDPRO
#sample 1946–2019
# ==========================================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Parsing and cleaning
df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series
y <- log(df$INDPRO)

#Decimal time (so 1973.3 / 2007.2 make sense)
yr  <- as.integer(format(df$observation_date, "%Y"))
mon <- as.integer(format(df$observation_date, "%m"))
time <- yr + (mon - 1) / 12

#Breakpoints 
b1 <- 1973.3
b2 <- 2007.2

#Hinge terms for slope breaks: (t-b)+
h1 <- pmax(0, time - b1)
h2 <- pmax(0, time - b2)

#Dummy for intercept break after 2007.2
D2 <- as.numeric(time >= b2)

#a) break in slope after 1973.3
m_a <- lm(y ~ time + h1)

#b) break in slope after 1973.3 and 2007.2
m_b <- lm(y ~ time + h1 + h2)

#c) (b) + break in intercept after 2007.2
m_c <- lm(y ~ time + h1 + h2 + D2)

fits <- list(fitted(m_a), fitted(m_b), fitted(m_c))

#Converting decimal-year breaks to Date for vertical lines
dec_year_to_date <- function(x) {
  yy <- floor(x)
  mm <- floor((x - yy) * 12) + 1
  as.Date(sprintf("%04d-%02d-01", yy, mm))
}

b1_date <- dec_year_to_date(b1)
b2_date <- dec_year_to_date(b2)

#Ploting
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), mgp = c(2.2, 0.7, 0))

#(a) break in slope after 1973.3, (b) break in slope after 1973.3, 2007.2
#(c) break in slope after 1973.3, 2007.2, break in intercept after 2007.2
panel_labels <- c("(a)", "(b)", "(c)")

for (i in 1:3) {
  plot(df$observation_date, y,
       type = "l",
       xlab = "", ylab = "",
       xaxt = "s", yaxt = "s")
  
  #fitted broken trend
  lines(df$observation_date, fits[[i]], col = "red", lwd = 2)
  
  #break lines
  abline(v = b1_date, col = "darkgreen", lwd = 2)
  if (i >= 2) abline(v = b2_date, col = "blue", lwd = 2)
  

  mtext("Time", side = 1, line = 2.2)
  mtext("log(INDPRO)", side = 2, line = 2.2)
  
 
  mtext(panel_labels[i], side = 4, line = 0.5, las = 0)
}


# ==========================================
#Figure 1.5: Deviations of quarterly log INDPRO from broken trends
#sample 1946–2019
# ==========================================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Parsing and cleaning
df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)
df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series
y <- log(df$INDPRO)

#Decimal time for breakpoints 1973.3 / 2007.2
yr  <- as.integer(format(df$observation_date, "%Y"))
mon <- as.integer(format(df$observation_date, "%m"))
time <- yr + (mon - 1) / 12

#Breakpoints
b1 <- 1973.3
b2 <- 2007.2

#Hinge terms for slope breaks
h1 <- pmax(0, time - b1)
h2 <- pmax(0, time - b2)

#Dummy for intercept break after 2007.2
D2 <- as.numeric(time >= b2)

#Models (a)-(c)
m_a <- lm(y ~ time + h1)               #slope break after 1973.3
m_b <- lm(y ~ time + h1 + h2)          #slope breaks after 1973.3 and 2007.2
m_c <- lm(y ~ time + h1 + h2 + D2)     #slope breaks after 1973.3 and 2007.2
                                       #plus intercept break after 2007.2

#Deviations (residuals)
u_a <- resid(m_a)
u_b <- resid(m_b)
u_c <- resid(m_c)

u_list <- list(u_a, u_b, u_c)

#Converting decimal-year breaks to Date for vertical lines
dec_year_to_date <- function(x) {
  yy <- floor(x)
  mm <- floor((x - yy) * 12) + 1
  as.Date(sprintf("%04d-%02d-01", yy, mm))
}

b1_date <- dec_year_to_date(b1)
b2_date <- dec_year_to_date(b2)

#Ploting
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1), mgp = c(2.2, 0.7, 0))

panel_labels <- c("(a)", "(b)", "(c)")

for (i in 1:3) {
  plot(df$observation_date, u_list[[i]],
       type = "p", pch = 3, cex = 0.6,   # pch=3 gives "+" 
       xlab = "", ylab = "", main = "")
  
  #zero line 
  abline(h = 0, lty = 2)
  
  #break lines
  abline(v = b1_date, col = "darkgreen", lwd = 2)
  if (i >= 2) abline(v = b2_date, col = "blue", lwd = 2)
  

  mtext("Time", side = 1, line = 2.2)
  mtext("Deviation", side = 2, line = 2.2)
  

  mtext(panel_labels[i], side = 4, line = 0.5)
}


# ==========================================
#Figure 1.6: Hodrick-Prescott filter applied to the quarterly log INDPRO
#(a) λ=50000000 (b) λ=1000000 (c) λ=500000 (d) λ=50000
# ==========================================

library(mFilter)

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

#Parsing and cleaning
df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

#Triming the sample
df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series
y <- log(df$INDPRO)

#Applying HP filter with different lambdas
hp1 <- hpfilter(y, freq = 50000000)   # λ = 50,000,000
hp2 <- hpfilter(y, freq = 1000000)   # λ = 1,000,000
hp3 <- hpfilter(y, freq = 500000)    # λ =   500,000
hp4 <- hpfilter(y, freq = 50000)     # λ =    50,000

trends <- list(
  hp1$trend,
  hp2$trend,
  hp3$trend,
  hp4$trend
)

#Ploting
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

panel_labels <- c("(a)  λ = 50,000,000",
                  "(b)  λ = 1,000,000",
                  "(c)  λ =   500,000",
                  "(d)  λ =    50,000")

for (i in 1:4) {
  
  plot(df$observation_date, y,
       type = "l",
       xlab = "Time",
       ylab = "log(INDPRO)",
       main = panel_labels[i])
  
  #HP trend
  lines(df$observation_date, trends[[i]],
        col = "red", lwd = 2)
}


#Large λ → very smooth trend (almost linear)
#Small λ → flexible trend (follows data closely)

# ==========================================
#Figure 1.7: Deviations of quarterly log INDPRO from HP trends
#(a) λ=50000000 (b) λ=1000000 (c) λ=500000 (d) λ=50000
# ==========================================

library(mFilter)

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series
y <- log(df$INDPRO)

#HP filters
hp1 <- hpfilter(y, freq = 5000000)
hp2 <- hpfilter(y, freq = 1000000)
hp3 <- hpfilter(y, freq = 500000)
hp4 <- hpfilter(y, freq = 50000)

cycles <- list(
  hp1$cycle,
  hp2$cycle,
  hp3$cycle,
  hp4$cycle)

#Ploting
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

panel_labels <- c("(a)  λ = 5,000,000",
                  "(b)  λ = 1,000,000",
                  "(c)  λ =   500,000",
                  "(d)  λ =    50,000")

for (i in 1:4) {
  
  plot(df$observation_date, cycles[[i]],
       type = "p",
       pch  = 16,       
       cex  = 0.6,
       xlab = "",
       ylab = "",
       main = panel_labels[i])
  
  abline(h = 0, lty = 2)
  
  mtext("Time", side = 1, line = 2.2)
  mtext("Deviation", side = 2, line = 2.2)
}


# ==========================================
#Figure 1.8: First differences of quarterly log INDPRO
#sample 1946–2019
# ==========================================

df <- read.csv("INDPRO.csv", stringsAsFactors = FALSE)

df$observation_date <- as.Date(df$observation_date)
df$INDPRO <- as.numeric(df$INDPRO)

df <- df[!is.na(df$observation_date) & is.finite(df$INDPRO) & df$INDPRO > 0, ]

df <- df[df$observation_date >= as.Date("1946-01-01") & df$observation_date <= as.Date("2019-12-31"), ]

#Log series
y <- log(df$INDPRO)

#First differences: Δlog(INDPRO)
dy <- diff(y)
dy_dates <- df$observation_date[-1]

#Ploting 
par(mfrow = c(1,1))

plot(dy_dates, dy,
     type = "o",          
     pch  = 16, cex = 0.6,
     xlab = "Time",
     ylab = expression(Delta*log(INDPRO)),
     main = "")
