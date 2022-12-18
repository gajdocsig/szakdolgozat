options(scipen = 999) # minden számot decimálisban ír ki
setwd("C:/Users/Dell/OneDrive/Kvantitatív pénzügy/Pénzügyi Matematika MSc/¤ Szakdoga/Számolások_R")

library(quantmod)

library(fitdistrplus)
# library(fBasics)
# library(fGarch)
library(ghyp) 
# library(sn)
# library(ald)
# library(VGAM)
# library(MASS)
library(rugarch)
# library(tseries)
# library(forecast)
# library(aTSA)
# library(MTS)
# library(stats)
library(VarianceGamma)
library(PerformanceAnalytics)


# start_time <- Sys.time()
# end_time <- Sys.time()
# end_time - start_time


### 0. Függvények ###

## FFT call price
FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
    m <- r - log(phi(-(0+1i)))
    phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
    psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * alpha + 1) * v)
    lambda <- (2 * pi)/(N * eta)
    b <- 1/2 * N * lambda
    ku <- -b + lambda * (0:(N - 1))
    v <- eta * (0:(N - 1))
    tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - ((1:N) - 1 == 0))/3
    ft <- fft(tmp)
    res <- exp(-alpha * ku) * ft/pi
    inter <- spline(ku, Re(res), xout = log(K/S0))
    return(inter$y * S0)
}

# VG kar fv
phiVG <- function(u) {
    omega <- (1/nu) * (log(1 - theta * nu - sigma^2 * nu/2))
    tmp <- 1 - (0+1i) * theta * nu * u + 0.5 * sigma^2 * u^2 * nu
    tmp <- tmp^(-1/nu)
    exp((0+1i) * u * log(S0) + u * (r + omega) * (0+1i)) * tmp
}

# BS call price
BlackScholes <- function(S, K, r, T, sig, type){
    
    if(type=="C"){
        d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
        d2 <- d1 - sig*sqrt(T)
        
        value <- S*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
        return(value)}
    
    if(type=="P"){
        d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
        d2 <- d1 - sig*sqrt(T)
        
        value <-  (K*exp(-r*T)*pnorm(-d2) - S*pnorm(-d1))
        return(value)}
}



# Put ár paritásból
put.price <- function(call.price, S0, K, r, T) {
    put.price <- call.price + exp(-r*T)*K - S0
    return(put.price)
}


### 1. Get Data ###

SP500 <- read.table("SP500.txt", header = FALSE)
SP500 <- as.data.frame(SP500)
SP500$V2 <- NULL
SP500$V1 <- as.Date.character(SP500$V1)
SP500_xts <- xts(SP500$V3, order.by = SP500$V1)



### 2. Loghozamok ###

reszvenyhozamok <- dailyReturn(SP500_xts, type = "log") 
loghozam_napi <- as.vector(reszvenyhozamok$daily.returns[-1])

reszvenyhozamok <- weeklyReturn(SP500_xts, type = "log") 
loghozam_heti <- as.vector(reszvenyhozamok$weekly.returns)

reszvenyhozamok <- monthlyReturn(SP500_xts, type = "log") 
loghozam_havi <- as.vector(reszvenyhozamok$monthly.returns)

reszvenyhozamok <- quarterlyReturn(SP500_xts, type = "log") 
loghozam_negyedevi <- as.vector(reszvenyhozamok$quarterly.returns)

reszvenyhozamok <- yearlyReturn(SP500_xts, type = "log") 
loghozam_évi <- as.vector(reszvenyhozamok$yearly.returns)


### 3. Leíró statisztika ###

summary(loghozam_napi)

adf_tesztek = tseries::adf.test(loghozam_napi, k = 25)
adf_tesztek

# Leverage effect

lh_10napi <- diff(loghozam_napi, k, arith = FALSE) - 1
lh_10napi10 <- rep(NA, 1491)
lh_10napi10[1491] <- 0


vol_napi <- rep(NA, 14909)
for (i in 10:14918) {
  vol_napi[i-9] <- sd(loghozam_napi[(i-10):i])
}
vol_10napi10 <- rep(NA, 1491)
vol_10napi10[1] <- 0

for (i in 1:14908) {
  if (i %% 10 == 0) {
    lh_10napi10[i/10] <- lh_10napi[i]
    vol_10napi10[i/10+1] <- vol_napi[i]
  }
}


cor(x = lh_10napi10, y = vol_10napi10)



### 4. ARIMA-GARCH modell ###

## 4.1. Modell specifikáció és illesztés


garchspec <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder=c(1, 1)),
                      mean.model = list(armaOrder=c(0, 0)), distribution.model = "ghyp")
garchfit <- ugarchfit(garchspec, data = loghozam_napi)
infocriteria(garchfit)

resfit <- residuals(garchfit, standardize = T)
resfit_num <- as.numeric(resfit)

acf(abs(resfit))
hist(resfit, breaks = 150)
plot(resfit, type = "l")

ghyptransform(mu = garchfit@fit$coef["mu"], sigma = 1, skew = garchfit@fit$coef["skew"], shape = garchfit@fit$coef["shape"], lambda = garchfit@fit$coef["ghlambda"])

## 4.2. Modell goodness of fit
# 
# # Q-Q plot ghyp
# probs<-(1:length(hozamok))/length(hozamok)
# qqplot(x=GeneralizedHyperbolic::qghyp(probs, mu = 0.2004764, alpha = 0.2312004, 
#                beta = -0.2047319, delta = 2.327656, lambda = -3.761949),
#        y=resfit_num, xlab='Theoretical quantiles', ylab='Empirical quantiles',
#        main='QQ plot for Goodness of Fit for Gamma distribution',
#        pch=16, col='blue')
# abline(a=0, b=1, col='red', lwd=2)
# 
# # Q-Q plot norm
# qqplot(x=qnorm(probs, mean = 0.0003104138, sd = 1),
#        y=resfit_num, xlab='Theoretical quantiles', ylab='Empirical quantiles',
#        main='QQ plot for Goodness of Fit for Gamma distribution',
#        pch=16, col='blue')
# abline(a=0, b=1, col='red', lwd=2)
# 



### 5. S&P 500 szimuláció ###

path <- 1000
day <- 3000
year <- 10
quarter <- 40
month <- 120

# # loghozam szimuláció
# R_matrix <- matrix(NA, day, path) # sorok a napi hozamok, oszlopok a szimulációk
# for (x in 1:path) {
#     sim_garch <- ugarchsim(garchfit, n.sim = day, startMethod = "sample", set.seed = x)
#     R_matrix[, x] <- as.numeric(sim_garch@simulation$seriesSim)
# }
# 
# # export CSV file-ba
# write.csv(R_matrix, "R_matrix.csv")



# ## árfolyam trajektóriák 2021.10.01-től
# S_matrix <- matrix(NA, day + 1, n)
# S_matrix[1,] <- SP500[14919, "V3"]
# 
# for (p in 1:path) {
#     for (d in 1:day) {
#         S_matrix[d+1,p] <- S_matrix[d,p] * exp(R_matrix[d,p])
#     }
# }
# 
# write.csv(S_matrix, "S_matrix.csv")

# # havi loghozamok
# R_M_matrix <- matrix(NA, nrow = 120, ncol = 1000) 
# for (p in 1:path) {
#     for (m in 1:month) {
#         R_M_matrix[m,p] <- sum(R_D_matrix[(21*m-20):(21*m),p])
#     }
# }
# write.csv(R_M_matrix, "R_M_matrix.csv")
# 
# 
# # évi loghozamok
# R_Y_matrix <- matrix(NA, nrow = 10, ncol = 1000) 
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_matrix[y,p] <- sum(R_D_matrix[(252*y-251):(252*y),p])
#     }
# }
# write.csv(R_Y_matrix, "R_Y_matrix.csv")
# 
# #10 éves loghozamok
# R_10Y_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#         R_10Y_vektor[p] <- sum(R_Y_matrix[,p])
# }
# write.csv(R_10Y_vektor, "R_10Y_vektor.csv")


# szimulált loghozamok import
R_D_matrix <- read.csv("R_D_matrix.csv")
R_D_matrix <- R_D_matrix[,-1]
R_M_matrix <- read.csv("R_M_matrix.csv")
R_M_matrix <- R_M_matrix[,-1]
R_Y_matrix <- read.csv("R_Y_matrix.csv")
R_Y_matrix <- R_Y_matrix[,-1]
R_10Y_vektor <- read.csv("R_10Y_vektor.csv")
R_10Y_vektor <- R_10Y_vektor[,-1]


# szimulált árfolyamok import
S_D_matrix <- read.csv("S_D_matrix.csv")
S_D_matrix <- S_D_matrix[,-1]



### 6. Portfóliók ###

year <- 10
quarter <- 40
month <- 120
C0 <- 1000000 # induló tőke

# # évzáró árfolyamadatok
# S_Y_matrix <- matrix(NA, year+2, path)
# for (p in 1:path) {
#     for (y in 1:(year+2)) {
#         S_Y_matrix[y,p] <- S_D_matrix[y*252-251,p]
#     }
# }
# 
# write.csv(S_Y_matrix, "S_Y_matrix.csv")

# # negyedévzáró árfolyamadatok
# S_Q_matrix <- matrix(NA, quarter+2, path)
# for (p in 1:path) {
#     for (q in 1:(quarter+2)) {
#         S_Q_matrix[q,p] <- S_D_matrix[q*63-62,p]
#     }
# }

# # hózáró árfolyamadatok
# S_M_matrix <- matrix(NA, quarter+2, path)
# for (p in 1:path) {
#     for (m in 1:(quarter+2)) {
#         S_M_matrix[m,p] <- S_D_matrix[m*21-20,p]
#     }
# }
# 
# write.csv(S_Y_matrix, "S_Y_matrix.csv")
# write.csv(S_Q_matrix, "S_Q_matrix.csv")
# write.csv(S_M_matrix, "S_M_matrix.csv")

S_Y_matrix <- read.csv("S_Y_matrix.csv") # import éves árfolyamadatok
S_Y_matrix <- S_Y_matrix[-1]
S_Q_matrix <- read.csv("S_Q_matrix.csv") # import negyedéves árfolyamadatok
S_Q_matrix <- S_Q_matrix[-1]
S_M_matrix <- read.csv("S_M_matrix.csv") # import havi árfolyamadatok
S_M_matrix <- S_M_matrix[-1]


r <- 0.02/252 # napi kamat
VG_napi <- vgFit(loghozam_napi)
theta <- VG_napi$param["theta"]
nu <- VG_napi$param["nu"]
sigma <- VG_napi$param["sigma"]

?vgFit


# ## 7.1. 1év/-0% portfólió (évi 1 opció, 0%-os put) ##
# T <- 252
# 
# C_Y0_matrix <- matrix(NA, year+2, path) # készpénz vektora
# C_Y0_matrix[1,] <- C0
# db_Y0_matrix <- matrix(NA, year+1, path) # részvények darabszáma
# call_Y0_matrix <- matrix(NA, year+1, path) # call árak
# put_Y0_matrix <- matrix(NA, year+1, path) # put árak
# 
# for (p in 1:path) {
#     for (y in 1:(year+1)) {
#         S0 <- S_Y_matrix[y,p]
#         call_Y0_matrix[y,p] <- FFTcall.price(phiVG, S0 = S_Y_matrix[y,p], K = S_Y_matrix[y,p], r = r, T = T)
#         put_Y0_matrix[y,p] <- put.price(call_Y0_matrix[y,p], S0 = S_Y_matrix[y,p], K = S_Y_matrix[y,p], r = r, T = T)
#         db_Y0_matrix[y,p] <- C_Y0_matrix[y,p] / (S_Y_matrix[y,p] + put_Y0_matrix[y,p])
#         C_Y0_matrix[y+1,p] <- ifelse(S_Y_matrix[y+1,p] > S_Y_matrix[y,p], db_Y0_matrix[y,p] * S_Y_matrix[y+1,p], db_Y0_matrix[y,p] * S_Y_matrix[y,p])
#     }
# }
# 
# C_Y0_matrix <- C_Y0_matrix[-12,]
# db_Y0_matrix <- db_Y0_matrix[-11,]
# call_Y0_matrix <- call_Y0_matrix[-11,]
# put_Y0_matrix <- put_Y0_matrix[-11,]
# 
# write.csv(C_Y0_matrix, "C_Y0_matrix.csv")
# write.csv(db_Y0_matrix, "db_Y0_matrix.csv")
# write.csv(call_Y0_matrix, "call_Y0_matrix.csv")
# write.csv(put_Y0_matrix, "put_Y0_matrix.csv")

C_Y0_matrix <- read.csv("C_Y0_matrix.csv")
C_Y0_matrix <- C_Y0_matrix[-1]
db_Y0_matrix <- read.csv("db_Y0_matrix.csv")
db_Y0_matrix <- db_Y0_matrix[-1]
call_Y0_matrix <- read.csv("call_Y0_matrix.csv")
call_Y0_matrix <- call_Y0_matrix[-1]
put_Y0_matrix <- read.csv("put_Y0_matrix.csv")
put_Y0_matrix <- put_Y0_matrix[-1]



# ## 7.2. 1év/-20% portfólió (évi 1 opció, 20%-os put) ##
# T <- 252
# 
# C_Y2_matrix <- matrix(NA, year+2, path) # készpénz vektora
# C_Y2_matrix[1,] <- C0
# db_Y2_matrix <- matrix(NA, year+1, path) # részvények darabszáma
# call_Y2_matrix <- matrix(NA, year+1, path) # call árak
# put_Y2_matrix <- matrix(NA, year+1, path) # put árak
# 
# for (p in 1:path) {
#     for (y in 1:(year+1)) {
#         S0 <- S_Y_matrix[y,p]
#         call_Y2_matrix[y,p] <- FFTcall.price(phiVG, S0 = S_Y_matrix[y,p], K = 0.8*S_Y_matrix[y,p], r = r, T = T)
#         put_Y2_matrix[y,p] <- put.price(call_Y2_matrix[y,p], S0 = S_Y_matrix[y,p], K = 0.8*S_Y_matrix[y,p], r = r, T = T)
#         db_Y2_matrix[y,p] <- C_Y2_matrix[y,p] / (S_Y_matrix[y,p] + put_Y2_matrix[y,p])
#         C_Y2_matrix[y+1,p] <- ifelse(S_Y_matrix[y+1,p] > 0.8 * S_Y_matrix[y,p], db_Y2_matrix[y,p] * S_Y_matrix[y+1,p], db_Y2_matrix[y,p] * 0.8 * S_Y_matrix[y,p])
#     }
# }
# 
# C_Y2_matrix <- C_Y2_matrix[-12,]
# db_Y2_matrix <- db_Y2_matrix[-11,]
# call_Y2_matrix <- call_Y2_matrix[-11,]
# put_Y2_matrix <- put_Y2_matrix[-11,]
# 
# write.csv(C_Y2_matrix, "C_Y2_matrix.csv")
# write.csv(db_Y2_matrix, "db_Y2_matrix.csv")
# write.csv(call_Y2_matrix, "call_Y2_matrix.csv")
# write.csv(put_Y2_matrix, "put_Y2_matrix.csv")

C_Y2_matrix <- read.csv("C_Y2_matrix.csv")
C_Y2_matrix <- C_Y2_matrix[-1]
db_Y2_matrix <- read.csv("db_Y2_matrix.csv")
db_Y2_matrix <- db_Y2_matrix[-1]
call_Y2_matrix <- read.csv("call_Y2_matrix.csv")
call_Y2_matrix <- call_Y2_matrix[-1]
put_Y2_matrix <- read.csv("put_Y2_matrix.csv")
put_Y2_matrix <- put_Y2_matrix[-1]



## 7.3. 1év/-40% portfólió (évi 1 opció, 40%-os put) ##

# T <- 252
# 
# C_Y4_matrix <- matrix(NA, year+2, path) # készpénz vektora
# C_Y4_matrix[1,] <- C0
# db_Y4_matrix <- matrix(NA, year+1, path) # részvények darabszáma
# call_Y4_matrix <- matrix(NA, year+1, path) # call árak
# put_Y4_matrix <- matrix(NA, year+1, path) # put árak
# 
# for (p in 1:path) {
#     for (y in 1:(year+1)) {
#         S0 <- S_Y_matrix[y,p]
#         call_Y4_matrix[y,p] <- FFTcall.price(phiVG, S0 = S_Y_matrix[y,p], K = 0.6*S_Y_matrix[y,p], r = r, T = T)
#         put_Y4_matrix[y,p] <- put.price(call_Y4_matrix[y,p], S0 = S_Y_matrix[y,p], K = 0.6*S_Y_matrix[y,p], r = r, T = T)
#         db_Y4_matrix[y,p] <- C_Y4_matrix[y,p] / (S_Y_matrix[y,p] + put_Y4_matrix[y,p])
#         C_Y4_matrix[y+1,p] <- ifelse(S_Y_matrix[y+1,p] > 0.6 * S_Y_matrix[y,p], db_Y4_matrix[y,p] * S_Y_matrix[y+1,p], db_Y4_matrix[y,p] * 0.6*S_Y_matrix[y,p])
#     }
# }
# 
# C_Y4_matrix <- C_Y4_matrix[-12,]
# db_Y4_matrix <- db_Y4_matrix[-11,]
# call_Y4_matrix <- call_Y4_matrix[-11,]
# put_Y4_matrix <- put_Y4_matrix[-11,]
# 
# write.csv(C_Y4_matrix, "C_Y4_matrix.csv")
# write.csv(db_Y4_matrix, "db_Y4_matrix.csv")
# write.csv(call_Y4_matrix, "call_Y4_matrix.csv")
# write.csv(put_Y4_matrix, "put_Y4_matrix.csv")

C_Y4_matrix <- read.csv("C_Y4_matrix.csv")
C_Y4_matrix <- C_Y4_matrix[-1]
db_Y4_matrix <- read.csv("db_Y4_matrix.csv")
db_Y4_matrix <- db_Y4_matrix[-1]
call_Y4_matrix <- read.csv("call_Y4_matrix.csv")
call_Y4_matrix <- call_Y4_matrix[-1]
put_Y4_matrix <- read.csv("put_Y4_matrix.csv")
put_Y4_matrix <- put_Y4_matrix[-1]



# 
# ## 7.4. 1negyedév/-0% portfólió (évi 1 opció, 0%-os put) ##
# 
# T <- 63
# 
# C_Q0_matrix <- matrix(NA, quarter+2, path) # készpénz vektora
# C_Q0_matrix[1,] <- C0
# db_Q0_matrix <- matrix(NA, quarter+1, path) # részvények darabszáma
# call_Q0_matrix <- matrix(NA, quarter+1, path) # call árak
# put_Q0_matrix <- matrix(NA, quarter+1, path) # put árak
# 
# for (p in 1:path) {
#     for (q in 1:(quarter+1)) {
#         S0 <- S_Q_matrix[q,p]
#         call_Q0_matrix[q,p] <- FFTcall.price(phiVG, S0 = S_Q_matrix[q,p], K = S_Q_matrix[q,p], r = r, T = T)
#         put_Q0_matrix[q,p] <- put.price(call_Q0_matrix[q,p], S0 = S_Q_matrix[q,p], K = S_Q_matrix[q,p], r = r, T = T)
#         db_Q0_matrix[q,p] <- C_Q0_matrix[q,p] / (S_Q_matrix[q,p] + put_Q0_matrix[q,p])
#         C_Q0_matrix[q+1,p] <- ifelse(S_Q_matrix[q+1,p] > S_Q_matrix[q,p], db_Q0_matrix[q,p] * S_Q_matrix[q+1,p], db_Q0_matrix[q,p] * S_Q_matrix[q,p])
#     }
# }
# 
# C_Q0_matrix <- C_Q0_matrix[-42,]
# db_Q0_matrix <- db_Q0_matrix[-41,]
# call_Q0_matrix <- call_Q0_matrix[-41,]
# put_Q0_matrix <- put_Q0_matrix[-41,]
# 
# write.csv(C_Q0_matrix, "C_Q0_matrix.csv")
# write.csv(db_Q0_matrix, "db_Q0_matrix.csv")
# write.csv(call_Q0_matrix, "call_Q0_matrix.csv")
# write.csv(put_Q0_matrix, "put_Q0_matrix.csv")

C_Q0_matrix <- read.csv("C_Q0_matrix.csv")
C_Q0_matrix <- C_Q0_matrix[-1]
db_Q0_matrix <- read.csv("db_Q0_matrix.csv")
db_Q0_matrix <- db_Q0_matrix[-1]
call_Q0_matrix <- read.csv("call_Q0_matrix.csv")
call_Q0_matrix <- call_Q0_matrix[-1]
put_Q0_matrix <- read.csv("put_Q0_matrix.csv")
put_Q0_matrix <- put_Q0_matrix[-1]




## 7.5. 1negyedév/-20% portfólió (évi 4 opció, 10%-os put) ##

# T <- 63
# 
# C_Q2_matrix <- matrix(NA, quarter+2, path) # készpénz vektora
# C_Q2_matrix[1,] <- C0
# db_Q2_matrix <- matrix(NA, quarter+1, path) # részvények darabszáma
# call_Q2_matrix <- matrix(NA, quarter+1, path) # call árak
# put_Q2_matrix <- matrix(NA, quarter+1, path) # put árak
# 
# for (p in 1:path) {
#     for (q in 1:(quarter+1)) {
#         S0 <- S_Q_matrix[q,p]
#         call_Q2_matrix[q,p] <- FFTcall.price(phiVG, S0 = S_Q_matrix[q,p], K = 0.9*S_Q_matrix[q,p], r = r, T = T)
#         put_Q2_matrix[q,p] <- put.price(call_Q2_matrix[q,p], S0 = S_Q_matrix[q,p], K = 0.9*S_Q_matrix[q,p], r = r, T = T)
#         db_Q2_matrix[q,p] <- C_Q2_matrix[q,p] / (S_Q_matrix[q,p] + put_Q2_matrix[q,p])
#         C_Q2_matrix[q+1,p] <- ifelse(S_Q_matrix[q+1,p] > 0.9 * S_Q_matrix[q,p], db_Q2_matrix[q,p] * S_Q_matrix[q+1,p], db_Q2_matrix[q,p] * 0.9*S_Q_matrix[q,p])
#     }
# }
# 
# C_Q2_matrix <- C_Q2_matrix[-42,]
# db_Q2_matrix <- db_Q2_matrix[-41,]
# call_Q2_matrix <- call_Q2_matrix[-41,]
# put_Q2_matrix <- put_Q2_matrix[-41,]
# 
# write.csv(C_Q2_matrix, "C_Q2_matrix.csv")
# write.csv(db_Q2_matrix, "db_Q2_matrix.csv")
# write.csv(call_Q2_matrix, "call_Q2_matrix.csv")
# write.csv(put_Q2_matrix, "put_Q2_matrix.csv")

C_Q2_matrix <- read.csv("C_Q2_matrix.csv")
C_Q2_matrix <- C_Q2_matrix[-1]
db_Q2_matrix <- read.csv("db_Q2_matrix.csv")
db_Q2_matrix <- db_Q2_matrix[-1]
call_Q2_matrix <- read.csv("call_Q2_matrix.csv")
call_Q2_matrix <- call_Q2_matrix[-1]
put_Q2_matrix <- read.csv("put_Q2_matrix.csv")
put_Q2_matrix <- put_Q2_matrix[-1]



## 7.6. 1negyedév/-40% portfólió (évi 4 opció, 20%-os put) ##

# T <- 63
# 
# C_Q4_matrix <- matrix(NA, quarter+2, path) # készpénz vektora
# C_Q4_matrix[1,] <- C0
# db_Q4_matrix <- matrix(NA, quarter+1, path) # részvények darabszáma
# call_Q4_matrix <- matrix(NA, quarter+1, path) # call árak
# put_Q4_matrix <- matrix(NA, quarter+1, path) # put árak
# 
# for (p in 1:path) {
#     for (q in 1:(quarter+1)) {
#         S0 <- S_Q_matrix[q,p]
#         call_Q4_matrix[q,p] <- FFTcall.price(phiVG, S0 = S_Q_matrix[q,p], K = 0.8*S_Q_matrix[q,p], r = r, T = T)
#         put_Q4_matrix[q,p] <- put.price(call_Q4_matrix[q,p], S0 = S_Q_matrix[q,p], K = 0.8*S_Q_matrix[q,p], r = r, T = T)
#         db_Q4_matrix[q,p] <- C_Q4_matrix[q,p] / (S_Q_matrix[q,p] + put_Q4_matrix[q,p])
#         C_Q4_matrix[q+1,p] <- ifelse(S_Q_matrix[q+1,p] > 0.8 * S_Q_matrix[q,p], db_Q4_matrix[q,p] * S_Q_matrix[q+1,p], db_Q4_matrix[q,p] * 0.8*S_Q_matrix[q,p])
#     }
# }
# 
# C_Q4_matrix <- C_Q4_matrix[-42,]
# db_Q4_matrix <- db_Q4_matrix[-41,]
# call_Q4_matrix <- call_Q4_matrix[-41,]
# put_Q4_matrix <- put_Q4_matrix[-41,]
# 
# write.csv(C_Q4_matrix, "C_Q4_matrix.csv")
# write.csv(db_Q4_matrix, "db_Q4_matrix.csv")
# write.csv(call_Q4_matrix, "call_Q4_matrix.csv")
# write.csv(put_Q4_matrix, "put_Q4_matrix.csv")

C_Q4_matrix <- read.csv("C_Q4_matrix.csv")
C_Q4_matrix <- C_Q4_matrix[-1]
db_Q4_matrix <- read.csv("db_Q4_matrix.csv")
db_Q4_matrix <- db_Q4_matrix[-1]
call_Q4_matrix <- read.csv("call_Q4_matrix.csv")
call_Q4_matrix <- call_Q4_matrix[-1]
put_Q4_matrix <- read.csv("put_Q4_matrix.csv")
put_Q4_matrix <- put_Q4_matrix[-1]




# ## 7.7. 1hónap/-0% portfólió (évi 12 opció, 0%-os put) ##
# 
# T <- 21
# 
# C_M0_matrix <- matrix(NA, month+2, path) # készpénz vektora
# C_M0_matrix[1,] <- C0
# db_M0_matrix <- matrix(NA, month+1, path) # részvények darabszáma
# call_M0_matrix <- matrix(NA, month+1, path) # call árak
# put_M0_matrix <- matrix(NA, month+1, path) # put árak
# 
# for (p in 1:path) {
#     for (m in 1:(month+1)) {
#         S0 <- S_M_matrix[m,p]
#         call_M0_matrix[m,p] <- FFTcall.price(phiVG, S0 = S_M_matrix[m,p], K = S_M_matrix[m,p], r = r, T = T)
#         put_M0_matrix[m,p] <- put.price(call_M0_matrix[m,p], S0 = S_M_matrix[m,p], K = S_M_matrix[m,p], r = r, T = T)
#         db_M0_matrix[m,p] <- C_M0_matrix[m,p] / (S_M_matrix[m,p] + put_M0_matrix[m,p])
#         C_M0_matrix[m+1,p] <- ifelse(S_M_matrix[m+1,p] > S_M_matrix[m,p], db_M0_matrix[m,p] * S_M_matrix[m+1,p], db_M0_matrix[m,p] * S_M_matrix[m,p])
#     }
# }
# 
# C_M0_matrix <- C_M0_matrix[-122,]
# db_M0_matrix <- db_M0_matrix[-121,]
# call_M0_matrix <- call_M0_matrix[-121,]
# put_M0_matrix <- put_M0_matrix[-121,]
# 
# write.csv(C_M0_matrix, "C_M0_matrix.csv")
# write.csv(db_M0_matrix, "db_M0_matrix.csv")
# write.csv(call_M0_matrix, "call_M0_matrix.csv")
# write.csv(put_M0_matrix, "put_M0_matrix.csv")

C_M0_matrix <- read.csv("C_M0_matrix.csv")
C_M0_matrix <- C_M0_matrix[-1]
db_M0_matrix <- read.csv("db_M0_matrix.csv")
db_M0_matrix <- db_M0_matrix[-1]
call_M0_matrix <- read.csv("call_M0_matrix.csv")
call_M0_matrix <- call_M0_matrix[-1]
put_M0_matrix <- read.csv("put_M0_matrix.csv")
put_M0_matrix <- put_M0_matrix[-1]




## 7.8. 1hónap/-20% portfólió (évi 12 opció, 5.77%-os put) ##

# T <- 21
# 
# C_M2_matrix <- matrix(NA, month+2, path) # készpénz vektora
# C_M2_matrix[1,] <- C0
# db_M2_matrix <- matrix(NA, month+1, path) # részvények darabszáma
# call_M2_matrix <- matrix(NA, month+1, path) # call árak
# put_M2_matrix <- matrix(NA, month+1, path) # put árak
# 
# for (p in 1:path) {
#     for (m in 1:(month+1)) {
#         S0 <- S_M_matrix[m,p]
#         call_M2_matrix[m,p] <- FFTcall.price(phiVG, S0 = S_M_matrix[m,p], K = 0.9423 * S_M_matrix[m,p], r = r, T = T)
#         put_M2_matrix[m,p] <- put.price(call_M2_matrix[m,p], S0 = S_M_matrix[m,p], K = 0.9423 * S_M_matrix[m,p], r = r, T = T)
#         db_M2_matrix[m,p] <- C_M2_matrix[m,p] / (S_M_matrix[m,p] + put_M2_matrix[m,p])
#         C_M2_matrix[m+1,p] <- ifelse(S_M_matrix[m+1,p] > 0.9423 * S_M_matrix[m,p], db_M2_matrix[m,p] * S_M_matrix[m+1,p], db_M2_matrix[m,p] * 0.9423 * S_M_matrix[m,p])
#     }
#     print(p)
# }
# 
# C_M2_matrix <- C_M2_matrix[-122,]
# db_M2_matrix <- db_M2_matrix[-121,]
# call_M2_matrix <- call_M2_matrix[-121,]
# put_M2_matrix <- put_M2_matrix[-121,]
# 
# write.csv(C_M2_matrix, "C_M2_matrix.csv")
# write.csv(db_M2_matrix, "db_M2_matrix.csv")
# write.csv(call_M2_matrix, "call_M2_matrix.csv")
# write.csv(put_M2_matrix, "put_M2_matrix.csv")

C_M2_matrix <- read.csv("C_M2_matrix.csv")
C_M2_matrix <- C_M2_matrix[-1]
db_M2_matrix <- read.csv("db_M2_matrix.csv")
db_M2_matrix <- db_M2_matrix[-1]
call_M2_matrix <- read.csv("call_M2_matrix.csv")
call_M2_matrix <- call_M2_matrix[-1]
put_M2_matrix <- read.csv("put_M2_matrix.csv")
put_M2_matrix <- put_M2_matrix[-1]





## 7.9. 1hónap/-40% portfólió (évi 12 opció, 11.55%-os put) ##

# T <- 21
# 
# C_M4_matrix <- matrix(NA, month+2, path) # készpénz vektora
# C_M4_matrix[1,] <- C0
# db_M4_matrix <- matrix(NA, month+1, path) # részvények darabszáma
# call_M4_matrix <- matrix(NA, month+1, path) # call árak
# put_M4_matrix <- matrix(NA, month+1, path) # put árak
# 
# for (p in 1:path) {
#     for (m in 1:(month+1)) {
#         S0 <- S_M_matrix[m,p]
#         call_M4_matrix[m,p] <- FFTcall.price(phiVG, S0 = S_M_matrix[m,p], K = 0.8845 * S_M_matrix[m,p], r = r, T = T)
#         put_M4_matrix[m,p] <- put.price(call_M4_matrix[m,p], S0 = S_M_matrix[m,p], K = 0.8845 * S_M_matrix[m,p], r = r, T = T)
#         db_M4_matrix[m,p] <- C_M4_matrix[m,p] / (S_M_matrix[m,p] + put_M4_matrix[m,p])
#         C_M4_matrix[m+1,p] <- ifelse(S_M_matrix[m+1,p] > 0.8845 * S_M_matrix[m,p], db_M4_matrix[m,p] * S_M_matrix[m+1,p], db_M4_matrix[m,p] * 0.8845 * S_M_matrix[m,p])
#     }
#     print(2*p)
# }
# 
# C_M4_matrix <- C_M4_matrix[-122,]
# db_M4_matrix <- db_M4_matrix[-121,]
# call_M4_matrix <- call_M4_matrix[-121,]
# put_M4_matrix <- put_M4_matrix[-121,]
# 
# write.csv(C_M4_matrix, "C_M4_matrix.csv")
# write.csv(db_M4_matrix, "db_M4_matrix.csv")
# write.csv(call_M4_matrix, "call_M4_matrix.csv")
# write.csv(put_M4_matrix, "put_M4_matrix.csv")

C_M4_matrix <- read.csv("C_M4_matrix.csv")
C_M4_matrix <- C_M4_matrix[-1]
db_M4_matrix <- read.csv("db_M4_matrix.csv")
db_M4_matrix <- db_M4_matrix[-1]
call_M4_matrix <- read.csv("call_M4_matrix.csv")
call_M4_matrix <- call_M4_matrix[-1]
put_M4_matrix <- read.csv("put_M4_matrix.csv")
put_M4_matrix <- put_M4_matrix[-1]





### 8. Stratégiák hozamai ###

# # havi loghozamok
# 
# R_M_M0_matrix <- matrix(NA, nrow = 120, ncol = 1000)
# for (p in 1:path) {
#     for (m in 1:(month)) {
#         R_M_M0_matrix[m,p] <- log(C_M0_matrix[m+1,p]/C_M0_matrix[m,p])
#     }
# }
# 
# R_M_M2_matrix <- matrix(NA, nrow = 120, ncol = 1000)
# for (p in 1:path) {
#     for (m in 1:(month)) {
#         R_M_M2_matrix[m,p] <- log(C_M2_matrix[m+1,p]/C_M2_matrix[m,p])
#     }
# }
# 
# R_M_M4_matrix <- matrix(NA, nrow = 120, ncol = 1000)
# for (p in 1:path) {
#     for (m in 1:(month)) {
#         R_M_M4_matrix[m,p] <- log(C_M4_matrix[m+1,p]/C_M4_matrix[m,p])
#     }
# }
# 
# # évi loghozamok
# 
# R_Y_M0_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_M0_matrix[y,p] <- sum(R_M_M0_matrix[(12*y-11):(12*y),p])
#     }
# }
# 
# R_Y_M2_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_M2_matrix[y,p] <- sum(R_M_M2_matrix[(12*y-11):(12*y),p])
#     }
# }
# 
# R_Y_M4_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_M4_matrix[y,p] <- sum(R_M_M4_matrix[(12*y-11):(12*y),p])
#     }
# }
# 
# R_Q_Q0_matrix <- matrix(NA, nrow = 40, ncol = 1000)
# for (p in 1:path) {
#     for (q in 1:(quarter)) {
#         R_Q_Q0_matrix[q,p] <- log(C_Q0_matrix[q+1,p]/C_Q0_matrix[q,p])
#     }
# }
# R_Y_Q0_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Q0_matrix[y,p] <- sum(R_Q_Q0_matrix[(4*y-3):(4*y),p])
#     }
# }
# 
# R_Q_Q2_matrix <- matrix(NA, nrow = 40, ncol = 1000)
# for (p in 1:path) {
#     for (q in 1:(quarter)) {
#         R_Q_Q2_matrix[q,p] <- log(C_Q2_matrix[q+1,p]/C_Q2_matrix[q,p])
#     }
# }
# R_Y_Q2_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Q2_matrix[y,p] <- sum(R_Q_Q2_matrix[(4*y-3):(4*y),p])
#     }
# }
# 
# R_Q_Q4_matrix <- matrix(NA, nrow = 40, ncol = 1000)
# for (p in 1:path) {
#     for (q in 1:(quarter)) {
#         R_Q_Q4_matrix[q,p] <- log(C_Q4_matrix[q+1,p]/C_Q4_matrix[q,p])
#     }
# }
# R_Y_Q4_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Q4_matrix[y,p] <- sum(R_Q_Q4_matrix[(4*y-3):(4*y),p])
#     }
# }
# 
# R_Y_Y0_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Y0_matrix[y,p] <- log(C_Y0_matrix[y+1,p]/C_Y0_matrix[y,p])
#     }
# }
# 
# R_Y_Y2_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Y2_matrix[y,p] <- log(C_Y2_matrix[y+1,p]/C_Y2_matrix[y,p])
#     }
# }
# 
# R_Y_Y4_matrix <- matrix(NA, nrow = 10, ncol = 1000)
# for (p in 1:path) {
#     for (y in 1:year) {
#         R_Y_Y4_matrix[y,p] <- log(C_Y4_matrix[y+1,p]/C_Y4_matrix[y,p])
#     }
# }
# 
# 
# # 10 éves loghozamok
# 
# R_10Y_M0_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#         R_10Y_M0_vektor[p] <- sum(R_Y_M0_matrix[,p])
# }
# 
# R_10Y_M2_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_M2_vektor[p] <- sum(R_Y_M2_matrix[,p])
# }
# 
# R_10Y_M4_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_M4_vektor[p] <- sum(R_Y_M4_matrix[,p])
# }
# 
# R_10Y_Q0_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Q0_vektor[p] <- sum(R_Y_Q0_matrix[,p])
# }
# 
# R_10Y_Q2_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Q2_vektor[p] <- sum(R_Y_Q2_matrix[,p])
# }
# 
# R_10Y_Q4_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Q4_vektor[p] <- sum(R_Y_Q4_matrix[,p])
# }
# 
# R_10Y_Y0_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Y0_vektor[p] <- sum(R_Y_Y0_matrix[,p])
# }
# 
# R_10Y_Y2_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Y2_vektor[p] <- sum(R_Y_Y2_matrix[,p])
# }
# 
# R_10Y_Y4_vektor <- rep(NA, 1000)
# for (p in 1:path) {
#     R_10Y_Y4_vektor[p] <- sum(R_Y_Y4_matrix[,p])
# }
# 
# 
# # adatok exportálása
# write.csv(R_M_M0_matrix, "R_M_M0_matrix.csv")
# write.csv(R_M_M2_matrix, "R_M_M2_matrix.csv")
# write.csv(R_M_M4_matrix, "R_M_M4_matrix.csv")
# write.csv(R_Y_M0_matrix, "R_Y_M0_matrix.csv")
# write.csv(R_Y_M2_matrix, "R_Y_M2_matrix.csv")
# write.csv(R_Y_M4_matrix, "R_Y_M4_matrix.csv")
# write.csv(R_Y_Q0_matrix, "R_Y_Q0_matrix.csv")
# write.csv(R_Y_Q2_matrix, "R_Y_Q2_matrix.csv")
# write.csv(R_Y_Q4_matrix, "R_Y_Q4_matrix.csv")
# write.csv(R_Y_Y0_matrix, "R_Y_Y0_matrix.csv")
# write.csv(R_Y_Y2_matrix, "R_Y_Y2_matrix.csv")
# write.csv(R_Y_Y4_matrix, "R_Y_Y4_matrix.csv")
# write.csv(R_10Y_M0_vektor, "R_10Y_M0_vektor.csv")
# write.csv(R_10Y_M2_vektor, "R_10Y_M2_vektor.csv")
# write.csv(R_10Y_M4_vektor, "R_10Y_M4_vektor.csv")
# write.csv(R_10Y_Q0_vektor, "R_10Y_Q0_vektor.csv")
# write.csv(R_10Y_Q2_vektor, "R_10Y_Q2_vektor.csv")
# write.csv(R_10Y_Q4_vektor, "R_10Y_Q4_vektor.csv")
# write.csv(R_10Y_Y0_vektor, "R_10Y_Y0_vektor.csv")
# write.csv(R_10Y_Y2_vektor, "R_10Y_Y2_vektor.csv")
# write.csv(R_10Y_Y4_vektor, "R_10Y_Y4_vektor.csv")

# adatok importálása

R_M_M0_matrix <- read.csv("R_M_M0_matrix.csv")
R_M_M2_matrix <- read.csv("R_M_M2_matrix.csv")
R_M_M4_matrix <- read.csv("R_M_M4_matrix.csv")
R_Y_M0_matrix <- read.csv("R_Y_M0_matrix.csv")
R_Y_M2_matrix <- read.csv("R_Y_M2_matrix.csv")
R_Y_M4_matrix <- read.csv("R_Y_M4_matrix.csv")
R_Y_Q0_matrix <- read.csv("R_Y_Q0_matrix.csv")
R_Y_Q2_matrix <- read.csv("R_Y_Q2_matrix.csv")
R_Y_Q4_matrix <- read.csv("R_Y_Q4_matrix.csv")
R_Y_Y0_matrix <- read.csv("R_Y_Y0_matrix.csv")
R_Y_Y2_matrix <- read.csv("R_Y_Y2_matrix.csv")
R_Y_Y4_matrix <- read.csv("R_Y_Y4_matrix.csv")
R_10Y_M0_vektor <- read.csv("R_10Y_M0_vektor.csv")
R_10Y_M2_vektor <- read.csv("R_10Y_M2_vektor.csv")
R_10Y_M4_vektor <- read.csv("R_10Y_M4_vektor.csv")
R_10Y_Q0_vektor <- read.csv("R_10Y_Q0_vektor.csv")
R_10Y_Q2_vektor <- read.csv("R_10Y_Q2_vektor.csv")
R_10Y_Q4_vektor <- read.csv("R_10Y_Q4_vektor.csv")
R_10Y_Y0_vektor <- read.csv("R_10Y_Y0_vektor.csv")
R_10Y_Y2_vektor <- read.csv("R_10Y_Y2_vektor.csv")
R_10Y_Y4_vektor <- read.csv("R_10Y_Y4_vektor.csv")
R_M_M0_matrix <- R_M_M0_matrix[,-1]
R_M_M2_matrix <- R_M_M2_matrix[,-1]
R_M_M4_matrix <- R_M_M4_matrix[,-1]
R_Y_M0_matrix <- R_Y_M0_matrix[,-1]
R_Y_M2_matrix <- R_Y_M2_matrix[,-1]
R_Y_M4_matrix <- R_Y_M4_matrix[,-1]
R_Y_Q0_matrix <- R_Y_Q0_matrix[,-1]
R_Y_Q2_matrix <- R_Y_Q2_matrix[,-1]
R_Y_Q4_matrix <- R_Y_Q4_matrix[,-1]
R_Y_Y0_matrix <- R_Y_Y0_matrix[,-1]
R_Y_Y2_matrix <- R_Y_Y2_matrix[,-1]
R_Y_Y4_matrix <- R_Y_Y4_matrix[,-1]
R_10Y_M0_vektor <- R_10Y_M0_vektor[,-1]
R_10Y_M2_vektor <- R_10Y_M2_vektor[,-1]
R_10Y_M4_vektor <- R_10Y_M4_vektor[,-1]
R_10Y_Q0_vektor <- R_10Y_Q0_vektor[,-1]
R_10Y_Q2_vektor <- R_10Y_Q2_vektor[,-1]
R_10Y_Q4_vektor <- R_10Y_Q4_vektor[,-1]
R_10Y_Y0_vektor <- R_10Y_Y0_vektor[,-1]
R_10Y_Y2_vektor <- R_10Y_Y2_vektor[,-1]
R_10Y_Y4_vektor <- R_10Y_Y4_vektor[,-1]


# vektor létrehozása

R_M_vektor <- as.vector(as.matrix(R_M_matrix))
R_Y_vektor <- as.vector(as.matrix(R_Y_matrix))
R_M_M0_vektor <- as.vector(as.matrix(R_M_M0_matrix))
R_M_M2_vektor <- as.vector(as.matrix(R_M_M2_matrix))
R_M_M4_vektor <- as.vector(as.matrix(R_M_M4_matrix))
R_Y_M0_vektor <- as.vector(as.matrix(R_Y_M0_matrix))
R_Y_M2_vektor <- as.vector(as.matrix(R_Y_M2_matrix))
R_Y_M4_vektor <- as.vector(as.matrix(R_Y_M4_matrix))
R_Y_Q0_vektor <- as.vector(as.matrix(R_Y_Q0_matrix))
R_Y_Q2_vektor <- as.vector(as.matrix(R_Y_Q2_matrix))
R_Y_Q4_vektor <- as.vector(as.matrix(R_Y_Q4_matrix))
R_Y_Y0_vektor <- as.vector(as.matrix(R_Y_Y0_matrix))
R_Y_Y2_vektor <- as.vector(as.matrix(R_Y_Y2_matrix))
R_Y_Y4_vektor <- as.vector(as.matrix(R_Y_Y4_matrix))


# ### 9. Tesztelés ### 
# #####################################
# 
# R_D_vektor <- as.vector(as.matrix(R_D_matrix))
# R_M_vektor <- as.vector(as.matrix(R_M_matrix))
# R_Y_vektor <- as.vector(as.matrix(R_Y_matrix))
# R_D_vektor <- R_D_vektor[1:2520000]
# R_D_vektor <- R_D_vektor[200001:304918]
# 
# 
# quantile(loghozam_napi, c(.01, 0.10, .25, .50, 0.75, 0.90, 0.99))
# quantile(R_D_vektor, c(.01, 0.10, .25, .50, 0.75, 0.90, 0.99))
# 
# 
# qqplot(loghozam_napi, R_D_vektor)
# 
# plot(quantile(R_D_vektor, seq(0.001, 0.999, by = 0.001)))
# plot(quantile(loghozam_napi, seq(0, 1, by = 0.001)))
# quantile(loghozam_napi, c(.001, .05, .10, 0.90, 0.95, 0.999))


# # Leverage effect
# k = 10
# lh_10napi <- diff(R_D_vektor, k, arith = FALSE) - 1
# lh_10napi10 <- rep(NA, 300000)
# lh_10napi10[300000] <- 0
# 
# vol_napi <- rep(NA, 2999991)
# for (i in 10:3000000) {
#   vol_napi[i-9] <- sd(R_D_vektor[(i-10):i])
# }
# vol_10napi10 <- rep(NA, 300000)
# vol_10napi10[1] <- 0
# 
# for (i in 1:2999991) {
#   if (i %% 10 == 0) {
#     lh_10napi10[i/10] <- lh_10napi[i]
#     vol_10napi10[i/10+1] <- vol_napi[i]
#   }
# }
# 
# 
# cor(x = lh_10napi10, y = vol_10napi10)


### 10. Kockázati mértékek (hozam) ### 
#####################################


# Hozam

mean(R_M_vektor)
mean(R_Y_vektor)
mean(R_10Y_vektor)

mean(R_M_M0_vektor)
mean(R_M_M2_vektor)
mean(R_M_M4_vektor)
mean(R_Y_M0_vektor)
mean(R_Y_M2_vektor)
mean(R_Y_M4_vektor)
mean(R_Y_Q0_vektor)
mean(R_Y_Q2_vektor)
mean(R_Y_Q4_vektor)
mean(R_Y_Y0_vektor)
mean(R_Y_Y2_vektor)
mean(R_Y_Y4_vektor)
mean(R_10Y_M0_vektor)
mean(R_10Y_M2_vektor)
mean(R_10Y_M4_vektor)
mean(R_10Y_Q0_vektor)
mean(R_10Y_Q2_vektor)
mean(R_10Y_Q4_vektor)
mean(R_10Y_Y0_vektor)
mean(R_10Y_Y2_vektor)
mean(R_10Y_Y4_vektor)


## Szórás

# sd(R_M_vektor)
# sd(R_Y_vektor)
# sd(R_10Y_vektor)
# 
# sd(R_M_M0_vektor) 
# sd(R_M_M2_vektor)
# sd(R_M_M4_vektor)
# sd(R_Y_M0_vektor)
# sd(R_Y_M2_vektor) 
# sd(R_Y_M4_vektor)
# sd(R_Y_Q0_vektor) 
# sd(R_Y_Q2_vektor) 
# sd(R_Y_Q4_vektor)
# sd(R_Y_Y0_vektor)
# sd(R_Y_Y2_vektor)
# sd(R_Y_Y4_vektor)
# sd(R_10Y_M0_vektor)
# sd(R_10Y_M2_vektor)
# sd(R_10Y_M4_vektor)
# sd(R_10Y_Q0_vektor)
# sd(R_10Y_Q2_vektor)
# sd(R_10Y_Q4_vektor)
# sd(R_10Y_Y0_vektor)
# sd(R_10Y_Y2_vektor)
# sd(R_10Y_Y4_vektor)


# VaR95

VaR(R_Y_vektor, p = 0.95, method = "historical")

VaR(R_Y_M0_vektor, p = 0.95, method = "historical")
VaR(R_Y_M2_vektor, p = 0.95, method = "historical")
VaR(R_Y_M4_vektor, p = 0.95, method = "historical")
VaR(R_Y_Q0_vektor, p = 0.95, method = "historical")
VaR(R_Y_Q2_vektor, p = 0.95, method = "historical")
VaR(R_Y_Q4_vektor, p = 0.95, method = "historical")
VaR(R_Y_Y0_vektor, p = 0.95, method = "historical")
VaR(R_Y_Y2_vektor, p = 0.95, method = "historical")
VaR(R_Y_Y4_vektor, p = 0.95, method = "historical")



# ES95

ES(R_Y_vektor, p = 0.95, method = "historical")

ES(R_Y_M0_vektor, p = 0.95, method = "historical")
ES(R_Y_M2_vektor, p = 0.95, method = "historical")
ES(R_Y_M4_vektor, p = 0.95, method = "historical")
ES(R_Y_Q0_vektor, p = 0.95, method = "historical")
ES(R_Y_Q2_vektor, p = 0.95, method = "historical")
ES(R_Y_Q4_vektor, p = 0.95, method = "historical")
ES(R_Y_Y0_vektor, p = 0.95, method = "historical")
ES(R_Y_Y2_vektor, p = 0.95, method = "historical")
ES(R_Y_Y4_vektor, p = 0.95, method = "historical")


# Sharpe-ráta

SharpeRatio(R_Y_vektor, Rf = 0.02, FUN = c("StdDev"))



# MDD

#függvények
drawdown <- function(ret) {
  cum.ret <- c(0, cumsum(ret))
  drawdown <- cum.ret - cummax(cum.ret)
  return(tail(drawdown, -1))
}

maxdrawdown <- function(ret) min(drawdown(ret))


mdd <- matrix(NA, nrow = 10, ncol = 1000)
for (j in 1:1000) {
    mdd[10,j] <- maxdrawdown(R_Y_Y4_matrix[,j])
}

for (i in 1:10) {
  print(mean(mdd[i,]))
}






# Maximum drawdown (elmúlt egy év)








### 11. Ábrák ###

dev.off()

# 1. S&P 500 plot


par(mfrow=c(1,1), mar=c(3,5,1,1))
plot(loghozam_napi, type = "l", ylab = "Napi loghozamok", col = "red4", xlab = NULL, ylim = c(-0.20,0.20)) 


# 2. S&P 500 normális illesztve
par(mfrow=c(1,1))
hist(loghozam_napi, breaks = 200, freq = FALSE, ylim = c(0,80),
     xlim = c(-0.30,0.20), col="red4", border = "red4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(loghozam_napi), max(loghozam_napi), length = 2000) 
yfit <- dnorm(xfit, mean = mean(loghozam_napi), sd = sd(loghozam_napi)) 
lines(xfit, yfit, col = "black", lwd = 2)

# 3. S&P 500 aggregált normalitás
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

hist(loghozam_napi, breaks = 250, freq = FALSE, ylim = c(0,60),
         xlim = c(-0.30,0.30), col="red4", border = "red4", main = NULL, 
     xlab = "Napi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(loghozam_napi), max(loghozam_napi), length = 2000) 
yfit <- dnorm(xfit, mean = mean(loghozam_napi), sd = sd(loghozam_napi)) 
lines(xfit, yfit, col = "black", lwd = 2)

hist(loghozam_heti, breaks = 100, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="red4", border = "red4", main = NULL, xlab = "Heti loghozamok", ylab = "Sűrűség")
xfit <- seq(min(loghozam_heti), max(loghozam_heti), length = 2000) 
yfit <- dnorm(xfit, mean = mean(loghozam_heti), sd = sd(loghozam_heti)) 
lines(xfit, yfit, col = "black", lwd = 2)

hist(loghozam_havi, breaks = 30, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="red4", border = "red4", main = NULL, xlab = "Havi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(loghozam_havi), max(loghozam_havi), length = 2000) 
yfit <- dnorm(xfit, mean = mean(loghozam_havi), sd = sd(loghozam_havi)) 
lines(xfit, yfit, col = "black", lwd = 2)


# 4. S&P 500 autokorreláció

par(mfrow=c(1,1), mar=c(5,5,1,1))
acf(loghozam_napi, lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "red4", main = NULL, lwd = 4)
acf(abs(loghozam_napi), lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "red4", main = NULL, lwd = 4)

acf((loghozam_napi)^2)   # nem azonos eloszlású, nem fehérzaj


# 5. GARCH-hatásoktól szűrt
par(mfrow=c(1,1), mar=c(5,5,1,1))
hist(resfit_simple, breaks = 200, freq = FALSE, col="red4", border = "red4", 
     main = NULL, xlab = "Sztenderdizált, szűrt napi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(resfit_simple), max(resfit_simple), length = 2000) 
yfit <- dnorm(xfit, mean = mean(resfit_simple), sd = sd(resfit_simple)) 
lines(xfit, yfit, col = "black", lwd = 2)


# 6. Szimulált S&P500 trajektóriák
matplot(S_D_matrix[1:2520,5:9], lwd = 2.5, lty = 1, type = "l",
        col = "blue4", ylab = "Árfolyam", xlab = "Kereskedési napok")


# 7. Tesztelés, összehasonlítás a historikussal
par(mfrow=c(4,2), mar=c(5,5,1,1))
hist(loghozam_napi, breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="red4", border = "red4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V1"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V2"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V3"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V4"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V5"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V7"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")
hist(R_D_matrix[,"V6"], breaks = 50, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.30,0.30), col="blue4", border = "blue4", main = NULL, xlab = "Napi loghozamok", ylab = "Sűrűség")


# 8. Tesztelés QQ plot

s = seq(0.01, 0.99, 0.01)
R_D_vektor <- as.vector(as.matrix(R_D_matrix))
q_emp <- as.vector(quantile(loghozam_napi, s))
q_sim <- as.vector(quantile(R_D_vektor, s))

qqplot(q_emp, q_sim, col = "red4", bg = "red4", lwd = 3, pch =20, 
       xlab = "Empirikus kvantilisek", ylab = "Szimulált kvantilisek")
abline(0,1, lwd = 2)


# 9. Tesztelés momentumok


R_D_matrix_new <- matrix(R_D_vektor, nrow = 14918, ncol = 201)
mean_vector <- rep(NA, 201)
sd_vector <- rep(NA, 201)
skew_vector <- rep(NA, 201)
kurt_vector <- rep(NA, 201)
for (i in 1:201) {
  mean_vector[i] <- mean(R_D_matrix_new[,i])
  sd_vector[i] <- sd(R_D_matrix_new[,i])
  skew_vector[i] <- skewness(R_D_matrix_new[,i])
  kurt_vector[i] <- kurtosis(R_D_matrix_new[,i])
}

par(mfrow=c(1,1))

hist(mean_vector, breaks = 35, freq = FALSE, 
     col="blue4", border = "blue4", main = NULL, 
     xlab = "Szórás", ylab = "Sűrűség")
abline(v=mean(loghozam_napi), col = "red4", lwd = 3)

par(mfrow=c(1,3))

hist(sd_vector, breaks = 20, freq = FALSE, 
     col="blue4", border = "blue4", main = NULL, 
     xlab = "Szórás", ylab = "Sűrűség")
abline(v=sd(loghozam_napi), col = "red4", lwd = 3)

hist(skew_vector, breaks = 35, freq = FALSE, 
     col="blue4", border = "blue4", main = NULL, 
     xlab = "Ferdeség", ylab = "Sűrűség")
abline(v=skewness(loghozam_napi), col = "red4", lwd = 3)

hist(kurt_vector, breaks = 20, freq = FALSE, 
     col="blue4", border = "blue4", main = NULL, 
     xlab = "Csúcsosság", ylab = "Sűrűség")
abline(v=kurtosis(loghozam_napi), col = "red4", lwd = 3)

par(mfrow=c(1,1))


# 10. S&P 500 aggregált normalitás
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

hist(R_D_vektor, breaks = 2000, freq = FALSE, ylim = c(0,60),
     xlim = c(-0.20,0.20), col="blue4", border = "blue4", main = NULL, 
     xlab = "Napi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(R_D_vektor), max(R_D_vektor), length = 2000) 
yfit <- dnorm(xfit, mean = mean(R_D_vektor), sd = sd(R_D_vektor)) 
lines(xfit, yfit, col = "black", lwd = 2)

hist(R_Y_vektor, breaks = 500, freq = FALSE,
     col="blue4", border = "blue4",  xlim = c(-2,2),
     main = NULL, xlab = "Éves loghozamok", ylab = "Sűrűség")
xfit <- seq(min(R_Y_vektor), max(R_Y_vektor), length = 2000) 
yfit <- dnorm(xfit, mean = mean(R_Y_vektor), sd = sd(R_Y_vektor)) 
lines(xfit, yfit, col = "black", lwd = 2)

hist(R_10Y_vektor, breaks = 100, freq = FALSE,
     col="blue4", border = "blue4",  xlim = c(-6,6),
     main = NULL, xlab = "10 éves loghozamok", ylab = "Sűrűség")
xfit <- seq(min(R_10Y_vektor), max(R_10Y_vektor), length = 2000) 
yfit <- dnorm(xfit, mean = mean(R_10Y_vektor), sd = sd(R_10Y_vektor)) 
lines(xfit, yfit, col = "black", lwd = 2)


# 11. Szimulált hozamok autokorrelációi

par(mfrow=c(2,2), mar=c(5,5,1,1))
R_D_vektor <- as.vector(as.matrix(R_D_matrix))

acf(loghozam_napi, lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "red4", main = NULL, lwd = 4)
acf(R_D_vektor, lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "blue4", main = NULL, lwd = 4)

acf(abs(loghozam_napi), lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "red4", main = NULL, lwd = 4)



acf(abs(R_D_vektor), lag.max = 25, xlab = "Késleltetés", ylab = "Autokorrelációs függvény", col = "blue4", main = NULL, lwd = 4)

acf((R_D_vektor)^2)   # nem azonos eloszlású, nem fehérzaj


# 12. Feltételes vastagfarkúság szimulált hozamokon

garchspec <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder=c(1, 1)),
                        mean.model = list(armaOrder=c(0, 0)), distribution.model = "ghyp")
garchfit <- ugarchfit(garchspec, data = R_D_vektor[1:14918])

resfit_simple <- residuals(garchfit, standardize = T)
resfit_num <- as.numeric(resfit)


par(mfrow=c(1,1), mar=c(5,5,1,1))
hist(resfit_simple, breaks = 200, freq = FALSE, col="blue4", border = "blue4", 
     main = NULL, xlab = "Sztenderdizált, szűrt napi loghozamok", ylab = "Sűrűség")
xfit <- seq(min(resfit_simple), max(resfit_simple), length = 2000) 
yfit <- dnorm(xfit, mean = mean(resfit_simple), sd = sd(resfit_simple)) 
lines(xfit, yfit, col = "black", lwd = 4)





# szimulált loghozam sfv-ek
par(mfrow = c(2,2), mar = c(3,3,2,1))
plot(density(R_D_matrix[1:2520,1]), main = "", xlab = "", xlim = c(-0.2,0.2), ylim = c(0,60))
polygon(density(R_D_matrix[1:2520,1]), col="royalblue", border="royalblue")
plot(density(R_D_matrix[1:2520,2]), main = "", xlab = "", xlim = c(-0.2,0.2), ylim = c(0,60))
polygon(density(R_D_matrix[1:2520,2]), col="navy", border="navy")
plot(density(R_D_matrix[1:2520,3]), main = "", xlab = "", xlim = c(-0.2,0.2), ylim = c(0,60))
polygon(density(R_D_matrix[1:2520,3]), col="steelblue", border="steelblue")
plot(density(R_D_matrix[1:2520,5]), main = "", xlab = "", xlim = c(-0.2,0.2), ylim = c(0,60))
polygon(density(R_D_matrix[1:2520,5]), col="mediumblue", border="mediumblue")