rm(list = ls())

time_start <- Sys.time()

#install.packages("writexl")
library("writexl")

# 0. Data import -------------------------------------------------------------

frec <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Frecuencia_PF.csv")
frec <- as.numeric(frec[,1])

sev <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Severidad_PF.csv")
sev <- as.numeric(sev[,1])

mort <- read.csv2("https://raw.githubusercontent.com/asiergs/csv_to_view/main/Mortalidad_PF.csv")

# 1. Data observation --------------------------------------------------------

## 1.1 Frecuency -------------------------------------------------------------

hist(frec,main = "Frequency histogram",xlab = "Reported accidents",
     ylab = "Count", right = F, xaxt = "n")
axis(1, at = seq(0, 8, 1), labels = seq(0, 8, 1))

## 1.2 Severity --------------------------------------------------------------

hist(sev, breaks = 50,main = "Severity histogram",xlab = "Cost",
     ylab = "Count", xaxt = "n")
axis(1, at = seq(0, 50000, 5000), labels = seq(0, 50000, 5000))

### 3.2.1 Kernel --------------------------------------------------------------

#### 3.2.1.1 Uniform kernel ----------------------------------------------------

f_unif_sample <- function(x,sample, h=0){
  if (h==0){
    h <-  0.9*min(sd(sample),(quantile(sample,0.75)-quantile(sample, 0.25))/1.34)/
      length(sample)^(1/5)}
  n <- length(x)
  out <- c(0)
  for (i in 1:length(x)){
    weights <- dunif((x[i]-sample)/h,0,2)
    out[i] <- sum(weights)/(length(sample)*h)
  }
  out
}

range <- seq(-10000,60000)
value <- f_unif_sample(range, sev)
min(value)
sum(value)

plot(sev, runif(length(sev),0,0.000005),col = "purple", pch = 4, lwd=0.02,
     cex=0.3,xlim=c(-10000,60000),ylim = c(0,0.00020), xlab = "Severity",
     ylab = "Probability", main = "Uniform Kernel")
curve(f_unif_sample(x, sev), from = -10000,
      to=60000, n=4000,col="blue", add = TRUE)

#### 3.2.1.2 Normal Kernel -----------------------------------------------------

f_norm_sample <- function(x,sample, h=0){
  if (h==0){h <- 1.06*sd(sample)/length(sample)^(1/5)}
  n <- length(x)
  out <- c(0)
  for (i in 1:length(x)){
    weights <- dnorm((x[i]-sample)/h,0,1)
    out[i] <- sum(weights)/(length(sample)*h)
  }
  out
}

range <- seq(-10000,60000)
value <- f_norm_sample(range, sev)
min(value)
sum(value)

plot(sev, runif(length(sev),0,0.0000035),col = "purple", pch = 4, lwd=0.02,
     cex=0.3, xlim=c(-10000,60000),ylim = c(0,0.00015), xlab = "Severity",
     ylab = "Probability", main = "Normal Kernel")
curve(f_norm_sample(x, sev), from = -5000,
      to=60000, n=4000,col="blue", add = TRUE)

# 2. Estimation method validation (DGP) --------------------------------------

## 2.1  Poisson --------------------------------------------------------------

### 2.1.1 Error functions ----------------------------------------------------

MM_pois <- function(param, sample){
  (mean(sample)-param)^2
}

PM_pois <- function(param, sample){
  (quantile(sample, 0.8)-qpois(0.8, param))^2
}

ML_pois <- function(param, sample){
  -sum(dpois(sample, param, log = TRUE))
}

### 2.1.2 Parameters calculation ---------------------------------------------

n <- length(frec)

param <- 25
sample_DGP <- rpois(n, param)
mean(sample_DGP)

optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)
optim(c(1), PM_pois, method = "L-BFGS-B", sample = sample_DGP)
optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)

# Since the distribution is discrete, the numerical methods do not properly 
# work. This leads us to the need of another numerical method to estimate the
# parameter based in the Percentile Matching. The following optimization
# function is proposed as a starting point.

optim_discrete <- function(param, FUN,sample, xi=0.1, xf = 20, dx = 1E-1,
                           yi = 0.1, yf = 0.99, dy = 1E-2){
  # for one parameter
  if (length(param)==1){
    range <- seq(xi,xf,dx)
    out <- c()
    for (i in 1:length(range)){
      out[i] <- FUN(range[i],sample)
    }
    sol <- range[out==min(out)]
    list(par = mean(sol), value = min(out))
  }
  # for two parameters
  else if (length(param)==2){
    range_1 <- seq(xi,xf,dx)
    range_2 <- seq(yi,yf,dy)
    out <- matrix(0,length(range_1),length(range_2),
                  dimnames = list(range_1, range_2))
    for (i in 1:length(range_1)){
      for (j in 1:length(range_2)){
        param <- c(range_1[i], range_2[j])
        out[i,j] <- FUN(param, sample)
      }
    }
    pos <- which(out == min(out), arr.ind = T)[1,]
    param <- c(range_1[pos[1]], range_2[pos[2]])
    list(par = param, value = min(out))
  }
}
optim_discrete(c(1), PM_pois, sample_DGP)

### 2.1.3 Bias calculation ---------------------------------------------------

sim_size <- 10000
n <- length(frec)

MM <- c()
PM <- c()
ML <- c()

bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rpois(n, param)
  MM[i] <- optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  PM[i] <- optim_discrete(c(1), PM_pois, sample_DGP, xi=5, xf=40, dx=1E-1)$par
  ML[i] <- optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}
bias_MM <- mean(MM)-param
bias_PM <- mean(PM)-param
bias_ML <- mean(ML)-param

results_pois <- data.frame(method = c("MM", "PM", "ML"),
                           bias = c(bias_MM, bias_PM, bias_ML))
results_pois

### 2.1.4 MSE calculation ----------------------------------------------------

MSE_MM <- var(MM)+bias_MM^2
MSE_PM <- var(PM)+bias_PM^2
MSE_ML <- var(ML)+bias_ML^2

results_pois <- cbind(results_pois, MSE = c(MSE_MM, MSE_PM, MSE_ML))

### 2.1.5 Consistency --------------------------------------------------------

size_min <- 50
size_max <- 10000
step <- 3
size <- seq(size_min, size_max, step)

MM <- c()
PM <- c()
ML <- c()
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rpois(size[i], param)
  MM[i] <- optim(c(1), MM_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  PM[i] <- optim_discrete(c(1), PM_pois, sample_DGP, xi=5, xf=40, dx=1E-1)$par
  ML[i] <- optim(c(1), ML_pois, method = "L-BFGS-B", sample = sample_DGP)$par
  setTxtProgressBar(bar, i)
}

MM_p <- MM-param
PM_p <- PM-param
ML_p <- ML-param

plot(size, MM_p, type = "l", main = "Consistency of MM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, PM_p, type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_p, type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

## 2.2  Negative Binomial ----------------------------------------------------

# Since the r parameter shall be a natural number, the method of moments is not
# valid since there is no way we can set the first moment and the second moment
# with the sample values

### 2.2.1 Error functions ----------------------------------------------------

PM_nbinom <- function(param,setparam, sample){
  r1 <- (quantile(sample, 0.8)-qnbinom(0.8, setparam, param))^2
  r2 <- (quantile(sample, 0.2)-qnbinom(0.2, setparam, param))^2
  r1 + r2
}

ML_nbinom <- function(param, setparam, sample){
  -sum(dnbinom(sample, setparam, param, log = TRUE))
}

### 2.2.2 Parameters calculation ---------------------------------------------

n <- length(frec)

param <- c(25,0.25)
sample_DGP <- rnbinom(n, param[1], param[2])
mean(sample_DGP)
var(sample_DGP)

setparam_range <- seq(23,27)

par_PM <- c()
val_PM <- c()
par_ML <- c()
val_ML <- c()

for (i in 1:length(setparam_range)){
  
  sol_PM <- optim(c(0.3), PM_nbinom, method = "L-BFGS-B",
                  setparam = setparam_range[i], sample = sample_DGP,
                  lower = c(1E-5), upper = (1-1E-5))
  par_PM[i] <- sol_PM$par
  val_PM[i] <- sol_PM$value
  
  sol_ML <- optim(c(0.3), ML_nbinom, method = "L-BFGS-B",
                  setparam = setparam_range[i], sample = sample_DGP,
                  lower = c(1E-5), upper = (1-1E-5))
  par_ML[i] <- sol_ML$par
  val_ML[i] <- sol_ML$value
  
}

pos <- which(val_PM==min(val_PM))[1]
c(setparam_range[pos],par_PM[pos])

pos <- which(val_ML==min(val_ML))
c(setparam_range[pos],par_ML[pos])

### 2.2.3 Bias calculation ---------------------------------------------------

sim_size <- 10000

PM <- matrix(0,sim_size,2,dimnames = list(c(),c("r", "p")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("r", "p")))

setparam_range <- seq(15,35)

par_PM <- c()
val_PM <- c()
par_ML <- c()
val_ML <- c()

bar <- txtProgressBar(0,sim_size,style=3)
for (j in 1:sim_size){
  sample_DGP <- rnbinom(length(frec), param[1], param[2])
  for (i in 1:length(setparam_range)){
    sol_PM <- optim(c(0.3), PM_nbinom, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(1E-5), upper = (1-1E-5))
    par_PM[i] <- sol_PM$par
    val_PM[i] <- sol_PM$value
    
    sol_ML <- optim(c(0.3), ML_nbinom, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(1E-5), upper = (1-1E-5))
    par_ML[i] <- sol_ML$par
    val_ML[i] <- sol_ML$value
    
  }
  
  pos <- which(val_PM==min(val_PM))[1]
  PM[j,] <- c(setparam_range[pos],par_PM[pos])
  
  pos <- which(val_ML==min(val_ML))
  ML[j,] <- c(setparam_range[pos],par_ML[pos])
  setTxtProgressBar(bar, j)
}

bias_r_PM <- mean(PM[,1])-param[1]
bias_r_ML <- mean(ML[,1])-param[1]

bias_p_PM <- mean(PM[,2])-param[2]
bias_p_ML <- mean(ML[,2])-param[2]

results_nbinom <- data.frame(method = c("PM", "ML"),
                             bias_r = c(bias_r_PM, bias_r_ML),
                             bias_p = c(bias_p_PM, bias_p_ML))

### 2.2.4 MSE calculation ----------------------------------------------------

MSE_r_PM <- var(PM[,1])+bias_r_PM^2
MSE_r_ML <- var(ML[,1])+bias_r_ML^2

MSE_p_PM <- var(PM[,2])+bias_p_PM^2
MSE_p_ML <- var(ML[,2])+bias_p_ML^2

results_nbinom <- cbind(results_nbinom,
                        MSE_r = c(MSE_r_PM, MSE_r_ML),
                        MSE_p = c(MSE_p_PM, MSE_p_ML))
results_nbinom

### 2.2.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 10000
step <- 3
size <- seq(size_min, size_max, step)

setparam_range <- seq(20,30)

PM <- matrix(0,length(size),2)
ML <- matrix(0,length(size),2)

par_PM <- c()
val_PM <- c()
par_ML <- c()
val_ML <- c()

bar <- txtProgressBar(0,length(size),style=3)
for (j in 1:length(size)){
  sample_DGP <- rnbinom(size[j], param[1],param[2])
  for (i in 1:length(setparam_range)){
    
    sol_PM <- optim(c(0.3), PM_nbinom, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(1E-5), upper = (1-1E-5))  
    par_PM[i] <- sol_PM$par
    val_PM[i] <- sol_PM$value
    
    sol_ML <- optim(c(0.3), ML_nbinom, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(0.0001), upper = (1-1E-5))
    par_ML[i] <- sol_ML$par
    val_ML[i] <- sol_ML$value
    
  }
  pos <- which(val_PM==min(val_PM))[1]
  PM[j,] <- c(setparam_range[pos],par_PM[pos])
  
  pos <- which(val_ML==min(val_ML))
  ML[j,] <- c(setparam_range[pos],par_ML[pos])
  setTxtProgressBar(bar, j)
}

#### 2.2.5.1 r parameter -----------------------------------------------------

PM_nb <- PM
ML_nb <- ML

PM_nb[,1] <- PM[,1]-param[1]
ML_nb[,1] <- ML[,1]-param[1]

plot(size, PM_nb[,1], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_nb[,1], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

#### 2.2.5.2 p parameter -----------------------------------------------------

PM_nb[,2] <- PM[,2]-param[2]
ML_nb[,2] <- ML[,2]-param[2]

plot(size, PM_nb[,2], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_nb[,2], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

## 2.3 Lognormal -------------------------------------------------------------

### 2.3.1 Error functions ----------------------------------------------------

MM_lognormal<- function(param, sample){
  r1 <- (mean(log(sample)) - param[1])^2
  r2 <- (var(log(sample)) - param[2]^2)^2
  r1 + r2
}

PM_lognormal <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qlnorm(0.8, param[1], param[2]))^2
  r2 <- (quantile(sample, 0.2)-qlnorm(0.2, param[1], param[2]))^2
  r1 + r2
}

ML_lognormal <- function(param, sample){
  -sum(dlnorm(sample, param[1], param[2], log = TRUE))
}

### 2.3.2 Parameters calculation ---------------------------------------------

n <- length(sev)
param <- c(6,1.5)
sample_DGP <- rlnorm(n, param[1], param[2])
mean(sample_DGP)
sd(sample_DGP)

optim(c(5,1.2), MM_lognormal, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(5,1.2), PM_lognormal, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(5,1.2), ML_lognormal, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))

### 2.3.3 Bias calculation ---------------------------------------------------

sim_size <- 10000
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("meanlog", "sdlog")))

n <- length(sev)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rlnorm(n, param[1], param[2])
  MM[i,] <- optim(c(5,1.2), MM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(5,1.2), PM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(5,1.2), ML_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}
bias_meanlog_MM <- mean(MM[,1])-param[1]
bias_meanlog_PM <- mean(PM[,1])-param[1]
bias_meanlog_ML <- mean(ML[,1])-param[1]

bias_sdlog_MM <- mean(MM[,2])-param[2]
bias_sdlog_PM <- mean(PM[,2])-param[2]
bias_sdlog_ML <- mean(ML[,2])-param[2]

results_lognormal <- data.frame(method = c("MM", "PM", "ML"),
                                bias_meanlog = c(bias_meanlog_MM, bias_meanlog_PM,
                                                 bias_meanlog_ML),
                                bias_sdlog = c(bias_sdlog_MM, bias_sdlog_PM,
                                               bias_sdlog_ML))

### 2.3.4 MSE calculation ----------------------------------------------------

MSE_meanlog_MM <- var(MM[,1])+bias_meanlog_MM^2
MSE_meanlog_PM <- var(PM[,1])+bias_meanlog_PM^2
MSE_meanlog_ML <- var(ML[,1])+bias_meanlog_ML^2

MSE_sdlog_MM <- var(MM[,2])+bias_sdlog_MM^2
MSE_sdlog_PM <- var(PM[,2])+bias_sdlog_PM^2
MSE_sdlog_ML <- var(ML[,2])+bias_sdlog_ML^2

results_lognormal <- cbind(results_lognormal,
                           MSE_meanlog = c(MSE_meanlog_MM, MSE_meanlog_PM,
                                           MSE_meanlog_ML),
                           MSE_sdlog = c(MSE_sdlog_MM, MSE_sdlog_PM,
                                         MSE_sdlog_ML))
results_lognormal

### 2.3.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 10000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("meanlog", "sdlog")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rlnorm(size[i], param[1], param[2])
  MM[i,] <- optim(c(5,1.2), MM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP)$par
  PM[i,] <- optim(c(5,1.2), PM_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, upper = c(20,2),
                  lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(5,1.2), ML_lognormal, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}

#### 2.3.5.1 meanlog parameter -----------------------------------------------

MM_ln <- MM
PM_ln <- PM
ML_ln <- ML

MM_ln[,1] <- MM[,1]-param[1]
PM_ln[,1] <- PM[,1]-param[1]
ML_ln[,1] <- ML[,1]-param[1]

plot(size, MM_ln[,1], type = "l", main = "Consistency of MM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, PM_ln[,1], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_ln[,1], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

#### 2.3.5.2 sdlog parameter -------------------------------------------------

MM_ln[,2] <- MM[,2]-param[2]
PM_ln[,2] <- PM[,2]-param[2]
ML_ln[,2] <- ML[,2]-param[2]

plot(size, MM_ln[,2], type = "l", main = "Consistency of MM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, PM_ln[,2], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_ln[,2], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")


## 2.4 Erlang ----------------------------------------------------------------

# Since the k parameter shall be a natural number, the method of moments is not
# valid since there is no way we can set the first moment and the second moment
# with the sample values

### 2.4.1 Error functions ----------------------------------------------------

PM_erlang <- function(param, setparam, sample){
  r1 <- (quantile(sample, 0.8)-qgamma(0.8, setparam,scale = param))^2
  r2 <- (quantile(sample, 0.2)-qgamma(0.2, setparam,scale = param))^2
  r1 + r2
}

ML_erlang <- function(param, setparam, sample){
  -sum(dgamma(sample, setparam, scale = param, log = TRUE))
}

### 2.4.2 Parameters calculation ---------------------------------------------

n <- length(sev)
param <- c(1,50)
sample_DGP <- rgamma(n, param[1],scale = param[2])
mean(sample_DGP)
sd(sample_DGP)

setparam_range <- seq(1,3)

# since the optimization function does not accept any constrain such as setting
# k as natural number, we need to "manually" evaluate the methods for each k
# value and select which of them fits it best

par_PM <- c()
val_PM <- c()
par_ML <- c()
val_ML <- c()

for (i in 1:length(setparam_range)){
  
  sol_PM <-optim(c(40), PM_erlang, method = "L-BFGS-B",
                 setparam=setparam_range[i],
                 sample = sample_DGP, lower = c(0.0001))
  par_PM[i] <- sol_PM$par
  val_PM[i] <- sol_PM$value
  
  sol_ML <- optim(c(40), ML_erlang, method = "L-BFGS-B",
                  setparam=setparam_range[i],
                  sample = sample_DGP, lower = c(0.0001))
  par_ML[i] <- sol_ML$par
  val_ML[i] <- sol_ML$value
  
}

pos <- which(val_PM==min(val_PM))
c(setparam_range[pos],par_PM[pos])

pos <- which(val_ML==min(val_ML))
c(setparam_range[pos],par_ML[pos])

### 2.4.3 Bias calculation ---------------------------------------------------

sim_size <- 10000

PM <- matrix(0,sim_size,2,dimnames = list(c(),c("k", "beta")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("k", "beta")))

setparam_range <- seq(1,3)

par_PM <- c()
val_PM <- c()
par_ML <- c()
val_ML <- c()

bar <- txtProgressBar(0,sim_size,style=3)
for (j in 1:sim_size){
  sample_DGP <- rgamma(length(sev), param[1], scale = param[2])
  for (i in 1:length(setparam_range)){
    sol_PM <- optim(c(40), PM_erlang, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(0.0001))
    par_PM[i] <- sol_PM$par
    val_PM[i] <- sol_PM$value
    
    sol_ML <- optim(c(40), ML_erlang, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(0.0001))
    par_ML[i] <- sol_ML$par
    val_ML[i] <- sol_ML$value
    
  }
  
  pos <- which(val_PM==min(val_PM))
  PM[j,] <-  c(setparam_range[pos],par_PM[pos])
  
  pos <- which(val_ML==min(val_ML))
  ML[j,] <- c(setparam_range[pos],par_ML[pos])
  setTxtProgressBar(bar, j)
}

bias_k_PM <-    mean(PM[,1])-param[1]
bias_beta_PM <- mean(PM[,2])-param[2]
bias_k_ML <-    mean(ML[,1])-param[1]
bias_beta_ML <- mean(ML[,2])-param[2]

results_erlang <- data.frame(method = c("PM", "ML"),
                             bias_k = c(bias_k_PM,bias_k_ML),
                             bias_beta = c(bias_beta_PM,bias_beta_ML))

### 2.4.4 MSE calculation ----------------------------------------------------

MSE_k_PM <-    var(PM[,1])+(mean(PM[,1]-param[1]))^2 
MSE_beta_PM <- var(PM[,2])+(mean(PM[,2]-param[2]))^2 
MSE_k_ML <-    var(ML[,1])+(mean(ML[,1]-param[1]))^2 
MSE_beta_ML <- var(ML[,2])+(mean(ML[,2]-param[2]))^2 

results_erlang <- cbind(results_erlang,
                        MSE_k = c(MSE_k_PM,MSE_k_ML),
                        MSE_beta = c(MSE_beta_PM,MSE_beta_ML))

results_erlang

### 2.4.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 10000
step <- 3
size <- seq(size_min, size_max, step)

setparam_range <- seq(1,3)

PM <- matrix(0,length(size),2)
ML <- matrix(0,length(size),2)

bar <- txtProgressBar(0,length(size),style=3)
for (j in 1:length(size)){
  sample_DGP <- rgamma(size[j], param[1], scale = param[2])
  for (i in 1:length(setparam_range)){
    
    sol_PM <- optim(c(40), PM_erlang, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(0.0001))
    par_PM[i] <- sol_PM$par
    val_PM[i] <- sol_PM$value
    
    sol_ML <- optim(c(40), ML_erlang, method = "L-BFGS-B",
                    setparam=setparam_range[i],
                    sample = sample_DGP, lower = c(0.0001))
    par_ML[i] <- sol_ML$par
    val_ML[i] <- sol_ML$value
    
  }
  
  pos <- which(val_PM==min(val_PM))       
  PM[j,] <- c(setparam_range[pos],par_PM[pos])
  
  pos <- which(val_ML==min(val_ML))
  ML[j,] <- c(setparam_range[pos],par_ML[pos]) 
  setTxtProgressBar(bar, j)
}

#### 2.4.5.1 k parameter -----------------------------------------------------

PM_e <- PM
ML_e <- ML

PM_e[,1] <- PM[,1]-param[1]
ML_e[,1] <- ML[,1]-param[1]

plot(size, PM_e[,1], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_e[,1], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

#### 2.4.5.2 beta parameter --------------------------------------------------

PM_e[,2] <- PM[,2]-param[2]
ML_e[,2] <- ML[,2]-param[2]

plot(size, PM_e[,2], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_e[,2], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

## 2.5 Gamma -----------------------------------------------------------------

### 2.5.1 Error functions ----------------------------------------------------

MM_gamma <- function(param, sample){
  r1 <- (mean(sample)-param[1]*param[2])^2
  r2 <- (var(sample)-param[1]*param[2]^2)^2
  r1+r2
}

PM_gamma <- function(param, sample){
  r1 <- (quantile(sample, 0.8)-qgamma(0.8, param[1], scale = param[2]))^2
  r2 <- (quantile(sample, 0.2)-qgamma(0.2, param[1], scale = param[2]))^2
  r1 + r2
}

ML_gamma <- function(param, sample){
  -sum(dgamma(sample, param[1], scale = param[2], log = TRUE))
}

### 2.5.2 Parameters calculation ---------------------------------------------

n <- length(sev)
param <- c(5,25)
sample_DGP <- rgamma(n, param[1], scale = param[2])
mean(sample_DGP)
sd(sample_DGP)

optim(c(4,20), MM_gamma, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(4,20), PM_gamma, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))
optim(c(4,20), ML_gamma, method = "L-BFGS-B", sample = sample_DGP,
      lower = c(0.001,0.001))

### 2.5.3 Bias calculation ---------------------------------------------------

sim_size <- 10000
MM <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))
PM <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))
ML <- matrix(0,sim_size,2,dimnames = list(c(),c("alpha", "theta")))

n <- length(sev)
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample_DGP <- rgamma(n, param[1], scale = param[2])
  MM[i,] <- optim(c(4,20), MM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  PM[i,] <- optim(c(4,20), PM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(4,20), ML_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}
bias_alpha_MM <- mean(MM[,1])-param[1]
bias_alpha_PM <- mean(PM[,1])-param[1]
bias_alpha_ML <- mean(ML[,1])-param[1]

bias_theta_MM <- mean(MM[,2])-param[2]
bias_theta_PM <- mean(PM[,2])-param[2]
bias_theta_ML <- mean(ML[,2])-param[2]

results_gamma <- data.frame(method = c("MM", "PM", "ML"),
                            bias_alpha = c(bias_alpha_MM,
                                           bias_alpha_PM,
                                           bias_alpha_ML),
                            bias_theta = c(bias_theta_MM, bias_theta_PM,
                                           bias_theta_ML))

### 2.5.4 MSE calculation ----------------------------------------------------

MSE_alpha_MM <- var(MM[,1])+bias_alpha_MM^2
MSE_alpha_PM <- var(PM[,1])+bias_alpha_PM^2
MSE_alpha_ML <- var(ML[,1])+bias_alpha_ML^2

MSE_theta_MM <- var(MM[,2])+bias_theta_MM^2
MSE_theta_PM <- var(PM[,2])+bias_theta_PM^2
MSE_theta_ML <- var(ML[,2])+bias_theta_ML^2

results_gamma <- cbind(results_gamma,
                       MSE_alpha = c(MSE_alpha_MM, MSE_alpha_PM,
                                     MSE_alpha_ML),
                       MSE_theta = c(MSE_theta_MM, MSE_theta_PM,
                                     MSE_theta_ML))
results_gamma

### 2.5.5 Consistency --------------------------------------------------------

size_min <- 150
size_max <- 10000
step <- 3
size <- seq(size_min, size_max, step)

MM <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
PM <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
ML <- matrix(0,length(size),2,dimnames = list(c(),c("alpha", "theta")))
bar <- txtProgressBar(0,length(size),style=3)
for (i in 1:length(size)){
  sample_DGP <- rgamma(size[i], param[1], scale = param[2])
  MM[i,] <- optim(c(4,20), MM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  PM[i,] <- optim(c(4,20), PM_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  ML[i,] <- optim(c(4,20), ML_gamma, method = "L-BFGS-B",
                  sample = sample_DGP, lower = c(0.001,0.001))$par
  setTxtProgressBar(bar, i)
}

#### 2.5.5.1 alpha parameter -------------------------------------------------

MM_g <- MM
PM_g <- PM
ML_g <- ML

MM_g[,1] <- MM[,1]-param[1]
PM_g[,1] <- PM[,1]-param[1]
ML_g[,1] <- ML[,1]-param[1]

plot(size, MM_g[,1], type = "l", main = "Consistency of MM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, PM_g[,1], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_g[,1], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

#### 2.5.5.2 theta parameter -------------------------------------------------

MM_g[,2] <- MM[,2]-param[2]
PM_g[,2] <- PM[,2]-param[2]
ML_g[,2] <- ML[,2]-param[2]

plot(size, MM_g[,2], type = "l", main = "Consistency of MM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red") 
plot(size, PM_g[,2], type = "l", main = "Consistency of PM method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")
plot(size, ML_g[,2], type = "l", main = "Consistency of ML method",
     xlab = "Sample size", ylab = "Bias")
abline(h = 0, col = "red")

# 3. Distributions estimation ------------------------------------------------

## 3.1 Frequency -------------------------------------------------------------

### 3.1.1 As poisson ---------------------------------------------------------

param_pois <- mean(frec)

#### 3.1.1.1 Graphical comparison --------------------------------------------

F_empiric_raw <- function(x, sample=sample_DGP){
  out <- c()
  for (i in 1:length(x)){
    out[i] <- mean(sample<=x[i])
  }
  out
}

dx <- 1E-5
range <- c(-1, -dx, dx)
for (i in 1:10) {
  range[2*i+2]=c(i-dx)
  range[2*i+3]=c(i+dx)
}
empiric <- F_empiric_raw(range,frec)
pois <- ppois(range, param_pois)

plot(range, empiric, type = "l", col = "blue", xlab = "Reported accidents",
     ylab = "Cumulated probability", main = "Poisson")
lines(range, pois, type = "l", col = "red")
legend(x = "bottomright",
       legend=c("Empiric", "Poisson Estimated           "),
       col=c("blue", "red"), lty=1)

qqplot(frec, rpois(1000000, param_pois),
       xlab = "Sample quantiles", ylab = "Theoretical quantiles",
       main = "Poisson")
lines(seq(0,8),seq(0,8), col="blue")

#### 3.1.1.2 Tests -----------------------------------------------------------

chisq_test <- function(sample, dFUN, param, range = seq(0,9)){
  observed <- c()
  for (i in 1:length(range)) observed[i] = sum(sample==range[i])
  if (length(param)==1){
    expected <- dFUN(range, param)*length(sample)
  }
  else{
    expected <- dFUN(range, param[1], param[2])*length(sample)
  }
  sum((observed-expected)^2/expected)
}

KS_test_discrete <- function(sample, pFUN, param, range = seq(0,12)){
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    pestimated <- pFUN(range, param)
  }
  else{
    pestimated <- pFUN(range, param[1], param[2])
  }
  max(abs(empiric-pestimated))
}

CVM_test_discrete <- function(sample, pFUN, dFUN, param, range = seq(0,12)){
  cvm <- c(0)
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    destimated <- dFUN(range, param)
    pestimated <- pFUN(range, param)
  }
  else{
    destimated <- dFUN(range, param[1], param[2])
    pestimated <- pFUN(range, param[1], param[2])
  }
  cvm <- (empiric-pestimated)^2*destimated
  sum(cvm)*length(sample)
}

AD_test_discrete <- function(sample, pFUN, dFUN, param, range = seq(0,12)){
  ad <- c(0)
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    destimated <- dFUN(range, param)
    pestimated <- pFUN(range, param)
  }
  else {
    destimated <- dFUN(range, param[1], param[2])
    pestimated <- pFUN(range, param[1], param[2])
  }
  ad <- (empiric-pestimated)^2/(pestimated*(1-pestimated))*destimated
  sum(ad)*length(sample)
}

sim_size <- 30000
chisq_pois <- c()
ks_pois <- c()
cvm_pois <- c()
ad_pois <- c()
for (i in 1:sim_size){
  sample <- rpois(length(frec),param_pois)
  chisq_pois[i] <- chisq_test(sample, dpois, param_pois, seq(0,9))
  ks_pois[i] <- KS_test_discrete(sample, ppois, param_pois, seq(0,12))
  cvm_pois[i] <- CVM_test_discrete(sample, ppois,dpois, param_pois, seq(0,12))
  ad_pois[i] <- AD_test_discrete(sample, ppois,dpois, param_pois, seq(0,12))
}

hist(chisq_pois)
quantile(chisq_pois, 0.95)
chisq_pois_frec <- chisq_test(frec, dpois, param_pois, seq(0,9))
chisq_pois_pvalue <- mean(chisq_pois>chisq_pois_frec)

# double check
chisq_pois_pvalue
1-pchisq(chisq_pois_frec,9)

hist(ks_pois)
quantile(ks_pois, 0.95)
ks_pois_frec <- KS_test_discrete(frec,ppois, param_pois, seq(0,12))
ks_pois_pvalue <- mean(ks_pois>ks_pois_frec)

hist(cvm_pois)
quantile(cvm_pois, 0.95)
cvm_pois_frec <- CVM_test_discrete(frec,ppois,dpois, param_pois, seq(0,12))
cvm_pois_pvalue <- mean(cvm_pois>cvm_pois_frec)

hist(ad_pois)
quantile(ad_pois, 0.95)
ad_pois_frec <- AD_test_discrete(frec,ppois,dpois, param_pois, seq(0,12))
ad_pois_pvalue <- mean(ad_pois>ad_pois_frec)

### 3.1.2 As negative binomial -----------------------------------------------

setparam_range <- seq(1,500)

par_ML <- c()
val_ML <- c()

for (i in 1:length(setparam_range)){
  
  sol_ML <- optim(c(0.3), ML_nbinom, method = "L-BFGS-B",
                  setparam = setparam_range[i], sample = frec,
                  lower = c(1E-5), upper = (1-1E-5))
  par_ML[i] <- sol_ML$par
  val_ML[i] <- sol_ML$value
  
}

pos <- which(val_ML==min(val_ML))
param_nbinom <- c(setparam_range[pos],par_ML[pos])

#### 3.1.2.1 Graphical comparison --------------------------------------------

dx <- 1E-5
range <- c(-1, -dx, dx)
for (i in 1:10) {
  range[2*i+2]=c(i-dx)
  range[2*i+3]=c(i+dx)
}
empiric <- F_empiric_raw(range,frec)
nbinom <- pnbinom(range, param_nbinom[1], param_nbinom[2])

plot(range, empiric, type = "l", col = "blue", xlab = "Reported accidents",
     ylab = "Cumulated probability", main = "Negative Binomial")
lines(range, nbinom, type = "l", col = "red")
legend(x = "bottomright",
       legend=c("Empiric","Negative Binomial Estimated                 "),
       col=c("blue", "red"), lty=1)

qqplot(frec, rnbinom(1000000, param_nbinom[1], param_nbinom[2]),
       xlab = "Sample quantiles", ylab = "Theoretical quantiles",
       main = "Negative Binomial")
lines(seq(0,8),seq(0,8), col="blue")

#### 3.1.2.2 Tests -----------------------------------------------------------

sim_size <- 30000
chisq_nbinom <- c()
ks_nbinom <- c()
cvm_nbinom <- c()
ad_nbinom <- c()
for (i in 1:sim_size){
  sample <- rnbinom(length(frec),param_nbinom[1], param_nbinom[2])
  chisq_nbinom[i] <- chisq_test(sample, dnbinom, param_nbinom, seq(0,9))
  ks_nbinom[i] <- KS_test_discrete(sample, pnbinom, param_nbinom, seq(0,12))
  cvm_nbinom[i] <- CVM_test_discrete(sample, pnbinom,dnbinom, param_nbinom,
                                     seq(0,12))
  ad_nbinom[i] <- AD_test_discrete(sample, pnbinom,dnbinom, param_nbinom,
                                   seq(0,12))
}

hist(chisq_nbinom)
quantile(chisq_nbinom, 0.95)
chisq_nbinom_frec <- chisq_test(frec, dnbinom, param_nbinom, seq(0,9))
chisq_nbinom_pvalue <- mean(chisq_nbinom>chisq_nbinom_frec)

# double check
chisq_nbinom_pvalue
1-pchisq(chisq_nbinom_frec,9)

hist(ks_nbinom)
quantile(ks_nbinom, 0.95)
ks_nbinom_frec <- KS_test_discrete(frec,pnbinom, param_nbinom, seq(0,12))
ks_nbinom_pvalue <- mean(ks_nbinom>ks_nbinom_frec)

hist(cvm_nbinom)
quantile(cvm_nbinom, 0.95)
cvm_nbinom_frec <- CVM_test_discrete(frec,pnbinom, dnbinom, param_nbinom)
cvm_nbinom_pvalue <- mean(cvm_nbinom>cvm_nbinom_frec)

hist(ad_nbinom)
quantile(ad_nbinom, 0.95)
ad_nbinom_frec <- AD_test_discrete(frec,pnbinom,dnbinom, param_nbinom)
ad_nbinom_pvalue <- mean(ad_nbinom>ad_nbinom_frec)

### 3.1.3 Cross validation --------------------------------------------------

n <- length(frec)

MSE_pois <- c()
MSE_nbinom <- c()

sim <- 1000
bar <- txtProgressBar(0,sim,style=3)
for (i in 1:sim){
  sample_train <- sample(frec,n)
  
  sample_test <- sample_train[1:ceiling(n/4)]
  sample_train <- sample_train[ceiling(n/4+1):n]
  
  param_pois_cv <- mean(sample_train)
  
  setparam_range <- seq(1,1200,1)
  
  par_ML <- c()
  val_ML <- c()
  
  for (k in 1:length(setparam_range)){
    
    sol_ML <- optim(c(0.7), ML_nbinom, method = "L-BFGS-B",
                    setparam = setparam_range[k], sample = sample_train,
                    lower = c(1E-5), upper = (1-1E-8))
    par_ML[k] <- sol_ML$par
    val_ML[k] <- sol_ML$value
    
  }
  
  pos <- which(val_ML==min(val_ML))
  param_nbinom_cv <- c(setparam_range[pos],par_ML[pos])
  
  range <- seq(0,8)
  observed <- c()
  for (j in 1:length(range)) observed[j] <- sum(sample_test==range[j])
  expected_pois <- dpois(range, param_pois_cv)*n/4
  expected_nbinom <- dnbinom(range, param_nbinom_cv[1], param_nbinom_cv[2])*n/4
  
  MSE_pois[i] <- sum((observed-expected_pois)^2)/length(range)
  MSE_nbinom[i] <- sum((observed-expected_nbinom)^2)/length(range)
  setTxtProgressBar(bar, i)
  
}

cv_expected_MSE_pois <- mean(MSE_pois)
cv_expected_MSE_nbinom <- mean(MSE_nbinom)

cv_sd_MSE_pois <- sd(MSE_pois)
cv_ad_MSE_nbinom <- sd(MSE_nbinom)


### 3.1.4 Results summary ----------------------------------------------------

results_frec <- data.frame(distribution = c("Poisson", "Negative binomial"),
                           chisq = c(chisq_pois_frec,chisq_nbinom_frec),
                           chisq_pvalue = c(chisq_pois_pvalue,chisq_nbinom_pvalue),
                           KS = c(ks_pois_frec,ks_nbinom_frec),
                           KS_pvalue = c(ks_pois_pvalue,ks_nbinom_pvalue),
                           CVM = c(cvm_pois_frec,cvm_nbinom_frec),
                           CVM_pvalue = c(cvm_pois_pvalue,cvm_nbinom_pvalue),
                           AD = c(ad_pois_frec,ad_nbinom_frec),
                           AD_pvalue = c(ad_pois_pvalue,ad_nbinom_pvalue),
                           CV_expect_MSE = c(cv_expected_MSE_pois,
                                             cv_expected_MSE_nbinom),
                           CV_sd_MSE = c(cv_sd_MSE_pois,cv_ad_MSE_nbinom))

results_frec

## 3.2 Severity --------------------------------------------------------------

### 3.2.3 As Lognormal -------------------------------------------------------

param_lnorm <- optim(c(8,1.3), ML_lognormal, method = "L-BFGS-B", sample = sev,
                     lower = c(0.001,0.001))$par

#### 3.2.3.1 Graphical comparison --------------------------------------------

range <- seq(-5000, 60000)
lnorm <- plnorm(range, param_lnorm[1],param_lnorm[2])

plot(ecdf(sev), col = "blue", xlab = "Severity",
     main = "Lognormal",ylab = "Cumulated probability")
lines(range, lnorm, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric",
                                   "Lognormal Estimated                                   "),
       col=c("blue", "red"), lty=1)

qqplot(sev, rlnorm(1000000, param_lnorm[1], param_lnorm[2]),
       xlab = "Sample quantiles", ylab = "Theoretical quantiles",
       main="Lognormal", ylim = c(0,250000))
lines(seq(0,50000),seq(0,50000), col="blue")

#### 3.2.3.2 Tests -----------------------------------------------------------

KS_test <- function(sample, pFUN, param, range = c(0,50000), n=100000){
  dx <- (range[2]-range[1])/n
  range <- seq(range[1], range[2], dx)
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    estimated <- pFUN(range,param)
  }
  if (length(param)==2){
    estimated <- pFUN(range,param[1], param[2])
  }
  max(abs(empiric-estimated))
}

CVM_test <- function(sample, pFUN, dFUN, param, range = c(0,50000), n=100000){
  dx <- (range[2]-range[1])/n
  range <- seq(range[1], range[2], dx)
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    pestimated <- pFUN(range,param)
    destimated <- dFUN(range,param)
  }
  if (length(param)==2){
    pestimated <- pFUN(range,param[1], param[2])
    destimated <- dFUN(range,param[1], param[2])
  }
  cvm <- (empiric-pestimated)^2*destimated*dx
  sum(cvm)*length(sample)
}

AD_test <- function(sample, pFUN, dFUN, param, range = c(0,50000), n=100000){
  dx <- (range[2]-range[1])/n
  range <- seq(range[1], range[2], dx)
  empiric <- c()
  for (i in 1:length(range)){
    empiric[i] <- mean(sample<=range[i])
  }
  if (length(param)==1){
    pestimated <- pFUN(range,param)
    destimated <- dFUN(range,param)
  }
  if (length(param)==2){
    pestimated <- pFUN(range,param[1], param[2])
    destimated <- dFUN(range,param[1], param[2])
  }
  cvm <- (empiric-pestimated)^2/(pestimated*(1-pestimated))*destimated*dx
  sum(cvm)*length(sample)
}

sim_size <- 30000
ks_lnorm <- c()
cvm_lnorm <- c()
ad_lnorm <- c()
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample <- rlnorm(length(sev),param_lnorm[1], param_lnorm[2])
  ks_lnorm[i] <- KS_test(sample, plnorm, param_lnorm,
                         c(0.001,50000), n=1000)
  cvm_lnorm[i] <- CVM_test(sample, plnorm, dlnorm, param_lnorm,
                           c(0.001,50000), n=1000)
  ad_lnorm[i] <- AD_test(sample, plnorm, dlnorm, param_lnorm,
                         c(0.001,50000), n=1000)
  
  setTxtProgressBar(bar, i)
}

hist(ks_lnorm)
quantile(ks_lnorm, 0.95)
ks_lnorm_sev <- KS_test(sev, plnorm, param_lnorm,c(0.001,50000), n=1000)
ks_lnorm_pvalue <- mean(ks_lnorm>ks_lnorm_sev)

hist(cvm_lnorm)
quantile(cvm_lnorm, 0.95)
cvm_lnorm_sev <- CVM_test(sev, plnorm, dlnorm, param_lnorm,
                          c(0.001,50000), n=1000)
cvm_lnorm_pvalue <- mean(cvm_lnorm>cvm_lnorm_sev)

hist(ad_lnorm)
quantile(ad_lnorm, 0.95)
ad_lnorm_sev <- AD_test(sev, plnorm, dlnorm, param_lnorm,
                        c(0.001,50000), n=1000)
ad_lnorm_pvalue <- mean(ad_lnorm>ad_lnorm_sev)

### 3.2.4 As Erlang -----------------------------------------------------------

setparam_range <- seq(1,3)

par_ML <- c()
val_ML <- c()
for (i in 1:length(setparam_range)){
  sol_ML <- try(optim(c(4000), ML_erlang, method = "L-BFGS-B",
                      setparam=setparam_range[i],
                      sample = sev, lower = c(0.0001)), T)
  if (typeof(sol_ML)!="list") sol_ML <- list(par = setparam_range[i], value= Inf)
  par_ML[i] <- sol_ML$par
  val_ML[i] <- sol_ML$value
}

pos <- which(val_ML==min(val_ML))
sol_ML <- optim(c(4000), ML_erlang, method = "L-BFGS-B",
                setparam=setparam_range[pos],
                sample = sev, lower = c(0.0001))         
param_erlang <- c(setparam_range[pos],sol_ML$par)

#### 3.2.4.1 Graphical comparison --------------------------------------------

range <- seq(-5000, 60000)
erlang <- pgamma(range, param_erlang[1],scale = param_erlang[2])

plot(ecdf(sev), col = "blue", xlab = "Severity",
     main = "Erlang",ylab = "Cumulated probability")
lines(range, erlang, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric",
                                   "Erlang Estimated                         "),
       col=c("blue", "red"), lty=1)

qqplot(sev, rgamma(100000, param_erlang[1],scale = param_erlang[2]),
       xlab = "Sample quantiles", ylab = "Theoretical quantiles", main="Erlang")
lines(seq(0,50000),seq(0,50000), col="blue")

#### 3.2.4.2 Tests -----------------------------------------------------------

sim_size <- 30000
ks_erlang <- c()
cvm_erlang <- c()
ad_erlang <- c()
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample <- rgamma(length(sev),param_erlang[1],scale = param_erlang[2])
  ks_erlang[i] <- KS_test(sample, pgamma,c(param_erlang[1],1/param_erlang[2]),
                          c(0.001,50000), n=1000)
  cvm_erlang[i] <- CVM_test(sample, pgamma, dgamma, c(param_erlang[1],
                                                      1/param_erlang[2]),
                            c(0.001,50000), n=1000)
  ad_erlang[i] <- AD_test(sample, pgamma, dgamma, c(param_erlang[1],
                                                    1/param_erlang[2]),
                          c(0.001,50000), n=1000)
  setTxtProgressBar(bar, i)
}

hist(ks_erlang)
quantile(ks_erlang, 0.95)
ks_erlang_sev <- KS_test(sev, pgamma,c(param_erlang[1],1/param_erlang[2]),
                         c(0.001,50000), n=1000)
ks_erlang_pvalue <- mean(ks_erlang>ks_erlang_sev)

hist(cvm_erlang)
quantile(cvm_erlang, 0.95)
cvm_erlang_sev <- CVM_test(sev,  pgamma, dgamma, c(param_erlang[1],
                                                   1/param_erlang[2]),
                           c(0.001,50000), n=1000)
cvm_erlang_pvalue <- mean(cvm_erlang>cvm_erlang_sev)

hist(ad_erlang)
quantile(ad_erlang, 0.95)
ad_erlang_sev <- AD_test(sev,  pgamma, dgamma, c(param_erlang[1],
                                                 1/param_erlang[2]),
                         c(0.001,50000), n=1000)
ad_erlang_pvalue <- mean(ad_erlang>ad_erlang_sev)

### 3.2.5 As Gamma -----------------------------------------------------------

param_gamma <- optim(c(0.8,100), ML_gamma, method = "L-BFGS-B", sample = sev,
                     lower = c(0.001,0.001))$par

#### 3.2.5.1 Graphical comparison --------------------------------------------

range <- seq(-5000, 60000)
gamma <- pgamma(range, param_gamma[1],scale = param_gamma[2])

plot(ecdf(sev), col = "blue", xlab = "Severity",
     main = "Gamma",ylab = "Cumulated probability")
lines(range, gamma, type = "l", col = "red")
legend(x = "bottomright", legend=c("Empiric",
                                   "Gamma Estimated                              "),
       col=c("blue", "red"), lty=1)

qqplot(sev, rgamma(1000000, param_gamma[1], scale = param_gamma[2]),
       xlab = "Sample quantiles", ylab = "Theoretical quantiles", main="Gamma")
lines(seq(0,50000),seq(0,50000), col="blue")

#### 3.2.5.2 Tests -----------------------------------------------------------

sim_size <- 30000
ks_gamma <- c()
cvm_gamma <- c()
ad_gamma <- c()
bar <- txtProgressBar(0,sim_size,style=3)
for (i in 1:sim_size){
  sample <- rgamma(length(sev),param_gamma[1], scale = param_erlang[2])
  ks_gamma[i] <- KS_test(sample, pgamma,c(param_gamma[1], 1/param_gamma[2]),
                         c(0.001,50000), n=1000)
  cvm_gamma[i] <- CVM_test(sample, pgamma, dgamma,
                           c(param_gamma[1], 1/param_gamma[2]),
                           c(0.001,50000), n=1000)
  ad_gamma[i] <- AD_test(sample, pgamma, dgamma,
                         c(param_gamma[1], 1/param_gamma[2]),
                         c(0.001,50000), n=1000)
  setTxtProgressBar(bar, i)
}

hist(ks_gamma)
quantile(ks_gamma, 0.95)
ks_gamma_sev <- KS_test(sample, pgamma,c(param_gamma[1], 1/param_gamma[2]),
                        c(0.001,50000), n=1000)
ks_gamma_pvalue <- mean(ks_gamma>ks_gamma_sev)

hist(cvm_gamma)
quantile(cvm_gamma, 0.95)
cvm_gamma_sev <- CVM_test(sample, pgamma, dgamma,
                          c(param_gamma[1], 1/param_gamma[2]),
                          c(0.001,50000), n=1000)
cvm_gamma_pvalue <- mean(cvm_gamma>cvm_gamma_sev)

hist(ad_gamma)
quantile(ad_gamma, 0.95)
ad_gamma_sev <- AD_test(sample, pgamma, dgamma,
                        c(param_gamma[1], 1/param_gamma[2]),
                        c(0.001,50000), n=1000)
ad_gamma_pvalue <- mean(ad_gamma>ad_gamma_sev)

### 3.2.6 Results summary -------------------------------------------------------

results_sev <- data.frame(distribution = c("log-normal", "Erlang", "Gamma"),
                          KS = c(ks_lnorm_sev,
                                 ks_erlang_sev,
                                 ks_gamma_sev),
                          KS_pvalue = c(ks_lnorm_pvalue,
                                        ks_erlang_pvalue,
                                        ks_gamma_pvalue),
                          CVM = c(cvm_lnorm_sev,
                                  cvm_erlang_sev,
                                  cvm_gamma_sev),
                          CVM_pvalue = c(cvm_lnorm_pvalue,
                                         cvm_erlang_pvalue,
                                         cvm_gamma_pvalue),
                          AD = c(ad_lnorm_sev,
                                 ad_erlang_sev,
                                 ad_gamma_sev),
                          AD_pvalue = c(ad_lnorm_pvalue,
                                        ad_erlang_pvalue,
                                        ad_gamma_pvalue))

results_sev     

# 4. Aggregated model for non-life -------------------------------------------

## 4.1 One person accident report --------------------------------------------

# in order to save computation time, we first create 1M scenarios we will
# sample to simulate our policies

sim <- 500000

n <- rpois(sim, param_pois)

cost_nl_single <- c()
for (i in 1:sim){
  cost_nl_single[i] <- sum(rgamma(n[i], param_erlang[1],
                                  scale = param_erlang[2]))
}

est_cost_nl_single <- mean(cost_nl_single)
sd_cost_nl_single <- sd(cost_nl_single)
VaR995_nl_single <- quantile(cost_nl_single,0.995)
TVaR995_nl_single <- mean(cost_nl_single[cost_nl_single>=VaR995_nl_single])


hist(cost_nl_single, breaks = 100, probability = TRUE, xlab = "Cost",
     main = "Single car insurance cost distribution", xaxt="n", right = TRUE,
     xlim = c(0,55000))
axis(1, at = seq(0, 1.2E5, 0.5E4), labels = seq(0, 1.2E5, 0.5E4))

## 4.2 All policies accident reports -----------------------------------------

policies <- 25234

sim <- 500000
cost_nl <- c()

bar <- txtProgressBar(0,sim,style=3)
for (i in 1:sim){
  n <- sum(rpois(policies, param_pois))
  cost_nl[i] <- sum(rgamma(n, param_erlang[1],scale = param_erlang[2]))
  setTxtProgressBar(bar, i)
}



hist(cost_nl, breaks = 100, probability = TRUE, xlab = "Cost (Millions)",
     main = "All car policies cost distribution", xaxt="n")
axis(1, at = seq(2.5E8, 2.8E8, 1E6), labels = seq(2.5E2, 2.8E2, 1))

x <- seq(250E6, 265E6, 1E4)
hist(cost_nl, breaks = 100, probability = TRUE, xlab = "Cost (Millions)",
     main = "All car policies cost distribution", xaxt="n")
lines(x, dnorm(x, mean(cost_nl), sd(cost_nl)), col="blue", n=4000, lwd = 2)
axis(1, at = seq(2.5E8, 2.8E8, 1E6), labels = seq(2.5E2, 2.8E2, 1))


est_cost_nl <- mean(cost_nl)
sd_cost_nl <- sd(cost_nl)
VaR995_nl <- quantile(cost_nl,0.995)
TVaR995_nl <- mean(cost_nl[cost_nl>=VaR995_nl])

# as check, based on the Central Limit theorem, we know that if all the 
# random variables come from the same distribution (which is the case), then
# the expected value shall be equal to n times the single variable expected
# value, also the variance shall be n times the single variable variance

est_cost_nl
est_cost_nl_single*policies

var(cost_nl)
var(cost_nl_single)*policies

# 5. Life model --------------------------------------------------------------

plot(mort$Age, mort$qx, type = "l", xlab = "Age", ylab = "Death probability",
     col="blue")
axis(1, at = seq(0,100,10))

mort$qx <- as.double(mort$qx)

sample_size <- 20000
deceased_sample <- matrix(0, length(mort$Age),sample_size,
                          dimnames = list(mort$Age))

for (i in 1:length(mort[,1])){
  deceased_sample[i,] <- rbinom(sample_size,size = mort$ni[i],prob = mort$qx[i])
}

# double check the death rate is as expected
mean(deceased_sample[81,])/mort$ni[81]
mort$qx[81]

sim <- 500000
deceased <- c()
bar <- txtProgressBar(0,sim,style=3)
for (i in 1:sim){
  total <- c()
  for (j in 1:length(mort[,1])){
    total[j] <- sample(deceased_sample[j,],1)
  }
  deceased[i] <- sum(total)
  setTxtProgressBar(bar, i)
}
amount <- 1E6
cost_l <- deceased*amount

hist(cost_l, breaks = 50, probability = TRUE, xlab = "Cost (Millions)",
     main = "All life policies cost distribution", xaxt="n")
axis(1, at = seq(2E7, 10E7, 5E6), labels = seq(20, 100, 5))

VaR995_deceased <- quantile(deceased, 0.995)
TVaR995_deceased <- mean(deceased[deceased>VaR995_deceased])

est_cost_l <- mean(cost_l)
sd_cost_l <- sd(cost_l)
VaR995_l <- VaR995_deceased*amount
TVaR995_l <- TVaR995_deceased*amount

# double check the expected deaths
sum(mort$qx*mort$ni)
mean(deceased)

sd(deceased)
sqrt(sum(mort$ni*mort$qx*(1-mort$qx)))

# 6. Life and Non-life aggregated model --------------------------------------

sim <- 10000000

# assuming independence between variables
cost_total <- sample(cost_l,sim,TRUE)+sample(cost_nl,sim,TRUE)


hist(cost_total, breaks = 50, probability = TRUE, xlab = "Cost (Millions)",
     main = "Total cost distribution", xaxt="n")
axis(1, at = seq(2E8, 4E8, 5E6), labels = seq(200, 400, 5))

est_cost_total <- mean(cost_total)
sd_cost_total <- sd(cost_total)
VaR995_total <- quantile(cost_total,0.995)
TVaR995_total <- mean(cost_total[cost_total>=VaR995_total])

# Summary

cost_summary <- data.frame(case = c("single_nl_cost", "total_nl_cost",
                                    "total_l_cost", "total_cost"),
                           expected_cost = c(est_cost_nl_single,
                                             est_cost_nl, est_cost_l,
                                             est_cost_total),
                           sd_cost = c(sd_cost_nl_single, sd_cost_nl,
                                       sd_cost_l, sd_cost_total),
                           VaR995 = c(VaR995_nl_single,VaR995_nl,
                                      VaR995_l, VaR995_total),
                           TVaR995 = c(TVaR995_nl_single,TVaR995_nl,
                                       TVaR995_l, TVaR995_total))
cost_summary

time_finish <- Sys.time()
run_time <- time_finish-time_start

# Results export -------------------------------------------------------------

dir.create("Results")
write_xlsx(results_pois,"Results/results_pois.xlsx")
write_xlsx(results_nbinom,"Results/results_nbinomial.xlsx")
write_xlsx(results_lognormal,"Results/results_lognormal.xlsx")
write_xlsx(results_erlang,"Results/results_erlang.xlsx")
write_xlsx(results_gamma,"Results/results_gamma.xlsx")
write_xlsx(results_frec,"Results/results_frec_test.xlsx")
write_xlsx(results_sev,"Results/results_sev.xlsx")
write_xlsx(cost_summary,"Results/cost_summary.xlsx")
