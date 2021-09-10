################################################################################################################## 
################################################################################################################## 
## ESTIMATE UNCERTAINTY IN THE MINIMUM POINT OF AN EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL
##
## EXAMPLE FOR MMT IN LONDON, 1993-2006 (MMT=19.1 [18.9,19.3])
## 
## Tobias A, Armstrong B, Gasparrini A. Investigating uncertainty in the minimum mortality temperature: 
## methods and application to 52 Spanish cities. Epidemiology 2017;28(1):72-76.
################################################################################################################## 
################################################################################################################## 

# LOAD LIBRARIES
library(mgcv); library(tsModel); library(dlnm); library(splines)

# LOAD FUNCTION TO FIND MMT
source("findmin.R")

# LOAD LONDON DATA FROM GitHub
data <- read.csv(url('https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/london.csv'))
head(data)

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"; vardegree <- NULL; varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21; lagnk <- 3

# DEFINE THE CROSSBASIS
argvar <- list(fun=varfun, knots=quantile(data$tmean, varper/100, na.rm=T))
arglag <- list(knots=logknots(lag,lagnk))
cb <- crossbasis(data$tmean, lag=lag, argvar=argvar, arglag=arglag)

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 8

# TIME-SERIES REGRESSION MODEL
formula <- death ~ cb + as.factor(dow) + ns(time, df=round(dfseas*length(date)/365.25))
model   <- glm(formula, family=quasipoisson, data=data, na.action="na.exclude")

# OBTAIN MMT   
mmt     <- findmin(cb, model)
simmt   <- findmin(cb, model, sim=T)
mmt_low <- quantile(simmt,c(2.5)/100)
mmt_upp <- quantile(simmt,c(97.5)/100)
cbind(mmt, mmt_low, mmt_upp)

# OBTAIN MMT RESTRICTIED TO 1st-99th PERCENTILES   
min <- quantile(data$tmean,c(1)/100,na.rm=T)
max <- quantile(data$tmean,c(99)/100,na.rm=T)
mmt2   <- findmin(cb,model,from=min,to=max)
simmt2 <- findmin(cb,model,sim=T,from=min,to=max)
mmt2_low <- quantile(simmt2,c(2.5)/100)
mmt2_upp <- quantile(simmt2,c(97.5)/100)
cbind(mmt, mmt2_low, mmt2_upp)

# OBTAIN PREDICTIONS AND PLOT EXPOSURE-RESPONSE WITH MMT
pred <- crosspred(cb, model, cen=mmt)
plot(pred, "overall", ylab="Relative Risk", xlab="Temperature (ºC)", lwd=1.5)
abline(v=mmt) 
abline(v=c(mmt_low, mmt_upp), lty=2)

################################################################################################################## 
################################################################################################################## 
## END OF EXAMPLE
################################################################################################################## 
################################################################################################################## 
