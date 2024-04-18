################################################################################
# SUPPLEMENTAL MATERIAL of the article:
#   "Health impact projections under climate change scenarios: a tutorial."
#   Ana M. Vicedo-Cabrera, Francesco Sera, Antonio Gasparrini.
#
# This code reproduces the analysis described in the tutorial.
#
# [date]
# * an updated version of this code compatible with future
#   versions of the software is available at the personal website of the
#   last author (https://github.com/gasparrini/....)
################################################################################

# This script is adapted by Eunice Lo for Regions in England
# Using observed ONS mortality & HadUK-Grid temp. data 
# To retrospectively estimate daily temperature-related deaths 
# For all years available, i.e. 1981-2022
# Without climate model projections
# Adapted from Ana M. Vicedo-Cabrera's:
#  - 01EstimationERassociation.r & 
#  - 02ProjTempMortalitySeries.r & 
#  - 04_05_06ExtCurveProjUncert.r 
# Updated 28/03/2023, then 30/12/2023

# LOAD THE PACKAGES
library(dlnm) ; library(splines) ; library(MASS)
library(plyr) ; library(climdex.pcic)

################################################################################
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
################################################################################

# Use the observed daily temperature-mortality series to estimate
#     the coefficients defining the exposure-response association.
# In this example, we estimate the temperature-related mortality association 
#     using data from England and Wales between 2005 and any year in 2014-2018.

# LOAD OBSERVED DATA - DAILY MORTALITY & TEMPERATURE-SERIES BETWEEN 1981-2020 in regions in England and Wales
regions <- c("East Midlands", "East of England", "London", "North East England",
              "North West England", "South East England", "South West England",
              "Wales", "West Midlands", "Yorkshire and Humber")
reg <- regions[menu(regions, title="Please select a region")]
reg1 <- gsub(" ", "", reg, fixed=TRUE)

# mortality and temperature data 
# 1981-2019 (2020 has covid in the earlier data)
obs1 <- read.csv(paste0("/home/bridge/yl17544/papers_repos/UKHSA_monitoring/data/ONS_HadUK-Grid_1981_2020_noleap_",reg1,".csv"),
                 colClasses=c("Date",NA,NA,NA,NA,NA,"NULL","NULL","NULL","NULL",NA,"NULL"))
#obs1 <- read.csv(paste0("/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/data/ONS_HadUK-Grid_25km_LSOApopweighted_1981_2019_noleap_",reg1,".csv"),colClasses=c("Date",NA,NA,NA,NA,NA,"NULL","NULL","NULL","NULL",NA,"NULL"))
obs1 <- subset(obs1, date<="2019-12-31")
# 2020-2022 (removed any deaths with covid on the certificate; ONS data)
obs2 <- read.csv(paste0("/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/data/ONS_HadUK-Grid_2020_2022_without_covidcert_noleap_",reg1,".csv"),colClasses=c("Date",NA,NA,NA,NA,NA,NA,"NULL"))
#obs2 <- read.csv(paste0("/home/bridge/yl17544/papers_repos/EnglandWales_climate_vs_covid_paper/data/ONS_HadUK-Grid_25km_LSOApopweighted_2020_2022_without_covidcert_noleap_",reg1,".csv"), colClasses=c("Date",NA,NA,NA,NA,NA,NA,"NULL"))
obs2 <- subset(obs2, date>="2020-01-01")
# combine obs in time
obs <- rbind(obs1, obs2)
obs$doy <- as.numeric(strftime(obs$date, format="%j"))

# DEFINITION OF THE CROSS-BASIS FOR TEMPERATURE
# - SPECIFICATION PARAMETERS OF THE EXPOSURE-RESPONSE DIMENSION OF THE CROSS-BASIS
# This bit is same as Ana's 01EstimationERassociation.R
# source: https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata/blob/master/01EstimationERassociation.R
  
# argvar: main model, cubic natural spline with three internal knots in 
#   the 10th, 75th, 90th percentiles of the year-round temperature distribution
argvar <- list(fun="ns", knots=quantile(obs$tmean,c(10,75,90)/100, na.rm=T),
  Bound=range(obs$tmean,na.rm=T)) 

# - SPECIFICATION PARAMETERS OF THE LAG-ASSOCIATION DIMENSION OF THE CROSS-BASIS
# Definition of the maximum lag, that is, 21 days
maxlag <- 21
# arglag: main model, it fits a cubic natural spline with three internal knots 
#   equally-spaced in the log-scale.
arglag <- list(fun="ns",knots=logknots(maxlag,nk=3))

# - CREATE CROSSBASIS OBJECTS
cb <- crossbasis(obs$tmean,lag=maxlag,argvar,arglag=arglag)  

# FIT THE MODEL
# Include in the model the crossbasis term of temperature, along with the 
#    indicator for day of the week (dow) and natural cubic spline of 
#    time with 8 df per year. 
m <- glm(death_allage ~ cb + dow + ns(date,df=round(8*length(date)/365)), 
  data=obs, family=quasipoisson)  

# - DEFINE PROVVISIONAL CENTERING POINT TO HAVE THE INITIAL PREDICTION
varcen <- 19 

# constrain centering point to certain percentiles
cenlim <- round(quantile(obs$tmean, c(0.02,0.98), na.rm=T),1)

# - ESTIMATE MMT FROM THE PREDICTED EXPOSURE-RESPONSE ASSOCIATION 
# MMT corresponds to the temperature of minimum mortality, which will be used as
#    as reference to estimate relative risks and as temperature threshold 
#    to differentiate the contribution of heat and cold to the total mortality 
#    attributable to non-optimal temperatures.
cp <- crosspred(cb,m,cen=varcen,by=0.1,from=cenlim[1],to=cenlim[2])
cen <- cp$predvar[which.min(cp$allRRfit)] 

# - RE-CENTER & GET PREDICTIONS FOR EACH MODEL CENTERING ON THE MMT 
pred <- crosspred(cb, m, cen=cen, by=1)   

##### plotting code
# x limits on plot, rounded to the nearest 1 deg C of obs
xmin <- round_any(min(obs$tmean), 1, f=ceiling)
xmax <- round_any(max(obs$tmean), 1, f=ceiling)
xlab <- expression(paste("Temperature (",degree,"C)"))

# plot settings
pdf(paste0("/home/bridge/yl17544/plots/EnglandWales_climate_vs_covid/ER1981-2022_yearround_21dayslag_MMT2to98_ONSdata_",reg1,".pdf"),
    height=6,width=6,pointsize=14)
# Trim off excess margin space (bottom, left, top, right)
par(mar=c(4.2, 4, 1.8, 0.8))

# actual plotting
plot(pred,"overall",type="n",ylim=c(0.9,2.5),xlim=c(xmin,xmax),axes=T,lab=c(6,5,7),xlab=xlab,
     ylab="RR",main=c(paste0(reg," 1981-2022"), paste0("MMT=",round(cen,digits=1),"°C")))
# red for hot
ind2 <- pred$predvar>=xmin & pred$predvar<=cen
ind3 <- pred$predvar>=cen & pred$predvar<=xmax
lines(pred$predvar[ind2],pred$allRRfit[ind2],col="blue",lwd=1.5)
lines(pred$predvar[ind3],pred$allRRfit[ind3],col="red",lwd=1.5)

# MMT line
abline(v=cen,lty=3)

dev.off()
cat("Saved ER curve!\n")

if(FALSE) {
################################################################################
# 02 CLIMATOLOGICALLY EXPECTED MORTALITY SERIES
################################################################################

# EXPECTED MORTALITY SERIES:
# It is computed as the average mortality for each day of the year 
#   from daily observed deaths, then repeated along the same projection period
#   of the modelled temperature series.

# first, year-round average deaths per day
deathdoy <- tapply(obs$death_allage,as.numeric(format(obs$date,"%j")),
  mean,na.rm=T)[seq(365)]
while(any(isna <- is.na(deathdoy)))
  deathdoy[isna] <- rowMeans(Lag(deathdoy,c(-1,1)),na.rm=T)[isna]
deathexp <- rep(deathdoy,length=nrow(obs))
  
################################################################################
# 03 BIAS CORRECTION - NOT NEEDED
################################################################################

################################################################################
# 04 EXTRAPOLATION OF THE EXPOSURE-RESPONSE CURVE
# 05 PROJECTION & QUANTIFICATION OF THE IMPACT
# 06 ENSEMBLE ESTIMATES & QUANTIFICATION OF THE UNCERTAINTY
################################################################################

# The three last steps of the analysis (extrapolation of the curve, 
#   impact projections and quantification of the uncertainty) can be performed 
#   sequentially using the following code.

# In brief, once we extrapolate the curve, we estimate the daily number of attributable 
#   deaths (AN). 
# By dividing between the total mortality, we estimate the corresponding 
#   attributable fractions (AFs).
# Uncertainty of the estimated impacts is expressed in terms of empirical 
#   confidence intervals, defined as the 2.5th and 97.5th percentiles of the 
#   empirical distribution of the impacts across coefficients samples. 
#   The distribution is obtained through Monte Carlo simulation of the coefficients.

# EXTRACT COEFFICIENTS AND VCOV FROM THE MODEL IN STAGE 1
# With the crossreduce function we can reduce the fit of the bidimensional DLNM
#   (of the stage 1 using the observed temperature-mortality series)
#   to summaries defined in the exposure-response dimension. We can then 
#   extract the coefficients and the covariance matrix defining the overall 
#   exposure-response association.

red <- crossreduce(cb,m,cen=cen)
coef <- coef(red)
vcov <- vcov(red)

# (1) DIMENSION - RANGE OF TEMPERATURES
temprange <- c("tot","cold","heat")

# (2) NUMBER OF ITERATION IN THE MONTE-CARLO SIMULATION 
nsim <- 100

# DEFINE THE DATAFRAME
ansim <- data.frame(matrix(NA, nrow=length(obs$date), ncol=nsim+4,
  dimnames=list(NULL,c("date","tmean","mmt","est",paste0("sim",seq(nsim))))))
        
# STORE DATES, TMEAN AND MMT 
ansim[,"date"] <- obs$date
ansim[,"tmean"] <- obs$tmean
ansim[,"mmt"] <- cen

# (4) EXTRAPOLATION OF THE CURVE: 
# - DERIVE THE CENTERED BASIS USING THE PROJECTED TEMPERATURE SERIES
#   AND EXTRACT PARAMETERS
bvar <- do.call(onebasis,c(list(x=obs$tmean),argvar))
cenvec <- do.call(onebasis,c(list(x=cen),argvar))
bvarcen <- scale(bvar,center=cenvec,scale=F)
    
# (5) IMPACT PROJECTIONS:
# - COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
an <- (1-exp(-bvarcen%*%coef(red)))*deathexp
    
# - STORE AN IN ARRAY BEFORE THE ITERATIONS
ansim[,"est"] <- an
    
# (6) ESTIMATE UNCERTAINTY OF THE PROJECTED AN:
# - SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
set.seed(13041975)
coefsim <- mvrnorm(nsim,coef,vcov)
    
# - LOOP ACROSS ITERATIONS
for(s in seq(nsim)) {
      
  # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
  an <- (1-exp(-bvarcen%*%coefsim[s,]))*deathexp
      
  # STORE THE ATTRIBUTABLE MORTALITY
  ansim[,s+4] <- an

}

################################################################################
# SAVE OUTPUT
################################################################################
options(digits = 3)
options(width = 5000)

fname <- paste0("outputs/daily_attributable_deaths_ER1981-2022_yearround_21dayslag_MMT2to98_nsim",nsim,
  "_1981-2022_ONSdata_",reg1,".csv")
write.csv(ansim,fname,row.names=FALSE) 
cat("Saved output \n")
}
