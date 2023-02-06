#################################################################
# Unified Code for FinStats Final Project    ####################
#################################################################

library(moments)
library(quantmod)
library(readxl)
library(quadprog)
library(PerformanceAnalytics)
library(boot)
library(sjstats)
library(fitdistrplus)
library(fGarch)
library(corrplot)
library(tseries)
library(copula)
library(ggfortify)
library(ggplot2)
library(xtable)

data = read.csv("FinStat_Datafile.csv", header = T)
n = dim(data)[1]


AAPL_ret = data$AAPL[-1] / data$AAPL[-n] - 1
GOOGL_ret = data$GOOGL[-1] / data$GOOGL[-n] - 1
AMZN_ret = data$AMZN[-1] / data$AMZN[-n] - 1
FB_ret = data$FB[-1] / data$FB[-n] - 1
NFLX_ret = data$NFLX[-1] / data$NFLX[-n] - 1
MSFT_ret = data$MSFT[-1] / data$MSFT[-n] - 1
JPM_ret = data$JPM[-1] / data$JPM[-n] - 1
BAC_ret = data$BAC[-1] / data$BAC[-n] - 1
GS_ret = data$GS[-1] / data$GS[-n] - 1
CSGN_ret = data$CSGN[-1] / data$CSGN[-n] - 1
TSLA_ret = data$TSLA[-1] / data$TSLA[-n] - 1
GM_ret = data$GM[-1] / data$GM[-n] - 1
PG_ret = data$PG[-1] / data$PG[-n] - 1
WMT_ret = data$WMT[-1] / data$WMT[-n] - 1
TGT_ret = data$TGT[-1] / data$TGT[-n] - 1
#Log Returns
AAPL_logret = log(data$AAPL[-1] / data$AAPL[-n])
GOOGL_logret = log(data$GOOGL[-1] / data$GOOGL[-n])
AMZN_logret = log(data$AMZN[-1] / data$AMZN[-n])
FB_logret = log(data$FB[-1] / data$FB[-n])
NFLX_logret = log(data$NFLX[-1] / data$NFLX[-n])
MSFT_logret = log(data$MSFT[-1] / data$MSFT[-n])
JPM_logret = log(data$JPM[-1] / data$JPM[-n])
BAC_logret = log(data$BAC[-1] / data$BAC[-n])
GS_logret = log(data$GS[-1] / data$GS[-n])
CSGN_logret = log(data$CSGN[-1] / data$CSGN[-n])
TSLA_logret = log(data$TSLA[-1] / data$TSLA[-n])
GM_logret = log(data$GM[-1] / data$GM[-n])
PG_logret = log(data$PG[-1] / data$PG[-n])
WMT_logret = log(data$WMT[-1] / data$WMT[-n])
TGT_logret = log(data$TGT[-1] / data$TGT[-n])


#Descriptive Statistics
#Mean 
mm = AAPL_logret
mean(mm)

#Skewness
skew = AAPL_logret
skewness(skew)

#Kurtosis
kur = AAPL_logret
kurtosis(kur)

#Shapiro-Wilk Normality Test (low p-value -> not normal)
shap = AAPL_logret
shapiro.test(shap)

#Cullen - Frey graph 
#Change the variable g for different assets
g = TGT_logret
graph = descdist(g, boot = 100)
plotdist(g, histo = TRUE, demp = TRUE)

#  QQ Normal Plot (Graphical Normality Check)
qq = AAPL_logret
qqnormPlot(qq)

# Betas for each asset
aux <- data.frame(AAPL_ret = data$AAPL[-1] / data$AAPL[-n] - 1,
                  GOOGL_ret = data$GOOGL[-1] / data$GOOGL[-n] - 1,
                  AMZN_ret = data$AMZN[-1] / data$AMZN[-n] - 1,
                  FB_ret = data$FB[-1] / data$FB[-n] - 1,
                  NFLX_ret = data$NFLX[-1] / data$NFLX[-n] - 1,
                  MSFT_ret = data$MSFT[-1] / data$MSFT[-n] - 1,
                  JPM_ret = data$JPM[-1] / data$JPM[-n] - 1,
                  BAC_ret = data$BAC[-1] / data$BAC[-n] - 1,
                  GS_ret = data$GS[-1] / data$GS[-n] - 1,
                  CSGN_ret = data$CSGN[-1] / data$CSGN[-n] - 1,
                  TSLA_ret = data$TSLA[-1] / data$TSLA[-n] - 1,
                  GM_ret = data$GM[-1] / data$GM[-n] - 1,
                  PG_ret = data$PG[-1] / data$PG[-n] - 1,
                  WMT_ret = data$WMT[-1] / data$WMT[-n] - 1,
                  TGT_ret = data$TGT[-1] / data$TGT[-n] - 1,
                  SP500 = data$SP500[-1] / data$SP500[-n] - 1)
trunc(100000*sapply(aux[, -ncol(aux)], mean))/100000
trunc(100000*sapply(aux[, -ncol(aux)], sd))/100000
sapply(aux[, -ncol(aux)], function(x) cor(x, aux$SP500)*(sd(x)/sd(aux$SP500)) )


x = g
n = length(x)
start = c(mean(x),sd(x),4)
loglik_t = function(beta) sum(-dt((x - beta[1])/beta[2],
                                  beta[3],log = TRUE) + log(beta[2]))
fit_t = optim(start,loglik_t,hessian = T, method = "L-BFGS-B",lower = c(-1,0.001,1))
fit_t$par
AIC_t = 2*fit_t$value + 2*3
BIC_t = 2 * fit_t$value + log(n) * 3
sd_t = sqrt(diag(solve(fit_sstd$hessian)))

#Fit Skewed T Distribution
y = g
n = length(y)
start_sstd = c(mean(y),sd(y),5,1)
loglik_sstd = function(beta) sum(-dsstd(y, mean =  beta[1],sd = beta[2],
                                        nu = beta[3],xi = beta[4],log = TRUE)) 
fit_sstd = optim(start_st,loglik_sstd,hessian = T, method = "L-BFGS-B",lower = c(0.1,0.01,2.1,-2))
fit_sstd$par
AIC_sstd = 2*fit_t$value + 2*4
BIC_sstd =2*fit_sstd$value + log(n) *4
sd_sstd = sqrt(diag(solve(fit_sstd$hessian)))

#Fit Normal Distribution
z = g
fit_norm = fitdist(z,"norm")
AIC_norm = fit_norm$aic
BIC_norm = fit_norm$bic

#Sharpe Ratio (Annualised Mean and SD used)
#Change variable s for different assets
s = TGT_logret
mean_ret = mean(s)*12
rf = 0.0045
sd_ret = sd(s)*sqrt(12)
Sharpe_Ratio = (mean_ret - rf) /sd_ret
Sharpe_Ratio

#Covariance Analysis
#Using package "corrplot"
mat = as.matrix(data[-1])
pairs(mat)
cov(mat)
cor(mat)
head(data[-1])
M = cor(data[-1])
head(round(M,2))
library(corrplot)
corrplot(M, method = "circle")

# Test for Stationarity 
# Graphical Approach : ACF (check if the acf decays to zero quickly : stationary)
library(tseries)
ss = AAPL_logret
acf(ss)
# Stationarity test : Using Hypothesis testing : Augmentented Dickety Fuller (ADF) test 
# low p- value (< 0.05) implies stationarity 
# p-value close to 0.05 will imply seasonality
adf.test(ss)

#Copulas fitted using Semi Parametric Pseudo-Maximum Liklihood
n = dim(data)[1]
cop_data = cbind(rank(AAPL_logret)/(n+1),rank(GOOGL_logret)/(n+1),rank(AMZN_logret)/(n+1),
                 rank(FB_logret)/(n+1),rank(NFLX_logret)/(n+1),rank(MSFT_logret)/(n+1), rank(JPM_logret)/(n+1)
                 , rank(BAC_logret)/(n+1), rank(GS_logret)/(n+1), rank(CSGN_logret)/(n+1),rank(TSLA_logret)/(n+1)
                 , rank(GM_logret)/(n+1), rank(PG_logret)/(n+1), rank(WMT_logret)/(n+1), rank(TGT_logret)/(n+1))

library(fGarch)
library(copula)
#Normal Copula
fnorm = fitCopula(copula = normalCopula(dim=15),data = cop_data,method = "ml")
AIC_fnorm = -2*(fnorm@loglik) + 2*length(fnorm@estimate)
AIC_fnorm

#T Copula
tCop = fitCopula(copula = tCopula(dim=15),data = cop_data,method = "ml")
AIC_tCop = -2*(tCop@loglik) + 2*length(tCop@estimate)
AIC_tCop

#Frank Copula
fr = fitCopula(copula = frankCopula(3,dim = 15),data = cop_data,method = "ml")
AIC_fr = -2*(fr@loglik) + 2*length(fr@estimate)
AIC_fr

#Clayton Copula
clay = fitCopula(copula = claytonCopula(1,dim = 15),data = cop_data,method = "ml")
AIC_clay = -2*(clay@loglik) + 2*length(clay@estimate)
AIC_clay

#Gumbel Copula
gumb = (fitCopula(copula = gumbelCopula(3,dim = 15),data = cop_data,method = "ml"))
AIC_gumb = -2*(gumb@loglik) + 2*length(gumb@estimate)
AIC_gumb





# VAR Calculation
S0=100000
alpha=0.05
var_data = TGT_logret
q = as.numeric(quantile(var_data, alpha))
VaR_nonp = -S0 * q
VaR_nonp
IEVaR = (var_data < q)
ES_nonp = -S0 * sum(var_data * IEVaR) / sum(IEVaR)
ES_nonp






table_summary <- read.csv("Table1.csv")

<<results=tex>>
xtable(table_summary)



#################################################
# Portolio Theory        ########################
#################################################

data_val = data[,!(names(data)%in%c('ï..Date'))]
returns_all = 100*(data.frame(diff(as.matrix(data_val))))/(head(data_val,-1))
returns = returns_all[,!(names(returns_all)%in%c('SP500'))]
mean_vect = colMeans(returns)
cov_mat = cov(returns)
sd_vect = sqrt(diag(cov_mat))


M = length(mean_vect)
muP = seq(.02,8,length=300)
sdP = muP
weights = matrix(0,nrow=300,ncol=M)
Amat = cbind(rep(1,M),mean_vect)

for (i in 1:length(muP))
{
  result =
    solve.QP(Dmat=cov_mat,dvec=rep(0,M), Amat=Amat,
             c(1,muP[i]), meq=2)
  sdP[i] = sqrt(2*result$value)
  weights[i,] = result$solution
}
plot(sdP,muP,type="l",xlim=c(-1,12),ylim=c(0,10))

mvpind <- sdP==min(sdP)
mvp_weights <- weights[mvpind]
muP[mvpind]

points(sdP[mvpind],muP[mvpind],cex=1,pch="*",col='blue')
text(sdP[mvpind],muP[mvpind],sprintf("(%.2f,%.2f)",sdP[mvpind],muP[mvpind]),cex=1,pos=1)

#VaR for MVP
investment = 100000
mvp_weights <- round(mvp_weights,5)
mvp_return <- rowSums(as.matrix(returns)%*%diag(mvp_weights))
q = quantile(0.01*mvp_return,0.05)
VaR_mvp = q*-investment # Value at risk
IEVaR = ((0.01*mvp_return)<q)
ESMVP = -investment*sum(0.01*mvp_return*IEVaR)/sum(IEVaR)

mufree = 0.45
points(0,mufree,cex=3,col="blue",pch="*")
sharpe =( muP-mufree)/sdP
tan_ind = (sharpe == max(sharpe))
lines(c(0,sdP[tan_ind]),c(mufree,muP[tan_ind]),col="red",lwd=3)
points(sdP[tan_ind],muP[tan_ind],col="blue",cex=3,pch="*")
ef_ind = (muP > muP[mvpind])
lines(sdP[ef_ind],muP[ef_ind],type="l",xlim=c(0,.25),ylim=c(0,.3),col="cyan",lwd=3)

asset_names <- colnames(returns)
for(i in 1:15) {
  text(sd_vect[i],mean_vect[i],asset_names[i], cex=0.5)
}

summary_mvp = data.frame(Metric = c(colnames(returns), "Anualized Mean", "Anualized Volatility", "VaR(0.5)", "ES(0.5)"),
                          Value = c(mvp_weights, 12*muP[mvpind], sqrt(12)*sdP[mvpind], VaR_mvp,  ESMVP ) )
<<results=tex>>
xtable(t(summary_mvp))
summary_tan = data.frame(Metric = c(colnames(returns), "Anualized Mean", "Anualized Volatility", "Sharpe"),
                         Value = c(weights[tan_ind,], 12*muP[tan_ind], sqrt(12)*sdP[tan_ind], sharpe[tan_ind] ) )
<<results=tex>>
xtable(t(summary_tan))


################################
#Short sells are not allowed  ##
################################

M = length(mean_vect)
muP = seq(.02,7,length=300)
sdP = muP
weights = matrix(0,nrow=300,ncol=M)
Amat = cbind(rep(1,M),mean_vect,diag(1,nrow=M),-diag(1,nrow=M))

for (i in 1:length(muP))
{
  result =
    solve.QP(Dmat=cov_mat,dvec=rep(0,M), Amat=Amat,
             c(1,muP[i],rep(0,M),rep(-8,M)), meq=2)
  sdP[i] = sqrt(2*result$value)
  weights[i,] = result$solution
}

plot(sdP,muP,type="l",xlim=c(0,20),ylim=c(0,10))

mvpind <- sdP==min(sdP)
mvp_weights <- weights[mvpind]
muP[mvpind]

points(sdP[mvpind],muP[mvpind],cex=3,pch="*",col='blue')
text(sdP[mvpind],muP[mvpind],sprintf("(%.2f,%.2f)",sdP[mvpind],muP[mvpind]),cex=1,pos=1)

#VaR for MVP
investment = 100000
mvp_weights <- round(mvp_weights,5)
mvp_return <- rowSums(as.matrix(returns)%*%diag(mvp_weights))
q = quantile(0.01*mvp_return,0.05)
VaR_mvp = q*-investment # Value at risk
IEVaR = ((0.01*mvp_return)<q)
ESMVP = -investment*sum(0.01*mvp_return*IEVaR)/sum(IEVaR)

mufree = 0.45
points(0,mufree,cex=3,col="blue",pch="*")
sharpe =( muP-mufree)/sdP
tan_ind = (sharpe == max(sharpe))
lines(c(0,sdP[tan_ind]),c(mufree,muP[tan_ind]),col="red",lwd=3)
points(sdP[tan_ind],muP[tan_ind],col="blue",cex=3,pch="*")
ef_ind = (muP > muP[mvpind])
lines(sdP[ef_ind],muP[ef_ind],type="l",xlim=c(0,.25),ylim=c(0,.3),col="cyan",lwd=3)

asset_names <- colnames(returns)
for(i in 1:15) {
  text(sd_vect[i],mean_vect[i],asset_names[i], cex=0.5)
}

#tangency portfolio
tan_weights = weights[tan_ind]
tan_return = muP[tan_ind]


#beta value of assets
returns_beta_calc = returns_all-mufree
summary(lm(returns_beta_calc$AAPL~returns_beta_calc$SP500))


summary_mvp = data.frame(Metric = c(colnames(returns), "Anualized Mean", "Anualized Volatility", "VaR(0.5)", "ES(0.5)"),
                         Value = c(mvp_weights, 12*muP[mvpind], sqrt(12)*sdP[mvpind], VaR_mvp,  ESMVP ) )
<<results=tex>>
xtable(summary_mvp)
summary_tan = data.frame(Metric = c(colnames(returns), "Anualized Mean", "Anualized Volatility", "Sharpe"),
                         Value = c(weights[tan_ind,], 12*muP[tan_ind], sqrt(12)*sdP[tan_ind], sharpe[tan_ind] ) )
<<results=tex>>
xtable(summary_tan)


#################################################
# Section - 4 Asset allocation     ##############
#################################################

# Identify the mean closest to 0.5:
target_ret_weights = weights[abs(muP-0.5) == min(abs(muP-0.5))]
# this results in 41.1% investment in CSGN and 58.9% in PG, note that this is
#not on the efficient frontier and more return can be obtained for the same risk
sdP[abs(muP-0.5) == min(abs(muP-0.5))] #risk

target_ret_weights <- round(target_ret_weights,5)
target_return <- rowSums(as.matrix(returns)%*%diag(target_ret_weights))
VaR_target = -100000*VaR(0.01*target_return, method="historical")
ES_target = -100000*ES(0.01*target_return, method="historical")

summary_5 = data.frame(Metric = c(colnames(returns), "Mean", "Volatility", "VaR(0.5)", "ES(0.5)"),
                         Value = c(target_ret_weights, 0.5, sdP[abs(muP-0.5) == min(abs(muP-0.5))],
                                   VaR_target[,1],  ES_target[,1] ) )
xtable(summary_5)

# with risk free asset

risk_weight = (0.5-mufree)/(tan_return-mufree) #98.3% invested in riskfree
risky_assets_allocate = risk_weight*tan_weights
risk_target = risk_weight*sdP[tan_ind]
#VaR and ES
risky_assets_allocate <- round(risky_assets_allocate,5)
target_return <- rowSums(as.matrix(returns)%*%diag(risky_assets_allocate))
VaR_target = -100000*VaR(0.01*target_return, method="historical")
ES_target = -100000*ES(0.01*target_return, method="historical")

risky_assets_allocate = round(risky_assets_allocate*1000)/1000
summary_5 = data.frame(Metric = c("Risk free", colnames(returns), "Mean", "Volatility", "VaR(0.5)", "ES(0.5)"),
                       Value = c((1-risk_weight), risky_assets_allocate, 0.5, risk_target,
                                 VaR_target[,1],  ES_target[,1] ) )
xtable(t(summary_5))

#################################################
# Section 6     #################################
#################################################

investment = 100000
returns_var = 0.01*returns


##################################
# parametric method        #######
##################################

# default alpha is 5%
VaRp = -investment*sapply(returns_var,function(x)VaR(x, method="gaussian"))

ESp = -investment*sapply(returns_var,function(x)ES(x, method="gaussian"))

VaRp == max(VaRp) # tesla has the highest
VaRp == min(VaRp) # msft has the least
ESp == max(ESp) # tesla has the highest
ESp == min(ESp) # msft has the least

varnp_fn <- function(d, i){
  d2 <- d[i,]
  result <- -100000*sapply(d2,function(x)VaR(x, method="gaussian"))
  return(result)
}

boot_res_var <- boot(returns_var, varnp_fn, R=100)
boot_sd_var<- NULL
boot_confidence_var <- NULL
temp <- NULL
for(i in 1:ncol(boot_res_var$t)) {      
  boot_sd_var[i] <- sd(boot_res_var$t[ , i])
  temp <- boot_ci(boot_res_var$t[,i],method = "quantile", ci.lvl = 0.95)
  boot_confidence_var$low[i] = temp$conf.low
  boot_confidence_var$VaR[i] = VaRp[i]
  boot_confidence_var$high[i] = temp$conf.high
}

xtable(t(data.frame(Asset=names(VaRnp), boot_confidence_var)))


esnp_fn <- function(d, i){
  d2 <- d[i,]
  result <- -100000*sapply(d2,function(x)ES(x, method="gaussian"))
  return(result)
}

boot_res_es <- boot(returns_var, esnp_fn, R=100)
boot_sd_es<- NULL
boot_confidence_es <- NULL
temp <- NULL
for(i in 1:ncol(boot_res_es$t)) {      
  boot_sd_es[i] <- sd(boot_res_es$t[ , i])
  temp <- boot_ci(boot_res_es$t[,i],method = "quantile", ci.lvl = 0.95)
  boot_confidence_es$low[i] = temp$conf.low
  boot_confidence_es$ES[i] = ESp[i]
  boot_confidence_es$high[i] = temp$conf.high
}

xtable(t(data.frame(Asset=names(ESnp), lapply(boot_confidence_es, round) ) ))


##################################
# non-parametric method    #######
##################################

VaRnp = -investment*sapply(returns_var,function(x)VaR(x, method="historical"))

ESnp = -investment*sapply(returns_var,function(x)ES(x, method="historical"))

VaRnp == max(VaRnp) # tesla has the highest
VaRnp == min(VaRnp) # wmt has the least
ESnp == max(ESnp) # tesla has the highest
ESnp == min(ESnp) # msft has the least

#bootstrap for standard error and confidence interval

varnp_fn <- function(d, i){
  d2 <- d[i,]
  result <- -100000*sapply(d2,function(x)VaR(x, method="historical"))
  return(result)
}

boot_res_var <- boot(returns_var, varnp_fn, R=100)
boot_sd_var<- NULL
boot_confidence_var <- NULL
temp <- NULL
for(i in 1:ncol(boot_res_var$t)) {      
  boot_sd_var[i] <- sd(boot_res_var$t[ , i])
  temp <- boot_ci(boot_res_var$t[,i],method = "quantile", ci.lvl = 0.95)
  boot_confidence_var$low[i] = temp$conf.low
  boot_confidence_var$VaR[i] = VaRnp[i]
  boot_confidence_var$high[i] = temp$conf.high
}

xtable(t( data.frame(Asset=names(VaRnp), lapply(boot_confidence_var, round) ) ))


esnp_fn <- function(d, i){
  d2 <- d[i,]
  result <- -100000*sapply(d2,function(x)ES(x, method="historical"))
  return(result)
}

boot_res_es <- boot(returns_var, esnp_fn, R=100)
boot_sd_es<- NULL
boot_confidence_es <- NULL
temp <- NULL
for(i in 1:ncol(boot_res_es$t)) {      
  boot_sd_es[i] <- sd(boot_res_es$t[ , i])
  temp <- boot_ci(boot_res_es$t[,i],method = "quantile", ci.lvl = 0.95)
  boot_confidence_es$low[i] = temp$conf.low
  boot_confidence_es$ES[i] = ESnp[i]
  boot_confidence_es$high[i] = temp$conf.high
}

xtable(t( data.frame(Asset=names(ESnp), lapply(boot_confidence_es, round) ) ))



#################################################
# plots         #################################
#################################################

# Price
par(mfrow=c(5,3))
for(k in 2:(ncol(data)-1)){
  print(colnames(data)[k])
  plot(ts( data[, k] ), col= "blue", xlab="Day", xaxt='n', main=colnames(data)[k],
       ylab="Price")
  
  indx <- c( 1, floor(nrow(data)/3), floor(2*nrow(data)/3), nrow(data) )
  axis(side=1, at=indx, labels=data[,1][indx], cex=0.7)
  
}
par(mfrow=c(1,1))
plot(ts( data[, ncol(data)] ), col= "blue", xlab="Day", xaxt='n',
     main=colnames(data)[ ncol(data)], lwd=2,
     ylab="Price")
indx <- c( 1, floor(nrow(data)/3), floor(2*nrow(data)/3), nrow(data) )
axis(side=1, at=indx, labels=data[,1][indx], cex=0.7)

# Returns
fechaux = data[-1,1]
par(mfrow=c(5,3))
for(k in 1:(ncol(aux)-1)){
  print(colnames(aux)[k])
  plot(ts( aux[, k] ), col= "blue", xlab="Day", xaxt='n',
       main=gsub("_ret", "", colnames(aux)[k] ),
       ylab="Monthly Return")
  
  indx <- c( 1, floor(nrow(aux)/3), floor(2*nrow(aux)/3), nrow(aux) )
  axis(side=1, at=indx, labels=fechaux[indx], cex=0.7)
  
}
par(mfrow=c(1,1))
plot(ts( aux[, ncol(aux)] ), col= "blue", xlab="Day", xaxt='n',
     main=gsub("_ret", "", colnames(aux)[ ncol(aux)]), lwd=2,
     ylab="Monthly Return")
indx <- c( 1, floor(nrow(aux)/3), floor(2*nrow(aux)/3), nrow(aux) )
axis(side=1, at=indx, labels=fechaux[indx], cex=0.7)

# Equity curves
par(mfrow=c(5,3))
for(k in 1:(ncol(aux)-1)){
  print(colnames(aux)[k])
  plot(ts( cumprod(1+aux[, k]) ), col= "blue", xlab="Day", xaxt='n',
       main=gsub("_ret", "", colnames(aux)[k] ),
       ylab="$1 Equity")
  
  indx <- c( 1, floor(nrow(aux)/3), floor(2*nrow(aux)/3), nrow(aux) )
  axis(side=1, at=indx, labels=fechaux[indx], cex=0.7)
  
}
par(mfrow=c(1,1))
plot(ts( cumprod(1+aux[, ncol(aux)]) ), col= "blue", xlab="Day", xaxt='n',
     main=gsub("_ret", "", colnames(aux)[ ncol(aux)]), lwd=2,
     ylab="$1 Equity")
indx <- c( 1, floor(nrow(aux)/3), floor(2*nrow(aux)/3), nrow(aux) )
axis(side=1, at=indx, labels=fechaux[indx], cex=0.7)
  

par(mfrow=c(5,3))
# Hist
for(k in 1:(ncol(aux)-1)){
  print(colnames(aux)[k])
  hist( aux[, k], xlab="Return", freq=F,
       main=gsub("_ret", "", colnames(aux)[k] ))
  lines(density(aux[, k]), lty=2)
}
# qqplot
par(mfrow=c(5,3))
for(k in 1:(ncol(aux)-1)){
  print(colnames(aux)[k])
  qqnormPlot( aux[, k], title=F, main=gsub("_ret", "", colnames(aux)[k] ))
}
# boxplot
par(mfrow=c(5,3))
for(k in 1:(ncol(aux)-1)){
  print(colnames(aux)[k])
  boxplot( aux[, k], main=gsub("_ret", "", colnames(aux)[k] ))
  
}
# Multiscatter
colnames(aux) = colnames(data)[-1]
pairs(aux[, -1], pch = 19)

