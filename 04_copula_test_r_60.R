# Preparing R -------------------------------------------------------------
rm(list = ls())

#setwd("/Users/konstantingobler/Desktop/Lebenseinkommen Projekt/New/DIW")

# Packages
source("packages.R")
source("functions.R")



### We should think about only considering nonzero obs in both variables

mydata=read.csv(file="test_60er.csv",head=TRUE,sep=";")

mydata <- select(mydata, one_of(c("schnittek_einzel_32", "par_inc_einzel")))

## check these family and parameter values in the copula package

# pseudo obs
var_a <- pobs(mydata)[,2]
var_b <- pobs(mydata)[,1]

selectedCopula60 <- BiCopSelect(var_a, var_b, familyset = NA)
selectedCopula60
selectedCopula60$family
selectedCopula60$par

# double check!
survgumbel <- surGumbelCopula(param = 1)
set.seed(500)
m <- pobs(as.matrix(mydata.2))
fit <- fitCopula(survgumbel,m,method="mpl")
coef(fit)
# muy buen!


# plotting it it looks like this:
persp(surGumbelCopula(1.17), dCopula)
#sampling from it we can do this easily:
u <- rCopula(3965,surGumbelCopula(1.17))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
pairs.panels(u)
pairs.panels(m)

mydata.1 <- mydata %>% 
  mutate(kidsincome = ifelse(schnittek_einzel_32==0,1, schnittek_einzel_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_einzel==0,1,par_inc_einzel))

###### parents income #######

descdist(mydata.1$parentsincome, discrete=FALSE, boot=5000)
# could again be Lognormal Weibull or Gamma
fit2_lognormal <- fitdist(mydata.1$parentsincome, "lnorm")
meanlog_a <- fit2_lognormal$estimate[1]
sdlog_a <- fit2_lognormal$estimate[2]

fit2_gamma <- fitdist(mydata.1$parentsincome, "gamma", method = "mme")
summary(fit2_gamma)
shape_a <- fit2_gamma$estimate[1]
rate_a <- fit2_gamma$estimate[2]

fit2_weibull <- fitdist(mydata.1$parentsincome, "weibull")
summary(fit2_weibull)
shapeweibull_a <- fit2_weibull$estimate[1]
scaleweibull_a <- fit2_weibull$estimate[2]

hist(mydata.1$parentsincome,breaks=80,main='parents income',density=30,col='cyan',freq=F)
lines(seq(1,141637,80),dlnorm(seq(1,141637,80),meanlog_a,sdlog_a),col='red',lwd=2)
lines(seq(1,141637,80),dgamma(seq(1,141637,80),shape_a, rate_a), col="green", lwd=2)
lines(seq(1,141637,80),dweibull(seq(1,141637,80), shapeweibull_a, scaleweibull_a), col="orange", lwd=2)
legend('topright',c('Fitted lognormal',"Fitted Gamma", "Fitted Weibull"),col=c('red',"green", "orange"),lwd=2)

# here choose gamma! 

####### kids income ######

descdist(mydata.1$kidsincome, discrete=FALSE, boot=5000)
# again looks more like a gamma mixture

fit2_lognormalkids <- fitdist(mydata.1$kidsincome, "lnorm")
meanlog_b <- fit2_lognormalkids$estimate[1]
sdlog_b <- fit2_lognormalkids$estimate[2]

fit2_gammakids <- fitdist(mydata.1$kidsincome, "gamma", method = "mme")
summary(fit2_gammakids)
shape_b <- fit2_gammakids$estimate[1]
rate_b <- fit2_gammakids$estimate[2]

fit2_weibullkids <- fitdist(mydata.1$kidsincome, "weibull")
summary(fit2_weibullkids)
shapeweibull_b <- fit2_weibullkids$estimate[1]
scaleweibull_b <- fit2_weibullkids$estimate[2]

hist(mydata.1$kidsincome,breaks=80,main='kids income',density=30,col='cyan',freq=F)
lines(seq(1,309455,80),dlnorm(seq(1,309455,80),meanlog_a,sdlog_a),col='red',lwd=2)
lines(seq(1,309455,80),dgamma(seq(1,309455,80),shape_a, rate_a), col="green", lwd=2)
lines(seq(1,309455,80),dweibull(seq(1,309455,80), shapeweibull_a, scaleweibull_a), col="orange", lwd=2)
legend('topright',c('Fitted lognormal',"Fitted Gamma", "Fitted Weibull"),col=c('red',"green", "orange"),lwd=2)

## here choose weibull!


### thus Copula: Survival Gumbel
### Marginal Parents income: Gamma
### Marginal Kids income: Weibull

my_dist_60 <- mvdc(surGumbelCopula(1.17), margins = c("gamma","weibull"), paramMargins = list(list(shape = shape_a, rate = rate_a), list(shape = shapeweibull_b, scale = scaleweibull_b)))

v <- rMvdc(50000, my_dist_60)
#write.csv(v, "v.csv")

# Compute the density
pdf_mvd <- dMvdc(v, my_dist_60)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist_60)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
#setting the range over childs income and parents income 
persp(my_dist_60, dMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Density")
contour(my_dist_60, dMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Contour plot")
persp(my_dist_60, pMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "CDF")
contour(my_dist_60, pMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Contour plot")

# Plot the data for a visual comparison
plot(mydata$par_inc_einzel, mydata$schnittek_einzel_32, main = 'Test dataset x and y', col = "blue")
points(v[,1], v[,2], col = '#FF000020', cex=0.5)
legend('bottomright', c('Observed', 'Simulated'), col = c('blue', 'red'), pch=21)

cor(mydata.1, method = "kendall")
cor(mydata.1, method = "spearman")

cor(v, method = "kendall")
cor(v, method = "spearman")

# looks good... 



