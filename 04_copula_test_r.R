# Preparing R -------------------------------------------------------------
rm(list = ls())

#setwd("/Users/konstantingobler/Desktop/Lebenseinkommen Projekt/New/DIW")

# Packages
source("packages.R")

### We should think about only considering nonzero obs in both variables

mydata=read.csv(file="test.csv",head=TRUE,sep=";")

#remove negative values in kids income
mydata <- mydata %>% 
  mutate(schnittek_einzel_32= ifelse(schnittek_einzel_32<0,0,schnittek_einzel_32))  

min(mydata$schnittek_einzel_32)

# fitting pseudo observations because copulas only accept values in the 
# unit interval

var_a <- pobs(mydata)[,2]
var_b <- pobs(mydata)[,1]


selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
selectedCopula
selectedCopula$family
selectedCopula$par

## checking some goodness of fit statistics
## here by basically doing the same thing backwards to see whether the parameter
## estimates are robust 
## BB7 copula is a special case of the typically single-parameter copulas
## subsumed by the so called Archimedean family: two-parameter Archimedean
## copulas such as the Joe-Clayton (BB7) copula's more flexible structure
## allows for different non-zero lower and upper tail dependence coefficien

# note that we can only test the first parameter:  
## building subsample
mydata.2 <- select(mydata, one_of(c("schnittek_einzel_32", "par_inc_einzel")))

BB7.cop <- joeBiCopula(param = 2)
set.seed(500)
m <- pobs(as.matrix(mydata.2))
fit <- fitCopula(BB7.cop,m,method="mpl")
coef(fit)
# 1.1 corresponds to the one we retrieved from BiCopSelect()
# muy bien

###################

# hence we now know that the copula (the joint rank transformed distribution)
# is the BB7 with parameters p1=1.05 and p2=0.31

# plotting it it looks like this:
persp(BB7Copula(c(1.05, 0.31)), dCopula)
#sampling from it we can do this easily:
u <- rCopula(3965,BB7Copula(c(1.05, 0.31)))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
pairs.panels(u)

# looks a lot like independece which is totally fine since the spearmanR is 
# only 0.22

#and compared to original data (as inteded the correlation structure prevails)
cor(mydata.2,method='spearman')
pairs.panels(mydata.2)
#compared to our true data (transformed to unit interval):
cor(m,method='spearman')
pairs.panels(m)

# it can also be nice to check the scatterplot fitting a linear line
plot(mydata.2$par_inc_einzel, mydata.2$schnittek_einzel_32, pch='.')
abline(lm(mydata.2$par_inc_einzel~mydata.2$schnittek_einzel_32),col='red',lwd=1)

## before generating the multivariate distribution using our copula BB7
## we should put some thought into our marginal distributions of par_inc_einzel
## and schnittek_einzel_32

# this is actually the crucial part:
#remove zeros in order to be able to check log- based distributions
mydata.3 <- mydata.2 %>% 
  mutate(kidsincome = ifelse(schnittek_einzel_32==0,1, schnittek_einzel_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_einzel==0,1,par_inc_einzel))

# evaluate marginal distribution via descdist
# start with parentsincome
descdist(mydata.3$parentsincome, discrete=FALSE, boot=5000)
#summary statistics
#------
#  min:  0   max:  362951 
#median:  39626 
#mean:  44014 
#estimated sd:  29676 
#estimated skewness:  2.7 
#estimated kurtosis:  19 

# this looks almost surely like a lognormal! 

fit2_lognormal <- fitdist(mydata.3$parentsincome, "lnorm")
summary(fit2_lognormal)
meanlog_a <- fit2_lognormal$estimate[1]
sdlog_a <- fit2_lognormal$estimate[2]

# also plot against Exponential

fit2_exp <- fitdist(mydata.3$parentsincome, "exp", method="mme")
summary(fit2_exp)
exp.rate_a <- fit2_exp$estimate[1]

# let us compare data histogram with the two densities:

hist(mydata.3$parentsincome,breaks=80,main='parents incomes',freq=F,density=30,col='cyan') 
lines(seq(1,362951,80),dlnorm(seq(1,362951,80),meanlog_a,sdlog_a),col='red',lwd=2)
lines(seq(1,362951,80),dexp(seq(1,362951,80), exp.rate_a), col="green", lwd=2)
legend('topright',c('Fitted lognormal',"Fitted Exponential"),col=c('red',"green"),lwd=2)

##### now kidsincome ####

descdist(mydata.3$kidsincome, discrete=FALSE, boot=5000)

#summary statistics
#------
#  min:  0   max:  309455 
#median:  37767 
#mean:  41401 
#estimated sd:  28449 
#estimated skewness:  1.4 
#estimated kurtosis:  8.6 

# this could be either gamma or lognormal or Weibull but the fit will not be very accurate
fit2_lognormalkids <- fitdist(mydata.3$kidsincome, "lnorm")
meanlog_b <- fit2_lognormalkids$estimate[1]
sdlog_b <- fit2_lognormalkids$estimate[2]

fit2_gammakids <- fitdist(mydata.3$kidsincome, "gamma", method = "mme")
summary(fit2_gammakids)
shape_b <- fit2_gammakids$estimate[1]
rate_b <- fit2_gammakids$estimate[2]

fit2_weibullkids <- fitdist(mydata.3$kidsincome, "weibull")
summary(fit2_weibullkids)
shapeweibull_b <- fit2_weibullkids$estimate[1]
scaleweibull_b <- fit2_weibullkids$estimate[2]


hist(mydata.3$kidsincome,breaks=80,main='childrens income',density=30,col='cyan',freq=F)
lines(seq(1,309455,80),dlnorm(seq(1,309455,80),meanlog_b,sdlog_b),col='red',lwd=2)
lines(seq(1,309455,80),dgamma(seq(1,309455,80),shape_b, rate_b), col="green", lwd=2)
lines(seq(1,309455,80),dweibull(seq(1,309455,80), shapeweibull_b, scaleweibull_b), col="orange", lwd=2)
legend('topright',c('Fitted lognormal',"Fitted Gamma", "Fitted Weibull"),col=c('red',"green", "orange"),lwd=2)

#we could try Weibull or Gamma


# Generate the multivariate distribution 
# here margins refers to the marginal distribution of the input variables which we choose 
# in a way that they represent the distributional structure of the marginals present

#my_dist <- mvdc(BB7Copula(param = c(1.05,0.31)), margins = c("lnorm","gamma"),   paramMargins = list(list(meanlog = meanlog_a, sdlog = sdlog_a), list(shape = shape_b, rate = rate_b)))
my_dist <- mvdc(BB7Copula(param = c(1.05,0.31)), margins = c("lnorm","weibull"), paramMargins = list(list(meanlog = meanlog_a, sdlog = sdlog_a), list(shape = shapeweibull_b, scale = scaleweibull_b)))

v <- rMvdc(50000, my_dist)
#write.csv(v, "v.csv")

# Compute the density
pdf_mvd <- dMvdc(v, my_dist)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
#setting the range over childs income and parents income 
persp(my_dist, dMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Density")
contour(my_dist, dMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Contour plot")
persp(my_dist, pMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "CDF")
contour(my_dist, pMvdc, xlim = c(0, 50000), ylim=c(0, 50000), main = "Contour plot")

# Plot the data for a visual comparison
plot(mydata$par_inc_einzel, mydata$schnittek_einzel_32, main = 'Test dataset x and y', col = "blue")
points(v[,1], v[,2], col = 'red')
legend('bottomright', c('Observed', 'Simulated'), col = c('blue', 'red'), pch=21)

cor(mydata.2, method = "kendall")
cor(mydata.2, method = "spearman")

cor(v, method = "kendall")
cor(v, method = "spearman")

# looks good... 

# we could play around with gamma and weibull..
pairs.panels(mydata.2)
pairs.panels(v)


#ggplot(data = mydata, mapping = aes(x = var_a, y = var_b)) +  geom_point() + geom_smooth()


