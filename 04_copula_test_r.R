# Preparing R -------------------------------------------------------------
rm(list = ls())

# Packages
source("packages.R")
source("functions.R")

mydata=read.csv(file="test.csv",head=TRUE,sep=";")

#remove negative values in kids income
mydata <- mydata %>% 
  mutate(schnittek_einzel_32= ifelse(schnittek_einzel_32<0,0,schnittek_einzel_32))  

min(mydata$schnittek_einzel_32)

mydata <- select(mydata, one_of(c("schnittek_einzel_32", "par_inc_einzel")))

selectedCopula <- pre.marginals.copula(data = mydata)
par <- selectedCopula[2]
par2 <- selectedCopula[3]


## checking some goodness of fit statistics
## here by basically doing the same thing backwards to see whether the parameter
## estimates are robust 
## BB7 copula is a special case of the typically single-parameter copulas
## subsumed by the so called Archimedean family: two-parameter Archimedean
## copulas such as the Joe-Clayton (BB7) copula's more flexible structure
## allows for different non-zero lower and upper tail dependence coefficien

survgumbel <- surGumbelCopula(param = 1)
set.seed(500)
m <- pobs(as.matrix(mydata))
fit <- fitCopula(survgumbel,m,method="mpl")
coef(fit)
# 1.2 corresponds to the one we retrieved from BiCopSelect()
# muy bien

bb7copula <- BB7Copula(param = c(1,1))
set.seed(500)
fit <- fitCopula(bb7copula,m,method="mpl")
param <- coef(fit)

overview <- BiCopEstList(m[,1], m[,2],rotations = T)

min(overview$summary$AIC) # BB7 copula!
min(overview$summary$BIC) # rotated Gumbel copula! (survGumbel)
max(overview$summary$logLik) #families 8,10


###################

# hence we now know that the copula (the joint rank transformed distribution)
# is the BB7 with parameters p1=1.05 and p2=0.31

# plotting it it looks like this:
persp(surGumbelCopula(par), dCopula)
#sampling from it we can do this easily:
u <- rCopula(3965,surGumbelCopula(par))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
pairs.panels(u)
pairs.panels(m)

# looks a lot like independece which is totally fine since the spearmanR is 
# only 0.22


# it can also be nice to check the scatterplot fitting a linear line
#plot(mydata.2$par_inc_einzel, mydata.2$schnittek_einzel_32, pch='.')
#abline(lm(mydata.2$par_inc_einzel~mydata.2$schnittek_einzel_32),col='red',lwd=1)

## before generating the multivariate distribution using our copula BB7
## we should put some thought into our marginal distributions of par_inc_einzel
## and schnittek_einzel_32

# this is actually the crucial part:
#remove zeros in order to be able to check log- based distributions
mydata.1 <- mydata %>% 
  mutate(kidsincome = ifelse(schnittek_einzel_32==0,1, schnittek_einzel_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_einzel==0,1,par_inc_einzel))


# evaluate marginal distribution via descdist
# start with parentsincome
descdist(mydata.1$parentsincome, discrete=FALSE, boot=5000)
# could again be Lognormal Weibull or Gamma
fit2_lognormal <- fitdist(mydata.1$parentsincome, "lnorm")

fit2_gamma <- fitdist(mydata.1$parentsincome, "gamma", method = "mme")

fit2_weibull <- fitdist(mydata.1$parentsincome, "weibull", method = "mge")

#ggplot of distributions

ggplot(data = mydata.1) +
  geom_histogram(data = as.data.frame(mydata.1$parentsincome), aes(x=mydata.1$parentsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=mydata.1$parentsincome, y=dgamma(mydata.1$parentsincome,fit2_gamma$estimate[1], fit2_gamma$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=mydata.1$parentsincome, y=dweibull(mydata.1$parentsincome,fit2_weibull$estimate[1], fit2_weibull$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=mydata.1$parentsincome, y=dlnorm(mydata.1$parentsincome,fit2_lognormal$estimate[1], fit2_lognormal$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Parents Income in €") +
  ggtitle("Comparing marginal densities for parents income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("parentsdist.pdf")

# here choose weibull!!!! 

####### kids income ######

descdist(mydata.1$kidsincome, discrete=FALSE, boot=5000)
# again looks more like a gamma mixture

fit2_lognormalkids <- fitdist(mydata.1$kidsincome, "lnorm")

fit2_gammakids <- fitdist(mydata.1$kidsincome, "gamma", method = "mme")

fit2_weibullkids <- fitdist(mydata.1$kidsincome, "weibull")

#ggplot of distributions

ggplot(data = mydata.1) +
  geom_histogram(data = as.data.frame(mydata.1$kidsincome), aes(x=mydata.1$kidsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=mydata.1$kidsincome, y=dgamma(mydata.1$kidsincome,fit2_gammakids$estimate[1], fit2_gammakids$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=mydata.1$kidsincome, y=dweibull(mydata.1$kidsincome,fit2_weibullkids$estimate[1], fit2_weibullkids$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=mydata.1$kidsincome, y=dlnorm(mydata.1$kidsincome,fit2_lognormalkids$estimate[1], fit2_lognormalkids$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Childrens Income in €") +
  ggtitle("Comparing marginal densities for childrens income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  coord_cartesian(ylim=c(0, 0.0000205))+
  theme_classic()

ggsave("kidsdist.pdf")
## here choose gamma!



# Generate the multivariate distribution 
# here margins refers to the marginal distribution of the input variables which we choose 
# in a way that they represent the distributional structure of the marginals present

my_dist <- mvdc(surGumbelCopula(par), margins = c("weibull","gamma"), paramMargins = list(list(shape = fit2_weibull$estimate[1], scale = fit2_weibull$estimate[2]), list(shape = fit2_gammakids$estimate[1], rate = fit2_gammakids$estimate[2])))

v <- rMvdc(50000, my_dist)
#write.csv(v, "v.csv")

# Compute the density
pdf_mvd <- dMvdc(v, my_dist)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")

den3d <- kde2d(v[,1],v[,2])
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface(  contours = list(
  z = list(
    show=TRUE,
    usecolormap=TRUE,
    highlightcolor="#ff0000",
    project=list(z=TRUE)
  )
)
) %>%
  layout(
    scene = list(
      camera=list(
        eye = list(x=1.87, y=0.88, z=-0.64)
      )
    )
  )



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
pairs.panels(mydata)
pairs.panels(v)


#ggplot(data = mydata, mapping = aes(x = var_a, y = var_b)) +  geom_point() + geom_smooth()

