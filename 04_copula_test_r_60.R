# Preparing R -------------------------------------------------------------
rm(list = ls())

# Packages
source("packages.R")
source("functions.R")

mydata=read.csv(file="test_60er.csv",head=TRUE,sep=";")

mydata <- select(mydata, one_of(c("schnittek_einzel_32", "par_inc_einzel")))

selectedCopula <- pre.marginals.copula(data=mydata)
par <- selectedCopula[2]
## check these family and parameter values in the copula package

# double check!
survgumbel <- surGumbelCopula(param = 1)
set.seed(500)
m <- pobs(as.matrix(mydata))
fit <- fitCopula(survgumbel,m,method="mpl")
coef(fit)
# muy buen!

# plotting it it looks like this:
persp(surGumbelCopula(par), dCopula)

#sampling from it we can do this easily:
u <- rCopula(3965,surGumbelCopula(par))
plot(u[,1],u[,2],pch='.',col='blue')
cor(u,method='spearman')
pairs.panels(u)
pairs.panels(m)

BiCopChiPlot(u[,1],u[,2], PLOT=TRUE, mode="NULL")
# λi measures a distance of a data point (ui1,ui2) to
# the center of the bivariate data set
# χi corresponds to a correlation coefficient
# (λi,χi) will tend to be located above zero for positively dependent margins

BiCopKPlot(u[,1],u[,2], PLOT=TRUE)

# recode zeros
mydata.1 <- mydata %>% 
  mutate(kidsincome = ifelse(schnittek_einzel_32==0,1, schnittek_einzel_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_einzel==0,1,par_inc_einzel))

###### parents income #######

descdist(mydata.1$parentsincome, discrete=FALSE, boot=5000)
# could again be Lognormal Weibull or Gamma
fit2_lognormal <- fitdist(mydata.1$parentsincome, "lnorm")

fit2_gamma <- fitdist(mydata.1$parentsincome, "gamma", method = "mme")

fit2_weibull <- fitdist(mydata.1$parentsincome, "weibull")

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

ggsave("parentsdist_60.pdf")

# here choose gamma! 

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

ggsave("kidsdist_60.pdf")
## here choose gamma as well!


### thus Copula: Survival Gumbel
### Marginal Parents income: Gamma
### Marginal Kids income: Weibull

my_dist_60 <- mvdc(surGumbelCopula(par), margins = c("gamma","gamma"), paramMargins = list(list(shape = fit2_gamma$estimate[1], rate = fit2_gamma$estimate[2]), list(shape = fit2_gammakids$estimate[1], rate = fit2_gammakids$estimate[2])))

v <- rMvdc(50000, my_dist_60)
#write.csv(v, "v.csv")

# Compute the density
pdf_mvd <- dMvdc(v, my_dist_60)
# Compute the CDF
cdf_mvd <- pMvdc(v, my_dist_60)

# 3D plain scatterplot of the generated bivariate distribution
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")

#the following two are equivalent!
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



#setting suitable range over childs income and parents income 
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



