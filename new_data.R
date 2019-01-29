# Preparing R -------------------------------------------------------------
rm(list = ls())
set.seed(500)

# Packages
source("packages.R")
source("functions.R")

mydata=read.csv(file="test.csv",head=TRUE,sep=";")

#einfach missing row drop
#Da hat sich ein negativer Konsum verirrt
mydata <- mydata %>% 
  mutate(schnittek_konsum_32 = ifelse(schnittek_konsum_32<0,NA,schnittek_konsum_32))

####### confirm previous copula models ########

baseline <- select(mydata, one_of(c("schnittek_einzel_32", "par_inc_einzel")))
#use separate fit in order to determine which copula model is adequate:
basepobs <- pobs(as.matrix(baseline))
baseoverview <- BiCopEstList(basepobs[,1], basepobs[,2],rotations = T)
basetest <- baseoverview$summary
arrange(basetest, AIC, BIC)
#seems to be appropriate to choose family=14 -> survGumbel!
survgumbel <- surGumbelCopula(param = 1)
basefit <- fitCopula(survgumbel,basepobs,method="mpl")
par <- coef(basefit)[[1]]
#plot
persp(surGumbelCopula(par), dCopula)
# and more beatiful:
persp(surGumbelCopula(par), dCopula,
      xlab = "transformed x income", ylab = "transformed y income",
      main = "Survival-Gumbel Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'survGumbelcopula.pdf')
dev.off()
#sampling from it:
basesample <- rCopula(3965,surGumbelCopula(par))
plot(basesample[,1],basesample[,2],pch='.',col='blue')
cor(basesample,method='spearman')
#simulated dependency structure:
pairs.panels(basesample)
#actual rank transformed dependency structure
pairs.panels(basepobs)


#### construct distribution: ####

#for lognormal
baseline1 <- baseline %>% 
  mutate(kidsincome = ifelse(schnittek_einzel_32==0,1, schnittek_einzel_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_einzel==0,1,par_inc_einzel))

# evaluate marginal distribution via descdist
# start with parentsincome
descdist(baseline1$parentsincome, discrete=FALSE, boot=5000)
# could again be Lognormal or Weibull
fit2_lognormal <- fitdist(baseline1$parentsincome, "lnorm")

fit2_gamma <- fitdist(baseline1$parentsincome, "gamma", method = "mme")

fit2_weibull <- fitdist(baseline1$parentsincome, "weibull", method = "mge")

#ggplot of distributions

ggplot(data = baseline1) +
  geom_histogram(data = as.data.frame(baseline1$parentsincome), aes(x=baseline1$parentsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=baseline1$parentsincome, y=dgamma(baseline1$parentsincome,fit2_gamma$estimate[1], fit2_gamma$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=baseline1$parentsincome, y=dweibull(baseline1$parentsincome,fit2_weibull$estimate[1], fit2_weibull$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=baseline1$parentsincome, y=dlnorm(baseline1$parentsincome,fit2_lognormal$estimate[1], fit2_lognormal$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Parents Income in €") +
  ggtitle("Comparing marginal densities for parents income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("baselineParentsdist.pdf")

# here choose weibull!!!! 

#now kids income

descdist(baseline1$kidsincome, discrete=FALSE, boot=5000)
# again looks more like a gamma mixture

fit2_lognormalkids <- fitdist(baseline1$kidsincome, "lnorm")

fit2_gammakids <- fitdist(baseline1$kidsincome, "gamma", method = "mme")

fit2_weibullkids <- fitdist(baseline1$kidsincome, "weibull")

#ggplot of distributions

ggplot(data = baseline1) +
  geom_histogram(data = as.data.frame(baseline1$kidsincome), aes(x=baseline1$kidsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=baseline1$kidsincome, y=dgamma(baseline1$kidsincome,fit2_gammakids$estimate[1], fit2_gammakids$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=baseline1$kidsincome, y=dweibull(baseline1$kidsincome,fit2_weibullkids$estimate[1], fit2_weibullkids$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=baseline1$kidsincome, y=dlnorm(baseline1$kidsincome,fit2_lognormalkids$estimate[1], fit2_lognormalkids$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Childrens Income in €") +
  ggtitle("Comparing marginal densities for childrens income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  coord_cartesian(ylim=c(0, 0.0000205))+
  theme_classic()

ggsave("baseKidsdist.pdf")
## here choose gamma!

# I don't understand the error... must have something to do with the new R version -> package conflictions
# maybe it works on windows? 
mybasedist <- mvdc(surGumbelCopula(par), margins = c("weibull","gamma"), paramMargins = list(list(shape = fit2_weibull$estimate[1], scale = fit2_weibull$estimate[2]), list(shape = fit2_gammakids$estimate[1], rate = fit2_gammakids$estimate[2])))
#anyways! it's the same copula as before! 


#####################################################################################
                        ########try netto copula#####
#####################################################################################

netto <- select(mydata, one_of(c("schnittek_netto_32", "par_inc_netto", "schnittek_einzel_32", "par_inc_einzel")))
any(is.na(netto$schnittek_netto_32)) # there are missings
sum(is.na(netto$schnittek_netto_32)) # 8 missings!

netto <- filter(netto, !is.na(schnittek_netto_32))

#use separate fit in order to determine which copula model is adequate:
nettopobs <- pobs(as.matrix(netto))
nettoverview <- BiCopEstList(nettopobs[,1], nettopobs[,2],rotations = T)
nettotest <- nettoverview$summary
arrange(nettotest, AIC, BIC)
#seems to be appropriate to choose family=2 -> Students t-Copula!

loglikenetto <- pre.marginals.copula(netto)

student <- BiCop(family = 2, par = loglikenetto[2], par2 = loglikenetto[3])

#plot
persp(tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]),dCopula)
# and more beatiful:
persp(tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]), dCopula,
      xlab = "transformed x net income", ylab = "transformed y net income",
      main = "Students t-Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'tcopula.pdf')
dev.off()


#note that the t-copula is symmetrical
#The t-copula emphasizes extreme results: it is usually good for modelling phenomena
#where there is high correlation in the extreme values (the tails of the distribution).

#sampling from it:
nettosample <- rCopula(3965,tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]))
plot(nettosample[,1],nettosample[,2],pch='.',col='blue')

cor(nettosample,method='spearman')
#simulated dependency structure:
pairs.panels(nettosample)
#actual rank transformed dependency structure
pairs.panels(nettopobs)
# we get a somewhat higher spearman's correlation in out copula model



# if we're only interested in the copula model this can be left out!

#parents netto income
descdist(netto$par_inc_netto, discrete=FALSE, boot=5000)

#for lognormal
netto1 <- netto %>% 
  mutate(kidsincome = ifelse(schnittek_netto_32==0,1, schnittek_netto_32)) %>% 
  mutate(parentsincome = ifelse(par_inc_netto==0,1,par_inc_netto))

fit2_lognormal_netto <- fitdist(netto1$parentsincome, "lnorm")

fit2_gamma_netto <- fitdist(netto1$parentsincome, "gamma", method = "mme")

fit2_weibull_netto <- fitdist(netto1$parentsincome, "weibull")


ggplot(data = netto1) +
  geom_histogram(data = as.data.frame(netto1$parentsincome), aes(x=netto1$parentsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=netto1$parentsincome, y=dgamma(netto1$parentsincome,fit2_gamma_netto$estimate[1], fit2_gamma_netto$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=netto1$parentsincome, y=dweibull(netto1$parentsincome,fit2_weibull_netto$estimate[1], fit2_weibull_netto$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=netto1$parentsincome, y=dlnorm(netto1$parentsincome,fit2_lognormal_netto$estimate[1], fit2_lognormal_netto$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Parents Net Income in €") +
  ggtitle("Comparing marginal densities for parents income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("nettoParentsdist.pdf")

#choose gamma 

fit2_lognormal_nettokids <- fitdist(netto1$kidsincome, "lnorm")

fit2_gamma_nettokids <- fitdist(netto1$kidsincome, "gamma", method = "mme")

fit2_weibull_nettokids <- fitdist(netto1$kidsincome, "weibull")

ggplot(data = netto1) +
  geom_histogram(data = as.data.frame(netto1$kidsincome), aes(x=netto1$kidsincome, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=netto1$kidsincome, y=dgamma(netto1$kidsincome,fit2_gamma_nettokids$estimate[1], fit2_gamma_nettokids$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=netto1$kidsincome, y=dweibull(netto1$kidsincome,fit2_weibull_nettokids$estimate[1], fit2_weibull_nettokids$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=netto1$kidsincome, y=dlnorm(netto1$kidsincome,fit2_lognormal_nettokids$estimate[1], fit2_lognormal_nettokids$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Kids Net Income in €") +
  ggtitle("Comparing marginal densities for kids income")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("nettoKidsdist.pdf")

#also choose gamma!

mynettodist <- mvdc(tCopula(dim=2,loglikenetto[2],df=round(loglikenetto[3])), margins = c("gamma","gamma"), paramMargins = list(list(shape = fit2_gamma_netto$estimate[1], rate = fit2_gamma_netto$estimate[2]), list(shape = fit2_gamma_nettokids$estimate[1], rate = fit2_gamma_nettokids$estimate[2])))
nettosim <- rMvdc(5000, mynettodist)
pdf_mvd <- dMvdc(nettosim, mynettodist)
# Compute the CDF
cdf_mvd <- pMvdc(nettosim, mynettodist)
par(mfrow = c(1, 2))
scatterplot3d(nettosim[,1],nettosim[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(nettosim[,1],nettosim[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")

den3d <- kde2d(nettosim[,1],nettosim[,2])
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

den3d.bb7 <- kde2d(x[,1],x[,2])
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


# Plot the data for a visual comparison
plot(netto$par_inc_netto, netto$schnittek_netto_32, main = 'Test dataset x and y', col = "blue")
points(nettosim[,1], nettosim[,2], col = '#FF000020', cex=0.5)
legend('bottomright', c('Observed', 'Simulated'), col = c('blue', 'red'), pch=21)

cor(netto, method = "kendall")
cor(netto, method = "spearman")

cor(nettosim, method = "kendall")
cor(nettosim, method = "spearman")

# looks good... 
pairs.panels(netto)
pairs.panels(nettosim)


#####################################################################################
                          ########try consumption copula#####
#####################################################################################

consumption <- select(mydata, one_of(c("schnittek_konsum_32", "par_inc_konsum", "schnittek_einzel_32", "par_inc_einzel")))
sum(is.na(consumption$schnittek_konsum_32)) # there are missings
sum(is.na(consumption$par_inc_konsum)) # 8 missings!

consumption <- filter(consumption, !((is.na(schnittek_konsum_32)) | (is.na(par_inc_konsum))))

#use separate fit in order to determine which copula model is adequate:
consumptionmat <- select(consumption,one_of(c("schnittek_konsum_32", "par_inc_konsum"))) 
consopobs <- pobs(as.matrix(consumptionmat))
consoverview <- BiCopEstList(consopobs[,1], consopobs[,2],rotations = T)
constest <- consoverview$summary
arrange(constest, AIC, BIC)
#seems to be appropriate to choose family=14 -> survGumbel!

loglikecons <- pre.marginals.copula(consumptionmat)

student <- BiCop(family = 2, par = loglikenetto[2], par2 = loglikenetto[3])

#plot
persp(tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]),dCopula)
# and more beatiful:
persp(tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]), dCopula,
      xlab = "transformed x net income", ylab = "transformed y net income",
      main = "Students t-Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'tcopula.pdf')
dev.off()


#note that the t-copula is symmetrical
#The t-copula emphasizes extreme results: it is usually good for modelling phenomena
#where there is high correlation in the extreme values (the tails of the distribution).

#sampling from it:
nettosample <- rCopula(3965,tCopula(dim=2,loglikenetto[2],df=loglikenetto[3]))
plot(nettosample[,1],nettosample[,2],pch='.',col='blue')

cor(nettosample,method='spearman')
#simulated dependency structure:
pairs.panels(nettosample)
#actual rank transformed dependency structure
pairs.panels(nettopobs)
# we get a somewhat higher spearman's correlation in out copula model





