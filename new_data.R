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
# or family <- BB7

survgumbel <- surGumbelCopula(param = 1)
basefit <- fitCopula(survgumbel,basepobs,method="mpl")
par <- coef(basefit)[[1]]


bb7 <- BB7Copula(param = c(1,1))
basefit2 <- fitCopula(bb7,basepobs,method="mpl")
param <- coef(basefit2)

#plot for  Survival Gumbel
persp(surGumbelCopula(par), dCopula)
# and more beatiful:
persp(surGumbelCopula(par), dCopula,
      xlab = "transformed x income", ylab = "transformed y income",
      main = "Survival-Gumbel Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'survGumbelcopula.pdf')
dev.off()


#plot for BB7
persp(BB7Copula(param), dCopula)
# and more beatiful:
persp(BB7Copula(param), dCopula,
      xlab = "transformed x income", ylab = "transformed y income",
      main = "Survival-Gumbel Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'bb7copula.pdf')
dev.off()


#sampling from it:
basesample <- rCopula(3965,BB7Copula(param)) #rCopula(3965,surGumbelCopula(par))
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
mybasedist <- mvdc(BB7Copula(param), margins = c("weibull","gamma"), paramMargins = list(list(shape = fit2_weibull$estimate[1], scale = fit2_weibull$estimate[2]), list(shape = fit2_gammakids$estimate[1], rate = fit2_gammakids$estimate[2])))

set.seed(1234)

v <- rMvdc(50000, mybasedist)
#x <- rMvdc(5000, my_dist_bb7)
write.csv(v, "v.csv")



#####################################################################################
                        ########try netto copula#####
#####################################################################################

netto <- select(mydata, one_of(c("schnittek_netto_32", "par_inc_netto")))
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


#perform distance test:
tobject <- tCopula(dim=2,loglikenetto[2],df=loglikenetto[3])

reps <- 500
draws <- 5000

x1 <- replicate(reps, copula.gen2(draws = 5000, copula = tobject))
x1 <- as.data.frame.matrix(x1)
x1.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[1]], y[[2]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))
repsvec <- c(1:reps)
x1.dist <- nullToNA(x1.dist)
x1.dist <- ldply(x1.dist, data.frame)
x1.dist <- bind_cols(x1.dist, reps = repsvec)

j <- ggplot(x1.dist, aes(reps, X..i.., group = 1)) +
  #geom_bar(stat = "identity") +
  geom_point() +
  geom_smooth() +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between gross and net income copula models") +
  scale_x_continuous(breaks = c(seq(from = 0, to = reps, by = 50)))
print(j)
ggsave("base_net_distance.pdf")



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
write.csv(nettosim, "netto.csv")



#####################################################################################
                          ######## consumption copula#####
#####################################################################################

consumption <- select(mydata, one_of(c("schnittek_konsum_32", "par_inc_konsum")))
sum(is.na(consumption$schnittek_konsum_32)) # there are 621 missings
sum(is.na(consumption$par_inc_konsum)) # 620 missings!

consumption <- filter(consumption, !((is.na(schnittek_konsum_32)) | (is.na(par_inc_konsum))))

#use separate fit in order to determine which copula model is adequate:
consopobs <- pobs(as.matrix(consumption))
consoverview <- BiCopEstList(consopobs[,1], consopobs[,2],rotations = T)
constest <- consoverview$summary
arrange(constest, AIC, BIC)
#seems to be appropriate to choose family=13 -> survClayton!
loglikecons <- pre.marginals.copula(consumption)

suclayton <- surClaytonCopula(param = loglikecons[2])
#indicates dependency particularly in the right tails!

#plot
#more beautiful:
persp(surClaytonCopula(param = loglikecons[2]), dCopula,
      xlab = "transformed x consumption", ylab = "transformed y consumption",
      main = "Survical Clayton Copula with ML paramter vector", phi = 20, theta = 30)
dev.copy(pdf,'surclaytoncopula.pdf')
dev.off()

#sampling from it:
#consample <- rCopula(3965,surClaytonCopula(param = loglikecons[2]))

x2 <- replicate(reps, copula.gen2(draws = 5000, copula = suclayton))
x2 <- as.data.frame.matrix(x2)
x2.dist <- sapply(x2, function(y) tryCatch({hellinger(y[[1]], y[[2]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))
repsvec <- c(1:reps)
x2.dist <- nullToNA(x2.dist)
x2.dist <- ldply(x2.dist, data.frame)
x2.dist <- bind_cols(x2.dist, reps = repsvec)

k <- ggplot(x2.dist, aes(reps, X..i.., group = 1)) +
  #geom_bar(stat = "identity") +
  geom_point() +
  geom_smooth() +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between gross income and consumtion copula models") +
  scale_x_continuous(breaks = c(seq(from = 0, to = reps, by = 50)))
print(k)
ggsave("base_consumption_distance.pdf")


# creating multivariate distribution:



#for lognormal
cons1 <- consumption %>% 
  mutate(kidsconsumption = ifelse(schnittek_konsum_32==0,1, schnittek_konsum_32)) %>% 
  mutate(parentsconsumption = ifelse(par_inc_konsum==0,1,par_inc_konsum))

#parents consumtion
descdist(cons1$parentsconsumption, discrete=FALSE, boot=5000)

fit2_lognormal_cons <- fitdist(cons1$parentsconsumption, "lnorm")

fit2_gamma_cons <- fitdist(cons1$parentsconsumption, "gamma", method = "mme")

fit2_weibull_cons <- fitdist(cons1$parentsconsumption, "weibull")


ggplot(data = cons1) +
  geom_histogram(data = as.data.frame(cons1$parentsconsumption), aes(x=cons1$parentsconsumption, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=cons1$parentsconsumption, y=dgamma(cons1$parentsconsumption,fit2_gamma_cons$estimate[1], fit2_gamma_cons$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=cons1$parentsconsumption, y=dweibull(cons1$parentsconsumption,fit2_weibull_cons$estimate[1], fit2_weibull_cons$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=cons1$parentsconsumption, y=dlnorm(cons1$parentsconsumption,fit2_lognormal_cons$estimate[1], fit2_lognormal_cons$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Parents consumption spendings in €") +
  ggtitle("Comparing marginal densities for parents consumtion")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("nettoParentsdist.pdf")

#choose gamma 

descdist(cons1$kidsconsumption, discrete=FALSE, boot=5000)


fit2_lognormal_conskids <- fitdist(cons1$kidsconsumption, "lnorm")

fit2_gamma_conskids <- fitdist(cons1$kidsconsumption, "gamma", method = "mme")

fit2_weibull_conskids <- fitdist(cons1$kidsconsumption, "weibull")

ggplot(data = cons1) +
  geom_histogram(data = as.data.frame(cons1$kidsconsumption), aes(x=cons1$kidsconsumption, y=..density.. , fill="histogram"), colour = NA ,alpha = 0.5) +
  geom_line(aes(x=cons1$kidsconsumption, y=dgamma(cons1$kidsconsumption,fit2_gamma_conskids$estimate[1], fit2_gamma_conskids$estimate[2]), col="gamma distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=cons1$kidsconsumption, y=dweibull(cons1$kidsconsumption,fit2_weibull_conskids$estimate[1], fit2_weibull_conskids$estimate[2]), col="weibull distribution"), alpha = 0.7, size = 1) + 
  geom_line(aes(x=cons1$kidsconsumption, y=dlnorm(cons1$kidsconsumption,fit2_lognormal_conskids$estimate[1], fit2_lognormal_conskids$estimate[2]), col="lognormal distribution"), alpha = 0.7, size = 1) + 
  xlab("Kids consumption spendings in €") +
  ggtitle("Comparing marginal densities for kids consumption")+
  guides(col=guide_legend(title="Distributional Families"))+
  guides(fill=guide_legend(title="Empirical Density"))+
  theme_classic()

ggsave("consKidsdist.pdf")

#also choose gamma!

myconsdist <- mvdc(suclayton, margins = c("gamma","gamma"), paramMargins = list(list(shape = fit2_gamma_cons$estimate[1], rate = fit2_gamma_cons$estimate[2]), list(shape = fit2_gamma_conskids$estimate[1], rate = fit2_gamma_conskids$estimate[2])))
conssim <- rMvdc(5000, mynettodist)
write.csv(conssim, "consumption.csv")









