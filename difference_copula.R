# Preparing R -------------------------------------------------------------
rm(list = ls())

# Packages
source("packages.R")
source("functions.R")
set.seed(500)

# overall copula:
load("pobsfull.RDA")
bb7copula <- BB7Copula(param = c(1,1))
fit <- fitCopula(bb7copula,m,method="mpl")
par.full <- coef(fit)
# 60s copula:
load("60sparameter.RDA")
par.60 <- par
rm(par)
# 70s copula:
load("pobs70.RDA")
bb7copula <- BB7Copula(param = c(1,1))
set.seed(500)
fit <- fitCopula(bb7copula,m,method="mpl")
par.70 <- coef(fit)
# 80s copula:
load("80sparameter.RDA")
par.80 <- par
rm(par)

#### simulating ####

reps <- 500

x1 <- replicate(reps, copula.generator(draws = 50000))
x1 <- as.data.frame.matrix(x1)

# computing hellinger distance between full cohorts and 60s cohort with 50 repetitions with  50.000 draws 
x1.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[1]], y[[2]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

# computing hellinger distance between full cohorts and 70s cohort with 50 repetitions with 50.000 
x2.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[1]], y[[3]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

# computing hellinger distance between full cohorts and 80s cohort with 50 repetitions with 50.000 
x3.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[1]], y[[4]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

repsvec <- c(1:500)
x1.dist <- nullToNA(x1.dist)
x1.dist <- ldply(x1.dist, data.frame)
x1.dist <- bind_cols(x1.dist, reps = repsvec)

j <- ggplot(x1.dist, aes(reps, X..i.., group = 1)) +
  #geom_bar(stat = "identity") +
  geom_point() +
  geom_smooth() +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between full data and 60's cohort") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 500, by = 50)))
print(j)
ggsave("60scopuladistance500.pdf")

mean60s.500reps <- mean(x1.dist$X..i.. , na.rm = T)
### now cohort 70s

x2.dist <- nullToNA(x2.dist)
x2.dist <- ldply(x2.dist, data.frame)
x2.dist <- bind_cols(x2.dist, reps = repsvec)

k <- ggplot(x2.dist, aes(reps, X..i.., group = 1)) +
  #geom_bar(stat = "identity") +
  geom_point() +
  geom_smooth() +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between full data and 70's cohort") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 500, by = 50)))
print(k)
ggsave("70scopuladistance500.pdf")

mean70s.500reps <- mean(x2.dist$X..i.. , na.rm = T)

### now cohort 80s

x3.dist <- nullToNA(x3.dist)
x3.dist <- ldply(x3.dist, data.frame)
x3.dist <- bind_cols(x3.dist, reps = repsvec)

l <- ggplot(x3.dist, aes(reps, X..i.., group = 1)) +
  #geom_bar(stat = "identity") +
  geom_point() +
  geom_smooth() +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between full data and 80's cohort") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 500, by = 50)))
print(l)
ggsave("80scopuladistance500.pdf")

mean80s.500reps <- mean(x3.dist$X..i.. , na.rm = T)

HellingerDistance <- c(mean60s.500reps,mean70s.500reps,mean80s.500reps)
Sample <- c("Full vs. 60s","Full vs. 70s","Full vs. 80s")

HellingerDistance <- cbind(Sample,HellingerDistance)

xtable(HellingerDistance, caption = "Comparison of Hellinger Distances between full model and cohort selections. Mean over 500 repititions of drawing random samples", auto = T)

# one graph with the fittet line:

i <- ggplot() +
  #geom_bar(stat = "identity") +
  #geom_point() +
  geom_smooth(data = x1.dist, aes(reps, X..i.., group = 1, color="60s distance")) +
  geom_smooth(data = x2.dist, aes(reps, X..i.., group = 1, color="70s distance")) +
  geom_smooth(data = x3.dist, aes(reps, X..i.., group = 1,color="80s distance")) +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between full data and cohorts compared") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 500, by = 50)))+
  guides(col=guide_legend(title="Cohorts"))+
  theme(legend.position = "top")
print(i)
ggsave("allcohortsfullcompared.pdf")

#---------------------------------------------------------------------------------#

###### testing Chetty's hypothesis of copula stability! #####

# computing hellinger distance between 60s cohort and 70s cohort with 50 repetitions with  50.000 draws 
y1.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[2]], y[[3]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

# computing hellinger distance between 60s cohort and 80s cohort with 50 repetitions with 50.000 
y2.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[2]], y[[4]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

# computing hellinger distance between 70s cohorts and 80s cohort with 50 repetitions with 50.000 
y3.dist <- sapply(x1, function(y) tryCatch({hellinger(y[[3]], y[[4]], lower = 0, upper = Inf, method = 1) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}))

y1.dist <- nullToNA(y1.dist)
y1.dist <- ldply(y1.dist, data.frame)
y1.dist <- bind_cols(y1.dist, reps = repsvec)

mean60s70s.500reps <- mean(y1.dist$X..i.. , na.rm = T)


y2.dist <- nullToNA(y2.dist)
y2.dist <- ldply(y2.dist, data.frame)
y2.dist <- bind_cols(y2.dist, reps = repsvec)

mean60s80s.500reps <- mean(y2.dist$X..i.. , na.rm = T)

y3.dist <- nullToNA(y3.dist)
y3.dist <- ldply(y3.dist, data.frame)
y3.dist <- bind_cols(y3.dist, reps = repsvec)

mean70s80s.500reps <- mean(y3.dist$X..i.. , na.rm = T)

h <- ggplot() +
  #geom_bar(stat = "identity") +
  #geom_point() +
  geom_smooth(data = y1.dist, aes(reps, X..i.., group = 1, color="60s and 70s")) +
  geom_smooth(data = y2.dist, aes(reps, X..i.., group = 1, color="60s and 80s")) +
  geom_smooth(data = y3.dist, aes(reps, X..i.., group = 1,color="70s and 80s")) +
  labs(x = " No. of Repetitions", 
       y = "Hellinger Distance", 
       title = "Hellinger distance between cohorts compared") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 500, by = 50)))+
  guides(col=guide_legend(title="Cohort Distance"))+
  theme(legend.position = "top")
print(h)
ggsave("chettycompared.pdf")

