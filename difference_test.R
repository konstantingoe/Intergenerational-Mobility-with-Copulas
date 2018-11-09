# Preparing R -------------------------------------------------------------
rm(list = ls())

# Packages
source("packages.R")
source("functions.R")

load("overallcop.RDA")
load("60ercop.RDA")
load("70ercop.RDA")
load("80ercop.RDA")

##### properly done now ######

reps <- 500
set.seed(0)

x1 <- replicate(reps, generator(draws = 50000))
x1 <- as.data.frame.matrix(x1)


# both work but sapply is faster:
#for (i in 1:50){
#  tryCatch({
#  x1.dist1[i] <- hellinger(x1[1,i][[1]], x1[2,i][[1]], lower = 0, upper = Inf, method = 1)
#  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }

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
ggsave("60splotdistance500.pdf")

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
ggsave("70splotdistance500.pdf")

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
ggsave("80splotdistance500.pdf")

mean80s.500reps <- mean(x3.dist$X..i.. , na.rm = T)

HellingerDistance <- c(mean60s.500reps,mean70s.500reps,mean80s.500reps)
Sample <- c("Full vs. 60s","Full vs. 70s","Full vs. 80s")

HellingerDistance <- cbind(Sample,HellingerDistance)

xtable(HellingerDistance, caption = "Comparison of Hellinger Distances between full model and cohort selections. Mean over 500 repititions of drawing random samples", auto = T)






