# Preparing R -------------------------------------------------------------
rm(list = ls())

# Packages
source("packages.R")
source("functions.R")

load("overallcop.RDA")
load("60ercop.RDA")
load("70ercop.RDA")
load("80ercop.RDA")

#ugly try

w1 <- rMvdc(500, my_dist_bb7)
x1 <- rMvdc(500, my_dist_60)
y1 <- rMvdc(500, my_dist_70_bb7)
z1 <- rMvdc(500, my_dist_80)

w2 <- rMvdc(1000, my_dist_bb7)
x2 <- rMvdc(1000, my_dist_60)
y2 <- rMvdc(1000, my_dist_70_bb7)
z2 <- rMvdc(1000, my_dist_80)

w3 <- rMvdc(50000, my_dist_bb7)
x3 <- rMvdc(50000, my_dist_60)
y3 <- rMvdc(50000, my_dist_70_bb7)
z3 <- rMvdc(50000, my_dist_80)

w4 <- rMvdc(100000, my_dist_bb7)
x4 <- rMvdc(100000, my_dist_60)
y4 <- rMvdc(100000, my_dist_70_bb7)
z4 <- rMvdc(100000, my_dist_80)

w5 <- rMvdc(500000, my_dist_bb7)
x5 <- rMvdc(500000, my_dist_60)
y5 <- rMvdc(500000, my_dist_70_bb7)
z5 <- rMvdc(500000, my_dist_80)



dist <- matrix(nrow = 5,ncol = 3)

dist[1,1] <- hellinger(w1, x1, lower = 0, upper = Inf, method = 1) 
dist[1,2] <- hellinger(w1, y1, lower = 0, upper = Inf, method = 1) 
dist[1,3] <- hellinger(w1, z1, lower = 0, upper = Inf, method = 1) 

dist[2,1] <- hellinger(w2, x2, lower = 0, upper = Inf, method = 1) 
dist[2,2] <- hellinger(w2, y2, lower = 0, upper = Inf, method = 1) 
dist[2,3] <- hellinger(w2, z2, lower = 0, upper = Inf, method = 1) 

dist[3,1] <- hellinger(w3, x3, lower = 0, upper = Inf, method = 1) 
dist[3,2] <- hellinger(w3, y3, lower = 0, upper = Inf, method = 1) 
dist[3,3] <- hellinger(w3, z3, lower = 0, upper = Inf, method = 1) 

dist[4,1] <- hellinger(w4, x4, lower = 0, upper = Inf, method = 1) 
dist[4,2] <- hellinger(w4, y4, lower = 0, upper = Inf, method = 1) 
dist[4,3] <- hellinger(w4, z4, lower = 0, upper = Inf, method = 1) 

dist[5,1] <- hellinger(w5, x5, lower = 0, upper = Inf, method = 1) 
dist[5,2] <- hellinger(w5, y5, lower = 0, upper = Inf, method = 1) 
dist[5,3] <- hellinger(w5, z5, lower = 0, upper = Inf, method = 1) 

save(dist, file = "tablediff.RDA")

colnames(dist, do.NULL = TRUE, prefix = "col")
colnames(dist) <- c("Hellinger-diff full vs. 60er","Hellinger-diff full vs. 70er","Hellinger-diff full vs. 80er")
rownames(dist, do.NULL = TRUE, prefix = "row")
rownames(dist) <- c("500 RD" ,"1.000 RD","50.000 RD", "100.000 RD" ,"500.000 RD")

stargazer(dist, type = "latex", title = "Comparison of Hellinger Distances between full model and cohort selections", out = "table.tex", style = "qje")

dist <- as.data.frame(dist, row.names = T)

diffplot <- ggplot(as.data.frame(w1), aes(w1[,1], w1[,2]))
diffplot + geom_bin2d()

