# beta with two informative components
bm <- mixbeta(inf = c(0.5, 10, 100), inf2 = c(0.5, 30, 80))

# robustified with mixcombine, i.e. a 10% uninformative part added
unif <- mixbeta(rob = c(1, 1, 1))
mixcombine(bm, unif, weight = c(9, 1))
