bmix <- mixbeta(inf1 = c(0.2, 8, 3), inf2 = c(0.8, 10, 2))
plot(bmix)

rbmix <- robustify(bmix, weight = 0.1, mean = 0.5)
rbmix
plot(rbmix)


gmnMix <- mixgamma(inf1 = c(0.2, 2, 3), inf2 = c(0.8, 2, 5), param = "mn")
plot(gmnMix)

rgmnMix <- robustify(gmnMix, weight = 0.1, mean = 2)
rgmnMix
plot(rgmnMix)


nm <- mixnorm(inf1 = c(0.2, 0.5, 0.7), inf2 = c(0.8, 2, 1), sigma = 2)
plot(nm)

rnMix <- robustify(nm, weight = 0.1, mean = 0, sigma = 2)
rnMix
plot(rnMix)
