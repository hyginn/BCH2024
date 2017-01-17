# 03_Models.R
#
# Purpose:  BCH2024 - Models
#
# Version: 1.0
#
# Date:    2017  01  16
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First final version
#
# TODO:
# Time-series, autocorrelation, Fourier transform
#
# ==============================================================================


# Initialize:
load("ygData.RData")
load("ygProfiles.RData")


# ==============================================================================
#      PART ONE: CORRELATION
# ==============================================================================

# In principle, correlation between two variables measures the degree to which
# the value of one variable can be predicted from the value of the other. In
# particular, we usually refer to "linear correlation", i.e. predicting the
# value of one variable as a solution to a linear equation with the other
# variable. (Note that correlation does not necessarily imply causality.) In R,
# the function cov() measures covariance and cor() measures the Pearson
# coefficient of correlation (a normalized measure of covariance). Pearson's
# coeffecient of correlation values range from -1 to 1, with 0 indicating no
# correlation.

# Lets explore intutively what correlation values mean for data:

set.seed(12357)
x <- rnorm(50) # 50 random deviates from a N(0,1)
y <- rnorm(50) # again 50 random values
plot(x, y)
cor(x, y)
# This is uncorrelated data, our variables x and y are drawn from a normal ("Gaussian)
# distribution. cor(x, y) has a small value.

y <- x
plot(x, y)
cor(x, y)
# This is perfectly correlated data, x is drawn from a normal ("Gaussian)
# distribution but y is just the same value. cor(x, y) is one.

# Let's explore this with a bit more variety: here is a function that has a value r as an argument, and a function. We compute y as values that are to r-parts a function of x, and to (1-r) parts random noise. Then we plot x and y, and compute cor(x,y)

plotCor <- function(x, r, f) {
    noise <- (1-r) * rnorm(length(x))
    y <- (r * f(x)) + noise
    plot(x, y)
    return(cor(x, y))
}

x <- rnorm(50)
plotCor(x, 0.99, function(x) { x })
plotCor(x, 0.90, function(x) { x })
plotCor(x, 0.80, function(x) { x })
plotCor(x, 0.40, function(x) { x })
plotCor(x, 0.01, function(x) { x })

# Correlations around 0.4 still correspond to a quite clearly visible trend. But
# note that cor() is not a good measure for non-linear correlations:


# periodic ...
plotCor(x, 0.9, function(x) { cos(x * pi) })

# polynomial ...
plotCor(x, 0.9, function(x) { x^2 })

# exponential
plotCor(x, 0.9, function(x) { exp(x) })

# circular ...
plotCor(cos(x*pi), 0.9, function(x) { sin(acos(x)) })

# In all of these cases, the clear functional relationship should have yielded a
# correlation of around 0.99 - but not a linear correlation, which turns out as
# much lower.

# ==============================================================================
#      PART TWO: REGRESSION
# ==============================================================================

# But let's stay with linear modelling for the moment, i.e. analyzinvariation
# under the assumption of a linear model. Then our task is, for a given data
# set, to infer what the parameters of it's linear model are.

# Here is a synthetic sample of observations that
# could come from measuring height and weight of
# a human cohort.

# We generate random heights in an interval, then
# calculate hypothetical weights according to a simple
# linear equation. Finally we add "errors" to the
# weights.

# The goal of our analysis is to recover the parameters of our synthetic data.
# Note: developing your analysis for synthetic data first is good practice: if
# you can't get the values for synthetic data correct, you can't expect to get
# correct values for real data.

n <- 50
set.seed(83)
hmin <- 1.5  # shortest proband
hmax <- 2.3  # tallest
HW <- data.frame(heights = numeric(n),
                 weights = numeric(n))
# generate a column of numbers in the interval
HW$heights <- runif(n, hmin, hmax)
# generate a column of "weights" with a linear model and add some noise
HW$weights <- 40 * HW$heights + 1 +  rnorm(n, 0, 15)

plot(HW$heights, HW$weights, xlab="Height (m)", ylab="Weight (kg)")

cor(HW$heights, HW$weights) # calculate correlation

# R provides lm() (linear model) for regression analysis:
# Remember: the true parameters were weight = 40 * height + 1
?lm
lm(HW$weights ~ HW$heights)

# What are these numbers?
# How do they relate to our question?
# Is the estimate any good?

# plot a regression line - abline() can take its coefficients directly from the
# output of lm()
abline(lm(HW$weights ~ HW$heights), col="firebrick", lwd=2)

# calculate residuals
res <- resid(lm(HW$weights ~ HW$heights))

# calculate idealized values
fit <- fitted(lm(HW$weights ~ HW$heights))

# plot differences
segments(HW$heights, HW$weights, HW$heights, fit, col="#AA000044")

# plot and analyze residuals. If the fit is good, the correlation of the
# residuals should be close to 0
plot(fit, res)
cor(fit, res)

# Calculate and plot prediction and confidence limits.
# PREDICTION limits give boundaries on future observations,
# they characterize how well the model is expected to
# accommodate new data.
#
# CONFIDENCE limits give boundaries on adequate models.
# They characterize how well we can trust our model
# parameters.

pp <- predict(lm(HW$weights ~ HW$heights), int="p")
pc <- predict(lm(HW$weights ~ HW$heights), int="c")
head(pc)

# Now plot pp and pc limits
# first sort on x
o <- order(HW$heights) # o is an index vector, sorted on x-values
HW2 <- HW[o, ]

# second, recompute pp, pc in sorted order
pc<-predict(lm(HW2$weights ~ HW2$heights), int="c")
pp<-predict(lm(HW2$weights ~ HW2$heights), int="p")

# Then plot
plot(HW2$heights, HW2$weights, xlab="Height (m)", ylab="Weight (kg)",
     ylim = range(HW2$weights, pp))
matlines(HW2$heights, pc, lty=c(1,2,2), col="slategrey")
matlines(HW2$heights, pp, lty=c(1,3,3), col="firebrick")

# This is the proper way to plot a linear regression: the inner boundaries (pc)
# show the possible range of our models - i.e. any regreesion line within these
# limits would explain our data; the outer boundaries show the possible range of
# our values, i.e. new values would be expected to fall within these boundaries.

# ==============================================================================
#      PART TWO: APPLYING REGRESSION TO GENE DISCOVERY
# ==============================================================================

# We can define a model for our expression profiles, and then search for genes that correlate with this model. Here is a model profile (with made-up parameters):

t <- seq(0, 120, by = 5)
myProf <- cos((t / 60) * 2 *pi)
plot(t, myProf)

myCor <- numeric(nrow(ygProfiles))
for (i in 1:nrow(ygProfiles)) {
    myCor[i] <- cor(ygProfiles[i, ], myProf)
}
# That's quite fast. What do we get?
hist(myCor)

# Some correlations are very high. Let's plot the 10 highest, and the 10 highest
# anticorrelations.

sel1 <- order(myCor, decreasing = TRUE)[1:10]
sel2 <- order(myCor, decreasing = FALSE)[1:10]


# ... list what genes these are ...
ygData[sel1, c("stdName", "alias")]
ygData[sel2, c("stdName", "alias")]

# ... and plot their expression profiles
plot(seq(0, 120, by = 5), rep(0, 25), type = "n",
     ylim = c(-1.5, 1.5),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA)   # G1

for (i in 1:10) {
    points(seq(0, 120, by = 5), ygProfiles[sel1[i], ], type = "b", col = "black")
    points(seq(0, 120, by = 5), ygProfiles[sel2[i], ], type = "b", col = "red")
}

# What do we learn? Model based correlations are easy to compute, and will of
# course find profiles that correlate with the models. But a linear correlation
# can be high even if the values are not high - as long as they vary smoothly
# with the model. In order to find genes of interest, we need to consider
# fitting the data to a model with differing parameters: we need non-linear
# curve fitting.


# ==============================================================================
#      PART THREE: NON-LINEAR REGRESSION
# ==============================================================================

# A cyclical expression model with parameters can take the following form with parameters for amplitude, frequency and phase

cycEx <- function(t, A, f, phi) { A * (cos((t - phi) / f) )}

t <- seq(0, 120, by = 5)

plot(t, cE(t, 0.15, 10, 10)) # example parameters.
# Lets consider a profile we saw before
iRow <- sel1[1]
plot(t, ygProfiles[iRow,])

# Our default parameters are not a very good fit:
points(0:120, cycEx(0:120, 0.15, 10, 10), type="l", col="#CC0000")

# Calculate parameters of non-linear fit
y <- ygProfiles[iRow,]

myFit <- nls(y ~ cycEx(t, A, f, phi),
                 start = list(A = 0.15,
                              f = 10,
                              phi = 10) )

points(0:120, cycEx(0:120,
                    coef(myFit)["A"],
                    coef(myFit)["f"],
                    coef(myFit)["phi"]),
       type="l", col="#00CC00")

cor(y, predict(myFit))

# With this, we can calculate nls-fits for all profiles, then select those that
# best match "interesting" models.

N <- nrow(ygProfiles)
nlsResults <- data.frame(A = numeric(N),
                         f = numeric(N),
                         phi = numeric(N),
                         cor = numeric(N))
for (i in 1:N) {
    pBar(i, N)
    y <- ygProfiles[i,]

    try(myFit <- nls(y ~ cycEx(t, A, f, phi),
                 start = list(A = 0.15,
                              f = 10,
                              phi = 10) ), silent = TRUE)
    if (length(myFit) > 0) {
        nlsResults$A[i] <- coef(myFit)["A"]
        nlsResults$f[i] <- coef(myFit)["f"]
        nlsResults$phi[i] <- coef(myFit)["phi"]
        nlsResults$cor[i] <- cor(y, predict(myFit))
    }
}

# What are some good fits? We'd like to see high amplitudes, and good
# correlations:
plot(nlsResults$A, nlsResults$cor)
( sel <- which(nlsResults$A > 0.4 & nlsResults$cor > 0.8) )

ygData[sel, c("stdName", "alias")]
# Interesting ... mostly histones.

plot(seq(0, 120, by = 5), rep(0, 25), type = "n",
     ylim = c(-1.5, 1.5),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA)   # G1

for (i in 1:length(sel)) {
    points(seq(0, 120, by = 5), ygProfiles[sel[i], ], type = "b", col = "black")
}

# Let's review what we have done here: rather than rely on any specific model of
# "cyclical expression", we have fitted cyclical models to the data, according
# to parameters of amplitude, frequency and phase, and recorded how good a fit
# we achieve. This has allowed us, for example, to discover genes with a high
# amplitude of differential expression that are well modelled by a cyclical
# model.

# You can explore different aspects of the fits. Are there genes with
# signifcantly higher frequency? Do we see genes with phase-shifts across the
# entire range of time-points?

# TASK
# How would you use this data to define genes that are, and genes that are not
# cyclically expressed? How would you draw the line? Write code to sort the
# genes into these two categories.



# [END]
