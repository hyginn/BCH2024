# 03_Models.R
#
# Purpose:  BCH2024 - Models
#
# Version: 1.0
#
# Date:    2017  01  17
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.1    Greatly expanded nls() section
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
myModel <- cos((t / 60) * 2 *pi)
plot(t, myModel, col="#CC0000", type = "l",
     xlab = "t (min.)", ylab = "expression log-ratio")

# We can easily calculate the correlation of a real expression profile with our
# synthetic profile - for example:
iRow <- 814
points(t, ygProfiles[iRow, ], type="b")
cor(ygProfiles[iRow, ], myModel)
# ... or
iRow <- 5571
points(t, ygProfiles[iRow, ], type="b", col="#00CC00")
cor(ygProfiles[iRow, ], myModel)

# Note that the _amplitude_ of our synthetic profile does not matter - the
# coefficient of correlation is high if one set of points can be transformed
# into the other with a linear equation - no matter what the coefficients of the
# equation are.

# Let's calculate correlations for all profiles ...

myCor <- numeric(nrow(ygProfiles))
for (iRow in 1:nrow(ygProfiles)) {
    myCor[iRow] <- cor(ygProfiles[iRow, ], myModel)
}
# That's quite fast. What do we get?
hist(myCor)

# Some correlations are very high. Let's plot the 10 highest correlations, and
# the 10 highest anticorrelations.

sel1 <- order(myCor, decreasing = TRUE)[1:10]
sel2 <- order(myCor, decreasing = FALSE)[1:10]


# ... list what genes these are ...
ygData[sel1, c("stdName", "alias")]
ygData[sel2, c("stdName", "alias")]

# ... and plot their expression profiles
plot(t, rep(0, 25), type = "n",
     ylim = c(-1.5, 1.5),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA)   # G1

for (i in 1:10) {
    points(t, ygProfiles[sel1[i], ], type = "b", col = "black")
    points(t, ygProfiles[sel2[i], ], type = "b", col = "red")
}

# What do we learn? Model based correlations are easy to compute, and will of
# course find profiles that correlate with the models. But a linear correlation
# can be high even if the absolute values are not high - as long as they vary
# linearly with the model. Moreover, we would not find genes that are
# phase-shifted, because these have near-zero correlations. Consider:

myShiftedModel <- cos(((t + 15)/ 60) * 2 *pi)  # shift by 15 minutes
plot(t, myModel, col="#CC0000", type = "l",
     xlab = "t (min.)", ylab = "expression log-ratio")
points(t, myShiftedModel, type="l", col="#AACCFF")

# Apply to an expression profile I picked:
points(t, ygProfiles[309, ], type = "b", col = "black")

# calculate correlations:
cor(ygProfiles[309, ], myModel)
cor(ygProfiles[309, ], myShiftedModel)

# Even though gene 309 (YBR071W) is cyclically expressed, it only has a measly
# 0.25 coefficient of correlation with our original model, whereas shifting
# the model profile by 15 minutes gives a high correlation of 0.88.

# In order to more generally find genes of interest, we need to consider
# fitting the data to a model with adjustable parameters: we need non-linear
# curve fitting.


# ==============================================================================
#      PART THREE: NON-LINEAR REGRESSION
# ==============================================================================

# A cyclical expression model with parameters can take the following form with parameters for amplitude, frequency and phase

cycEx <- function(t, A, phi, f) {
    # cosine function with amplitude A, phase phi (in minutes), and
    # frequency 1/f, scaled for one full cycle corresponding to 60 min.
    A * (cos((((t - phi) * 2 * pi) / 60) / f) )
}

# time, as usual for our profiles, is in minutes
t <- seq(0, 120, by = 5)

# What does this function look like? Let's write a small function to
# conveniently explore parameters:
plotModel <- function(t, A, phi, f, thisCol = "#CC0000", plt = TRUE) {

    ex <- cycEx(t, A, phi, f)
    if (plt) {
        plot(t, ex, col = thisCol, type = "l",
             xlab = "t (min.)", ylab = "expression log-ratio",
             main = "Model",
             sub = sprintf("A: %5.3f, f: %5.3f, phi: %5.3f", A, f, phi)
        )
        abline(h =  0, col = "#DDEEFF")
        abline(v = 60, col = "#DDEEFF")
    } else {
        points(t, ex, col = thisCol, type = "l")
    }
}

# Let's explore a few parameters for cycEx():

# Varying A
plotModel(t, A =  1.0, phi = 0, f = 1.0)
plotModel(t, A =  0.5, phi = 0, f = 1.0, thisCol = "#DD99CC", plt = FALSE)
plotModel(t, A = -1.0, phi = 0, f = 1.0, thisCol = "#FFDDEE", plt = FALSE)

# Varying 1/f
plotModel(t, A = 1.0, phi = 0, f = 1.0)
plotModel(t, A = 1.0, phi = 0, f = 2.0, thisCol = "#DD99CC", plt = FALSE)
plotModel(t, A = 1.0, phi = 0, f = 4.0, thisCol = "#FFDDEE", plt = FALSE)
plotModel(t, A = 1.0, phi = 0, f = 0.5, thisCol = "#990000", plt = FALSE)

# Varying phi
plotModel(t, A = 1.0, phi =  0, f = 1.0)
plotModel(t, A = 1.0, phi =  5, f = 1.0, thisCol = "#DD99CC", plt = FALSE)
plotModel(t, A = 1.0, phi = 15, f = 1.0, thisCol = "#EEBBDD", plt = FALSE)
plotModel(t, A = 1.0, phi = 30, f = 1.0, thisCol = "#FFDDEE", plt = FALSE)


# Let's consider a profile we found in our linear regression analysis: gene 814
# (YDL095W, PMT1), a gene that is required for cell wall rigidity.
iRow <- 814
plot(t, ygProfiles[iRow, ], col = "black", type = "b",
     ylim = c(-0.2, 0.2),
     xlab = "t (min.)", ylab = "expression log-ratio",
     main = sprintf("%s (%s)", ygData$sysName[iRow], ygData$stdName[iRow]))
abline(h =  0, col = "#DDEEFF")
abline(v = 60, col = "#DDEEFF")


# Our default parameters are not bad - after all, we discovered the gene using
# this model (with fixed parameters). (I'm arbitrarily using A = 0.2 here, for
# the optics):
plotModel(t, A =  0.2, phi = 0, f = 1.0, thisCol = "#DD99CC", plt = FALSE)
cor(ygProfiles[iRow, ], cycEx(t, A =  0.2, phi = 0, f = 1.0)) # 0.911

# But let's see if we can improve the fit:

#  1: assign the data to a variable
y <- ygProfiles[iRow,]

#  2: Use nls() to calculate a non-linear least squares fit. While linear
#  least-squares fits have an analytical solution, non-linear fits need to be
#  optimized numerically. This is typically done by varying parameters while
#  calculating the fit, then reporting the results once the best possible choice
#  of parameters has been found. As a numerical (not analytic) procedure, this
#  is subject to problems of convergence (no solution can be found in reasonable
#  time), and to getting stuck in local minima ... as we'll see later.
myFit <- nls(y ~ cycEx(t, A, phi, f),
                 start = list(A = 0.2,
                              phi = 0,
                              f = 1.0) )

myFit

# ... and we can plot the fitted function
plotModel(t, A = coef(myFit)["A"],
             phi = coef(myFit)["phi"],
             f = coef(myFit)["f"],
          thisCol = "#CC0000", plt = FALSE)

cor(ygProfiles[iRow, ], predict(myFit)) # 0.921

# You can see that the curve is closer to the data points, and that the already
# good correlation has improved a bit more.

# Here is a function to plot a profile, and its fitted curve

checkFit <- function(iRow, fit) {
    t <- seq(0, 120, by = 5)
    plot(t, ygProfiles[iRow, ], col = "black", type = "b",
         xlab = "t (min.)", ylab = "expression log-ratio",
         main = sprintf("%s (%s)", ygData$sysName[iRow], ygData$stdName[iRow]))
    abline(h =  0, col = "#DDEEFF")
    abline(v = 60, col = "#DDEEFF")
    mtext(sprintf("Parameters: cor: %5.3f, %s",
                  cor(y, predict(fit)),
                  paste(names(coef(fit)),
                        sprintf("%5.3f", coef(fit))
                        , sep = ": ", collapse = ", ")),
          col = "#DD99CC", side = 1, line = 4)
    points(t, predict(fit), col = "#DD99CC", type = "l")
}


# Using nls(), we can calculate curve-fits of our model for all expression
# profiles, then select those that best match "interesting" parameters.

N <- nrow(ygProfiles)
nlsResults <- data.frame(A = numeric(N),
                         phi = numeric(N),
                         f = numeric(N),
                         cor = numeric(N))
for (i in 1:N) {
    pBar(i, N)  # print a progress bar (function in .utilities.R)
    y <- ygProfiles[i,]

    try(myFit <- nls(y ~ cycEx(t, A, phi, f),
                 start = list(A = 0.15,
                              phi = 0.1,
                              f = 1.0) ), silent = TRUE)
    if (length(myFit) > 0) {
        nlsResults$A[i] <- coef(myFit)["A"]
        nlsResults$phi[i] <- coef(myFit)["phi"]
        nlsResults$f[i] <- coef(myFit)["f"]
        nlsResults$cor[i] <- cor(y, predict(myFit))
    }
}

# What are some good fits? For example, we could look for high amplitudes, and
# good correlations:
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
# model. As you can see in the plot, some of these would have had very poor
# correlations with our original model, because they are phase shifted.
# Consider:
iRow <- sel[1]
points(seq(0, 120, by = 5), ygProfiles[iRow, ], type = "b", col = "green")
cor(ygProfiles[iRow, ], myModel)
# vs:
iRow <- sel[6]
points(seq(0, 120, by = 5), ygProfiles[iRow, ], type = "b", col = "red")
cor(ygProfiles[iRow, ], myModel)


# You can explore different aspects of the fits. Are there genes with
# signifcantly higher frequency? Do we see genes with phase-shifts across the
# entire range of time-points?

# TASK
# How would you use this data to define genes that are, and genes that are not
# cyclically expressed? How would you draw the line?
#

# == REVIEWING THE MODEL
#

# Quite a few of our profiles have rather poor correlations. We could just
# discard them and say: these are not cyclically expressed genes. But it's of
# course better to spend some time and ask what is really going on there. For
# examle, here are genes with high amplitude and poor correlations:

sel <- which(nlsResults$A > 0.4 & nlsResults$cor < 0.5)
nlsResults[sel,]

# Let's write a nice function to plot the profiles and the fitted parameters so
# we can explore them more easily:
#

plotFit <- function(iRow, myA, myPhi, myF) {
    # Without parameters, we just plot the current fit with parameters
    # from nlsResults.
    # With parameters, we try a new fit with the given starting values.
    t <- seq(0, 120, by = 5)
    y <- ygProfiles[iRow,]
    origCor <- nlsResults$cor[iRow]
    origA <- nlsResults$A[iRow]
    origPhi <- nlsResults$phi[iRow]
    origF <- nlsResults$f[iRow]
    plot(t, y, type="b",
         xlab = "", ylab = "log-ratio",
         main = sprintf("%d: %s (%s)",
                        iRow,
                        ygData$sysName[iRow],
                        ygData$stdName[iRow]))
    mtext(sprintf("Original fit:  cor: %5.3f, A: %5.3f, phi: %5.3f, f: %5.3f",
                  origCor,
                  origA,
                  origPhi,
                  origF),
          col = "#AA0000", side = 1, line = 3)
    points(0:120, cycEx(0:120, origA, origPhi, origF),
           type="l", col="#AA0000")
    if (! missing(myA)) { # Try a new fit with these parameters
        myFit <- nls(y ~ cycEx(t, A, phi, f),
                     start = list(A = myA,
                                  phi = myPhi,
                                  f = myF),
                     control = nls.control(maxiter = 200))
        points(0:120, cycEx(0:120,
                            coef(myFit)["A"],
                            coef(myFit)["phi"],
                            coef(myFit)["f"]),
               type="l", col="#00DD88")
        mtext(sprintf("New fit:  cor: %5.3f, A: %5.3f, phi: %5.3f, f: %5.3f",
                      cor(y, predict(myFit)),
                      coef(myFit)["A"],
                      coef(myFit)["phi"],
                      coef(myFit)["f"]),
              col = "#00DD88", side = 1, line = 4)
    }
}

iRow <- 1199  # sel[1], when I ran the code
plotFit(iRow)
# Clearly, this fit has not converged. But if we guess possible parameters ...
plotFit(iRow, 0.05, 10, 1)
# ... we get a pretty decent fit with _much_ better correlation.


iRow <- 5058  # sel[2], when I ran the code
plotFit(iRow)
# Another case of non-convergence - I couldn't find parameters that worked for
# this fit :-(
plotFit(iRow, 0.04, 10, 1.5)


iRow <- 5059  # sel[3], when I ran the code
plotFit(iRow)
# Again, non-convergence.
plotFit(iRow, 0.02, 30, 0.5)
# Converges with these parameters, but the result shows that the function can't
# model the data well.


# == IMPROVING OUR FITTING STRATEGY:
#    (I) Try different parameters and select the best result

# This is pretty trivial - we'll just write a function that tries starting our
# parameter search from a few different options, then checks which one has the
# best correlation and returns these values. The function could be:

bestFit <- function(y) {
    # Tries different parameter settings for nls() and returns the best
    # fit object.
    nlsFits <- list()
    nlsCors <- numeric()
    myPars <- data.frame(A   = c(0.1, 0.1, 0.1, 0.1,  0.03,  0.03),
                         phi = c(0.1,  10,  30,  40,   0.1,   0.1),
                         f   = c(1.0, 1.0, 1.0, 1.0, 0.618, 1.618))
    for (i in 1:nrow(myPars)) {
        try(myFit <- nls(y ~ cycEx(t, A, phi, f),
                         start = list(A = myPars$A[i],
                                      phi = myPars$phi[i],
                                      f = myPars$f[i]) ),
            silent = TRUE)
        if (length(myFit) > 0) {
            nlsFits[[i]] <- myFit
            nlsCors[i] <- cor(y, predict(myFit))
        }
    }
    best <- which(nlsCors == max(abs(nlsCors)))[1]
    return(nlsFits[[best]])
}

# Let's try the procedure for our three problem profiles, and plot the results
iRow <- 1199
( newFit <- bestFit(ygProfiles[iRow, ]) )
checkFit(1199, newFit)

iRow <- 5058
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- 5059
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

# No magic here, but we get much more reasonable results than we had before. We
# also increase our processing requirements by a factor of six, because we are
# trying six parameter combinations for every fit. Let's see if this is worth
# the trouble:

N <- nrow(ygProfiles)
nlsBestFitResults <- data.frame(A = numeric(N),
                                phi = numeric(N),
                                f = numeric(N),
                                cor = numeric(N))
for (i in 1:N) {
    pBar(i, N)  # print a progress bar (function in .utilities.R)
    y <- ygProfiles[i,]
    myFit <- bestFit(y)
    if (length(myFit) > 0) {
        nlsBestFitResults$A[i]   <- coef(myFit)["A"]
        nlsBestFitResults$phi[i] <- coef(myFit)["phi"]
        nlsBestFitResults$f[i]   <- coef(myFit)["f"]
        nlsBestFitResults$cor[i] <- cor(y, predict(myFit))
    }
}

# Let's plot correlation / Amplitude side by side. For ease of comparison, we'll flip all nlsResult amplitudes to be negative, and we'll flip all nlsBestFitResults to be positive. And we'll plot them as solid dots, with high transparency to better visualize the density.

plot(0, 0, type = "l",
     xlim = c(-0.5, 0.5), ylim = c(-1, 1),
     xlab = "A", ylab = "cor")
points(-abs(nlsResults$A), nlsResults$cor, pch=19, col="#CC555505")
points(abs(nlsBestFitResults$A), nlsBestFitResults$cor, col="#33DD3F07")

# Qualitatively, we see great improvement, and quantitatively, eg. considering
# the number of gene expression profiles we have fit with a coefficient of
# correlation better than 0.8 ...

sum(nlsResults$cor > 0.8)
sum(nlsBestFitResults$cor > 0.8)

# For good measure, let's inspect some of the profiles of the "worst of the best", eg. ranked at position 400:

sel <- order(nlsBestFitResults$cor, decreasing = TRUE)[400:405]
sel  # [1] 1207 2578 4428  281 2501  439

nlsBestFitResults[sel, ]  # What do the parameters mean?

# Inspect the profiles and fits ...
iRow <- sel[1]
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- sel[2]
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- sel[3]
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- sel[4]
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- sel[5]
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

# This could hardly have turned out better: the last two rows show an important
# limitation on our fit: our model is constrained to be symmetric about 0, and
# if the data is not symmetric, but shifted, as in rows 2501 an 439, we can't
# get a good fit. Which leads us directly to:

# == IMPROVING OUR FITTING STRATEGY:
#    (II) Add parameters to our model, for flexibility

# A model can have too many parameters, and that will give rise to "overfitting"
# - satisfying the mathematics of the model, but not the physics or the biology
# of the data. But in our case, with 25 data points and only three parameters to
# fit, adding one or two parameters to the model should be fine. I would like to
# add two parameters: (I) a vertical offset, to release the constraint of our
# fitted model to be symmetric around 0, and (II) a damping function, to handle
# attenuation of expression - after all we saw a large component of global
# attenuation in our first Principal Component.

# Let's see what mathematical form such parameters can take.
# This was our original function ...
cycEx <- function(t, A, phi, f) {
    # cosine function with amplitude A, phase phi (in minutes), and
    # frequency 1/f, scaled for one full cycle corresponding to 60 min.
    A * (cos((((t - phi) * 2 * pi) / 60) / f) )
}

# ... and here we add a vertical offset B and an exponential damping term,
# exp(-k * t)...
cycEx2 <- function(t, A, phi, f, k, B) {
    # cosine function with amplitude A, phase phi (in minutes), and
    # frequency 1/f, scaled for one full cycle corresponding to 60 min,
    # with damping term exp(-k *t) and vertical offset B
    (exp(-k *t) * (A * (cos((((t - phi) * 2 * pi) / 60) / f) ))) + B
}

# Let's overwrite plotModel()  to
# conveniently explore parameters:
plotModel <- function(t, A = 1.0, phi = 0, f = 1.0, k = 0, B = 0,
                      thisCol = "#CC0000", plt = TRUE) {

    ex <- cycEx2(t, A, phi, f, k, B)
    if (plt) {
        plot(t, ex, col = thisCol, type = "l",
             ylim = c(min(ex) * 1.2, max(ex) * 1.2),
             xlab = "t (min.)", ylab = "expression log-ratio",
             main = "Model",
             sub = sprintf("A: %5.3f, f: %5.3f, phi: %5.3f, k: %5.3f, B: %5.3f",
                           A, f, phi, k, B)
        )
        abline(h =  0, col = "#DDEEFF")
        abline(v = 60, col = "#DDEEFF")
    } else {
        points(t, ex, col = thisCol, type = "l")
    }
}

# Varying B .. trivial
plotModel(t, B = 0)
plotModel(t, B = 0.2,  thisCol = "#DD99CC", plt = FALSE)
plotModel(t, B = -0.2, thisCol = "#FFDDEE", plt = FALSE)
plotModel(t, A = 0.5, B = -0.5, thisCol = "#CC99DD", plt = FALSE)

# Varying k
plotModel(t, k = 0)
plotModel(t, k = 0.01, thisCol = "#DD99CC", plt = FALSE)
plotModel(t, k = 0.02, thisCol = "#FFDDEE", plt = FALSE)
plotModel(t, A = 0.5, k = -0.008, thisCol = "#22EE66", plt = FALSE)


# Ok ... but does it fit?
# Let's update our bestFit function - and let's see if we can get away with not
# having to try additional parameters.

bestFit2 <- function(y) {
    # Tries different parameter settings for nls() and returns the best
    # fit object.
    nlsFits <- list()
    nlsCors <- numeric()
    myPars <- data.frame(A   = c(0.1, 0.1, 0.1, 0.1,  0.03,  0.03),
                         phi = c(0.1,  10,  30,  40,   0.1,   0.1),
                         f   = c(1.0, 1.0, 1.0, 1.0, 0.618, 1.618))
    for (i in 1:nrow(myPars)) {
        try(myFit <- nls(y ~ cycEx2(t, A, phi, f, k, B),
                         start = list(A = myPars$A[i],
                                      phi = myPars$phi[i],
                                      f = myPars$f[i],
                                      k = 0.01,
                                      B = 0.01)),
            silent = TRUE)
        if (length(myFit) > 0) {
            nlsFits[[i]] <- myFit
            nlsCors[i] <- cor(y, predict(myFit))
        }
    }
    best <- which(nlsCors == max(abs(nlsCors)))[1]
    return(nlsFits[[best]])
}

iRow <- sel[1] # 1207
( newFit <- bestFit2(ygProfiles[iRow, ]) )
checkFit(iRow, newFit)
# the old coefficient of correlation was 0.858

iRow <- sel[2] # 2578
checkFit(iRow, bestFit2(ygProfiles[iRow, ]))

iRow <- sel[3] # 4428
checkFit(iRow, bestFit2(ygProfiles[iRow, ]))

iRow <- sel[4] # 281
checkFit(iRow, bestFit2(ygProfiles[iRow, ]))
# This was our problem fit, with the vertical offset ... Compare to the old fit
checkFit(iRow, bestFit(ygProfiles[iRow, ]))

iRow <- sel[5]  # 2501
checkFit(iRow, bestFit2(ygProfiles[iRow, ]))
# This was the second problem fit, MUCH better now ... Compare to the old fit
checkFit(iRow, bestFit(ygProfiles[iRow, ]))
# Once again ...
checkFit(iRow, bestFit2(ygProfiles[iRow, ]))

# Obviously, you could now calculate all the fits again, and mine the results
# for genes with similar parameters (coexpressed?), phase shifted (causally
# related), or other interesting parameter combinations. There is no substitute
# for playing with the data and exploring it, to develop intuitions that will
# allow you to make discoveries.
#

# In the end, what have we learned through this?

# Non-linear modelling gives us a flexible way to query data for internal
# structure. We can now easily find expression profiles that correspond to
# "interesting" models, given our understanding of amplitude, phase shift and
# attenuation of the expression. Where we first were looking for structure in 25
# time-points each, and later in a handful of principal components, we can now
# query five parameters, e.g. to find genes with significant, or similar
# expression profiles.

# But to actually find "the most interesting" genes cannot be automated. Our
# tools help us view the data, they do not interpret the results. As always:
# data does not interpret itself. We have constructed some nicely sophisticated
# tools, but in a sense that has only shifted our task: from looking at raw
# data, to looking at parameter values. I would argue that much is gained by
# being able to query the data in more principled ways than just visual
# appearance, but still, the problem of biological relevance does not solve
# itself.


# ==============================================================================
#      PART FOUR: MASTERY
# ==============================================================================

# TASK:
# Here is a nice and useful challenge, to help you master this important
# topic. Write code for nls() fits of either of the below:

# - sigmoidal data (such as a protein denaturation curve) - see:
#   https://en.wikipedia.org/wiki/Logistic_function

# - a gene impulse expression model (see Chechik and Koller, 2009 -
#   https://www.ncbi.nlm.nih.gov/pubmed/19193146; see also
#   Sander et al. 2016 https://www.ncbi.nlm.nih.gov/pubmed/27797772)

# - any other non-linear regression problem that may be useful in the
#   context of your thesis.

# Contact me (or the list) if you need guidance how to tackle such an
# undertaking. In principle, you should simply be able to adapt code you have
# encountered above, begin by defining some synthetic data, define the function
# you wish to use, plot a few parameter settings, then set up the fit with
# reasonable starting values and plot the results. Work in small steps, and make
# sure to validate at every step that what you are doing is correct.

# Enjoy!


# [END]
