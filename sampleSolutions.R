# sampleSolutions.R
#
# Purpose:  BCH2024 - Sample solutions to tasks in the various scripts
#
# Version: 1.0
#
# Date:    2017  01  13
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    Data
#
# TODO:
#          Add analysis of row-wise and column-wise data for Data unit
#
# ==============================================================================

# Do not look at this code if you have not tried solving the task yourself.
#
# There are often many ways to achieve the same goal. My coding style for
#    teaching code is slow, defensive, and explicit. Your code may be faster,
#    and more compact. Develop an intution as to what style suits you well. But
#    if your style leads to slower code, take note of what I do differently.




# ==============================================================================
#      PART ONE: TASKS FOR 01_Data
# ==============================================================================

# Analyzing gset ...
#  Are values present for each row?
#     We could recast this to ask: what is the maximum number of missing values
#     in a row. We can use an apply() call, or a for-loop. Here's a for-loop:

GSE3635names <- featureNames(gset)
nNAs <- numeric(length(GSE3635names))
for (i in 1:length(GSE3635names)) {
    nNAs[i] <- sum(is.na(exprs(gset)[i, ]))
}
# What's  the maximum?
max(nNAs)
# confirm
( maxNA <- which(nNAs == max(nNAs)) )
GSE3635names[maxNA]
exprs(gset)[maxNA, ]

    rows <-
x <- rows[-(grep("^Y[A-Z][LR]", rows))]
x <- rows[duplicated(rows)]
# What's the fraction of spots that have at least one value missing?
sum(nNAs > 0) / length(GSE3635names) # about 2%

#  Are all rows genes?
#  What identifiers are being used?
#     (cf. http://www.yeastgenome.org/help/community/nomenclature-conventions)

( x <- GSE3635names[-(grep("^Y[A-Z][LR]", GSE3635names))] )

#  Are all rows/genes unique?
GSE3635names[duplicated(GSE3635names)]

#  Are all yeast genes accounted for?
#  ... this involves reading and processing a features file from SGD:

# Note! Reading this file as-is crashes RStudio!
# There is an unbalanced quotation mark.
# I have removed it in the source by removing the alias that contained it (B").

SGD_features <- read.delim(file = "data/SGD_features.tab",
                           header = FALSE,
                           stringsAsFactors = FALSE)
SGD_features <- SGD_features[ , c(1, 4, 5, 6, 16)]
nrow(SGD_features)
SGD_features <- SGD_features[ SGD_features$V4 != "", ]
nrow(SGD_features)
SGD_features <- SGD_features[grep("^Y[A-Z][LR]", SGD_features$V4), ]
nrow(SGD_features)
length(unique(SGD_features$V4))

x <- SGD_features$V4[-(grep("^Y[A-Z][LR]", SGD_features$V4))]

y1 <- setdiff(featureNames(gset), SGD_features$V4)
y2 <- setdiff(SGD_features$V4, featureNames(gset))

# TBC ...
# === Column-wise analyses
# Simple data descriptors
x <- exprs(gset)[ , 1]
mean(x)     # This needs fixing !
median(x)
IQR(x)
var(x)
sd(x)
summary(x)


# TBC ...
# === Row-wise analyses
# Expression level plot for selected genes:
#
# Task: Plot expression for
# Mbp1, Swi6, Swi4, Nrm1, Cln1, Clb6, Act1, and Alg9



# ==============================================================================
#      PART TWO: TASKS FOR 02_Features
# ==============================================================================


#   - Write code to plot the genes with the five highest and lowest loadings
#     of PC1, highest in black, and lowest in red.

# get the five highest values from pcaYG$x[ ,1] ...
sel1 <- order(pcaYG$x[ ,2], decreasing = TRUE)[1:5]
# and the five lowest values from pcaYG$x[ ,1] ...
sel2 <- order(pcaYG$x[ ,2], decreasing = FALSE)[1:5]

plot(seq(0, 120, by = 5), rep(0, 25), type = "n",
     ylim = c(-1.5, 1.5),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA)   # G1

for (i in 1:5) {
    points(seq(0, 120, by = 5), ygProfiles[sel1[i], ], type = "b", col = "black")
    points(seq(0, 120, by = 5), ygProfiles[sel2[i], ], type = "b", col = "red")
}

# ==============================================================================
#      PART THREE: TASKS FOR 03_Models
# ==============================================================================

# Review: Curve fit of gene expression: measuring GFP fluorescence in one-minute
# intervals yields the following data.
myDat <- data.frame( fl = c(
    710.1697,  689.2170,  714.3733,  688.9666,  653.9151,  628.6372,
    600.6179,  558.6891,  522.2780,  507.8960,  573.6212,  539.6243,
    525.4207,  539.3329,  524.4395,  501.4929,  492.4278,  480.4025,
    511.6782,  458.4573,  445.0892,  443.1762,  387.9263,  412.9964,
    392.0517,  364.3173,  371.6392,  392.1216,  364.8162,  375.5662,
    379.2305,  341.3690,  285.1799,  309.8994,  303.8815,  334.6760,
    342.3248,  283.7887,  272.3799,  234.8944,  294.6621,  278.1166,
    297.8715,  238.9789,  262.9021,  216.6459,  245.6330,  235.6556,
    241.1547,  267.7355,  200.9985,  173.5843,  223.2996,  181.2174,
    191.2544,  175.2843,  207.8503,  177.1721,  187.3568,  190.0217,
    121.8007,  136.6225,  144.2138,  190.5096,  150.6899,  145.2135,
    124.6761,  115.2702,  138.4352,  96.1188,   107.3551,  103.1700,
    141.4267,  103.9856,  137.1076,  88.6598,   106.0739),
    t = 1:77)

plot(myDat$t, myDat$fl)
# N_t = N_0 * exp(-lambda * t)
#
# 1/2 = exp(-lambda * t_half)
# log(2)/t_half = lambda

expDecay <- function(t, N_0, t_half) {
    lambda <- log(2)/t_half
    return(N_0 * exp(-lambda * t))
}

points(0:80, expDecay(0:80, 729, 28), type="l", col="#00CC00")

# How can we model this?
y <- myDat$fl
myFit <- nls(y ~ expDecay(myDat$t, N_0, t_half),
             start = list(N_0 = 729, t_half = 28))

points(myDat$t, expDecay(myDat$t, coef(myFit)["N_0"], coef(myFit)["t_half"]),
       type="l", col="#00CC00")

# What are the relevant parameters?

myResid <- resid(myFit)
plot(myDat$t, myResid )

cor(myDat$t, myResid )  # small. Good.
mean(myResid)
# ==============================================================================

# Impulse function:
# (Code courtesy of Scott :-)

# Model for impulse data
impEx <- function(h0, h1, h2, t1, t2, B, x) {
    S1 <- 1 / (1 + exp(1)^(-B * (x - t1)))
    S2 <- 1 / (1 + exp(1)^(B * (x - t2)))

    s1 <- h0 + (h1 - h0) * S1
    s2 <- h2 + (h1 - h2) * S2

    f <- ((1 / h1) * s1 * s2)

    return(f)
}

# Adding error to the impulse data calcualtion
impEx2 <- function(h0, h1, h2, t1, t2, B, x, err = 1) {
    noise <- (1-err) * rnorm(length(x))

    S1 <- 1 / (1 + exp(1)^(-B * (x - t1)))
    S2 <- 1 / (1 + exp(1)^(B * (x - t2)))

    s1 <- h0 + (h1 - h0) * S1
    s2 <- h2 + (h1 - h2) * S2

    f <- ((1 / h1) * s1 * s2) + noise

    return(f)
}

# Function to plot the model
plotModel <- function(h0, h1, h2, t1, t2, B, x, thisCol = "#CC0000", plt = TRUE) {

    ex <- impEx(h0, h1, h2, t1, t2, B, x)
    if (plt) {
        plot(x, ex, col = thisCol, type = "b",
             xlab = "t (min.)", ylab = "expression log-ratio",
             main = "Model",
             sub = sprintf("ho: %5.3f, h1: %5.3f, h2: %5.3f, t1: %5.3f, t2: %5.3f, B: %5.3f", h0, h1, h2, t1, t2, B)
        )
        abline(h =  0, col = "#DDEEFF")
        abline(v = 60, col = "#DDEEFF")
    } else {
        points(x, ex, col = thisCol, type = "l")
    }
}

# Plot the model
x <- 1:100
plotModel(1, 10, 5, 20, 30, 1, x)

x <- 1:100
plotModel(4.8, -8, 10.8, 19.8, 30, 1, x)



# Generate and plot some synthetic data with errors
y <- impEx2(1, 10, 5, 20, 30, 1, x, 0.5)
plot(x, y, type = "b")



# Perform a fit on the synthetic data
myFit <- nls(y ~ impEx(h0, h1, h2, t1, t2, B, x),
             start = list(h0 = 6,
                          h1 = -8,
                          h2 = 0,
                          t1 = 18,
                          t2 = 32,
                          B = 1.3 ))

# Plot the fit over the synthetic data
plotModel(x = x, h0 = coef(myFit)["h0"],
          h1 = coef(myFit)["h1"],
          h2 = coef(myFit)["h2"],
          t1 = coef(myFit)["t1"],
          t2 = coef(myFit)["t2"],
          B = coef(myFit)["B"],
          thisCol = "#CC0000", plt = FALSE)

# Try this with real data

GSE59784 <- read.delim("GSE59784_rna.norm.txt")



# [END]
