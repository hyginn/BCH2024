# 02_Features.R
#
# Purpose:  BCH2024 - Features
#
# Version: 1.2
#
# Date:    2017  01  14
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.2    Add ygData table to hold gene-information;
#          Rename combinedProfiles to ygProfiles.
# V 1.1    Merge GSE3635 and GSE4987 expression data sets.
# V 1.0    First final version
#
# TODO:
#
#
# ==============================================================================

# ==============================================================================
#      PART ONE: MERGING GSE3635 AND GSE4987
# ==============================================================================

# 1. Construct a datamodel
# 2. Make an object to hold the data
# 3. Merge the data
#
# == DATAMODEL

# Expression data is typical 2D data that can be organized in a
# spreasdsheet-like format: every row represents values for one gene/ probe/
# spot, and every column represents values obtained from one experiment.

# We downloaded two datasets from GEO: the low-resolution dataset GSE3635
# (sampling at 10' intervals), and the higher-resolution dataset GSE4987
# (sampling at 5' intervals and including a dye-swap replicate). We can store
# these datasets in a common matrix - but first we need to validate whether
# featureNames() retrieves the identical row names. I expect them to, after all
# the experiments are both associeted with the GPL1914 platform, but we need to
# check to make sure. First we load the two datasets from the .RData files I
# have uploaded to the GitHub repository:
#
load("GSE3635.RData")
load("GSE4987.RData")

GSE3635names <- featureNames(GSE3635)
GSE4987names <- featureNames(GSE4987)

# confirm:
length(GSE3635names) == length(GSE4987names) # TRUE
identical(GSE3635names, GSE4987names)  # FALSE


# FALSE ???
# That's not good. Why?
# Let's check if we can find a difference
sum(GSE3635names != GSE4987names)

# There seems to be a non-identical element ...
which(GSE3635names != GSE4987names)
GSE3635names[1671]
GSE4987names[1671]

# Well - that's surprising. But we can confirm that these rows are meant to
# reference identical genes, even though the spelling of the identifier is not
# consistent.

# confirm this:
identical(toupper(GSE3635names), toupper(GSE4987names))

# Good. That solves this question.

# This means we can iterate over the rows, extract data from the expression
# sets, and combine them into a single matrix. Note that we _can_ use a matrix,
# and don't need a data frame, because all our data points will be numeric
# values - we don't need to store strings in this object. We will label the
# rownames with the systematic IDs, and the columns with time. There are just
# some small issues to consider:

# A: Rownames for the 11 control-rows:
GSE3635names[1:15]
# Rownames of a matrix or dataframe must be unique, and they can't include
# blanks or special characters. There is a function that takes care of this ...
make.names(GSE3635names[1:15])

# B: Does the dye-swap reveal a bias?
# Lets randomly pick 500 rows from GSE4987, collect the values for forward and
# reverse measurements, and plot them:
set.seed (112358)
selRows <- sample(1:length(GSE4987names), 1000)
mySet <- exprs(GSE4987)[selRows, ]
# Forward experiments are in columns 1:25, reverse experiments are in columns
# 26:50. Remember that our boxplots had shown that the values range from
# -2 to 2.

# We can collect the values into two vectors:
forward <- numeric()
reverse <- numeric()
for (i in 1:25) {
    forward <- c(forward, mySet[ , i])
    reverse <- c(reverse, mySet[ , i + 25])
}
reverse = -reverse # dye-swap changes the expression sign!

length(forward)  # 25,000 values

plot(forward, reverse)
abline(h = 0, col = "#CCCCFF")
abline(v = 0, col = "#CCCCFF")

# if the samples were identical (no bias) we would expect them to lie
# approximately on a diagonal:
abline(0, 1, col = "#CC0000")

# I think this looks close enough to identical that we will not need to worry
# about correcting for bias. There are a number of spots that have close to 0
# values in one direction but are saturated in the other direction, and there is
# a slight tendency for down-regulated genes in the forward direction to be
# measured a bit less repressed in the reverse direction, but overall we won't
# be making huge mistakes if we simply average the values. However, the plot
# clearly shows us one thing: values of 2 or -2 appear to be measurement
# artefacts and should be excluded (i.e set to NA).

# We can get and set the values with the exprs() function - eg.  for column 1:
sum(exprs(GSE4987)[ , 1] ==  2, na.rm = TRUE)
sum(exprs(GSE4987)[ , 1] == -2, na.rm = TRUE)

# loop over all columns ...
for (i in 1:length(sampleNames(GSE4987))) {
    sel <- ! is.na(exprs(GSE4987))[ , i] &       # not NA, and
           ( exprs(GSE4987)[ , i] ==  2 |        # either 2, or
             exprs(GSE4987)[ , i] == -2 )        # -2
    exprs(GSE4987)[sel, i] <- NA                 # replace with NA
}

# Confirm - repeat the subset and plot:
#
set.seed (112358)
mySet <- exprs(GSE4987)[sample(1:length(GSE4987names), 1000), ]
forward <- numeric()
reverse <- numeric()
for (i in 1:25) {
    forward <- c(forward, mySet[ , i])
    reverse <- c(reverse, mySet[ , i + 25])
}
reverse = -reverse # dye-swap changes the expression sign!
plot(forward, reverse, xlim = c(-2, 2), ylim = c(-2, 2))
abline(h = 0, col = "#CCCCFF")
abline(v = 0, col = "#CCCCFF")
abline(0, 1, col = "#CC0000")

# Wait ... why are there still -2 values in "reverse"? Or are there?
min(reverse, na.rm = TRUE)

# You should appreciate that testing for equality with digital numbers can be
# tricky, due to the finite accurracy of the internal representation. If this is
# ever an issue, it is better to test for near-equality, eg. like
# eps <- 0.001
# abs(a - b) < eps

# let's not forget to remove the saturated values from the other dataset.
for (i in 1:length(sampleNames(GSE3635))) {
    sel <- ! is.na(exprs(GSE3635))[ , i] &       # not NA, and
           ( exprs(GSE3635)[ , i] ==  2 |        # either 2, or
             exprs(GSE3635)[ , i] == -2 )        # -2
    exprs(GSE3635)[sel, i] <- NA                 # replace with NA
}

# Before merging the sets, we should still test whether the high-res and the
# low-res experiments give comparable values: remember that the low-res set has
# samples from 0 to 120 minutes in 10 minute intervals. The high-res set has the
# same in 5 minute intervals. We can define a logical vector that pulls
# corresponding values from the high-res dataset - perhaps like so:
mask <- unlist(strsplit("10101010101010101010101010000000000000000000000000", ""))
( mask <- as.logical(as.numeric(mask)) )
# ... or we could have written
# mask <- c(rep(c(TRUE, FALSE), 13), rep(FALSE, 24))
# ... or defined the vector "by hand".

# collect values ...
set.seed (112358)
sel <- sample(1:length(GSE4987names), 1000)
myHighSet <- exprs(GSE4987)[sel, mask]
myLowSet  <- exprs(GSE3635)[sel, ]
high <- numeric()
low <- numeric()
for (i in 1:13) {
    high <- c(high, myHighSet[ , i])
    low  <- c(low,   myLowSet[ , i])
}
plot(high, low, xlim = c(-2, 2), ylim = c(-2, 2))
abline(h = 0, col = "#CCCCFF")
abline(v = 0, col = "#CCCCFF")
abline(0, 1, col = "#CC0000")

# Hm. What do you think?

# I think that the "low-res" samples are measured in the same sense that we had
# called "reverse" in the high-res samples. So which is it? Lets look at a
# cell-cycle gene that we know is regulated in a cyclical fashion - Swi4
# (YER111C), in row 1726 of our expression sets. According to Figure 1A of
# Pramilla et al (2006), Swi4 expression should peak at 20, and 80 minutes
# (columns 3 and 9), just before the G1/S transition:

plot(exprs(GSE3635)[1726, ], type = "b", ylim = c(-1/2, 1/2))
points(exprs(GSE4987)[1726, mask], type = "b", col = "#CC0000")
abline(h = 0, col = "#CCCCFF")

# It becomes clear that the low-res samples are low where expression is expected
# to be high, and if we want to combine the values, we need to invert them, just
# as we did with the "reverse" samples from the high-res dataset.

# But should we combine samples at all? It's worthwhile to plot a few genes
# explicitly to check. We'll write a function that pulls out corresponding
# values for a plot.
#
# === Digression: yeast gene data

# Let's get a little annoyance out of the way first: in order to look up an
# expression profile by standard name we need to find its systematic name in the
# SGD_features table. Then we need to search for the systematic name in the
# vector that featureNames(GSE3635) returns. Consider:

( x <- which(SGD_features$V5 == toupper("Swi4")) )
( y <- which(featureNames(GSE3635) == SGD_features$V4[x]) )

# This is complex, and error prone - we could easily forget for example that the
# "x" above references a row in SGD_features, not a row in the expression data.
# It would be better to have a data frame in which we can store all data for
# yeast genes for which we have expression profiles, and to make sure that the
# rows in that data match the rows of the expression profiles. At the same time,
# we can fix the pesky problem of upper/lower case in systematic names by
# setting everything to uppercase. Lets call this table "ygData" (yg for "yeast
# gene").

( N <- length(featureNames(GSE3635)) )

# define the data frame
ygData <- data.frame(SGD = character(N),
                     sysName = character(N),
                     stdName = character(N),
                     alias = character(N),
                     description = character(N),
                     stringsAsFactors = FALSE)

# fill in the data (this takes about half a minute ...)
for (i in 1:N) {
    # get the systematic name
    sn <- toupper(featureNames(GSE3635)[i])
    # find the row in SGD_features that contains this systematic name
    sgdRow <- which(toupper(SGD_features$V4) == sn)
    # if found
    if (length(sgdRow) == 1) {    #  collect the data and write to ygData
        ygData$SGD[i] <- SGD_features$V1[sgdRow]
        ygData$sysName[i] <- SGD_features$V4[sgdRow]
        ygData$stdName[i] <- SGD_features$V5[sgdRow]
        ygData$alias[i] <- SGD_features$V6[sgdRow]
        ygData$description[i] <- SGD_features$V16[sgdRow]
    }
}
row.names(ygData) <- featureNames(GSE4987)
save(ygData, file = "ygData.RData")

# confirm:
SGD_features[x, ]
ygData[y,]

#
# === Back to merging the two expression data sets
#
# We were planning to write a function that pulls out corresponding
# values from the two expression data sets, for a plot.
#
plotProfiles <- function(name) {
    # Plots expression profiles, given a standard name.
    # Datasets GSE3635 and GSE4987 must be present, as well as
    # the ygData table.
    # Non-observed time-points in GSE3635 are extrapolated
    # from adjacent measurements.
    thisRow <- which(toupper(name) == ygData$stdName)
    y <- -exprs(GSE3635)[thisRow, ]
    y1 <- numeric()
    y2 <-  exprs(GSE4987)[thisRow, 1:25]
    y3 <- -exprs(GSE4987)[thisRow, 26:50]
    y4 <- numeric()
    for (i in 1:25) {
        if (i %% 2) { # i is odd
            y1[i] <- y[floor(i / 2) + 1]
        } else {      # i is even
            y1[i] <- mean(c(y[i / 2], y[(i / 2) + 1]), na.rm = TRUE)  # extrapolate
        }
        y4[i] <- mean(c(y1[i], y2[i], y3[i]), na.rm = TRUE) # average
    }
    sampleTime <- seq(0, 120, by = 5)
    yMin <- min(c(y1, y2, y3), na.rm = TRUE)
    yMax <- max(c(y1, y2, y3), na.rm = TRUE)
    plot(sampleTime, y1,
         ylim = 1.2 * c(yMin, yMax),
         xlab = "time (min.)", ylab = "logRatio", main = name,
         col = "#999999", type = "b")
    points(sampleTime, y2, type = "b", col = "#CC0000")
    points(sampleTime, y3, type = "b", col = "#00CC00")
    points(sampleTime, y4, type = "b", lwd = 2, col = "#000000")
    abline(h =  0, col = "#CCCCFF")
    abline(v = 60, col = "#CCCCFF")
}

plotProfiles("Swi4")
plotProfiles("Tos4")
plotProfiles("Mbp1")

# What do you think ? I think the averaging procedure appears successful in that
# it seems to smooth out consecutive values in the curves, as we would expect
# from time-series data.

# Thus, to put everything together, we can use code like what we wrote above
# to write the averaged values into a matrix, row by row:
nRows <- length(featureNames(GSE4987))
nCols <- 25
ygProfiles <- matrix(numeric(nRows * nCols), nrow = nRows, ncol = nCols)
for (i in 1:nRows) {
    y <- -exprs(GSE3635)[i, ]
    y1 <- numeric()
    y2 <-  exprs(GSE4987)[i, 1:25]
    y3 <- -exprs(GSE4987)[i, 26:50]
    for (j in 1:nCols) {
        if (j %% 2) { # j is odd
            y1[j] <- y[floor(j / 2) + 1]
        } else {      # j is even
            y1[j] <- mean(c(y[j / 2], y[(j / 2) + 1]), na.rm = TRUE)  # extrapolate
        }
        ygProfiles[i, j] <- mean(c(y1[j], y2[j], y3[j]), na.rm = TRUE)
    }
}
row.names(ygProfiles) <- featureNames(GSE4987)
colnames(ygProfiles) <- sprintf("t.%d", seq(0, 120, by = 5))
save(ygProfiles, file = "ygProfiles.RData")


# ==============================================================================
#      PART TWO: USE limma TO CALCULATE DIFFERENTIAL EXPRESSION
# ==============================================================================

# 1. Based on the plots from the last unit, define two groups.
# 2. Adapt the limma code (from GEO2R) - copy the code below and
#      edit it in myScript.R. To get this code, I labeled samples
#      from the forward-dye experiment of GSE4987 as "Mid" (25, 30, 35 and
#      85, 90, 95 min.) and as "End" (55, 60, 65 and 115, 120 min.) and
#      calculated the Top 250 differentially expressed genes.
#    ... adapt it to work with your defined groups, and work on
#    the merged expression sets.

source("https://bioconductor.org/biocLite.R")

if (!require(limma, quietly=TRUE)) {
    biocLite("limma")
    library(limma)
}


# === original GEO2R code begins here ==========================================

fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "XXXXX000XXX111XXX000XXX11XXXXXXXXXXXXXXXXXXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

# Note: Careful! Is this correct? Are the values not log-transformed already?
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
# === original GEO2R code ends here ============================================

#
# What genes does this discover?



# ==============================================================================
#      PART THREE: PRINCIPAL COMPONENT ANALYSIS
# ==============================================================================

# == PCA INTRODUCTION

# Synthetic data example
# 500 normally distributed samples each: uncorrelated
set.seed(2707)
x1 <- rnorm(500,0,1)
y1 <- rnorm(500,0,1)

# generate y2 corrleated with (dependent on) x1
y2 <- 2*x1 + y1
mean(y2)
y2 <- y2-mean(y2)
mean(y2)
sd(y2)
y2 <- y2 / sd(y2)
sd(y2)
print(sd(y2), digits=22)


# Create a lattice plot with two rows and two columns
oPar <- par(mfrow = c(2,2)) # set new and save old graphics state
par(mfrow = c(2,2))
# four plots ...
hist(x1)
hist(y2)
plot(x1, y1)
plot(x1, y2)
par(oPar) # restore graphics state parameters


# calculate a PCA of x1 and y2
pcaSample <- prcomp(cbind(x1,y2))

# here are the information items from the returned list of results
pcaSample
pcaSample$sdev
pcaSample$rotation
summary(pcaSample)
head(pcaSample$x)
plot(pcaSample$x, xlim=c(-5,5), ylim=c(-5,5))

# Compare the histograms before and after the rotation:
oPar <- par(mfrow = c(2,2))
hist(x1, xlim=c(-4,4), ylim=c(0,150), main="")
hist(y2, xlim=c(-4,4), ylim=c(0,150), main="")
hist(pcaSample$x[,1], xlim=c(-4,4), ylim=c(0,150),
     main="", col=rgb(0.86,0,0,0.5))
hist(pcaSample$x[,2], xlim=c(-4,4), ylim=c(0,150),
     main="", col=rgb(0.31, 0.5, 0.74, 0.5))
par(oPar) # restore graphics state parameters

# Plot the sample along the Principal Components as axes
plot(pcaSample$x[,1],pcaSample$x[,2], xlim=c(-4,4), ylim=c(-4,4))

typeInfo(pcaSample)

?prcomp

# ==================================================
# EDA with PCA
# The relative importance of PCs
# ==================================================

# load one of the sample data sets in the R distribution

library(MASS)
data(crabs)

head(crabs)
# Two types: blue and orange
# Two genders: female and male
# FL frontal lobe size (mm)
# RW rear width (mm)
# CL carapace length (mm)
# CW carapace width (mm)
# BD body depth (mm)

# annotate...
fac <- as.factor(paste(crabs[, 1], crabs[, 2],sep="."))
head(fac)
c(fac[1], fac[51], fac[101], fac[151])
as.numeric(c(fac[1], fac[51], fac[101], fac[151]))

plot(crabs[, 4:8], pch=as.numeric(fac))
plot(crabs[, 4:5], pch=as.numeric(fac))
plot(crabs[, 5:6], pch=as.numeric(fac))

# Apply principal components analysis to the five measured dimensions
head(crabs)
pcaCrabs <- prcomp(crabs[, 4:8])

plot(pcaCrabs)
summary(pcaCrabs)
str(pcaCrabs)

# Plot projections along the components into a scatterplot.
# Axes for points are scaled as values, for vectors as variance
# Default for biplot() is the first and second component.

biplot(pcaCrabs, xlabs=as.numeric(fac))
legend(81, -63,c("1: B.F", "2: B.M", "3: O.F", "4: O.M"), box.col=1, bg="lightgrey")

# Plot the first against the third principal component
biplot(pcaCrabs, xlabs=as.numeric(fac), choices = c(1, 3))
legend(84, -63,
       c("1: B.F", "2: B.M", "3: O.F", "4: O.M"),
       box.col=1,
       bg="lightgrey")

# Plot the second against the third principal component
biplot(pcaCrabs, xlabs=as.numeric(fac), choices = c(2, 3))
legend(-14.8, 16.2,
       c("1: B.F", "2: B.M", "3: O.F", "4: O.M"),
       box.col=1,
       bg="lightgrey")

# ===================================================
# Task: plot the last plot (without vectors) with plotting
# symbols that correspond to the gender and type of crab:
# orange and blue circles for females and triangles for males.

# Advanced: also make symbol-size depend on carapace length or
# the mean of all five measurements.
# ===================================================

plot(pcaCrabs$x[,2], pcaCrabs$x[,2], type ="n")
points(pcaCrabs$x[  1: 50,2], pcaCrabs$x[  1: 50,3], pch=17, col="blue")
points(pcaCrabs$x[ 51:100,2], pcaCrabs$x[ 51:100,3], pch=19, col="blue")
points(pcaCrabs$x[101:150,2], pcaCrabs$x[101:150,3], pch=17, col="orange")
points(pcaCrabs$x[151:200,2], pcaCrabs$x[151:200,3], pch=19, col="orange")


# My solution.
# 1. write a function to color/scale points

crabPoints <- function(n) {
    if (crabs[n,"sp"] == "B") {
        sp = "#0066FF"
    } else {
        sp = "#FF9900"
    }
    if (crabs[n,"sex"] == "M") {
        sex = 24  # triangle
    } else {
        sex = 21  # circle
    }
    points(pcaCrabs$x[n,2],
           pcaCrabs$x[n,3],
           col = sp, bg = sp,
           pch = sex,
           cex = scalePC(pcaCrabs$x[n,1]))
}

# scale by first PC
scalePC <- function(x) {
    x <- x +30
    x <- (x / 55) * 2.5
    x <- x + 0.5
    return(x)
}

# 2. Create only the plotting frame
plot(pcaCrabs$x[,2], pcaCrabs$x[,3], type = "n")

# 3. Plot the points individually

for (i in 1:nrow(crabs)) {
    crabPoints(i)
}

# ==============================================================================
#      PART FOUR: PCA OF THE EXPRESSION DATA SET
# ==============================================================================

# To be updated - do this at home  ...



# ==============================================================================
#      PART FIVE: READING GO AND GOA DATA
# ==============================================================================

# To be updated - do this at home  ...





# [END]
