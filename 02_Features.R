# 02_Features.R
#
# Purpose:  BCH2024 - Features
#
# Version: 1.4.1
#
# Date:    2017  01  18
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.4.1  Typo
# V 1.4    Add PCA of expression analysis
# V 1.3.1  Typos and accidentally deleted code in the GEO2R section.
# V 1.3    Add section to adapt GEO2R code.
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

# === ADAPTING LIMMA CODE

# The GEO2R limma code is virtually uncommented, and has not been written for
# clarity, but for being easily produced by a code-generator that is triggered
# by the parameters on the GEO Web-site. Most of it is actually dispensable for
# our work with the merged datasets. But to figure that out requires a bit of
# reverse-engineering. To get this code from the GEO2R Website, I labeled
# samples from the forward-dye experiment of GSE4987 as "Mid" (25, 30, 35 and
# 85, 90, 95 min.) and as "End" (55, 60, 65 and 115, 120 min.) and calculated
# the Top 250 differentially expressed genes. Here we will adapt it to work with
# our merged dataset.


# In principle, the code goes through three steps:
# 1. Prepare the data
# 2. Define groups that are to be contrasted
# 3. Find genes that are significantly differentially expressed across groups
# 4. Format results.


source("https://bioconductor.org/biocLite.R")

if (!require(limma, quietly=TRUE)) {
    biocLite("limma")
    library(limma)
}

# The GEO dataset that this code was produced for is GSE4987
load("GSE4987.RData")
# The script uses the object name "gset" for the data
gset <- GSE4987


# === original GEO2R code begins here ==========================================

# With the first step, names within the gset object are made conformant to
# R's row- and column name requirements.

# GEO2R > fvarLabels(gset) <- make.names(fvarLabels(gset))

# We won't be needing this since we are handling our own annotations in ygData.
# We could use the anotations in gset, but they are now more than ten years old.

# The next step defines groups through a string that was constructed by the GEO
# server. It contains an "X" for all columns that should be ignored, a number
# for all columns that should be analyzed, where the number corresponds to the
# group. Here the groups are 0 and 1,

# GEO2R > gsms <- "XXXXX000XXX111XXX000XXX11XXXXXXXXXXXXXXXXXXXXXXXXX"
# GEO2R > sml <- c()
# GEO2R > for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# GEO2R >  eliminate samples marked as "X"
# GEO2R > sel <- which(sml != "X")
# GEO2R > sml <- sml[sel]
# GEO2R > gset <- gset[ ,sel]

# We won't be needing this either, but let's rememeber that "sml" is now a vector
# of numbers that label each column with a group it should belong to.

# The next part - calculating log values - may not be correct anyway: the data
# seem to contain log-ratio values already, since they are symmetric about 0.

# GEO2R > # log2 transform
# GEO2R > ex <- exprs(gset)
# GEO2R > qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

# GEO2R > LogC <- (qx[5] > 100) ||
# GEO2R >     (qx[6]-qx[1] > 50 && qx[2] > 0) ||
# GEO2R >     (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
# GEO2R > if (LogC) { ex[which(ex <= 0)] <- NaN
# GEO2R > exprs(gset) <- log2(ex) }

# In the next step, the numbers are prefixed with a letter G (G0, G1 ...) and
# declared to be "factors" in the statistical design. They are attached to gset.
# We will have to do something equivalent for our ygProfiles set ...

# GEO2R > # set up the data and proceed with analysis
# GEO2R > sml <- paste("G", sml, sep="")    # set group names
# GEO2R > fl <- as.factor(sml)
# GEO2R > gset$description <- fl

# The limma functions and objects depend on being of a certain class. It's not
# trivial to make a valid object from raw numbers. But we can just copy the
# GSE4987 dataset, remove all but 11 columns, and overwrite the data with
# our combined data values.  Our orginal dataset grouped columns 6,7,8 and 18,
# 19, 20 into group 0, columns 12, 13, 14 and 24, 25 into group 1.
sel <- c(6:8, 18:20, 12:14, 24, 25)
mySet <- GSE4987[ , sel]
colnames(mySet) <- colnames(ygProfiles)[sel]
exprs(mySet) <- ygProfiles[ , sel]

# GEO2R > design <- model.matrix(~ description + 0, gset)
# GEO2R > colnames(design) <- levels(fl)
#
# We define a vector of group labels and build the "design Matrix" for the
# statistical test:
( myGroups <- as.factor(paste("G", c(rep(0, 6), rep(1, 5)), sep="")) )
myDesign <- model.matrix(~ myGroups + 0, mySet)
colnames(myDesign) <- levels(myGroups)
myDesign

# GEO2R > fit <- lmFit(gset, design)
#
# Now we can calculate the fit to a linear model, depending on the two factors
myFit <- lmFit(mySet, myDesign)


# GEO2R > cont.matrix <- makeContrasts(G1-G0, levels=design)
# GEO2R > fit2 <- contrasts.fit(fit, cont.matrix)
# GEO2R > fit2 <- eBayes(fit2, 0.01)
# GEO2R > tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# GEO2R > tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
# GEO2R > write.table(tT, file=stdout(), row.names=F, sep="\t")

myCont.matrix <- makeContrasts(G1-G0, levels = myDesign)
myFit2 <- contrasts.fit(myFit, myCont.matrix)
myFit2 <- eBayes(myFit2, 0.01)
myTable <- topTable(myFit2, adjust="fdr", sort.by="B", number = 50)

mytT <- subset(myTable, select=c("ID","P.Value","B","Gene.symbol","Gene.title"))
write.table(mytT, file=stdout(), row.names=F, sep="\t")

# Careful - these IDs are now factors, which will mess up our use to extract
# values according to gene IDs. Not going into the details here but:
str(mytT)
mytT$ID <- as.character(mytT$ID)
str(mytT)

# Let's see what we got: let's plot the expression profiles for the top ten
# genes (note that we are taking the values we are plotting from ygProfiles, our
# mySet was only needed to discover the genes...):

plot(seq(0, 120, by = 5), rep(0, 25), type = "n",
     ylim = c(-0.6, 0.6),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA) # G1

for (i in 1:10) {
    thisID <- mytT$ID[i]
    points(seq(0, 120, by = 5), ygProfiles[thisID, ], type = "b")
}

# ... and let's also plot them according to the values we fed to limma:

plot(1:11, rep(0, 11), type = "n",
     ylim = c(-0.5, 0.5),
     xlab = "sample", ylab = "log-ratio expression")
rect(0.5, -2,  6.5, 2, col = "#dfeaf4", border = NA)   # G0
rect(6.5, -2, 11.5, 2, col = "#f4dfdf", border = NA)   # G1

sel <- c(6:8, 18:20, 12:14, 24, 25)  # selection of groups ...

for (i in 1:10) {
    thisID <- mytT$ID[i]
    points(1:11, ygProfiles[thisID, sel], type = "b")
}



# What do we learn? limma will return to us genes that are significantly
# different in expression according to the groups we define. This is not
# necessarily the same as a cyclically varying gene, nor does it necessarily
# find the genes whose expression levels are _most_ different, i.e. if the
# variance of a highly, differentially expressed gene within a group is large,
# it may not be very significant. Also, we are not exploiting the fact that
# these values are time series. Nevertheless, we find genes for which we see a
# change in expression levels along two cell-cycles.


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

# Principal component analysis can determine structure in our data. Let's apply
# it to our yGProfiles.

pcaYG <- prcomp(ygProfiles)
summary(pcaYG)

# Lets plot the importance of the components, with a color gradient that
# reflects their proportion of variance. Using color appropriately is a very
# important aspect of data analysis, and crucial for the interpretation of
# plots. I have uploaded a special script on colour - Colour.R - and I encourage
# you to work through this before continuing here.

myCol <- colorRampPalette(c("#be0a23",
                            "#666669", "#99999b", "#bbbbbb","#eeeeec"))

barplot(pcaYG$sdev, col=myCol(length(pcaYG$sdev)))

# We see that the first component explains more than half of the variance. But
# what does this component correspond to? The contribution of each dimension to
# this principal component is stored in the $rotation values. These are
# sometimes called the "eigengenes" for expression analysis - in analogy to the
# "eigenvalue" of a matrix equation.

plot(seq(0, 120, by = 5), pcaYG$rotation[,1], type = "b",
     ylim = c(-0.5, 0.5), col = myCol(length(pcaYG$sdev))[1],
     xlab = "time", ylab = "rotations")
# This tells us that the most common variation of the data corresponds to a
# global response in expression, right after the start of the experiment.

# What about the second component?
points(seq(0, 120, by = 5), pcaYG$rotation[,2], type = "b",
       col = myCol(length(pcaYG$sdev))[2])
# Here we see a cyclical response.

# Let's plot the first twelve components
opar <- par(mfrow = c(3,4), mar = c(0.5,0.5,0.5,0.5))
for (i in 1:12) {
    plot(seq(0, 120, by = 5), pcaYG$rotation[,i], type = "b",
         col = myCol(length(pcaYG$sdev))[i],
         xlab = "", ylab = "",
         axes = FALSE, frame.plot = TRUE)
    abline(h =  0, col = "#CCDDFF")
    abline(v = 60, col = "#CCDDFF")
}
par(opar)

# This is pretty interesting: we see global decaying response (PC1), cyclical
# responses that correspond to two cycles (PC2, 3, 4, and 6) - which are
# phase-shifted, another more global change (PC5), and components that may
# correspond to noise. Now we can ask: which genes correspond to these
# prototypes?

# The actual values (sometimes called "loadings") of the PCs for each gene are
# stored in the pcaYG$x matrix. For example, we can retrieve the ten genes that
# are most correlated with PC2. I use the order() function, which returns
# indices of values in a vector. This is a very useful function. Consider:
#
set.seed(123)
( x <- sample(1:30, 10) )
order(x, decreasing = TRUE)

# This means: the highest value is in position 9, the second highest in 5, the
# third highest in 4, etc. So I can use this, e.g. to pick the four highest
# values from x
  order(x, decreasing = TRUE)[1:4]
x[order(x, decreasing = TRUE)[1:4]]

# order() is similar to sort(), but different in that it returns indices, which
# we can use for selections.

# get the ten highest values from pcaYG$x[,2] ...
( sel <- order(pcaYG$x[ ,2], decreasing = TRUE)[1:10] )
pcaYG$x[sel, 2]

# ... list what genes these are ...
ygData[sel, c("stdName", "alias")]

# ... and plot their expression profiles
plot(seq(0, 120, by = 5), rep(0, 25), type = "n",
     ylim = c(-1.5, 1.5),
     xlab = "time", ylab = "log-ratio expression")
rect( 22.5, -2,  37.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 82.5, -2,  97.5, 2, col = "#dfeaf4", border = NA)   # G0
rect( 52.5, -2,  67.5, 2, col = "#f4dfdf", border = NA)   # G1
rect(112.5, -2, 122.5, 2, col = "#f4dfdf", border = NA)   # G1

for (i in 1:10) {
    points(seq(0, 120, by = 5), ygProfiles[sel[i], ], type = "b")
}

# Note that these are well-defined, cyclically expressed genes - however we
# would likely NOT have discovered them with our differential expression
# analysis, because they are NOT CONSTANT in the regions we have defined for
# grouping.

# This is an important difference: while we used intuition to define groups for
# discovery of differential expression, PCA gives us a way to look at the data
# with less bias, giving a greater voice to the dynamics within the data itself.
# We will explore next how we can apply regression analysis for alternative
# approaches to discover interesting genes.

# == TASK
#  -  Explore this some more and plot additional selections - e.g. the ten
#     highest and lowest values from other PCs.
#
# Hand in:
#   - Write code to plot the genes with the five highest and lowest loadings
#     of PC1, highest in black, and lowest in red.





# ==============================================================================
#      PART FIVE: READING GO AND GOA DATA
# ==============================================================================

# To be updated - later.





# [END]
