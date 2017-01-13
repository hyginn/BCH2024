# 02_Features.R
#
# Purpose:  BCH2024 - Features
#
# Version: 1.0
#
# Date:    2017  01  12
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
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
# Does the dye-swap reveal a bias?
# (If it does - we should remove the bias ...)
#
# Impute missing values. How?
#


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
