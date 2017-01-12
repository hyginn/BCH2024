# 01_Data.R
#
# Purpose:  BCH2024 - Data
#
# Version: 1.0
#
# Date:    2016  12  30
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First final version
#
# TODO:
#
#
# ==============================================================================


# Subset genes
# Plot expression profiles
# What is a cell-cycle gene?
# Where can we find such anotations?



# ==============================================================================
#      PART ONE: SETUP
# ==============================================================================

# Load the bioconductor package installer
source("https://bioconductor.org/biocLite.R")


if (!require(Biobase, quietly=TRUE)) {
    biocLite("Biobase")
    library(Biobase)
}


if (!require(GEOquery, quietly=TRUE)) {
    biocLite("GEOquery")
    library(GEOquery)
}


# ==============================================================================
#      PART TWO: THE DATA
# ==============================================================================

# The R code below is adapted from the GEO2R scripts produced for GSE3635
# and GSE4987.

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Wed Jan 11 17:39:46 EST 2017


# Load series and platform data from GEO
gset <- getGEO("GSE3635", GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) {
    idx <- grep("GPL1914", attr(gset, "names"))
} else {
    idx <- 1
}

gset <- gset[[idx]]
GSE3635 <- gset
# save(GSE3635, file="GSE3635.RData")

# This is an "Expression Set" - cf.
# https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf

# What does this contain?
help("ExpressionSet-class")

# Print gset
gset

# Access contents via methods:
featureNames(gset)[1:20]   # rows
sampleNames(gset)[1:10]    # columns

# Access contents by subsetting:
( tmp <- gset[12:17, 1:6] )

# Access data
exprs(tmp)

# What are the data:
#  ... in each cell
#  ... in each column
#  ... in each row
#  Are values present for all columns?
#  What do the columns mean?

#  Are values present for each row?
#  Are all rows genes?
#  What identifiers are being used?
#     (cf. http://www.yeastgenome.org/help/community/nomenclature-conventions)
#  Are all rows/genes unique?
#  Are all yeast genes accounted for?

#     To answer this question, I have provided the file "SGD_features.tab" in
#     the data directory of this project. I have downloaded it from SGD
#     (http://www.yeastgenome.org/download-data/curation), for a description of
#     its contents see here:
file.edit("data/SGD_features.README.txt")
#     Note: the file as-is crashed RStudio due to an unbalanced quotation
#           mark. I have removed the alias (B") that contained this abomination,
#           by editing the file. Your version in the data directory can be
#           read.

#     Task: - read this into a data frame "SGD_features"
#           - remove unneeded columns
#           - make meaningful column names
#           - check that systematic names are unique
#           - use systematic names as row names
#           - save as RData file
#           - confirm: are all rows of the expression data set represented in
#                      the feature table?
#           - confirm: how many / which genes in the feature table do not
#                      have expression data?

#  How do we handle rows/columns that are missing or not unique?


# Overview of the data
boxplot(exprs(gset))

#  Boxplot for selected GEO samples
#  GEO2R version
#  Study this - what's going on?
#  group names for all samples in a series
gsms <- "0123450123450"
sml <- c()
for (i in 1:nchar(gsms)) {
    sml[i] <- substr(gsms,i,i)
}
sml <- paste("G", sml, sep="")  # set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("t0","t10","t20","t30","t40","t50")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5","#dff4e4","#f4dff4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE3635", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

# === Column-wise analyses
# Simple data descriptors
x <- exprs(gset)[ , 1]
mean(x)     # This needs fixing !
median(x)
IQR(x)
var(x)
sd(x)
summary(x)

# === Row-wise analyses
# Expression level plot for selected genes:
#
# Task: Plot expression for
# Mbp1, Swi6, Swi4, Nrm1, Cln1, Clb6, Act1, and Alg9



# == GSE4987: HIGHER RESOLUTION DATA, TECHNICAL REPLICATE
# https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE4987

GSE4987 <- getGEO("GSE4987", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(GSE4987) > 1) {
    idx <- grep("GPL1914", attr(GSE4987, "names"))
} else {
    idx <- 1
}
GSE4987 <- GSE4987[[idx]]
# save(GSE4987, file="GSE4987.RData")

# What is the relationship of these columns to the previous set?
# Can we/should we merge the two sets?



# ==============================================================================
#      PART TWO: SUBSETTING REVIEW
# ==============================================================================

# ===== Sample data ==================================================

# Let's start with a small datframe of synthetic data to go through the main
# principles of subsetting. The same principles apply to matrices and vetors -
# however, data frames are more flexible because their columns can contain data
# of different types (character, numeric, logical ...). Values in vectors and
# matrices must always have the same type.

# Imagine you are a naturalist that has collected some living things and keeps observations in a table ...

set.seed(112358)
N <- 10

dat <- data.frame(name = sample(LETTERS, N, replace = TRUE),
                  legs = sample(c(2 * (0:5), 100), N, replace = TRUE),
                  type = character(N),
                  matrix(rnorm(5 * N), ncol = 5),
                  stringsAsFactors=FALSE)

# Some auxiliary data ...
dict <- c("fish", "bird", "beast", "bug", "spider", "crab", "centipede")
names(dict) <- c(2 * (0:5), 100)
# ... to populate the >>type<< column:
dat$type <- dict[as.character(dat$legs)]

# If you already understand the expression above, you're doing pretty well with
# the topic of this tutorial. If you don't, don't worry - by the end of the
# tutorial you will. Now let's see what we have:

head(dat)
str(dat)

# Note that we have given names to some columns, but R made names for the five
# columns of random values that were created as a matrix. Let us look at the
# basic ways to subset such objects. Basically, all these methods work with the
# subsetting operator "[".

?"["


# ===== Subsetting by index ==========================================

# Elements can be uniquely identified by indices in the range of their length
# (for vectors), or their rows and columns (in dataframes and matrices). The
# order is row, column.

dat[2,3]   # one element
dat[2, ]   # empty columns: use all of them
dat[ , 3]  # empty rows, use all of them

# If you want a particular set of row and columns, pass a vector of positive
# integers.
dat[c(2, 3), c(1, 2, 3)]

# Any function that returns a vector of integers can be used. Most frequently we
# use the range operator ":" . Retrieving ranges of rows and/or columns from a
# matrix or data frame is also called "slicing".

dat[1:4, 1:3]
dat[4:1, 1:3]   # same in reverse order
dat[seq(2, N, by=2), ]   # even rows

# But we can do more interesting things, since the indices don't have to be
# unique, or in any order:

dat[c(1, 1, 1, 2, 2, 3), 1:3]

# In particular we can select random subsets...
dat[sample(1:N, 3), 1:3]
dat[sample(1:N, 3), 1:3]
dat[sample(1:N, 3), 1:3]

# ... or sort the dataframe. Sorting requires the order() function, not sort().

sort(dat[ , 2])    # ... gives us the sorted values

order(dat[ , 2])   # ... tells us in which row the sotrted values are
dat[order(dat[ , 2]), 1:3]  # ordered by number of legs
dat[order(dat[ , 1]), 1:3]  # ordered by lexical order of names

# Note: I am indenting expressions so you can first evaluate the expressions
# individually, then see how they fit into the brackets to subset the data.


# ==== Negative indices

# If you specify a negative index, that element is excluded.

dat[-1, 1:3]   # not the first row
dat[-N, 1:3]   # not the last row

dat[-1:-3, 1:3]
dat[-(1:3), 1:3]  # same effect



# ===== Subsetting by boolean ========================================

# Instead of indices, we can specify sets of rows or columns by boolean values
# (type: logical): TRUE or FALSE. If we place a vector of logicals into the
# square brackets, only the rows resp. columns for which the expression is TRUE
# are returned.

dat[1:3, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)]

# You need to take care that the number of elements exactly matches the number
# of rows or columns, otherwise "vector recycling" will apply and this is
# probably unintended. Thus explicitly specifying a boolean selection like above
# is not all that useful. However, many R functions are "vectorized" and
# applying a logical expression or function to an entire column gives a vector
# of TRUE and FALSE values of the same length. If we place this vector into the
# square brackets, only the rows resp. columns for which the expression is TRUE
# are returned.

dat[ , 2]
dat[ , 2] > 4          # See how this creates a vector
dat[dat[ , 2] > 4, 1:3]

# Expressions can be combined with the "&" (and) and the "|" (or) operators.

dat[ , 4] > 0
dat[ , 5] < 0
dat[ , 4] > 0 & dat[ , 5] < 0
dat[dat[ , 4] > 0 & dat[ , 5] < 0, ]

# In this context, the any() and all() functions may be useful. But take care -
# you can't simply apply them to a range of columns: that would apply the
# condition to all elements of a selection at once. You need to use the apply()
# function to first return a vector. apply()'s second argument switches between
# row-wise and column-wise evaluation. Here, 1 means operate on rows.

apply(dat[ , 4:8], 1, max)           # row-wise, fetch the maximum
apply(dat[ , 4:8], 1, max) > 1       # max() > 1 ?
dat[apply(dat[ , 4:8], 1, max) > 1, ]

# To use any() and all(), we define our own function.

myF <- function(x){any(x > 1.5)}
myF(dat[3, 4:8])

apply(dat[ , 4:8], 1, myF)
dat[apply(dat[ , 4:8], 1, myF), ]

# But we can also write the definition "in place"...
apply(dat[ , 4:8], 1, function(x){all(x < 0.5)})
#-------------------------
dat[apply(dat[ , 4:8], 1, function(x){all(x < 0.5)}), ]

# ====== String matching expressions

# The function grep(), and the %in% operator can be used to subset via string
# matching:

grep("r", dat[ , 3])          # types that contain "r"
dat[grep("r", dat[ , 3]), 1:3]

grep("^c", dat[ , 3])         # types that begin with "c"
dat[grep("^c", dat[ , 3]), 1:3]


scary <- c("spider", "centipede")
dat[ , 3] %in% scary
dat[dat[ , 3] %in% scary, 1:3]



# ===== Subsetting by name ===========================================

# If rownames and/or columnnames have been defined, we can use these for
# selection. If not defined, they default to the row/column numbers as character
# strings(!).

rownames(dat)  # the row numbers, but note that they are strings!
colnames(dat)  # the ones we have defined

# If we place a string or a vector of strings into the brackets, R matches the
# corresponding row/ column names:

dat[1:5, "name"]
dat[1:5, c("name", "legs")]
dat[1:5, "eyes"]   # error: that name does not exist

# We can combine the techniques e.g. to flexibly select columns. Here we select
# the X1 to X5 columns:

colnames(dat)
grep("^X", colnames(dat))
colnames(dat)[grep("^X", colnames(dat))]
dat[1:3, colnames(dat)[grep("^X", colnames(dat))]]

# This is very useful when the exact position of columns may have changed during
# the analysis. Actually, rows and columns should really never be selected by
# number even though we have done so above. Such numbers are "magic numbers" and
# code that relies on such magic numbers is heard to read and very hard to
# maintain. It is always better to expose the logic with which your columns are
# selected and to make the selection explicit and robust. An exception may be
# when you need a slice of the data for testing purposes, but even then it may
# be preferrable to use the head() or tail() functions.

# ===== The "$" operator

# The "$" operator returns a single column as a vector. It is not strictly
# necessary - the column can just as well be named in quotation marks within the
# brackets - but I think it makes for more readable code.
dat[1:3, "legs"]
dat$legs[1:3]    # same result. This is the preferred version.
dat$"legs"[1:3]  # works, but isn't necessary
dat[1:3, legs]   # this returns an error - hopefully; but if for any
# reason the object DOES exist, you'll get an un-
# expected result. Know when to quote!


# Three more functions that I use all the time for data manipulation:
?which
?unique
?duplicated



# ==============================================================================
#      PART THREE: ANNOTATION AND FEATURE DATAMODEL
# ==============================================================================

# TASKS:
# GOA features:
#    - Find a source for yeast gene GOAs
#    - Download the data
#    - Design a datamodel to store it.
#        (How? What is the purpose ...)
#    -
# Gene features:
#    - design a datamodel to integrate expression data, functional
#      annotations, and relationships.


# [END]
