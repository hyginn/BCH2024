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




# [END]
