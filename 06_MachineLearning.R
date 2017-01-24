# 02_MachineLearning.R
#
# Purpose:  BCH2024 - Machine Learning
#
# Version: 1.1
#
# Date:    2017  01  12
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.1    Added h2o for crabs example
# V 1.0    First version
#
# TODO:
#    https://rstudio.github.io/tensorflow/
#    https://www.r-bloggers.com/what-are-the-best-machine-learning-packages-in-r/
#    https://cran.r-project.org/web/views/MachineLearning.html
#    https://cran.r-project.org/web/packages/h2o/index.html

#
# ==============================================================================




# ==============================================================================
#      PART ONE: A SIMPLE EXAMPLE WITH CRABS DATA USING THE caret PACKAGE.
# ==============================================================================


library(MASS)
data(crabs)

head(crabs)
str(crabs)
# Two species: blue and orange
# Two sexes:   female and male
#
# Rows   1: 50  blue males
# Rows  51:100  blue females
# Rows 101:150  orange males
# Rows 151:200  orange females
#
# FL frontal lobe size (mm)
# RW rear width (mm)
# CL carapace length (mm)
# CW carapace width (mm)
# BD body depth (mm)

# Lets make a combined factor column for species/sex and change it to numbers,
# as class labels.

crabs$spsx <- as.factor(paste(crabs[, 1], crabs[, 2],sep="."))

# ... and drop the first three columns:
crabs <- crabs[,-(1:3)]
head(crabs)

# It is common for machine learning datasets to have the target categories
# (class labels) in the last column. Commonly, these are assumed to be factors.

# Can machine learning distinguish the crabs?
#
# We will use the caret package, which includes functions for cross-validation
# and tools for fitting. We load it with the (non-default) option of also
# loading "suggested" packages, which loads many, many packages that are useful
# for statistical learning and analysis.

if (! require(caret, quietly = TRUE)) {
    install.packages("caret", dependencies = c("Depends", "Suggests"))
    library(caret)
}

# Patience ...

# We will randomly remove 20% of the data into a separate "validation" dataset.
# Machine learning operates on training- and test-data to optimize its
# parameters - but after we have built our models, we would still like to
# validate whether our prediction also works on completely "unknown" data.

set.seed(112358)
N <- nrow(crabs)
sel <- sample(1:N, round(N * 0.2))
crabsVal <- crabs[sel, ]
crabs <- crabs[-sel, ]
str(crabs)


# Define control parameters:
# 10-fold cross validation
myControl <- trainControl(method="cv", number=10)

# Accuracy: this is our target metric - correctly predicted instances vs. total
# instances in the test set, in %.
myMetric <- "Accuracy"


# Try a number of "typical" ML algorithms

# === linear algorithms
#     lda (linear discriminant analysis)
set.seed(112358)
fit.lda <- train(spsx~., data=crabs, method="lda",
                 metric=myMetric, trControl=myControl)

# === nonlinear algorithms
# CART (Classification And Regression Trees)
set.seed(112358)
fit.cart <- train(spsx~., data=crabs, method="rpart",
                  metric=myMetric, trControl=myControl)

# kNN (k-Nearest Neighbours)
set.seed(112358)
fit.knn <- train(spsx~., data=crabs, method="knn",
                 metric=myMetric, trControl=myControl)

# === other algorithms
# SVM (Support Vector Machine)
set.seed(112358)
fit.svm <- train(spsx~., data=crabs, method="svmRadial",
                 metric=myMetric, trControl=myControl)

# Random Forest (often the "general purpose" method of first choice in ML)
set.seed(112358)
fit.rf <- train(spsx~., data=crabs, method="rf",
                metric=myMetric, trControl=myControl)


# == Evaluate
# summarize accuracy of models
myMLresults <- resamples(list(lda =  fit.lda,
                              cart = fit.cart,
                              knn =  fit.knn,
                              svm =  fit.svm,
                              rf =   fit.rf))
summary(myMLresults)

# The kappa statistic compares observed accuracy with expected accuracy, and
# thus takes into account that random chance may also give correct
# classifications. For a gentle, plain-english explanation see:
# http://stats.stackexchange.com/questions/82162/kappa-statistic-in-plain-english


dotplot(myMLresults)

# Linear discriminant analysis performed the best, both regarding accuracy and
# kappa statistic. Which method performs the best is obviously highly dependent
# on the data - and all of the methods allow optimizations ... this is a pretty
# big topic overall.

# How well did we do?
print(fit.lda)

# How can we use the classifier for predictions on unknown data? Remember that
# we had "unknown" data in our validation dataset. The functions called by caret
# are set up in similar ways as lm(), nls() or other modeling functions -
# specifically, they have a predict() method that allows to make predictions on
# new data with the larned parameters. Since we know the correct category labels
# in our validation set, we can easily check how often our prediction was right,
# or wrong: first we make predictions ...

myPredictions <- predict(fit.lda, crabsVal)

# ... and then we analyze them in a "confusion matrix": predictions vs. known
# class labels.
confusionMatrix(myPredictions, crabsVal$spsx)


# ==============================================================================
#      PART TWO: "INDUSTRY STRENGTH" ML WTIH h2o
# ==============================================================================

# h2o is a large, open platform for data science written in Java. After
# installing the package, an instance of h2o will run as a server for analysis
# and allow the R h2o package functions to interact with it. Installation from
# CRAN should be straightforward - even though the CRAN package has no actual
# h2o code, the required java "jar" file will be downloaded when the h2o.init()
# is called for the first time.
#

if (! require(h2o, quietly = TRUE)) {
    install.packages("h2o")
    library(h2o)
}

H2O <- h2o.init()

# Prepare data again, to be sure ...
data(crabs)
crabs$spsx <- as.factor(paste(crabs[, 1], crabs[, 2],sep="."))
crabs <- crabs[,-(1:3)]
set.seed(112358)
N <- nrow(crabs)
sel <- sample(1:N, round(N * 0.2))
crabsVal <- crabs[sel, ]
crabs <- crabs[-sel, ]

# Prepare our dataset for h2o in h2o's .hex (hexadecimal) format:
crabs.hex <- as.h2o(crabs)
str(crabs.hex)

# Let's run a "Deep Neural Network" model (cf.
# https://en.wikipedia.org/wiki/Deep_learning for the concepts and vocabulary,
# also see http://docs.h2o.ai/h2o/latest-stable/h2o-docs/glossary.html for an
# h2o glossary) out of the box with all-default parameters:
( fit.h2o <- h2o.deeplearning(x = 1:5, y = 6, training_frame = crabs.hex) )

( myH2OPred <- h2o.predict(fit.h2o, as.h2o(crabsVal)) )
h2o.confusionMatrix(fit.h2o, as.h2o(crabsVal))

# This result is not very impressive - 21/40 errors! cf. our lda result:
confusionMatrix(myPredictions, crabsVal$spsx)
# ... with only 4/40 errors.

# Try to improve these results by tuning the parameters: ten-fold
# cross-validation (default is none), three hidden layers of smaller size
# (default is c(200, 200)), set activation function to tanh (sigmoidal) (default
# is "Rectifier"), and use LOTS of iterations (default is only 10):
fit.h2o.2 <- h2o.deeplearning(
    x = 1:5,
    y = 6,
    training_frame = crabs.hex,
    hidden = c(8, 8, 8),
    activation = "Tanh",
    epochs = 5000,
    nfolds = 10,
    fold_assignment = "AUTO"
)
h2o.confusionMatrix(fit.h2o.2, as.h2o(crabsVal))

# Impressively, the result is now almost perfect (only 1/40 errors). It's
# worthwhile to play around with the parameters aned see how they influence
# processing time and accuracy. You will find that not always will more layers,
# more nodes, more iterations lead to better results. You can also get a sense
# that you can burn A LOT of processing power building these models - but the
# results can also be very, very good. And that's the whole point, after all.

# Before you leave, don't forget to shut down the h2o server instance again, or
# it will keep on running in the background!
h2o.shutdown(prompt=FALSE)


# ==============================================================================
#      PART THREE: ML OF EXPRESSION DATA
# ==============================================================================

# Try to distinguish Mbp1 target genes from Swi4 target genes - is that
# possible?

# ... TBC
#
#



# [END]
