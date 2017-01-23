# quiz03-04.R
#
# Purpose:  BCH2024 - Final quiz: clustering
#
# Version: 1.0
#
# Date:    2017  01  23
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    Final
#
# TODO:
#
#
# ==============================================================================

# You have 30 minutes time.

# In the "Features" unit, we used PCA to distinguish crabs of different species and sexes, according to  morphometric data. Could we have distinguished the crabs by clustering?

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

# Lets make a combined factor column for species/sex ... we can conveniently use
# these values as class labels.

crabs$spsx <- as.numeric(as.factor(paste(crabs[, 1], crabs[, 2],sep=".")))



# QUIZ QUESTION :

# Use hierarchical clustering to separate the crabs into four groups. Note that
# you MUST subset the dataframe when you cluster, to select the appropriate
# columns. Use the morphometric columns only (obviously you should not cluster
# on class-labels! That would be pointless.)
#
# Hand in your code (minimal working code only)

# For each of the four groups, tabulate it's membership in terms of the four
# catagory labels we created for crabs$spsx. Hint: subset the rows by class,
# then use table() on the column of interest.
#
# Copy/paste the resulting four table() outputs.

# Was hierarchical clustering able to correctly categorize the crabs?
# (Write one sentence.)


# When you are done, copy your code and deliverables, paste it into an e-mail as
# plain text (no attached documents please) and send it to me (NOT to the
# mailing list :-)



# [END]
