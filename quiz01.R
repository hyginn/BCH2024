# quiz01.R
#
# Purpose:  BCH2024 - First Quiz - subsetting
#
# Version: 1.0
#
# Date:    2017  01  13
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    Final
#
# TODO:
#
#
# ==============================================================================

# You have 20 minutes time.

# Let's create a small subset dataframe of expression values.

load("GSE3635.RData")
set.seed(1357)
selRows <- sample(1:length(featureNames(GSE3635)), 100)

mySet <- data.frame(exprs(GSE3635)[selRows, ])
rownames(mySet) <- featureNames(GSE3635)[selRows]
str(mySet)

# QUIZ QUESTION 1 (2 marks):

# How many genes in this selection come from yeast Chromosome 2?
# Write an expression to count them.
# Hint: remember that the chromosome number is encoded in the second letter
#       of the systematic name ...






# QUIZ QUESTION 2 (2 marks):
# write code to find all rows of mySet that contain an NA value. Assign the row
# numbers to the variable iNAs







# QUIZ QUESTION 3 (6 marks):

# write a loop that iterates over iNAs to replace all NA values in a row with a
# randomly chosen value from the existing values. Use the function sample(), and
# execute:
set.seed(12345)
# ... before calling sample for the first time.
# Be careful: there may be more than one NA value in a row.
# Hint: use your subsetting skills to create a vector that holds the
#       indices of elements that are NA, and another vector that holds
#       indices of elements that are not NA. Iterate over the former, and
#       sample from the latter.
# Hint: remember that you can replace a value in a dataframe like so:
#       mySet[i, j] <- x









# When you are done, copy your code, paste it into an e-mail and send it to me
# (NOT to the mailing list :-)



# [END]
