# 05_Clustering.R
#
# Purpose:  BCH2024 - Clustering
#
# Version: 1.0
#
# Date:    2017  01  20
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First final version
#
# TODO:
#
#
# ==============================================================================

# Clustering cell-cycle genes


# ==============================================================================
#      PART ONE: ....
# ==============================================================================

load("ygData.RData")
load("ygProfiles.RData")
load("nlsParams.RData")

init()

# subset a selection of 260 genes with highest A and cor and f within 30% of 1
selCyc <- which(nlsParams$cor > 0.90  &
                    nlsParams$A   > 0.2 &
                    nlsParams$f > 0.75 &
                    nlsParams$f < 1.25)

# Let's check our choices ...
myKey <- ""
while (myKey != "x") {
    iRow <- selCyc[sample(1:length(selCyc), 1)]
    checkFit(iRow, bestFit2(ygProfiles[iRow, ]))
    myKey <- readline(prompt = "<enter> for next profile, <x> to abort >")
}

# subset for convenience

dat <- ygProfiles[selCyc, ]

# Clustering attempts to find which profiles are ore similar to each other than
# to other groups of profiles.
#
# ==================================================
# First explorations: Heatmap
# ==================================================
#
# Heatmaps are a staple of gene expression analysis.
# You can tweak many of the parameters, but for a first look
# we'll just heatmap the data with default parameters.

# This is a standard view that can be applied to all manners
# of multidimensional data, not just genes.
heatmap(dat)

# Just for illustration and readability let's map only
# every fifth gene
heatmap(dat[seq(1, nrow(dat), by=5), ])

# What's the actual range of values?
range(dat[,1])

# Study the heatmap, and consider what it tells you.
# For example, there seem to be genes that are high at t15 and t10
# but low at t45, t50 ...
set1 <- c("YER190W", "YIL177C", "YNL209W")
# ... and there are genes for which the inverse is true:
set2 <- c("YML033W", "YCL012W", "YMR215W")

# We can use a "parallel coordinates" plot - matplot()
# to look at the actual expression levels. Note that
# matplot expects the values column-wise ordered, thus
# we have to transpose - t() - the data!
matplot(t(dat[set1,]),
        type="l", lwd=2, col="skyblue", lty=1,
        ylim = c(-0.4, 0.4),
        xlab="time", ylab="expression")

# Then we can use lines() to superimpose the genes for set2.
# No transpose here :-)
for (i in 1:length(set2)) {
    lines(dat[set2[i], ], type="l", lwd=2, col="firebrick")
}

# Indeed, these genes - visibly different in the heatmap
# are mutualy similar in their expression profiles and different

# ==================================================
# Hierarchical clustering
# ==================================================
#
# Hierarchical clustering is probably the most basic technique.
# The dendrograms on the rows and columns of the heatmap
# were created by hierarchical clustering.

# For hierarchical clustering, first we need to produce
# a distance table. There are many ways to define distances
# let's just go with the default: "Euclidian distance".
distDat <-dist(dat)

# Then we use the clustering distance matrix to produce a
# dendrogram in which the most similar genes are connected, and then
# similar genes or connected groups are added. There are
# several ways to define "most-similar", lets just go with the
# default for now: "complete linkage" hierarchical clustering
hc <- hclust(distDat)

plot(hc)

# Not bad. But do note that both distance as well as clustering
# method matter, and there is not really a "best" way that
# works for all data. You'll need to explore: what you are looking for
# is a distance metric that gives the clearest block structure.

if (!require(minerva, quietly=TRUE)) {
    install.packages("minerva")
    library(minerva)
}
?mine


dEu <- function(x) dist(x, method="euclidian")
heatmap(dat, distfun = dEu)

dCan <- function(x) dist(x, method="canberra")
heatmap(dat, distfun = dCan)

dMax <- function(x) dist(x, method="maximum")
heatmap(dat, distfun = dMax)

dMink <- function(x) dist(x, method="minkowski")
heatmap(dat, distfun = dMink)

# You are not confined to the default distance functions, it
# is quite straightforward to define your own, for example
# using correlation properties. Here is a distance function
# defined as 1 - abs(pearson correlation)...

dCor <- function(x) as.dist(1 - abs(cor(t(x))))
heatmap(dat, distfun = dCor)

# ... and a similar one with the maximum information
# coefficient (MIC) as implemented in the minerva
# package
dMIC <- function(x) as.dist(mine(t(x))$MIC)
heatmap(dat, distfun = dMIC)

# I would choose dMax, but the differences are not
# really obvious ...

hc <- hclust(dMax(dat))

# Back to our original dendrogram ...
plot(hc)

# To get clusters from a dendrogram, we need to "cut" it at some
# level. The tree then falls apart into sub-trees and each of these
# is one "cluster"...

# Draw rectangles at different cut-levels, to give the desired number
# of clusters.
rect.hclust(hc,k=2)
rect.hclust(hc,k=5)
rect.hclust(hc,k=10)
rect.hclust(hc,k=20)
rect.hclust(hc,k=50)

# Now retrieve the actual indices and use them to generate
# parallel coordinate plots.

class <-cutree(hc, k = 20)

# Explain the output...
class

# The table() function allows us to count the number of
# occurences in each class ...
table(class)
sort(table(class))

# Let's plot the four largest classes (in parallel, into the same window)
# Look at this carefully. See how the selection statement on class
# generates a logical vector: TRUE in all rows for which the statement is true,
# and how this is used to select the rows of dat that we want to plot ...

oPar <- par(mfrow=c(2,2))
matplot(t(dat[class==10,]),type="l", xlab="time",ylab="log-ratio expression")
matplot(t(dat[class==2,]),type="l", xlab="time",ylab="log-ratio expression")
matplot(t(dat[class==11,]),type="l", xlab="time",ylab="log-ratio expression")
matplot(t(dat[class==4,]),type="l", xlab="time",ylab="log-ratio expression")
par(oPar)


# As an alternative, try Wards- linkage clustering (and read up on the
# options: single-, complete- and average-linkage clustering)
hc.ward <-hclust(distDat, method = "ward", members=NULL)

plot(hc.ward)

# draw rectangles
rect.hclust(hc.ward,k=9)

# This looks reasonable ...
# Now retrieve the actual indices and use them to generate
# paralell coordinate plots.

class.ward<-cutree(hc.ward, k = 9)
sort(table(class.ward))

# get some nice colors
if (!require(RColorBrewer, quietly=TRUE)) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
}

# what spectra are there in the package .. ?
display.brewer.all()

niceCols <- brewer.pal(9, "Paired")

oPar <- par(mfrow=c(3,3))
for (i in 1:9) {
    matplot(t(dat[class.ward == i,]),
            type="l", col=niceCols[i],
            xlab="time",ylab="log expression value")
}
par(oPar)

# trying the same with 1-cor(dat)

dCor <- as.dist((1 - cor(t(dat)))/2)
hCor <- hclust(dCor, method = "ward")
plot(hCor)


classCor <-cutree(hCor, k = 9)
sort(table(classCor))

niceCols <- brewer.pal(9, "Paired")

oPar <- par(mfrow=c(3,3))
for (i in 1:9) {
    matplot(t(dat[classCor == i,]),
            type="l", col=niceCols[i],
            xlab="time",ylab="log - ratio expression")
}
par(oPar)


# While this may be aesthetically somewhat satisfactory, it is clear that
# the clusters are not homogenous as we might need them for biological
# interpretation. This is a general problem with clustering methods that
# fix the number of cluster centres either directly as in Kmeans (see
# below), or indirectly by cutting trees at a fixed level. It is also
# a problem with the data, where differences in absolute values might
# override separation into clusters that might better be defined in terms
# of relative values.

# Here is a package that adresses the dynamic range problem.
# Read about it here: http://cran.r-project.org/web/packages/dynamicTreeCut/dynamicTreeCut.pdf
if (!require(dynamicTreeCut, quietly=TRUE)) {
    install.packages("dynamicTreeCut")
    library(dynamicTreeCut)
}

class.dynamic <- cutreeDynamic(dendro = hc.ward, distM = as.matrix(distDat),
                               cutHeight=100)

niceCols <- brewer.pal(8, "Spectral")

oPar <- par(mfrow=c(3,3))
for (i in 1:8) {
    matplot(t(dat[class.dynamic == i,]),
            type="l",
            col=niceCols[i],
            xlab="time",
            ylab="log expression value")
}
par(oPar)

# One thing our clustering algorithms do is to pull apart profiles
# that have similar shape, but different absolute levels. This is
# because we have not normalized our data. Let's thus try
# clustering merely based on profile shape, i.e.
# relative expression levels, by scaling all rows between zero
# and one.

datNorm <- t(apply(dat, 1, function(x)(x-min(x))/(max(x)-min(x))))
distDatNorm <-dist(datNorm)

hc.Norm <-hclust(distDatNorm)

class.dynamic <- cutreeDynamic(dendro = hc.Norm, distM = as.matrix(distDatNorm), cutHeight=15)

niceCols <- brewer.pal(6, "Spectral")

oPar <- par(mfrow=c(3,2))
for (i in 1:6) {
    matplot(t(datNorm[class.dynamic == i,]),
            type="l",
            col=niceCols[i],
            xlab="time",
            ylab="log expression value")
}
par(oPar)

# With hierarchical clustering, this is probably as good
# as we can get - the clusters are of reasonable size -
# but from a biological point of view one would argue
# that several of them are not really different.
# ==================================================
# Partitioning clustering
# ==================================================

# === K-means ======================================

# K-means clusters by assigning elements to a fixed
# number of cluster centres, so that similarity
# within a cluster is maximized.

?kmeans

k <- 4
cl<-kmeans(dat, k)

niceCols <- brewer.pal(k, "Spectral")

plot(dat[,"t.30"],dat[,"t.60"],col=niceCols[cl$cluster])
points(cl$centers, col = niceCols[1:k], pch = 8, cex=2)

# But: be aware ...
# ... K-means does not guarantee a globally optimal solution,
# merely a locally converged one.

# === K-medoids ======================================

# load library "cluster" for K-medoid partitioning
if (!require(cluster, quietly=TRUE)) {
    install.packages("cluster")
    library(cluster)
}
set.seed(112358)
k <- 4
cl<-pam(dat, 4)
plot(dat[,"t.30"],dat[,"t.60"], col=niceCols[cl$cluster])
# plot(cl) # shows boundary and silhouette plots

# ==================================================
# Affinity propagation clustering
# ==================================================
# Based on B. J. Frey and D. Dueck. Clustering by
# passing messages between data points.
# Science, 315(5814):972–976, 2007

if (!require(apcluster, quietly=TRUE)) {
    install.packages("apcluster")
    library(apcluster)
}

apRes <- apcluster(negDistMat(r=2), dat)
apRes

heatmap(apRes)

# try this on the normalized data
apRes <- apcluster(negDistMat(r=2), datNorm)
heatmap(apRes)

# The clear and pronounced block structure shows that this
# is a successful clustering...

length(apRes)
cutree(apRes)

oPar <- par(mfrow=c(2,2))
matplot(t(datNorm[unlist(apRes[2]),]),type="l",xlab="time",
        ylab="log-ratio expression")
matplot(t(datNorm[unlist(apRes[3]),]),type="l",xlab="time",
        ylab="log-ratio expression")
matplot(t(datNorm[unlist(apRes[8]),]),type="l",xlab="time",
        ylab="log-ratio expression")
matplot(t(datNorm[unlist(apRes[12]),]),type="l",xlab="time",
        ylab="log-ratio expression")
par(oPar)

# ==================================================
# Cluster quality metrics
# ==================================================
# So .. which should we use?
# I don't think that there is an obvious biological
# criterium to decide. Typically we should take some
# orthogonal information (e.g. shared transcription
# factor binding site), so if functionally related
# genes end up in similar clusters, and judge our
# cluster success based on its biological predictive
# value.
#
# But we can certainly say something about whether
# our clusters look good in a mathematical sense.
if (!require(clValid, quietly=TRUE)) {
    install.packages("clValid")
    library(clValid)
}
if (!require(kohonen, quietly=TRUE)) {
    install.packages("kohonen")
    library(kohonen)
}
if (!require(mclust, quietly=TRUE)) {
    install.packages("mclust")
    library(mclust)
}
?clValid

# This is pretty nice: we _can_ use biological
# knowledge for validation...
# But for our example, we'll try internal validation
# on all available methods.

valClust <- clValid(dat,
                    nClust = 2:9,
                    clMethods = c("hierarchical",
                                  "kmeans",
                                  "diana",
                                  "fanny",
                                  "som",
                                  "pam",
                                  "sota",
                                  "clara",
                                  "model"),
                    validation = "internal")
summary(valClust)
plot(valClust)
vignette("clValid")

# Task:
# 1 - What appears to be the best clustering method?
# 2 - How can you actually apply it to the data?

# ==================================================
# t-stochastic Neighbour embedding
# ==================================================

# tsne - pioneered by Geoff Hinton - is an embedding
# procedure that guarantees that points in a high-
# dimensional feature space that are close together
# remain close together when projected into a low-
# dimensional space. It can give astoundingly good
# results that help figuring out the internal structure
# of a dataset.

if (!require(tsne, quietly=TRUE)) {
    install.packages("tsne")
    library(tsne)
}

# The clustering method uses a "callback" to  execute
# a plotting routine as it improves its embedding
# through many cycles.

# define a plotting routine
plotProgress <- function(x){
    plot(x, type='n');
    #	text(x, labels = rownames(dat), cex=0.5)
    points(x, pch=21, col="#6677FF", bg="firebrick")
}

set.seed(11235)
tsneDat <- tsne(dat, epoch_callback = plotProgress, perplexity = 10)

# As with the examples above: points that are close in projection space are
# similar in feature space.


# [END]
