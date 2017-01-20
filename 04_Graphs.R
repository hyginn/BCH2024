# 04_Graphs.R
#
# Purpose:  BCH2024 - Graphs
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

# Building graphs for analysis requires to define edges that connect nodes. In
# our domain, nodes are usually genes, and the edges we are interested in can
# cover many aspect of functional relationships: protein-protein interactions,
# collaboration in systems, regulation ... and much more.




# GO - semantic similarity
# collect GO features from semantic similarity to reference genes
#
# STRING
# Co-expression
# Graph clustering
# Annotation propagation

# Reload data we have created previously:
load("ygData.RData")
load("ygProfiles.RData")



# ==============================================================================
#      PART ONE: SEMANTIC SIMILARITY
# ==============================================================================

# Many functional aspects of proteins are documented in the Gene Ontology. This
# will give us an idea of what a gene is or does, but importantly we can also
# use it to determine how similar two genes are - this is called "semantic
# similarity".
#

# GOSemSim is an R-package in the bioconductor project. It is not installed via
# the usual install.packages() comand (via CRAN) but via an installation script
# that is run from the bioconductor Website.

source("http://bioconductor.org/biocLite.R")
biocLite("GOSemSim")

library(GOSemSim)

# This loads the library and starts the Bioconductor environment.
# You can get an overview of functions by executing ...
browseVignettes()
# ... which will open a listing in your Web browser. Open the
# introduction to GOSemSim PDF. As the introduction suggests,
# now is a good time to execute ...
?GOSemSim

# Before we do anything with the function, we need to access an appropriate
# BioConductor OrgDb object ...
biocLite("AnnotationHub")
library(AnnotationHub)

# Which keytypes are available in the Saccharomyces annotation database?
keytypes(org.Sc.sgd.db)

# "COMMON", "GENENAME" and "SGD" look like they could be what we have in our
# ygData object ... let's see some samples:
#
head(keys(org.Sc.sgd.db, keytype = "COMMON"))
head(keys(org.Sc.sgd.db, keytype = "GENENAME"))
head(keys(org.Sc.sgd.db, keytype = "SGD"))

# It seems the best mapping will be via SGD IDs (we don't have "COMMON" standard
# names for all of our genes.) So let's load and prepare the GO-data for
# gosemsim:
scGObp <- godata("org.Sc.sgd.db", keytype = "SGD", ont="BP")

# The simplest function is to measure the semantic similarity of two GO terms.
# For example, we could consider GO:0035019 (somatic stem cell maintenance),
# GO:0045454 (cell redox homeostasis), and GO:0009786 (regulation of asymmetric
# cell division). Lets calculate these similarities.

goSim("GO:0035019", "GO:0009786", semData = scGObp, measure="Wang")
goSim("GO:0035019", "GO:0045454", semData = scGObp, measure="Wang")
goSim("GO:0009786", "GO:0045454", semData = scGObp, measure="Wang")

# Fair enough. Three numbers. Clearly we would appreciate an idea of the values
# that high similarity and low similarity can take. But in any case -
# we are really less interested in the similarity of GO terms - these
# are a function of how the Ontology was constructed. We are more
# interested in the functional similarity of our genes, and these
# have a number of GO terms associated with them.

# GOSemSim provides the functions ...
?geneSim()
?mgeneSim()
# ... to compute these values. Refer to the vignette for details, in
# particular, consider how multiple GO terms are combined, and how to
# keep/drop evidence codes.
# Here is a pairwise similarity example: ...

ygData[281, ]
ygData[2501, ]

geneSim("S000000247", "S000003527",
        semData = scGObp, measure="Wang", combine="BMA")

# Another number. And the list of GO terms that were considered.

# ... TBC



# [END]
