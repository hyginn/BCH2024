# ================================================== #
# Boris Steipe <boris.steipe@utoronto.ca>            #
#                                                    #
#                                                    #
# Plotting reference                                 #
# Colour                                             #
#                                                    #
# (c) Boris Steipe  2011-2017                        #
# ================================================== #

# ==================================================
# 2 - Colors
# ==================================================

# Colors can be specified by number, by name, as hex-triplets
# as rgb or hsv values, and through color palettes.

# Colors by number =================================

# You can specify colors by using a small integer. What this color will
# translate to depends on the current value of palette(). A new session of R and
# RStudio has the palette:
# palette(c("black", "red", "green3", "blue", "cyan", "magenta",
#            "yellow", "gray"))
# and these are really ugly color combinations. If you have
# previously run GEO2R code, this may hace changed the palette, which may now be
# palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5",
#           "#dff4e4","#f4dff4","#AABBCC"))
# i.e. the pastel colors used at the GEO2R Website...

# any case - we can visualize the colors ...
barplot(rep(1,9), col = 0:8, axes = FALSE)

# It's usualy a good idea to specify color values explicitly.



# Colors by name ===================================
# You may have noticed that "red", "green", and "blue"
# work for the col=... parameter, but you probably would
# not have imagined that "peachpuff", "firebrick" and
# "goldenrod" are valid as well. In fact, there are
# 657 named colors in R. Access them all by typing:
colors()

pie(c(1,1,2,3,5,8),
    col=c("peachpuff", "firebrick", "goldenrod", "papayawhip",
          "whitesmoke", "moccasin")
    )

# Read more about named colors (and related topics) at
# http://research.stowers-institute.org/efg/R/Color/Chart/

# Colors as hex-triplets ===========================
# Hex triplets in R work exactly as in HTML: a triplet of
# RGB values in two-digit hexadecimal representation. The
# first two digits specify the red value, the second two
# are for green, then blue. R accepts a fourth pair of
# digits to optionally specify the transparency, the
# semantics of the code is thus "#RRGGBB" or "#RRGGBBAA".
# Read more e.g. at http://en.wikipedia.org/wiki/Web_colors

# The function col2rgb() converts color names to rgb values ...
col2rgb("violetred")

# ... and rgb() converts rgb values to hex-code:
rgb(1, 0.5, 0.23)

# Unfortunately the output of col2rgb does not quite match
# rgb(). col2rgb creates rows with values between 0 and 255,
# and rgb by default expects columns with intensities from
# 0 to 1, you have to transpose and divide.
rgb(t(col2rgb("red"))/255)        # "#FF0000"
rgb(t(col2rgb("peachpuff"))/255)  # "#FFDAB9"

# but rgb() is very useful, it takes numeric values from 0 to 1
# thus it can be used to compute color ramps and palettes -
# although we'll meet a more sophisticated function for that below,


# There are many tools on the Web that help to generate
# pleasing palettes.

# Here is an example -"Creative Cloud"- taken from
#    https://kuler.adobe.com/

CC <- c("#011640", "#024059", "#F2F0D0", "#BE6C5C", "#8C3037" )
hist(rnorm(1000), breaks=20 , col=CC)

# Since colors are specified as quartets - (R,G,B,Transparency) - we can adjust
# the transparency with values other than "FF" (the default) for very crowded
# plots, or for creating overlays.

x <- rnorm(2000)
y <- x^3 * 0.25 + rnorm(2000, 0, 0.75)
# compare:
plot(x,y, pch=19, col="#EE3A8C")
plot(x,y, pch=19, col="#EE3A8C12") # Alpha at ~ 10%

# or with multiple overlays of varying size ...
plot(x,y, pch=16, col="#AA330009")
points(x,y, pch=19, cex=2, col="#44558803")
points(x,y, pch=20, cex=0.5, col="#EE3A8C08")

# Color palettes ===================================
# R has several inbuilt color palettes, or you can build your own.

# Inbuilt palettes =================================
?rainbow
# view the palettes
opar <- par(mfrow=c(3,2))
n <- 20
sq <- rep(1, n)
barplot(sq, col=rainbow(n),        axes=F, main="rainbow(n)")
barplot(sq, col=cm.colors(n),      axes=F, main="cm.colors(n)")
barplot(sq, col=topo.colors(n),    axes=F, main="topo.colors(n)")
barplot(sq, col=terrain.colors(n), axes=F, main="terrain.colors(n)")
barplot(sq, col=heat.colors(n),    axes=F, main="heat.colors(n)")
par(opar)


# Perceptually equally spaced color ================
# If you consider the rainbow spectrum ...
n <- 20
barplot(rep(1, n), col=rainbow(n),        axes=F, main="rainbow(n)")
# ... you will notice that some of the bars appear very similar.
# This may be individually different, but I feel that the four central
# greens, the two last blues and the two first reds are almost
# impossible to distinguish. The solution to this problem is to spread
# the colors non-linearly, so that adjacent colors are perceptually
# equally different. At the same time, the perceived lightness of
# the colors should be constant. Here is an attempt:

eqPal <- colorRampPalette(c("#f2003c",  # red
                            "#F0A200",  # orange
                            "#f0ea00",  # yellow
                            "#62C923",  # green
                            "#0A9A9B",  # blue
                            "#1958C3",  # indigo
                            "#8000D3",  # violet
                            "#D0007F"), # red
                          space="Lab",
                          interpolate="linear")


# Such a perceptually tuned spectrum is quite a bit more pleasing
# than one that is computed from extrapolating between rgb values.
# Adjacent colors are easier to distinguish, in particular hues that
# are close to the primary reds, greens and blues, and the lightness
# of the colors is more balanced, if you compare e.g. blue and yellow.

n <- 15
sq <- rep(1, n)
oPar <- par(mfrow = c(2,1))
barplot(sq, col = rainbow(n), axes=F, main="rainbow")
barplot(sq, col = eqPal(n),   axes=F, main="equalized palette")
par(oPar)

# plotting many bars shows that the range of colors we can distinguish
# is now more homogenous.

n <- 50
sq <- rep(1, n)
oPar <- par(mfrow=c(2,1))
barplot(sq, col=rainbow(n), axes=F, border= NA, main="rainbow")
barplot(sq, col=eqPal(n),   axes=F, border= NA, main="equalized palette")
par(oPar)


# Examples: coloring points by value ==================

# Coloring points by value can identify qualitative and quantitative
# differences. For quantitative differences, we employ gradients.

# In this example, we generate random points and calculate a "density"
# at each point. The we plot each point and color it according
# to its density.

n <- 1000
x <- rnorm(n)
y <- x^3 * 0.25 + rnorm(n, sd=0.75)
z <- rep(0, n) # initialize z
for (i in 1:n) { # calculate a density from the proximity of other points
	dx <- x-x[i]; dx <- dx * dx # square of the distance in x
	dy <- y-y[i]; dy <- dy * dy # square of the distance in y
	d <- dx + dy                # square of the sum
	d <- d[-i]                  # remove the self-element
	z[i] <- 1/sum(d)            # let density decay with 1/r^2
}
z <- z/max(z) - 0.00001 # normalize, but keep the max below 1.0
# now map each of the values of z into an interval of the palette
n <- 20                   # number of intervals
z <- floor(z * n) + 1     # convert to integers
pal <- rainbow(n)         # get a vector of colors from the palette
cz <- pal[z]              # apply the color for each density value
plot(x,y, col=cz, pch=16) # plot

# use a different palette
pal <- heat.colors(n)
cz <- pal[z]
plot(x,y, col=cz, pch=16)

# Custom palettes =================================
# "Cold" values should really be black, not red. Lets define a
# custom palette: colorRampPalette() is a function that returns
# a function. The returned function can be used to calculate
# a palette, a trajectory along a number of waypoints in
# colorspace.

# And let's use the actual values of glowing steel - as taken
# from here for example:
#     http://i.imgur.com/qEk0QJf.jpg
# and pick colors with this tool:
#     http://imagecolorpicker.com/en

fc <- colorRampPalette(c("#301C25", "#F53A43", "#FFC754", "#FEFDDE"), bias=0.5)
# assigns a function to fc

fc # look at the function

fc(n) # use the function to get n values
pal <- fc(n)
cz <- pal[z] # assign the colors according to z
plot(x,y, col=cz, pch=16)

# Useful palettes have also been described specifically
# for cartography. http://colorbrewer2.org/ has palettes for seqential
# and qualitative diferences, and options for colorblind-safe and
# photocopy friendly palettes. You can use them via an R package:

if (! require(RColorBrewer, quietly = TRUE)) {
    install.packages("RColorBrewer")
    library(RColorBrewer)
}
display.brewer.all()

# In the following example, we apply a Brewer palette to a Voronoi tesselation
# of a point set.

if (! require(deldir, quietly = TRUE)) {
    install.packages("deldir")
    library(deldir)
}

# Create a point set along a logarithmic spiral, with a bit
# of added noise.
li <- 0.1
n <- 45
dl <- 1.06
ncirc <- 13
da <- (2*pi)/ncirc
fnoise <-0.13

# create a matrix of points
x <- numeric(n)
x <- cbind(x, numeric(n))

set.seed(16180)
for (i in 1:n) {
	l <- li * (dl^(i-1))
	x[i,1] <- (l+(rnorm(1)*fnoise*l)) * cos((i-1)*da)
	x[i,2] <- (l+(rnorm(1)*fnoise*l)) * sin((i-1)*da)
}
plot(x[,1], x[,2])
ts <- deldir(x[,1], x[,2])       # calculate the tesselation
tl <- tile.list(ts)      # calculate the list of tiles
plot.tile.list(tl)       # plot it

# Let's color the cells by distance from a defined point
# using a Brewer palette
points(x[25,1], x[25,2], pch=20, col="red")   # pick a point

vec <- c(x[25,1], x[25,2]) # define a point

# define a function for Euclidian distance
vDist <- function(x,v) { sqrt(sum((x-v)^2)) }  # calculates Euclidian distance
d <- apply(x,1,vDist, v=vec)                   # apply this to the point set

dCol <- floor(((d-min(d))/(max(d)-min(d)) * 10)) + 1 # map d into 10 intervals
dCol[which(dCol>10)] <- 10                           # demote the largest one

pal <- brewer.pal(10, "RdGy")   # create the palette

# plot the tesselation, color by palette
plot.tile.list(tl, fillcol = pal[dCol], cex=0.8, pch=20, col.pts="slategrey")


# ==============================================================================
#    COLOUR PALETTE RESOURCE
# ==============================================================================

# Below is a large number of different palettes, as a resource that you can
# choose from, for all your sequential, divergent, or categorical colour needs:



# Color palettes can be extracted from images, for example like here:
#    http://www.pictaculous.com/

# Or you can use the "eye-dropper" tool of standard image editing
# software to access to color value of a particular pixel, or online here ...
#    http://html-color-codes.info/colors-from-image/
#    http://imagecolorpicker.com/

# Try your own - here are some examples.
# Glowing steel: a value going from small to large
fCol <- colorRampPalette(c("#301C25", "#F53A43", "#FFC754", "#FEFDDE"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# Fire: a value going from calm to urgent
fCol <- colorRampPalette(c("#004466", "#000000", "#9F2704", "#DE4F03", "#FFFF0F"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# Red/Blue: A value going from pos to negative: centre on hue, 85% sat;
# down to 50%b/95%s towards the ends, then 33%s/100%b towards the middle
fCol <- colorRampPalette(c("#802906", "#ed5d23", "#ffc3ab", "#FFFFFF", "#ace5ff", "#23b1f3", "#005780"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# subdued Red/Blue:
fCol <- colorRampPalette(c("#fa4c4b", "#c55a67", "#a16578", "#653e5e", "#71708f", "#4a799f"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# Grey/Orange: A value going from neutral to highlight
fCol <- colorRampPalette(c("#333333", "#4d4d4d", "#666666", "#808080", "#999999", "#b3b3b3", "#cccccc", "#fc9973", "#ff0d35"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# light/dark/cardinal red: A value going from neutral to highlight
fCol <- colorRampPalette(c("#eeeeec", "#bbbbbb", "#99999b", "#666669", "#44444c","#ac3747", "#be0a23"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")



# Active Lava: Colours taken from video stills of Lava flow - for
# highlighting extreme values
fCol <- colorRampPalette(c("#B3BAB7", "#A4ABA9", "#959C9B", "#868D8D", "#787E80",
                           "#696F72", "#5A6064", "#4C5157", "#3D4249", "#2E333B",
                           "#20242E", "#5e4157", "#c2575e", "#f8451c", "#ffb57d"))
n <- 100
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# Scales based on Color Brewer scales
# Greens
fCol <- colorRampPalette(c(
"#16371f",
"#1a5b2d",
"#00823d",
"#20b874",
"#71c69e",
"#acd98e",
"#fcf6a0",
"#fdfbe3",
"#ffffff"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# Two-tone purples
fCol <- colorRampPalette(c(
"#ce1366",
"#c3bfe1"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# Blues
fCol <- colorRampPalette(c(
"#303e9d",
"#9e5daa",
"#5960b0",
"#009ab4",
"#98baea",
"#cfe3fa",
"#ffffff"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# Reds
fCol <- colorRampPalette(c(
"#a2191a",
"#f0423e",
"#f77d81",
"#f89f9e",
"#fca037",
"#febf7e",
"#ffe2c2",
"#ffe0e0",
"#ffffff"))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# German Banknotes D-Mark hues @ (approx.) S:50, B:67
fCol <- colorRampPalette(c(
"#ab9e55",   #   50
"#8eab55",   #    5
"#55ab95",   #   20
"#5580ab",   #  100
"#ab55a7",   #   10
"#ab556b",   #  500
"#a35d48",   # 1000
"#a81e13"    #  200
))
n <- 32
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# GG-Plot spectrum hues (approx. equalized spectrum)
fCol <- colorRampPalette(c(
"#f9766d",   #   50
"#d99100",   #    5
"#a4a600",   #   20
"#39b700",   #  100
"#00b1f7",   #   10
"#9691ff",   #  500
"#e86bf4",   # 1000
"#ff62bd"   #  200
))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# modified from "tableau" spectrum hues (approx. equalized spectrum)
fCol <- colorRampPalette(c(
"#867ac4",   #
"#9976b9",   #
"#cd64b0",   #
"#d56283",   #
"#e66844",   #
"#ff981c",   #
"#fdb716",   #
"#fec657",   #
"#d5c63f",   #
"#a4bb49",   #
"#43ad68",   #
"#33abb4",   #
"#3f94be",
"#5565c6"
))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# petrol to dark
fCol <- colorRampPalette(c(
  "#39b6bc",   #
  "#a1d5d4",   #
  "#839493",   #
  "#283268"    #
))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# grey/egg/slate/dark
fCol <- colorRampPalette(c(
  "#eae7df",   #
  "#f5ef9f",   #
  "#6c7e8f",   #
  "#333f48"    #
))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# grey/dark/orange
fCol <- colorRampPalette(c(
  "#eef2f5",   #
  "#afbfcc",   #
  "#7891a7",   #
  "#4e4e4e",   #
  "#ff9674"    #
), bias = 0.7)
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# deepsea blue to bright yellow
fCol <- colorRampPalette(c(
"#4e196b",   #
"#492b7d",   #
"#413d86",   #
"#355a8f",   #
"#2a6f90",   #
"#1c848e",   #
"#10998a",   #
"#23b17a",   #
"#4fc663",   #
"#87d53c",   #
"#c6e000",   #
"#fbf54e"    #
))
n <- 20
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# bright rainbow - inspired by
# http://www.datapointed.net/2010/09/men-women-color-names/
fCol <- colorRampPalette(c(
  "#6627d9",   #
  "#8f26b3",   #
  "#bc09e6",   #
  "#e719b8",   #
  "#ee2c5b",   #
  "#fa0f3f",   #
  "#db1d1c",   #
  "#f84c15",   #
  "#ff7800",   #
  "#ffbd00",   #
  "#ffff00",   #
  "#e9f502",   #
  "#d4f502",   #
  "#b9f106",   #
  "#2af108",   #
  "#3ee192",   #
  "#29e3cb",   #
  "#00d0f2",   #
  "#3b8fe6",   #
  "#2936cc",   #
  "#2121a8",   #
  "#3b0bd9"    #
))
n <- 30
barplot(rep(1, n), col=fCol(n), axes=F, main="")

# ==== The Brewer scales ===============================
# see http://colorbrewer2.org/

# Sequential, single hue

fCol <- colorRampPalette(
c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#f7fcf5","#e5f5e0","#c7e9c0","#a1d99b","#74c476","#41ab5d","#238b45","#006d2c","#00441b"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525","#000000"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fcfbfd","#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



# sequential, multi-hue

fCol <- colorRampPalette(
c("#f7fcfd","#e0ecf4","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#f7fcf0","#e0f3db","#ccebc5","#a8ddb5","#7bccc4","#4eb3d3","#2b8cbe","#0868ac","#084081"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff7fb","#ece2f0","#d0d1e6","#a6bddb","#67a9cf","#3690c0","#02818a","#016c59","#014636"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#f7f4f9","#e7e1ef","#d4b9da","#c994c7","#df65b0","#e7298a","#ce1256","#980043","#67001f"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177","#49006a"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#ffffe5","#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d","#238443","#006837","#004529"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



# Diverging

fCol <- colorRampPalette(
c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#762a83","#9970ab","#c2a5cf","#e7d4e8","#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#b2182b","#d6604d","#f4a582","#fddbc7","#ffffff","#e0e0e0","#bababa","#878787","#4d4d4d"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")



fCol <- colorRampPalette(
c("#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd"))
n <- 10
barplot(rep(1, n), col=fCol(n), axes=F, main="")


# ToDo ... make various preview styles: histogram, scatterplot, stacked barchart illustrating sequential, divergent and categorical color choices.


# ====================================
# More information:
# http://en.wikipedia.org/wiki/Color_difference
# http://stackoverflow.com/questions/9018016/how-to-compare-two-colors
# http://www.scribblelive.com/blog/2013/03/25/seeing-color-through-infographics-and-data-visualizations/
# http://www.informationisbeautifulawards.com/
#

?convertColor
?hcl   # hue/chroma/luminance (equiperceptual) with demo...
?hsv
hsvCol <- matrix(seq(1,0,by= -0.02), ncol=1)
hsvCol <- cbind(hsvCol, 1)
hsvCol <- cbind(hsvCol, 1)
cols <- apply(hsvCol, 1, function(x) hsv(x[1], x[2], x[3]))
barplot(rep(1, length(cols)), col=cols)   # These colors have different brightness...
