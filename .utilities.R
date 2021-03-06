# .utilities.R
#
# Miscellaneous R code to suppport the project
#
# Version: 1.0
# Date:    2016 12
# Author:  Boris Steipe
#
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================

objectInfo <- function(x) {
    # Function to combine various information items about R objects
    #
    # Input: an R object
    # Value: none - prints information as side-effect

    cat("object contents:")
    print(x, digits = 22)  # print value at maximal precision

    cat("\nstructure of object:\n")
    str(x)

    if (! is.list(x)) { # Don't use cat() if x is a list. cat() can't handle lists.
        cat("\nmode:   ", mode(x), "\n")
        cat("typeof: ", typeof(x), "\n")
        cat("class:  ", class(x), "\n")
    }

    # if the object has attributes, print them too
    if (! is.null(attributes(x))) {
        cat("\nattributes:\n")
        attributes(x)
    }
    # Done
}

pBar <- function(i, l, nCh = 50) {
    # Draw a progress bar in the console
    # i: the current iteration
    # l: the total number of iterations
    # nCh: width of the progress bar
    ticks <- round(seq(1, l-1, length.out = nCh))
    if (i < l) {
        if (any(i == ticks)) {
            p <- which(i == ticks)
            p1 <- paste(rep("#", p), collapse = "")
            p2 <- paste(rep("-", nCh - p), collapse = "")
            cat(sprintf("\r|%s%s|", p1, p2))
            flush.console()
        }
    }
    else { # done
        cat("\n")
    }
}

checkFit <- function(iRow, fit) {
    t <- seq(0, 120, by = 5)
    y <- ygProfiles[iRow, ]
    plot(t, ygProfiles[iRow, ], col = "black", type = "b",
         xlab = "t (min.)", ylab = "expression log-ratio",
         main = sprintf("%d: %s (%s)",
                        iRow,
                        ygData$sysName[iRow],
                        ygData$stdName[iRow]))
    abline(h =  0, col = "#DDEEFF")
    abline(v = 60, col = "#DDEEFF")
    mtext(sprintf("Parameters: cor: %5.3f, %s",
                  cor(y, predict(fit)),
                  paste(names(coef(fit)),
                        sprintf("%5.3f", coef(fit))
                        , sep = ": ", collapse = ", ")),
          col = "#DD99CC", side = 1, line = 4)
    points(t, predict(fit), col = "#DD99CC", type = "l")
}



# [END]
