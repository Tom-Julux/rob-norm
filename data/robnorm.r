# Reference RobNorm implementation in R
# Copied from: https://github.com/mwgrassgreen/RobNorm
# License:
#                   GNU LESSER GENERAL PUBLIC LICENSE
#                       Version 3, 29 June 2007
#
# Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
# Everyone is permitted to copy and distribute verbatim copies
# of this license document, but changing it is not allowed.
#
#
#  This version of the GNU Lesser General Public License incorporates
# the terms and conditions of version 3 of the GNU General Public
# License, supplemented by the additional permissions listed below.
#
#  0. Additional Definitions.
#
#  As used herein, "this License" refers to version 3 of the GNU Lesser
# General Public License, and the "GNU GPL" refers to version 3 of the GNU
# General Public License.
#
#  "The Library" refers to a covered work governed by this License,
# other than an Application or a Combined Work as defined below.
#
#  An "Application" is any work that makes use of an interface provided
# by the Library, but which is not otherwise based on the Library.
# Defining a subclass of a class defined by the Library is deemed a mode
# of using an interface provided by the Library.
#
#  A "Combined Work" is a work produced by combining or linking an
# Application with the Library.  The particular version of the Library
# with which the Combined Work was made is also called the "Linked
# Version".
#
#  The "Minimal Corresponding Source" for a Combined Work means the
# Corresponding Source for the Combined Work, excluding any source code
# for portions of the Combined Work that, considered in isolation, are
# based on the Application, and not on the Linked Version.
#
#  The "Corresponding Application Code" for a Combined Work means the
# object code and/or source code for the Application, including any data
# and utility programs needed for reproducing the Combined Work from the
# Application, but excluding the System Libraries of the Combined Work.
#
#  1. Exception to Section 3 of the GNU GPL.
#
#  You may convey a covered work under sections 3 and 4 of this License
# without being bound by section 3 of the GNU GPL.
#
#  2. Conveying Modified Versions.
#
#  If you modify a copy of the Library, and, in your modifications, a
# facility refers to a function or data to be supplied by an Application
# that uses the facility (other than as an argument passed when the
# facility is invoked), then you may convey a copy of the modified
# version:
#
#   a) under this License, provided that you make a good faith effort to
#   ensure that, in the event an Application does not supply the
#   function or data, the facility still operates, and performs
#   whatever part of its purpose remains meaningful, or
#
#   b) under the GNU GPL, with none of the additional permissions of
#   this License applicable to that copy.
#
#  3. Object Code Incorporating Material from Library Header Files.
#
#  The object code form of an Application may incorporate material from
# a header file that is part of the Library.  You may convey such object
# code under terms of your choice, provided that, if the incorporated
# material is not limited to numerical parameters, data structure
# layouts and accessors, or small macros, inline functions and templates
# (ten or fewer lines in length), you do both of the following:
#
#   a) Give prominent notice with each copy of the object code that the
#   Library is used in it and that the Library and its use are
#   covered by this License.
#
#   b) Accompany the object code with a copy of the GNU GPL and this license
#   document.
#
#  4. Combined Works.
#
#  You may convey a Combined Work under terms of your choice that,
# taken together, effectively do not restrict modification of the
# portions of the Library contained in the Combined Work and reverse
# engineering for debugging such modifications, if you also do each of
# the following:
#
#   a) Give prominent notice with each copy of the Combined Work that
#   the Library is used in it and that the Library and its use are
#   covered by this License.
#
#   b) Accompany the Combined Work with a copy of the GNU GPL and this license
#   document.
#
#   c) For a Combined Work that displays copyright notices during
#   execution, include the copyright notice for the Library among
#   these notices, as well as a reference directing the user to the
#   copies of the GNU GPL and this license document.
#
#   d) Do one of the following:
#
#       0) Convey the Minimal Corresponding Source under the terms of this
#       License, and the Corresponding Application Code in a form
#       suitable for, and under terms that permit, the user to
#       recombine or relink the Application with a modified version
#       of the Linked Version to produce a modified Combined Work, in the
#       manner specified by section 6 of the GNU GPL for conveying
#       Corresponding Source.
#
#       1) Use a suitable shared library mechanism for linking with the
#       Library.  A suitable mechanism is one that (a) uses at run time
#       a copy of the Library already present on the user's computer
#       system, and (b) will operate properly with a modified version
#       of the Library that is interface-compatible with the Linked
#       Version.
#
#   e) Provide Installation Information, but only if you would otherwise
#   be required to provide such information under section 6 of the
#   GNU GPL, and only to the extent that such information is
#   necessary to install and execute a modified version of the
#   Combined Work produced by recombining or relinking the
#   Application with a modified version of the Linked Version. (If
#   you use option 4d0, the Installation Information must accompany
#   the Minimal Corresponding Source and Corresponding Application
#   Code. If you use option 4d1, you must provide the Installation
#   Information in the manner specified by section 6 of the GNU GPL
#   for conveying Corresponding Source.)
#
#  5. Combined Libraries.
#
#  You may place library facilities that are a work based on the
# Library side by side in a single library together with other library
# facilities that are not Applications and are not covered by this
# License, and convey such a combined library under terms of your
# choice, if you do both of the following:
#
#   a) Accompany the combined library with a copy of the same work based
#   on the Library, uncombined with any other library facilities,
#   conveyed under the terms of this License.
#
#   b) Give prominent notice with the combined library that part of it
#   is a work based on the Library, and explaining where to find the
#   accompanying uncombined form of the same work.
#
#  6. Revised Versions of the GNU Lesser General Public License.
#
#  The Free Software Foundation may publish revised and/or new versions
#  of the GNU Lesser General Public License from time to time. Such new
#  versions will be similar in spirit to the present version, but may
#  differ in detail to address new problems or concerns.
#
#  Each version is given a distinguishing version number. If the
# Library as you received it specifies that a certain numbered version
# of the GNU Lesser General Public License "or any later version"
# applies to it, you have the option of following the terms and
# conditions either of that published version or of any later version
# published by the Free Software Foundation. If the Library as you
# received it does not specify a version number of the GNU Lesser
# General Public License, you may choose any version of the GNU Lesser
# General Public License ever published by the Free Software Foundation.
#
#  If the Library as you received it specifies that a proxy can decide
# whether future versions of the GNU Lesser General Public License shall
# apply, that proxy's public statement of acceptance of any version is
# permanent authorization for you to choose that version for the
# Library.

#' @title sim.dat.fn
#' @description To simulation an expression matrix
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param row.frac Outlier fraction in the rows of the expression matrix.
#' @param col.frac Outlier fraction in the columns of the expression matrix.
#' @param mu.up The up-shifted mean of the outliers.
#' @param mu.down The down-shifted mean of the outliers.
#' @param n Number of rows (genes).
#' @param m Number of colunms (samples).
#' @param nu.fix Logic value indicating underying nu=0 (default: TRUE).
#' @return Simulated data and its parameters.
#' @import invgamma
#' @export


sim.dat.fn = function(row.frac, col.frac, mu.up, mu.down, n, m, nu.fix=TRUE) {
	mu.00 = rnorm(n, 0, 1)
	var.00 = rinvgamma(n, 5, 2)
    X.0 = matrix(rnorm(n*m, outer(mu.00, rep(1, m)), sqrt(outer(var.00, rep(1, m)))), n, m) # the null matrix
    
    if (nu.fix) {
    	nu.00 = rep(0,m)
    } else {
       nu.00 = rnorm(m, 0, 1)
	   nu.00[1:round(m*0.2)] = rnorm(round(m * 0.2), 1, 1)	
    }
    B = outer(rep(1, n), nu.00) # the sample effect matrix

	S = matrix(0, n, m) 
	if (row.frac*col.frac > 0) {
	   	   bk.nm = round(n*row.frac * m*col.frac)
		   a = rbinom(1, bk.nm, 0.8)
		   S[ 1:round(n*row.frac), 1:(m*col.frac)] = sample(c(rep(mu.up, a), rep(0, bk.nm-a)), bk.nm) # the shifted mean of the signal mx
		   a = rbinom(1, bk.nm, 0.8)
		   S[ (n-round(n*row.frac)+1):n, (m-m*col.frac+1):m] = sample(c(rep(mu.down, a), rep(0, bk.nm-a)), bk.nm) # the signal mx of shifted mean 		  
	}
	
	X = X.0 + B + S
	S.ind = matrix(0, nrow(X), ncol(X))
	S.ind[S != 0] = 1
	
	rownames(X) = paste("prt", 1:nrow(X), sep=".")
	colnames(X) = paste("s", 1:ncol(X), sep=".")
    return(list(dat=X, mu.00=mu.00, var.00=var.00, nu.00=nu.00, sig.ind=S.ind))	
}

#' @title RobNorm
#' @description To robustly normalize expression data.
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param X.0 The expression matrix in log scale.
#' @param gamma.0 The density exponent parameter gamma, in practice, taking gamma.0 = 0.5 or 1.
#' @param tol The tolerance for interations (default: 10^(-4)).
#' @param step The step limit (default: 50).
#' @return Normalized expression data
#' @export


RobNorm = function(X.0, gamma.0 = 0.5, tol = 10^(-4), step = 200) {
    # to select the rows that are nonmissing in at least half of the samples
    id.norm = rownames(X.0)[rowSums(!is.na(X.0)) >= (ncol(X.0) / 2)]
    X.1 = X.0[id.norm, ]

    # to construct the standard samples from row medians
    x.stand = apply(X.1, 1, median)
    names(x.stand) = rownames(X.1)

    X = cbind(x.stand, X.1)
    I = nrow(X)
    J = ncol(X)

    # to initialize parameters
    Nu.0 = apply(X, 2, function(x) median(x - x.stand, na.rm = TRUE))
    Mu.0 = apply(X - outer(rep(1, I), Nu.0), 1, mean, na.rm = TRUE)
    sigma2.0 = apply((X - outer(rep(1, I), Nu.0) - outer(Mu.0, rep(1, J)))^2, 1, mean, na.rm = TRUE)

    mean.0.mx = Mu.0 %*% t(rep(1, J)) + rep(1, I) %*% t(Nu.0)
    var.0.mx = sigma2.0 %*% t(rep(1, J))
    den.0 = dnorm(X, mean.0.mx, sqrt(var.0.mx))
    dum = den.0^gamma.0
    dum[is.na(den.0)] = NA
    w.mx = dum / (rowSums(dum, na.rm = TRUE) %*% t(rep(1, J)))
    w.size = rowSums(dum, na.rm = TRUE)

    # to obtain the initial divergence
    divergence.0 = sum(w.size * ((gamma.0 + 1) * rowSums(w.mx * (X - Mu.0 %*% t(rep(1, J)) - rep(1, I) %*% t(Nu.0))^2, na.rm = TRUE) / sigma2.0 + log(2 * pi * sigma2.0 / (gamma.0 + 1))))

    # iterations
    para.diff.int = c()
    divergence.int = divergence.0
    int = 0
    flag = FALSE
    while (flag == FALSE) {
        # to update mu and sigma-squared
        Mu.new = rowSums((X - rep(1, I) %*% t(Nu.0)) * w.mx, na.rm = TRUE)
        sigma2.new = (rowSums((X - rep(1, I) %*% t(Nu.0) - Mu.new %*% t(rep(1, J)))^2 * w.mx, na.rm = TRUE)) * (1 + gamma.0)

        # to break the iterations if some sigma-squared is too small (locally trapped)
        if (min(sigma2.new, na.rm = TRUE) < 10^(-10)) {
            Mu.0 = NULL
            Nu.0 = NULL
            sigma2.0 = NULL
            divergence.int = NULL
            para.diff.int = NULL
            stop("the pre-chosen gamma is large, to choose a smaller nonnegative gamma")
        }

        # to udpate nu
        weight.nu = w.mx * ((w.size / (sigma2.new / (1 + gamma.0))) %*% t(rep(1, J)))
        Nu.new = colSums(weight.nu * (X - Mu.new %*% t(rep(1, J))), na.rm = TRUE) / colSums(weight.nu, na.rm = TRUE)
        Mu.new = Mu.new + Nu.new[1]
        Nu.new = Nu.new - Nu.new[1]

        para.diff.new = sum(abs(Mu.new - Mu.0), na.rm = TRUE) + sum(abs(Nu.new - Nu.0), na.rm = TRUE) + sum(abs(sigma2.new - sigma2.0), na.rm = TRUE)
        para.diff.int = c(para.diff.int, para.diff.new)

        # to update the divergence
        divergence.new = sum(w.size * ((gamma.0 + 1) * rowSums(w.mx * (X - Mu.new %*% t(rep(1, J)) - rep(1, I) %*% t(Nu.new))^2, na.rm = TRUE) / sigma2.new + log(2 * pi * sigma2.new / (gamma.0 + 1))), na.rm = TRUE)
        divergence.int = c(divergence.int, divergence.new)

        int = int + 1
        flag = (para.diff.new < tol) | (int >= step)

        if (flag == TRUE) break

        Mu.0 = Mu.new
        Nu.0 = Nu.new
        sigma2.0 = sigma2.new
        mean.0.mx = Mu.0 %*% t(rep(1, J)) + rep(1, I) %*% t(Nu.0)
        var.0.mx = sigma2.0 %*% t(rep(1, J))
        den.0 = dnorm(X, mean.0.mx, sqrt(var.0.mx))

        dum = den.0^gamma.0
        dum[is.na(den.0)] = NA
        w.mx = dum / (rowSums(dum, na.rm = TRUE) %*% t(rep(1, J)))
        w.size = rowSums(dum, na.rm = TRUE)
    }

    nu.rob.est = Nu.0[-1]
    X.0.rob = X.0 - outer(rep(1, nrow(X.0)), nu.rob.est)

    return(list(norm.data = X.0.rob, id.norm = id.norm, stand.s = x.stand, nu.est = nu.rob.est, mu.est = Mu.0, sigma2.est = sigma2.0, w.mx = w.mx, divergence = divergence.int, para.diff = para.diff.int))
}


library(invgamma)

sim.result = sim.dat.fn(row.frac = 0.2, col.frac = 0.2, mu.up = 3, mu.down = -3, n = 500, m = 20, nu.fix = TRUE)
X.0 = sim.result$dat

# write the input to a text file for comparison with the Python implementation
#write.table(X.0, file = "./data/simulated_measurements.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
start.time <- Sys.time()

norm.result = RobNorm(X.0, gamma.0 = 0.25)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
X.0.norm = norm.result$norm.data

# write the normalized output to a text file for comparison with the Python implementation
#write.table(X.0.norm, file = "./data/simulated_measurements_normalized.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
