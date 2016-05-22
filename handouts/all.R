
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig" <- function(x, Y, cov.function = "stationary.cov", 
    lambda = NA, df = NA, GCV = FALSE, Z = NULL, cost = 1, knots = NA, 
    weights = NULL, m = 2, nstep.cv = 200, scale.type = "user", 
    x.center = rep(0, ncol(x)), x.scale = rep(1, ncol(x)), rho = NA, 
    sigma2 = NA, method = "REML", verbose = FALSE, mean.obj = NA, 
    sd.obj = NA, null.function = "Krig.null.function", wght.function = NULL, 
    offset = 0,  na.rm = TRUE, cov.args = NULL, 
    chol.args = NULL, null.args = NULL, wght.args = NULL, W = NULL, 
    give.warnings = TRUE, ...) # the verbose switch prints many intermediate steps as an aid in debugging.
#
{ 
    #
    # create output list
    out <- list()
    ###########################################################
    #  First series of steps simply store pieces of the passed
    #    information to output list (i.e. the Krig object)
    ##########################################################
        out$call <- match.call()
    #   turn off warning based on options
        if( options()$warn < 0 ){
        	give.warnings<- FALSE
        }
    #
    # save covariance function as its name
    #
    if( !is.character( cov.function)){
    out$cov.function.name <- as.character(substitute(cov.function))
    }
    else{ 
    	out$cov.function.name<-cov.function
    	} 
    #
    # save null space function as its name
    #
    out$null.function.name <- as.character(substitute(null.function))
    #
    # save weight  function as its name if it is not a NULL
    #
    if (is.null(wght.function)) {
        out$wght.function.name <- NULL
    }
    else {
        out$wght.function.name <- as.character(substitute(wght.function))
    }
    out$W <- W
    if (verbose) {
        print(out$cov.function.name)
        print(out$null.function.name)
        print(out$wght.function.name)
    }
    #
    # logical to indicate if the 'C' argument is present in cov.function
    #
    C.arg.missing <- all(names(formals(get(out$cov.function.name))) != 
        "C")
    if (C.arg.missing) 
        stop("Need to have C argument in covariance function\nsee Exp.cov.simple as an example")
    #
    # save parameters values possibly passed to the covariance function
    # also those added to call are assumed to be covariance arguments.
    if (!is.null(cov.args)) 
        out$args <- c(cov.args, list(...))
    else out$args <- list(...)
    #
    # default values for null space function
    out$null.args <- null.args
    #
    #       set degree of polynomial null space if this is default
    #       mkpoly is used so often is it helpful to include m argument
    #       by default in Krig call.
    if (out$null.function.name == "Krig.null.function") {
        out$null.args <- list(m = m)
        out$m <- m
    }
    #
    # default values for Cholesky decomposition, these are important
    # for sparse matrix decompositions used in Krig.engine.fixed.
    if (is.null(chol.args)) {
        out$chol.args <- list(pivot = FALSE)
    }
    else {
        out$chol.args <- chol.args
    }
    # additional arguments for weight matrix.
    out$wght.args <- wght.args
    #
    # the offset is the effective number of parameters used in the GCV
    # calculations -- unless this is part of an additive model this
    # is likely zero
    out$offset <- offset
    #
    # the cost is the multiplier applied to the GCV eff.df
    # lambda and df are two ways of parameterizing the smoothness
    # and are related by a monotonic function that unfortunately
    # depends on the locations of the data.
    # lambda can be used directly in the linear algebra, df
    # must be transformed to lambda numerically using the monotonic trransformation
    # sigma2 is the error variance and rho the multiplier for the covariance
    # method is how to determine lambda
    # the GCV logical forces the code to do the more elaborate decompositions
    # that faclitate estimating lambda -- even if a specific lambda value is
    # given.
    out$cost <- cost
    out$lambda <- lambda
    out$eff.df <- df
    out$sigma2 <- sigma2
    out$rho <- rho
    out$method <- method
    out$GCV <- GCV
    #
    # correlation model information
    #
    out$mean.obj <- mean.obj
    out$sd.obj <- sd.obj
    out$correlation.model <- !(is.na(mean.obj[1]) & is.na(sd.obj[1]))
    #
    # transformation info
    out$scale.type <- scale.type
    out$x.center <- x.center
    out$x.scale <- x.scale
    #
    # verbose block
    if (verbose) {
        cat("  Cov function arguments in call  ", fill = TRUE)
        print(out$args)
        cat(" covariance function used is : ", fill = TRUE)
        print(out$cov.function.name)
    }
    ###############################################################
    # Begin modifications and transformations of input information
    # note that many of these manipulations follow a strategy
    # of passing the Krig object (out) to a function and
    # then appending the information from this function to
    # the Krig object. In this way the Krig object  is built up
    # in steps and it is hoped easier to follow.
    ###############################################################
    # various checks on x and  Y including removal of NAs in Y
    # Here is an instance of adding to the Krig object
    # in this case some onerous bookkeeping making sure arguments are consistent
    out2 <- Krig.check.xY(x, Y, Z, weights, na.rm, verbose = verbose)
    out <- c(out, out2)
    # transform to correlation model (if appropriate)
    # find replicates and collapse to means and pool variances.
    # Transform unique x locations and knots.
    if (out$correlation.model) {
        out$y <- Krig.cor.Y(out, verbose = verbose)
    }
    out2 <- Krig.transform.xY(out, knots, verbose = verbose)
    out <- c(out, out2)
    # NOTE: knots have been transformed after this step
    #############################################################
    #  Figure out what to do
    #############################################################
    #
    # this functions works through the logic of
    # what has been supplied for lambda
    out2 <- Krig.which.lambda(out)
    out[names(out2)] <- out2  
    # Make weight matrix for observations
    #    ( this is proportional to the inverse square root of obs covariance)
    #     if a weight function or W has not been passed then this is
    #     diag( out$weightsM) for W
    #     The checks represent a limitation of this model to
    #     the  WBW type decoposition and no replicate observations.
    out$nondiag.W <- (!is.null(wght.function)) | (!is.null(W))
    # Do not continue if there there is a nondiagonal weight matrix
    # and replicate observations.
    if (out$nondiag.W) {
        if (out$knot.model | out$fixed.model) {
            stop("Non diagonal weight matrix for observations not supported\nwith knots or fixed lambda.")
        }
        if (!is.na(out$shat.pure.error)) {
            stop("Non diagonal weight matrix not implemented with replicate\nlocations")
        }
    }
    #  make weight matrix and its square root having passed checks
    out <- c(out, Krig.make.W(out, verbose = verbose))
    ########################################################
    #  You have reached the Engines!
    ########################################################
    #   Do the intensive linear algebra to find the solutions
    #   this is where all the heavy lifting happens.
    #
    #   Note that all the information is passed as a list
    #   including arguments to the cholesky decomposition
    #   used within Krig.engine.fixed
    #
    # The results are saved in the component matrices
    #
    # if method=='user' then just evaluate at single lambda
    #  fixed here means a fixed lambda
    #
    # For fixed lambda the decompositions with and without knots
    # are surprisingly similar and so are in one engine.
    ###########################################################
    if (out$fixed.model) {
        out$matrices <- Krig.engine.fixed(out, verbose = verbose)
        # can't find the trace of A matrix in fixed lambda case so set this to NA.
        out$eff.df <- NA
    }
    #
    # alternative are
    # matrix decompositions suitable for
    # evaluation at many lambdas to facilitate GCV/REML estimates  etc.
    #
    if (!out$fixed.model) {
        if (out$knot.model) {
            # the knot model engine
            out$matrices <- Krig.engine.knots(out, verbose = verbose)
            out$pure.ss <- out$matrices$pure.ss
        }
        else {
            # standard engine following the basic computations for thin plate splines
            out$matrices <- Krig.engine.default(out, verbose = verbose)
        }
    }
    #
    # store basic information about decompositions
    out$nt <- out$matrices$nt
    out$np <- out$matrices$np
    out$decomp <- out$matrices$decomp
    #
    # Now determine a logical vector indices for coefficients tied to  the
    # the 'spatial drift' i.e. the fixed part of the model
    # that is not due to the Z covariates.
    # NOTE that the spatial drift coefficients must be the first columns of the
    # M matrix
    if (is.null(out$Z)) {
        out$ind.drift <- rep(TRUE, out$nt)
    }
    else {
        
        mZ <- ncol(out$ZM)
        out$ind.drift <- c(rep(TRUE, out$nt - mZ), rep(FALSE, 
            mZ))
    }
    if (verbose) {
        cat("null df: ", out$nt, "drift df: ", sum(out$ind.drift), 
            fill = TRUE)
    }
    #########################
    # End of engine block
    #########################
    #################################################
    # Do GCV and REML search over lambda if not fixed or if GCV variable is TRUE
    #################################################
    if (!out$fixed.model | out$GCV) {
        if (verbose) {
            cat("call to gcv.Krig", fill = TRUE)
        }
        gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset, give.warnings=FALSE)
        out$gcv.grid <- gcv.out$gcv.grid
        #  a handy summary table of the search results
        out$lambda.est <- gcv.out$lambda.est
        out$warningTable<- gcv.out$warningTable
        if( verbose){
        	cat("summaries from grid search/optimization", fill=TRUE)
        	print(out$lambda.est)
        	print(out$warningTable)
        }
        if( give.warnings){
        	#NOTE: only print out grid search warning forthe method of interest.
        	printGCVWarnings( gcv.out$warningTable, method=method)
        }
          # assign the preferred lambda either from GCV/REML/MSE or the user value
        # NOTE: gcv/reml can be done but the estimate is
        # still evaluted at the passed user values of lambda (or df)
        # If df is passed need to calculate the implied lambda value
        if (out$method != "user") {
            out$lambda <- gcv.out$lambda.est[out$method, 1]
            out$eff.df <- out$lambda.est[out$method, 2]
        }
        else {
            if (!is.na(out$eff.df)) {
                out$lambda <- Krig.df.to.lambda(out$eff.df, out$matrices$D)
            }
            else {
                out$eff.df <- Krig.ftrace(out$lambda, out$matrices$D)
            }
        }
    }
    ##########################
    # end GCV/REML block
    ##########################
    #
    # Now we clean up what has happened and stuff 
    # information into output object.
    #
    ##########################################
    # find coefficients at prefered lambda
    # and evaluate the solution at observations
    ##########################################
    #   pass replicate group means -- no need to recalculate these.

    out2 <- Krig.coef(out, yM = out$yM, verbose = verbose)
    out <- c(out, out2)
    #######################################################################
    # fitted values and residuals and predicted values for full model and
    # also on the null space (fixed
    # effects). But be sure to do this at the nonmissing x's.
    ##################################################################
    out$fitted.values <- predict.Krig(out, x = out$x, Z = out$Z, 
        eval.correlation.model = FALSE)
    out$residuals <- out$y - out$fitted.values
    #
    # this is just M%*%d  note use of do.call using function name
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$x, Z = out$Z)))
    out$fitted.values.null <- as.matrix(Tmatrix) %*% out$d
    #
    # verbose block
    if (verbose) {
        cat("residuals", out$residuals, fill = TRUE)
    }
    #
    # find various estimates of sigma and rho
    out2 <- Krig.parameters(out)
    out <- c(out, out2)
    ################################################
    # assign the 'best' model as a default choice
    # either use the user supplied values or the results from
    # optimization
    ################################################
    passed.sigma2 <- (!is.na(out$sigma2))
    if (out$method == "user" & passed.sigma2) {
        out$best.model <- c(out$lambda, out$sigma2, out$rho)
    }
    else {
        # in this case lambda is from opt. or supplied by user
        out$best.model <- c(out$lambda, out$shat.MLE^2, out$rhohat)
    }
    # Note: values in best.model are used in subsquent functions as the choice
    # for these parameters!
    # set class
    class(out) <- c("Krig")
    return(out)
}


# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

Krig.check.xY <- function(x, Y, Z, weights, na.rm, 
    verbose = FALSE) {
    #
    # check for missing values in Y or X.
    #
    # save logical indicating where there are NA's
    # and check for NA's
    #
    ind <- is.na(Y)
    if (any(ind) & !na.rm) {
        stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }
    #
    # coerce x to be a matrix
    x <- as.matrix(x)
    #
    # coerce Y to be a vector
    #
    Y <- as.matrix(Y)
    if (ncol(Y) != 1) {
        stop("Krig can handle matrix Y data")
    }
    #
    #default weights ( reciprocal variance of errors).
    #
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    #
    # check that dimensions agree
    #
    if (length(Y) != nrow(x)) {
        stop(" length of y and number of rows of x differ")
    }
    if (length(Y) != length(weights)) {
        stop(" length of y and weights differ")
    }
    #  if Z is not NULL coerce to be  a matrix
    # and check  # of rows
    if (verbose) {
        print(Z)
    }
    if (!is.null(Z)) {
        if (!is.matrix(Z)) {
            Z <- matrix(c(Z), ncol = 1)
        }
        if (length(Y) != nrow(Z)) {
            stop(" length of y and number of rows of Z differ")
        }
    }
    # if NAs can be removed then remove them and warn the user
    if (na.rm) {
        ind <- is.na(Y)
        if(all(ind)){
        	stop("All Y values are missing!")
        }
        if (any(ind)) {
            Y <- Y[!ind]
            x <- as.matrix(x[!ind, ])
            if (!is.null(Z)) {
                Z <- Z[!ind, ]
            }
            weights <- weights[!ind]
        }
    }
    #
    # check for NA's in x matrix -- there should not be any !
    if (any(c(is.na(x)))) {
        stop(" NA's in x matrix")
    }
    #
    # check for NA's in Z matrix
    if (!is.null(Z)) {
        if (any(c(is.na(Z)))) {
            stop(" NA's in Z matrix")
        }
    }
    #
    # verbose block
    if (verbose) {
        cat("Y:", fill = TRUE)
        print(Y)
        cat("x:", fill = TRUE)
        print(x)
        cat("weights:", fill = TRUE)
        cat(weights, fill = TRUE)
    }
    #
    # save x, weights  and Y w/o NAs
    N <- length(Y)
    return(list(N = N, y = Y, x = x, weights = weights, Z = Z, 
        NA.ind = ind))
}

"Krig.coef" <- function(out, lambda = out$lambda, 
    y = NULL, yM = NULL, verbose = FALSE) {
    #
    # NOTE default value of lambda used from Krig object.
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
        #
        #old code   beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D))%d*%u)
        #
        ind <- (nt + 1):np
        D2 <- out$matrices$D[ind]
        #
        # note use of efficient diagonal multiply in next line
        temp2 <- (D2/(1 + lambda * D2)) %d*% u[ind, ]
        beta2 <- out$matrices$V %*% temp2
        temp.c <- rbind(matrix(0, nrow = nt, ncol = ndata), beta2)
        temp.c <- qr.qy(out$matrices$qr.T, temp.c)
        temp.c <- out$W2 %d*% temp.c
        temp <- temp.yM - do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots, C = temp.c)))
        temp <- out$W2 %d*% temp
        temp.d <- qr.coef(out$matrices$qr.T, temp)
    }
    #
    # case with knots
    # any lambda
    #
    if (out$decomp == "DR") {
        # X is the monster matrix ...  X = [ M | K]
        X <- cbind(do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM))), do.call(call.name, 
            c(out$args, list(x1 = out$xM, x2 = out$knots))))
        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM %d*% 
            temp.yM)
        beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) %d*% 
            u)
        temp.d <- beta[1:nt, ]
        temp.c <- beta[(nt + 1):np, ]
        temp <- X %*% out$matrices$G %*% u
        temp <- sum(out$weightsM * (temp.yM - temp)^2)
        #### ????
        out2$pure.ss <- temp + out2$pure.ss
    }
    #
    # fixed lambda knots == unique x's
    #
    if (out$decomp == "cholesky") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        temp.d <- qr.coef(out$matrices$qr.VT, forwardsolve(out$matrices$Mc, 
            transpose = TRUE, temp.yM, upper.tri = TRUE))
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp.yM - Tmatrix %*% temp.d, upper.tri = TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    #
    # fixed lambda with knots
    #
    if (out$decomp == "cholesky.knots") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        # form K matrix
        K <- do.call(call.name, c(out$args, list(x1 = out$xM, 
            x2 = out$knots)))
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM)))
        wY <- out$weightsM * temp.yM
        temp0 <- t(K) %*% (out$weightsM * Tmatrix)
        temp1 <- forwardsolve(out$matrices$Mc, temp0, transpose = TRUE, 
            upper.tri = TRUE)
        qr.Treg <- qr(t(Tmatrix) %*% (out$weightsM * Tmatrix) - 
            t(temp1) %*% temp1)
        temp0 <- t(K) %*% wY
        temp3 <- t(Tmatrix) %*% wY - t(temp1) %*% forwardsolve(out$matrices$Mc, 
            temp0, transpose = TRUE, upper.tri = TRUE)
        temp.d <- qr.coef(qr.Treg, temp3)
        temp1 <- t(K) %*% (wY - out$weightsM * (Tmatrix) %*% 
            temp.d)
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp1, upper.tri = TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    return(list(c = temp.c, d = temp.d, shat.rep = out2$shat.rep, 
        shat.pure.error = out2$shat.pure.error, pure.ss = out2$pure.ss))
}

Krig.cor.Y <- function(obj, verbose = FALSE) {
    # subtract mean
    if (!is.na(obj$mean.obj[1])) {
        Y <- obj$y - predict(obj$mean.obj, obj$x)
    }
    # divide by sd
    if (!is.na(obj$sd.obj[1])) {
        Y <- Y/predict(obj$sd.obj, obj$x)
    }
    Y
}

Krig.Amatrix <- function(object, x0 = object$x, lambda = NULL, 
    eval.correlation.model = FALSE, ...) {
    if (is.null(lambda)) {
        lambda <- object$lambda
    }
    M <- nrow(object$xM)
    N <- nrow(x0)
    # create output matrix
    out <- matrix(NA, N, M)
    #
    # loop through unique data locations predicting response
    # using unit vector
    # NOTE that the y vector has already been collapsed onto means.
    #
    for (k in 1:M) {
        ytemp <- rep(0, M)
        ytemp[k] <- 1
        out[, k] <- predict(object, x = x0, yM = ytemp, lambda = lambda, 
            eval.correlation.model = eval.correlation.model, 
            ...)
    }
    return(out)
}
"Krig.df.to.lambda" <- function(df, D, guess = 1, 
    tol = 1e-05) {
    if (is.list(D)) {
        D <- D$matrices$D
    }
    if (is.na(df)) 
        return(NA)
    if (df < sum(D == 0)) {
        warning("df too small to match with a lambda value")
        return(NA)
    }
    if (df > length(D)) {
        warning(" df too large to match a lambda value")
        return(NA)
    }
    l1 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l1 * D))
        if (tr <= df) 
            break
        l1 <- l1 * 4
    }
    l2 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l2 * D))
        if (tr >= df) 
            break
        l2 <- l2/4
    }
    info <- list(D = D, df = df, N = length(D))
    out <- bisection.search(log(l1), log(l2), Krig.fdf, tol = tol, 
        f.extra = info)$x
    +exp(out)
}

"Krig.engine.default" <- function(out, verbose = FALSE) {
    #
    # matrix decompositions for computing estimate
    #
    # Computational outline:( '.' is used for subscript)
    #
    # The form of the estimate is
    #    fhat(x) = sum phi.j(x) d.j  + sum psi.k(x) c.k
    #
    # the {phi.j} are the fixed part of the model usually low order polynomials
    # and is also referred to as spatial drift.
    #
    # the {psi.k} are the covariance functions evaluated at the unique observation
    # locations or 'knots'.  If xM.k is the kth unique location psi.k(x)= k(x, xM.k)
    # xM is also out$knots in the code below.
    #
    # the goal is find decompositions that facilitate rapid solution for
    # the vectors d and c. The eigen approach below was identified by
    # Wahba, Bates Wendelberger and is stable even for near colinear covariance
    # matrices.
    # This function does the main computations leading to the matrix decompositions.
    # With these decompositions the coefficients of the solution are found in
    # Krig.coef and the GCV and REML functions in Krig.gcv.
    #
    #  First is an outline calculations with equal weights
    #  T the fixed effects regression matrix  T.ij = phi.j(xM.i)
    #  K the covariance matrix for the unique locations
    # From the spline literature the solution solves the well known system
    # of two eqautions:
    #    -K( yM - Td - Kc) + lambda *Kc = 0
    #                 -T^t ( yM-Td -Kc) = 0
    #
    # Mulitple through by K inverse and substitute, these are equivalent to
    #
    #  -1-   -( yM- Td - Kc) + lambda c = 0
    #  -2-                        T^t c = 0
    #
    #
    #  A QR decomposition is done for   T= (Q.1,Q.2)R
    #   by definition  Q.2^T T =0
    #
    #  equation  -2- can be thought of as a constraint
    # with  c= Q.2 beta2
    # substitute in  -1-  and multiply through by Q.2^T
    #
    #      -Q.2^T yM  + Q.2^T K Q.2 beta2  + lambda beta2 = 0
    #
    #   Solving
    #   beta2 = {Q.2^T K Q.2 + lambda I )^ {-1} Q.2^T yM
    #
    # and so one sloves this linear system for beta2 and then uses
    #     c= Q.2 beta2
    #   to determine c.
    #
    #  eigenvalues and eigenvectors are found for M= Q.2^T K Q.2
    #     M = V diag(eta) V^T
    #  and these facilitate solving this system efficiently for
    #  many different values of lambda.
    #  create eigenvectors, D = (0, 1/eta)
    #  and G= ( 0,0) %*% diag(D)
    #         ( 0,V)
    # so that
    #
    #          beta2 = G%*% ( 1/( 1+ lambda D)) %*% u
    # with
    #
    #          u = (0, V Q.2^T W2 yM)
    #
    # Throughout keep in mind that M has smaller dimension than G due to
    # handling the null space.
    #
    # Now solve for d.
    #
    # From -1-  Td = yM - Kc - lambda c
    #      (Q.1^T) Td =  (Q.1^T) ( yM- Kc)
    #
    #   ( lambda c is zero by -2-)
    #
    #   so Rd = (Q.1^T) ( yM- Kc)
    # use qr functions to solve triangular system in R to find d.
    #
    #----------------------------------------------------------------------
    # What about errors with a general precision matrix, W?
    #
    # This is an important case because with replicated observations the
    # problem will simplify into a smoothing problem with the replicate group
    # means and unequal measurement error variances.
    #
    # the equations to solve are
    #     -KW( yM - Td - Kc) + lambda *Kc = 0
    #     -T^t W( yM-Td -Kc) =0
    #
    # Multiple through by K inverse and substitute, these are equivalent to
    #
    #  -1b-      -W( yM- Td - Kc) + lambda c = 0
    #  -2b-      (WT)^t c = 0
    #
    # Let W2 be the symmetric square root of W,  W= W2%*% W2
    # and W2.i be the inverse of W2.
    #
    #  -1c-      -(  W2 yM - W2 T d - (W2 K W2) W2.ic) + lambda W2.i c = 0
    #  -2c-      (W2T)^t  W2c = 0
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    if (verbose) {
        cat(" Model Matrix: spatial drift and Z", fill = TRUE)
        print(Tmatrix)
    }
    # Tmatrix premultiplied by sqrt of wieghts
    Tmatrix <- out$W2 %d*% Tmatrix
    qr.T <- qr(Tmatrix)
    if( qr.T$rank < ncol( Tmatrix)){
      stop("Regression matrix for fixed part of model is colinear")}
    #
    #verbose block
    if (verbose) {
        cat("first 5 rows of qr.T$qr", fill = TRUE)
        print(qr.T$qr[1:5, ])
    }
    #
    # find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
    #
    tempM <- t(out$W2 %d*% do.call(out$cov.function.name, c(out$args, 
        list(x1 = out$knots, x2 = out$knots))))
    tempM <- out$W2 %d*% tempM
    tempM <- qr.yq2(qr.T, tempM)
    tempM <- qr.q2ty(qr.T, tempM)
    np <- nrow(out$knots)
    nt <- (qr.T$rank)
    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)
    }
    #
    # Full set of decompositions for
    # estimator for nonzero lambda
    tempM <- eigen(tempM, symmetric = TRUE)
    D <- c(rep(0, nt), 1/tempM$values)
    #
    # verbose block
    if (verbose) {
        cat("eigen values:", fill = TRUE)
        print(D)
    }
    #
    # Find the transformed data vector used to
    # evaluate the solution, GCV, REML  at different lambdas
    #
    
    u <- c(rep(0, nt), t(tempM$vectors) %*% qr.q2ty(qr.T, c(out$W2 %d*% 
        out$yM)))
    if (verbose) {
        cat("u vector:", fill = TRUE)
        print(u)
    }
    #
    #
    return(list(D = D, qr.T = qr.T, decomp = "WBW", V = tempM$vectors, 
        u = u, nt = nt, np = np))
}

"Krig.engine.fixed" <- function(out, verbose = FALSE, 
    lambda = NA) {
    #
    # Model:
    #     Y_k=  f_k + e_k
    #  var( e_k) = sigma^2/W_k
    #
    #   f= Td + h
    #    T is often a low order polynomial
    #   E(h)=0    cov( h)= rho *K
    #
    # let M = (lambda W^{-1} + K)
    # the surface estimate  depends on coefficient vectors d and c
    #    The implementation in Krig/fields is that K are the
    #    cross covariances among the observation locations and the knot locations
    #    H is the covariance among the knot locations.
    #    Thus if knot locs == obs locs we have the obvious collapse to
    #    the simpler form for M above.
    #
    #   With M in hand ...
    #
    #   set
    #   d =  [(T)^t M^{-1} (T)]^{-1} (T)^t M^{-1} Y
    #  this is just the generalized LS estimate for d
    #
    #   lambda= sigma**2/rho
    #  the estimate for c is
    #   c=  M^{-1}(y - Td)
    #
    # This particular numerical strategy takes advantage of
    # fast Cholesky factorizations for positive definite matrices
    # and also provides a seamless framework for sparse matrix implementations
    #
    if (is.na(lambda)) 
        lambda <- out$lambda
    call.name <- out$cov.function.name
    if (!out$knot.model) {
        ####################################################
        # case of knot locs == obs locs  out$knots == out$xM
        ####################################################
        # create T matrix
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        if (verbose) {
            cat("Tmatrix:", fill = TRUE)
            print(Tmatrix)
        }
        np <- nrow(out$knots)
        nt <- ncol(Tmatrix)
        # form K
        tempM <- do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots)))
        # form M
        diag(tempM) <- (lambda/out$weightsM) + diag(tempM)
        #
        # find cholesky factor
        #  tempM = t(Mc)%*% Mc
        #  V=  Mc^{-T}
        # call cholesky but also add in the args supplied in Krig object.
        Mc <- do.call("chol", c(list(x = tempM), out$chol.args))
        VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE, 
            upper.tri = TRUE)
        qr.VT <- qr(VT)
        # find GLS covariance matrix of null space parameters.
        Rinv <- solve(qr.R(qr.VT))
        Omega <- Rinv %*% t(Rinv)
        #
        # now do generalized least squares for d
        # and then find c.
        d.coef <- qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
            out$yM, upper.tri = TRUE))
        if (verbose) {
            print(d.coef)
        }
        c.coef <- forwardsolve(Mc, transpose = TRUE, out$yM - 
            Tmatrix %*% d.coef, upper.tri = TRUE)
        c.coef <- backsolve(Mc, c.coef)
        # return all the goodies,  include lambda as a check because
        # results are meaningless for other values of lambda
        return(list(qr.VT = qr.VT, d = c(d.coef), c = c(c.coef), 
            Mc = Mc, decomp = "cholesky", nt = nt, np = np, lambda.fixed = lambda, 
            Omega = Omega))
    }
    else {
        ####################################################
        # case of knot locs != obs locs
        ####################################################
        # create weighted T matrix
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM)))
        nt <- ncol(Tmatrix)
        np <- nrow(out$knots) + nt
        # form H
        H <- do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots)))
        # form K matrix
        K <- do.call(call.name, c(out$args, list(x1 = out$xM, 
            x2 = out$knots)))
        #
        Mc <- do.call("chol", c(list(x = t(K) %*% (out$weightsM * 
            K) + lambda * H), out$chol.args))
        # weighted Y
        wY <- out$weightsM * out$yM
        temp0 <- t(K) %*% (out$weightsM * Tmatrix)
        temp1 <- forwardsolve(Mc, temp0, transpose = TRUE, upper.tri = TRUE)
        qr.Treg <- qr(t(Tmatrix) %*% (out$weightsM * Tmatrix) - 
            t(temp1) %*% temp1)
        temp0 <- t(K) %*% wY
        temp3 <- t(Tmatrix) %*% wY - t(temp1) %*% forwardsolve(Mc, 
            temp0, transpose = TRUE, upper.tri = TRUE)
        d.coef <- qr.coef(qr.Treg, temp3)
        temp1 <- t(K) %*% (wY - out$weightsM * (Tmatrix) %*% 
            d.coef)
        c.coef <- forwardsolve(Mc, transpose = TRUE, temp1, upper.tri = TRUE)
        c.coef <- backsolve(Mc, c.coef)
        list(qr.Treg = qr.Treg, d = c(d.coef), c = c(c.coef), 
            Mc = Mc, decomp = "cholesky.knots", nt = nt, np = np, 
            lambda.fixed = lambda, Omega = NA)
    }
    #
    # should not get here.
    #
}

"Krig.engine.knots" <- function(out, verbose = FALSE) {
    #
    # matrix decompostions for computing estimate when
    # knots are present
    # QR decomposition of null space regression matrix
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    qr.T <- qr(c(sqrt(out$weightsM)) * Tmatrix)
    nt <- ncol(Tmatrix)
    np <- nrow(out$knots) + nt
    if (verbose) {
        cat(nt, np, fill = TRUE)
    }
    # H is the penalty matrix in the ridge regression format
    # first part is zero because no penalty on part of estimator
    # spanned by T matrix
    H <- matrix(0, ncol = np, nrow = np)
    H[(nt + 1):np, (nt + 1):np] <- do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$knots, x2 = out$knots)))
    # X is the monster ...
    X <- cbind(do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM))), do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$xM, x2 = out$knots))))
    if (verbose) {
        cat("first lines of X", fill = TRUE)
        print(X[1:5, ])
    }
    #    sqrt(weightsM) * X
    XTwX <- t(X * out$weightsM) %*% X
    #
    # then  B= G(I-D)G^T
    # New version of diagonalize may be more stable
    out2 <- fields.diagonalize2((XTwX), H)
    D <- out2$D
    if (verbose) {
        cat("D;", fill = TRUE)
        cat(out2$D, fill = TRUE)
    }
    #
    #  G should satisfy:
    #     t(G) %*% XTwX %*%G = I and  t(G)%*%H%*%G = D
    #
    #     and
    #      solve( XtwX + lambda H) =  G%*%diag( 1/(1+ lambda*D))%*%t(G)
    #
    
    #  save XG to avoid an extra multiplication.
    XG <- X %*% out2$G
    
    u <- t(XG) %*% (out$weightsM * out$yM)
    #
    # adjust pure sum of squares to be that due to replicates
    # plus that due to fitting all the basis functions without
    # any smoothing. This will be the part of the RSS that does not
    # change as lambda is varied ( see e.g. gcv.Krig)
    #
    pure.ss <- sum(out$weightsM * (out$yM - XG %*% u)^2) + out$pure.ss
    if (verbose) {
        cat("total pure.ss from reps, reps + knots ", fill = TRUE)
        print(out$pure.ss)
        print(pure.ss)
    }
    
    #
    # in this form  the solution is (d,c)= G( I + lambda D)^-1 u
    # fitted.values = X ( d,c)
    #
    # output list
    # last D eigenvalues are zero due to null space of penalty
    # OLD code:    D[(np - nt + 1):np] <- 0
    # this should be enforced to machine precision from diagonalization.
    
    
    list(u = u, D = D, G = out2$G, qr.T = qr.T, decomp = "DR", 
        nt = nt, np = np, pure.ss = pure.ss)
}

"Krig.fdf" <- function(llam, info) {
    sum(1/(1 + exp(llam) * info$D)) - info$df
}

"Krig.fgcv" <- function(lam, obj) {
    #
    # GCV that is leave-one-group out
    #
    lD <- obj$matrices$D * lam
    RSS <- sum(((obj$matrices$u * lD)/(1 + lD))^2)
    MSE <- RSS/length(lD)
    if ((obj$N - length(lD)) > 0) {
        MSE <- MSE + obj$pure.ss/(obj$N - length(lD))
    }
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, MSE/den^2, 1e20)
}

"Krig.fgcv.model" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    MSE <- sum(((obj$matrices$u * lD)/(1 + lD))^2)/length(lD)
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    ifelse(den > 0, obj$shat.pure.error^2 + MSE/den^2, 1e20)
}

"Krig.fgcv.one" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    trA <- sum(1/(1 + lD))
    den <- 1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/obj$N
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, (RSS/obj$N)/den^2, 1e+20)
}

"Krig.flplike" <- function(lambda, obj) {
    #  - log profile likelihood for lambda
    # See section 3.4 from Nychka  Spatial Processes as Smoothers paper.
    # for equation and derivation
    D2 <- obj$matrices$D[obj$matrices$D > 0]
    u2 <- obj$matrices$u[obj$matrices$D > 0]
    lD <- D2 * lambda
    N2 <- length(D2)
    # MLE estimate of rho for fixed lambda
    rho.MLE <- (sum((D2 * (u2)^2)/(1 + lD)))/N2
    #
    # ln determinant of    K + lambda*WI
    lnDetCov <- -sum(log(D2/(1 + lD)))
    
    -1 * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(rho.MLE) - 
        (1/2) * lnDetCov)
      
    
}

"Krig.fs2hat" <- function(lam, obj) {
    lD  <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    den <- obj$N - (sum(1/(1 + lD)) + obj$offset)
    if (den < 0) {
        return(NA)
    }
    else {
        RSS/(den)
    }
}

"Krig.ftrace" <- function(lam, D) {
    sum(1/(1 + lam * D))
}

"Krig.make.W" <- function(out, verbose = FALSE) {
    if (verbose) {
        cat("W", fill = TRUE)
        print(out$W)
    }
    if (out$nondiag.W) {
        #
        # create W from scratch or grab it from passed object
        if (is.null(out$W)) {
            if (verbose) {
                print(out$wght.function.name)
            }
            W <- do.call(out$wght.function.name, c(list(x = out$xM), 
                out$wght.args))
            #       adjust W based on diagonal weight terms
            #
            W <- sqrt(out$weightsM) * t(sqrt(out$weightsM) * 
                W)
        }
        else {
            W <- out$W
        }
        #
        # symmetric square root
        temp <- eigen(W, symmetric = TRUE)
        W2 <- temp$vectors %*% diag(sqrt(temp$values)) %*% t(temp$vectors)
        return(list(W = W, W2 = W2))
    }
    else {
        #
        #  These are created only for use with default method to stay
        #   consistent with nondiagonal elements.
        if (out$fixed.model) {
            return(list(W = NULL, W2 = NULL))
        }
        else {
            return(list(W = out$weightsM, W2 = sqrt(out$weightsM)))
        }
    }
}

"Krig.make.Wi" <- function(out, verbose = FALSE) {
    #
    # If a weight matrix has been passed use it.
    #
    # Note that in either case the weight matrix assumes that
    # replicate observations have been collapses to the means.
    #
    if (out$nondiag.W) {
        temp <- eigen(out$W, symmetric = TRUE)
        Wi <- temp$vectors %*% diag(1/(temp$values)) %*% t(temp$vectors)
        W2i <- temp$vectors %*% diag(1/sqrt(temp$values)) %*% 
            t(temp$vectors)
        return(list(Wi = Wi, W2i = W2i))
    }
    else {
        #
        #  These are created only for use with default method to stay
        # consistent with nondiagonal elements.
        return(list(Wi = 1/out$weightsM, W2i = 1/sqrt(out$weightsM)))
    }
}

"Krig.make.u" <- function(out, y = NULL, yM = NULL, 
    verbose = FALSE) {
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
    }
    #
    # case with knots
    # any lambda
    #
    if (out$decomp == "DR") {
        # X is the monster matrix ...  X = [ M | K]
        X <- cbind(do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM))), do.call(call.name, 
            c(out$args, list(x1 = out$xM, x2 = out$knots))))
        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM %d*% 
            temp.yM)
    }
    return(list(u = u, shat.rep = out2$shat.rep, shat.pure.error = out2$shat.pure.error, 
        pure.ss = out2$pure.ss))
}

Krig.null.function <- function(x, Z = NULL, drop.Z = FALSE, 
    m) {
    # default function to create matrix for fixed part of model
    #  x, Z, and drop.Z are required
    #  Note that the degree of the polynomial is by convention (m-1)
    # returned matrix must have the columns from Z last!
    #
    if (is.null(Z) | drop.Z) {
        return(fields.mkpoly(x, m = m))
    }
    else {
        return(cbind(fields.mkpoly(x, m = m), Z))
    }
}

"Krig.parameters" <- function(obj, mle.calc = obj$mle.calc) {
    # if nondiag W is supplied then use it.
    # otherwise assume a diagonal set of weights.
    #
    # NOTE: calculation of  shat involves full set of obs
    # not those colllapsed to the mean.
    if (obj$nondiag.W) {
        shat.GCV <- sqrt(sum((obj$W2 %d*% obj$residuals)^2)/(length(obj$y) - 
            obj$eff.df))
    }
    else {
        shat.GCV <- sqrt(sum((obj$weights * obj$residuals^2)/(length(obj$y) - 
            obj$eff.df)))
    }
    if (mle.calc) {
        rho.MLE <- sum(c(obj$c) * c(obj$yM))/obj$N
        # set rho estimate to zero if negtive. Typically this
        # is an issue of machine precision and very small negative value.
        rho.MLE <- ifelse(rho.MLE < 0, 0, rho.MLE)
        
        #    commented out code for debugging ...
        #      if( rho.MLE< 0) {
        #        stop('problems computing rho.MLE')}
        # commented out is the REML estimate -- lose null space df because of
        # the restiction to orthogonal subspace of T.
        # rhohat<- rho.MLE <- sum(obj$c * obj$yM)/(obj$N - obj$nt)
        # .
        rhohat <- rho.MLE
        shat.MLE <- sqrt(rho.MLE * obj$lambda)
    }
    else {
        rhohat <- rho.MLE <- shat.MLE <- NA
    }
    list(shat.GCV = shat.GCV, rho.MLE = rho.MLE, shat.MLE = shat.MLE, 
        rhohat = rhohat)
}

"Krig.replicates" <- function(out=NULL, x,y, Z=NULL, weights=rep( 1, length(y)),
                               verbose = FALSE) {
    if( is.null(out)){
      out<- list( x=x, y=y, N= length(y), Z=Z, weights=weights)
    }
    rep.info <- cat.matrix(out$x)
    if (verbose) {
        cat("replication info", fill = TRUE)
        print(rep.info)
    }
    # If no replicates are found then reset output list to reflect this condition
    uniquerows <- !duplicated(rep.info)
    if (sum(uniquerows) == out$N) {
        shat.rep <- NA
        shat.pure.error <- NA
        pure.ss <- 0
        # coerce 'y' data vector as a single column matrix
        yM <- as.matrix(out$y)
        weightsM <- out$weights
        xM <- as.matrix(out$x[uniquerows, ])
        # coerce ZM to matrix
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z)
        }
        else {
            ZM <- NULL
        }
    }
    # collapse over spatial replicates
    else {
        rep.info.aov <- fast.1way(rep.info, out$y, out$weights)
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        # copy  replicate means as a single column matrix
        yM <- as.matrix(rep.info.aov$means)
        weightsM <- rep.info.aov$w.means
        xM <- as.matrix(out$x[uniquerows, ])
        # choose some Z's for replicate group means
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z[uniquerows, ])
        }
        else {
            ZM <- NULL
        }
        pure.ss <- rep.info.aov$SSE
        if (verbose) 
            print(rep.info.aov)
    }
    return(list(yM = yM, xM = xM, ZM = ZM, weightsM = weightsM, 
        uniquerows = uniquerows, shat.rep = shat.rep, shat.pure.error = shat.pure.error, 
        pure.ss = pure.ss, rep.info = rep.info))
}

Krig.transform.xY <- function(obj, knots, verbose = FALSE) {
    # find all replcates and  collapse to unique locations and mean response
    # and pooled variances and weights.
    out <- Krig.replicates(obj, verbose = verbose)
    if (verbose) {
        cat("yM from Krig.transform.xY", fill = TRUE)
        print(out$yM)
    }
    #
    # save information about knots.
    if (is.na(knots[1])) {
        out$knots <- out$xM
        out$mle.calc <- TRUE
        out$knot.model <- FALSE
    }
    else {
        out$mle.calc <- FALSE
        out$knot.model <- TRUE
        out$knots <- knots
    }
    #
    # scale x, knot locations and  save transformation info
    #
    out$xM <- transformx(out$xM, obj$scale.type, obj$x.center, 
        obj$x.scale)
    out$transform <- attributes(out$xM)
    out$knots <- scale(out$knots, center = out$transform$x.center, 
        scale = out$transform$x.scale)
    #
    #
    #verbose block
    #
    if (verbose) {
        cat("transform", fill = TRUE)
        print(out$transform)
    }
    if (verbose) {
        cat("knots in transformed scale", fill = TRUE)
        print(knots)
    }
    return(out)
}

"Krig.updateY" <- function(out, Y, verbose = FALSE, 
    yM = NA) {
    #given new Y values but keeping everything else the same finds the
    #new u vector and pure error SS associated with the Kriging estimate
    # the steps are
    # 1) standardize if neccesary
    # 2) find means, in the case of replicates
    # 3) based on the decomposition, multiply a weighted version of yM
    #    with a large matrix extracted from teh Krig object out.
    #
    # The out object may be large. This function is written so that out is # #not changed with the hope that it is not copied locally  in this  #function .
    # All of the output is accumulated in the list out2
    #STEP 1
    #
    # transform Y by mean and sd if needed
    #
    if (out$correlation.model) {
        Y <- (Y - predict(out$mean.obj, out$x))/predict(out$sd.obj, 
            out$x)
        if (verbose) 
            print(Y)
    }
    #
    #STEP 2
    if (is.na(yM[1])) {
        out2 <- Krig.ynew(out, Y)
    }
    else {
        out2 <- list(yM = yM, shat.rep = NA, shat.pure.error = NA, 
            pure.ss = NA)
    }
    if (verbose) {
        print(out2)
    }
    #
    #STEP3
    #
    # Note how matrices are grabbed from the Krig object
    #
    if (verbose) 
        cat("Type of decomposition", out$decomp, fill = TRUE)
    if (out$decomp == "DR") {
        #
        #
        u <- t(out$matrices$G) %*% t(out$matrices$X) %*% (out$weightsM * 
            out2$yM)
        #
        # find the pure error sums of sqaures.
        #
        temp <- out$matrices$X %*% out$matrices$G %*% u
        temp <- sum((out$W2 %d*% (out2$yM - temp))^2)
        out2$pure.ss <- temp + out2$pure.ss
        if (verbose) {
            cat("pure.ss", fill = TRUE)
            print(temp)
            print(out2$pure.ss)
        }
    }
    #####
    ##### end DR decomposition block
    #####
    ####
    #### begin WBW decomposition block
    ####
    if (out$decomp == "WBW") {
        #### decomposition of Q2TKQ2
        u <- c(rep(0, out$nt), t(out$matrices$V) %*% qr.q2ty(out$matrices$qr.T, 
            out$W2 %d*% out2$yM))
        if (verbose) 
            cat("u", u, fill = TRUE)
        #
        # pure error in this case from 1way ANOVA
        #
        if (verbose) {
            cat("pure.ss", fill = TRUE)
            print(out2$pure.ss)
        }
    }
    #####
    ##### end WBW block
    #####
    out2$u <- u
    out2
}
Krig.which.lambda <- function(out) {
    #
    # determine the method for finding lambda
    #  Note order
    # default is to do 'gcv/REML'
    out2 <- list()
    # copy all all parameters to out2 just to make this
    # easier to read.
    out2$method <- out$method
    out2$lambda.est <- NA
    out2$lambda <- out$lambda
    out2$eff.df <- out$eff.df
    out2$rho <- out$rho
    out2$sigma2 <- out$sigma2
    if (!is.na(out2$lambda) | !is.na(out2$eff.df)) {
        #
        # this indicates lambda has been supplied and leads to
        # the cholesky type computational approaches
        #        -- but only if GCV is FALSE
        #
        out2$method <- "user"
    }
    out2$GCV <- out$GCV
    if (!is.na(out2$eff.df)) {
        #
        # this indicates df has been supplied and needs
        # GCV to be true to compute the lambda
        # that matches the df
        #
        out2$GCV <- TRUE
    }
    if (!is.na(out2$rho) & !is.na(out2$sigma2)) {
        out2$method <- "user"
        out2$lambda <- out2$sigma2/out2$rho
    }
    #
    # NOTE: method='user' means that a value of lambda has been supplied
    #        and so GCV etc to determine lambda is not needed.
    #  gcv TRUE means that the decompositions will be done to
    #    evaluate the estimate at arbitrary lambda (and also be
    #    able to compute the effective degrees of freedom).
    #
    #    The fixed lambda calculations are very efficient but
    #    do not make it feasible for GCV/REML  or effective degrees of
    #    freedom calculations.
    #
    out2$fixed.model <- (out2$method == "user") & (!out2$GCV)
    #
    return(out2)
}

"Krig.ynew" <- function(out, y = NULL, yM = NULL) {
    #
    # calculates the collapsed y (weighted) mean vector based on the
    # X matrix and weights from the out object.
    # or just passes through the collapsed mean data if passed.
    #
    #
    # If there are no replicated obs. then return the full vector
    # pure error ss is zero
    #
    shat.rep <- NA
    shat.pure.error <- NA
    pure.ss <- 0
    # if no y's are given then it is assumed that one should use the
    # yM from the original data used to create the Krig object
    if (is.null(yM) & is.null(y)) {
        yM <- out$yM
    }
    #
    # case when yM is passed no calculations are needed
    #
    if (!is.null(yM)) {
        return(list(yM = as.matrix(yM), shat.rep = NA, shat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    # no reps case
    #
    if (length(unique(out$rep.info)) == out$N) {
        return(list(yM = as.matrix(y), shat.rep = NA, shat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    #  check that y is the right length
    #
    if (length(y) != out$N) {
        stop(" the new y vector is the wrong length!")
    }
    #
    # case when full y data is passed and replicate means need to be found
    #
    if (length(unique(out$rep.info)) < out$N) {
        #
        # calculate means by pooling Replicated obseravations but use the
        # the right weighting.
        #
        rep.info.aov <- fast.1way(out$rep.info, y, out$weights)[c("means", 
            "MSE", "SSE")]
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        return(list(yM = rep.info.aov$means, shat.rep = shat.rep, 
            shat.pure.error = shat.pure.error, pure.ss = rep.info.aov$SSE))
    }
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"as.image" <- function(Z, ind = NULL, grid = NULL, 
    x = NULL,  weights = rep(1, length(Z)), na.rm = FALSE, 
    nx = 64, ny = 64, boundary.grid = FALSE, nrow = NULL, ncol = NULL,
    FUN=NULL) {
    # NOTE that throughout ind is a two column integer matrix of
    # discretized locations in the image matrix.
    # Thanks to J. Rougier for fixing bugs in this function.
    # set some default values for arguments
    #
    # coerce Z to a vector
    Z <- c(Z)
    if( !is.null(ind)){
      x<- ind
     }
    if( !is.null(nrow)&!is.null(ncol)){
      nx<- nrow
      ny<- ncol
    }
    #
    # check for x or weights having missing values
    # we do not like these ...
    if (any(is.na(weights)) | any(is.na(c(x)))) {
        stop("missing values in weights or x")
    }
    # discretize locations to grid boxes
    # this function will also create a default grid based on range of
    # locations if grid is NULL
    #
    temp <- discretize.image(x, m = nx, n = ny, grid = grid, 
        boundary.grid = boundary.grid)
    grid <- temp$grid
    # index is a two component list that indexes  the x and y grid points.
    # points outside of grid are assigned as NA
    #
    # empty image matrices to hold weights and  weighted means
     w<- z <- matrix( NA, nrow=temp$m, ncol=temp$n)
     # find stats
     tempw<- tapply( weights, temp$index, sum, na.rm=FALSE)
     if( is.null(FUN)){
# usual weighted means case:     
     tempz<- tapply( Z*weights, temp$index,sum, na.rm=FALSE )
     tempz<- tempz/ tempw
     }
     else{
# just apply FUN to values in the grid box -- no weighting!     	
     	tempz<- tapply( Z, temp$index,FUN, na.rm=FALSE )
     	}
     # these are the indices that are represented by the locations
     # they may not include the entire set ( 1:nx and 1:ny)
     # so define what they do have.
  
     # insert the tabled values into the right rows and columns.
      z[ temp$ix, temp$iy] <- tempz
      w[ temp$ix, temp$iy] <- tempw
     # save call
     # xd created because it is a  pain to do otherwise and handy to have
    call <- match.call()
    list(x = grid$x, y = grid$y, z = z, call = call, ind = cbind(temp$index[[1]], temp$index[[2]]) , 
        weights = w, xd = cbind(grid$x[temp$index[[1]]], grid$y[temp$index[[2]]] ), 
        call = match.call(), FUN = FUN )
}
mKrig.MLE <- function(x, y, weights = rep(1, nrow(x)), 
    Z = NULL, ..., par.grid = NULL, lambda = NULL, lambda.profile = TRUE, 
    verbose = FALSE, relative.tolerance = 1e-04) {
    # these are all the arguments needed to call mKrig except lambda and those in par.grid.
    mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z), 
        list(...))
    lnProfileLike.max <- -1e+20
    # find NG --  number of parameters to try
    par.grid <- data.frame(par.grid)
    if (nrow(par.grid) == 0) {
        if (is.null(lambda)) {
            NG <- 1
        }
        else {
            NG <- length(lambda)
        }
    }
    else {
        NG <- nrow(par.grid)
    }
    # output matrix to summarize results
    summary <- matrix(NA, nrow = NG, ncol = 8)
    dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", 
        "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE", "counts eval", 
        "counts grad"))
    lambda.best <- NA
    # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
    # this is controlled by NAs for lambda starting values.
    if (is.null(lambda)) {
        lambda <- rep(NA, NG)
    }
    # default starting value for lambda is 1 or log lambda is 0
    llambda.opt <- 0
    optim.counts <- c(NA, NA)
    lnLike.eval <- list()
    # Define the objective function as a tricksy call to mKrig
    # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
    temp.fn <- function(x) {
        # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
        # assign to hold only a few components returned by mKrig
        hold <- do.call("mKrig", c(mKrig.args, list(find.trA = FALSE), 
            list(lambda = exp(x)), cov.args.temp))[c("lambda.fixed", 
            "rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
        # add this evalution to an  object (i.e. here a matrix) in the calling frame
        temp.eval <- get("capture.evaluations")
        assign("capture.evaluations", rbind(temp.eval, unlist(hold)), 
            envir = capture.env)
        return(hold$lnProfileLike.FULL)
    }
    #
    # begin loop over covariance arguments
    for (k in 1:NG) {
        llambda.start <- ifelse(is.na(lambda[k]), llambda.opt, 
            log(lambda[k]))
        # list of covariance arguments from par.grid with right names (some R arcania!)
        # note that this only works because temp.fn will search in this frame for this object
        # par.grid has been coerced to a data frame so one has a concept of a row subscript.
        cov.args.temp <- as.list(par.grid[k, ])
        names(cov.args.temp) <- names(par.grid)
        if (lambda.profile) {
            # set up matrix to store evaluations from within optim
            capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                dimnames = list(NULL, c("lambda", "rho.MLE", 
                  "sigma.MLE", "lnProfileLike.FULL")))
            capture.env <- environment()
            # call to optim
            look <- optim(llambda.start, temp.fn, method = "BFGS", 
                control = list(fnscale = -1, parscale = 0.1, 
                  ndeps = 0.05, reltol = relative.tolerance))
            llambda.opt <- look$par
            optim.counts <- look$counts
            # call to 1-d search
            #            opt.summary     <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
            #            llambda.opt <- opt.summary$maximum
            #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
            # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
            lnLike.eval <- c(lnLike.eval, list(capture.evaluations[-1, 
                ]))
        }
        else {
            # no refinement for lambda so just save the the 'start' value as final one.
            llambda.opt <- llambda.start
        }
        # final fit at optimal value (or starting value if not refinement/maximization for lambda)
        obj <- do.call("mKrig", c(mKrig.args, cov.args.temp, 
            list(lambda = exp(llambda.opt))))
        if (obj$lnProfileLike.FULL > lnProfileLike.max) {
            lnProfileLike.max <- obj$lnProfileLike.FULL
            cov.args.MLE <- cov.args.temp
            lambda.best <- exp(llambda.opt)
        }
        # save results of the kth covariance model evaluation
        summary[k, 1:8] <- c(obj$eff.df, obj$lnProfileLike.FULL, 
            obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, llambda.opt, 
            optim.counts)
        if (verbose) {
            cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
        }
    }
    return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
        mKrig.args = list(...), lambda.best = lambda.best, lambda.MLE = lambda.best, 
        call = match.call(), lnLike.eval = lnLike.eval))
}

# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

mKrig <- function(x, y, weights = rep(1, nrow(x)), Z = NULL, lambda = 0, cov.function = "stationary.cov", 
	m = 2, chol.args = NULL, cov.args = NULL, find.trA = TRUE, NtrA = 20, iseed = 123, llambda = NULL, 
	...) {
	# grab other arguments for covariance function
	cov.args <- c(cov.args, list(...))
	#
	if (!is.null(llambda)) {
		lambda <- exp(llambda)
	}
	# see comments in Krig.engine.fixed for algorithmic commentary
	#
# check for duplicate x's.
# stop if there are any
if (any(duplicated(cat.matrix(x)))) 
		stop("locations are not unique see help(mKrig) ")
	if (any(is.na(y))) 
		stop("Missing values in y should be removed")
	if (!is.null(Z)) {
		Z <- as.matrix(Z)
	}

	# create fixed part of model as m-1 order polynomial
	Tmatrix <- cbind(fields.mkpoly(x, m), Z)
	# set some dimensions
	np <- nrow(x)
	nt <- ncol(Tmatrix)
	nZ <- ifelse(is.null(Z), 0, ncol(Z))
	ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ))
	# as a place holder for reduced rank Kriging, distinguish between
	# observations locations and  the locations to evaluate covariance.
# (this is will also allow predict.mKrig to handle a Krig object)
knots <- x
	# covariance matrix at observation locations
	# NOTE: if cov.function is a sparse constuct then tempM will be sparse.
# see e.g. wendland.cov
tempM <- do.call(cov.function, c(cov.args, list(x1 = x, x2 = x)))
	#
	# decide how to handle the pivoting.
# one wants to do pivoting if the matrix is sparse.
# if tempM is not a matrix assume that it is in sparse format.
#
  sparse.flag <- !is.matrix(tempM)
# quantify sparsity of tempM for the mKrig object
  nzero <- ifelse(sparse.flag, length(tempM@entries), np^2)
#

# add diagonal matrix that is the observation error Variance
# NOTE: diag must be a overloaded function to handle  sparse format.
  if (lambda != 0) {
		diag(tempM) <- (lambda/weights) + diag(tempM)
	}
# At this point tempM is proportional to the covariance matrix of the
# observation vector, y.
#
# cholesky decoposition of tempM
# do.call used to supply other arguments to the function
# especially for sparse applications.
# If chol.args is NULL then this is the same as
#              Mc<-chol(tempM)
# set arguments that are passed to cholesky
  if (is.null(chol.args)) {
		chol.args <- list(pivot = sparse.flag)
	} 
  Mc <- do.call("chol", c(list(x = tempM), chol.args))
# Efficent way to multply inverse of Mc times the Tmatrix
	VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE,
	      upper.tri = TRUE)
	qr.VT <- qr(VT)
# start linear algebra to find solution
# Note that all these expressions make sense if y is a matrix
# of several data sets and one is solving for the coefficients
# of all of these at once. In this case d.coef and c.coef are matrices
#
# now do generalized least squares for d
d.coef <- as.matrix(qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, y, upper.tri = TRUE)))
	# and then find c.
	# find the coefficents for the spatial part.
c.coef <- as.matrix(forwardsolve(Mc, transpose = TRUE, y - Tmatrix %*% d.coef, upper.tri = TRUE))
	# save intermediate result this is   t(y- T d.coef)( M^{-1}) ( y- T d.coef)
	quad.form <- c(colSums(as.matrix(c.coef^2)))
	# find c coefficients
	c.coef <- as.matrix(backsolve(Mc, c.coef))
	# GLS covariance matrix for fixed part.
	Rinv <- solve(qr.R(qr.VT))
	Omega <- Rinv %*% t(Rinv)
	# MLE estimate of rho and sigma
	#    rhohat <- c(colSums(as.matrix(c.coef * y)))/(np - nt)
# NOTE if y is a matrix then each of these are vectors of parameters.
    rho.MLE <- quad.form/np
	rhohat <- c(colSums(as.matrix(c.coef * y)))/np
	shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
# the  log profile likehood with  rhohat  and  dhat substituted
# leaving a profile for just lambda.
# NOTE if y is a matrix then each of this is a vector of log profile
# likelihood values.
    lnDetCov <- 2 * sum(log(diag(Mc)))
    lnProfileLike <- (-np/2 - log(2 * pi) * (np/2) - (np/2) *
                         log(rho.MLE) - (1/2) * lnDetCov)
	rho.MLE.FULL <- mean(rho.MLE)
	sigma.MLE.FULL <- sqrt(lambda * rho.MLE.FULL)
	# if y is a matrix then compute the combined likelihood
	# under the assumption that the columns of y are replicated
# fields
    lnProfileLike.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - (np/2) *
                          log(rho.MLE.FULL) - (1/2) * lnDetCov))
#
# return coefficients and   include lambda as a check because
# results are meaningless for other values of lambda
# returned list is an 'object' of class mKrig (micro Krig)
# also save the matrix decompositions so coefficients can be
# recalculated for new y values.
out <- list(
                  d = (d.coef),
                  c = (c.coef),
                 nt = nt,
                 np = np,
       lambda.fixed = lambda,
                  x = x, 
		          y = y, 
	        weights = weights,
        	  knots = knots,
  cov.function.name = cov.function, args = cov.args, 
	              m = m,
	      chol.args = chol.args,
	           call = match.call(),
    nonzero.entries = nzero,
           shat.MLE = sigma.MLE, 
		  sigma.MLE = sigma.MLE,
	  	    rho.MLE = rho.MLE,
	      	 rhohat = rho.MLE,
      lnProfileLike = lnProfileLike, 
	   rho.MLE.FULL = rho.MLE.FULL,
     sigma.MLE.FULL = sigma.MLE.FULL,
 lnProfileLike.FULL = lnProfileLike.FULL, 
		   lnDetCov = lnDetCov,
		  quad.form = quad.form,
		      Omega = Omega,
		      qr.VT = qr.VT,
		         Mc = Mc, 
		    Tmatrix = Tmatrix,
		  ind.drift = ind.drift,
		         nZ = nZ)
	#
	# find the residuals directly from solution
# to avoid a call to predict
        out$residuals <- lambda * c.coef/weights
	out$fitted.values <- y - out$residuals
	# estimate effective degrees of freedom using Monte Carlo trace method.
	if (find.trA) {
		out2<-  mKrig.trace(out, iseed, NtrA)
		out<- c( out, out2)
	} else {
		out$eff.df <- NA
		out$trA.info <- NA
		out$GCV <- NA
	}
	class(out) <- "mKrig"
	return(out)
}# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

mKrig.trace <- function(object, iseed, NtrA) {
    set.seed(iseed)
    # if more tests that number of data points just
    # find A exactly by np predicts.
    np<- object$np
    if (NtrA >= object$np) {
        Ey <- diag(1, np)
        NtrA <- np
        hold <- diag(predict.mKrig(object, ynew = Ey))
        trA.info<- NA
        trA.est <- sum(hold)
    }
    else {
        # if fewer tests then use random trace method
        # find fitted.values  for iid N(0,1) 'data' to calculate the
        # the Monte Carlo estimate of tr A(lambda)
        # basically repeat the steps above but take some
        # short cuts because we only need fitted.values
        # create random normal 'data'
        Ey <- matrix(rnorm(np * NtrA), nrow = np, 
            ncol = NtrA)
        trA.info <- colSums(Ey * (predict.mKrig(object, ynew = Ey)))
        trA.est <- mean(trA.info)
    }
    if (NtrA < np) {
     MSE<-(sum(object$residuals^2)/np) 
     GCV <-       MSE/(1 - trA.est /np)^2
     GCV.info <- MSE/( 1 - trA.info/np)^2
    }
    else{
    	GCV<- NA
    	GCV.info <- NA
    }	
    return(
    list(trA.info = trA.info, eff.df = trA.est,
             GCV= GCV, GCV.info=GCV.info)
    )
}

mKrig.coef <- function(object, y) {
    # given new data y and the matrix decompositions in the
    # mKrig object find coefficients d and c.
    # d are the coefficients for the fixed part
    # in this case hard coded for a low order polynomial
    # c are coefficients for the basis functions derived from the
    # covariance function.
    #
    # see mKrig itself for more comments on the linear algebra
    #
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # generalized least squares for d
    if( any(is.na(y))){
    	stop("mKrig can not omit missing values in observation vecotor")
    }
    d.coef <- as.matrix(qr.coef(object$qr.VT, forwardsolve(object$Mc, 
        transpose = TRUE, y, upper.tri = TRUE)))
    #  residuals from subtracting off fixed part
    #  of model as m-1 order polynomial
    resid <- y - object$Tmatrix %*% d.coef
    # and now find c.
    c.coef <- forwardsolve(object$Mc, transpose = TRUE, resid, 
        upper.tri = TRUE)
    c.coef <- as.matrix(backsolve(object$Mc, c.coef))
    out <- list(d = (d.coef), c = (c.coef))
    return(out)
}
print.mKrig <- function(x, digits = 4, ...) {
    
    if (is.matrix(x$residuals)) {
        n <- nrow(x$residuals)
        NData <- ncol(x$residuals)
    }
    else {
        n <- length(x$residuals)
        NData <- 1
    }
    
    c1 <- "Number of Observations:"
    c2 <- n
    
    if (NData > 1) {
        c1 <- c(c1, "Number of data sets fit:")
        c2 <- c(c2, NData)
    }
    
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    c2 <- c(c2, x$m - 1)
    c1 <- c(c1, "Total number of parameters in base model")
    c2 <- c(c2, x$nt)
    if (x$nZ > 0) {
        c1 <- c(c1, "Number of additional covariates (Z)")
        c2 <- c(c2, x$nZ)
    }
    if (!is.na(x$eff.df)) {
        c1 <- c(c1, " Eff. degrees of freedom")
        c2 <- c(c2, signif(x$eff.df, digits))
        if (length(x$trA.info) < x$np) {
            c1 <- c(c1, "   Standard Error of estimate: ")
            c2 <- c(c2, signif(sd(x$trA.info)/sqrt(length(x$trA.info)), 
                digits))
        }
    }
    c1 <- c(c1, "Smoothing parameter")
    c2 <- c(c2, signif(x$lambda.fixed, digits))
    
    if (NData == 1) {
        c1 <- c(c1, "MLE sigma ")
        c2 <- c(c2, signif(x$shat.MLE, digits))
        c1 <- c(c1, "MLE rho")
        c2 <- c(c2, signif(x$rho.MLE, digits))
    }
    
    c1 <- c(c1, "Nonzero entries in covariance")
    c2 <- c(c2, x$nonzero.entries)
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat(" ", fill = TRUE)
    cat(" Covariance Model:", x$cov.function, fill = TRUE)
    if (x$cov.function == "stationary.cov") {
        cat("   Covariance function:  ", ifelse(is.null(x$args$Covariance), 
            "Exponential", x$args$Covariance), fill = TRUE)
    }
    if (!is.null(x$args)) {
        cat("   Non-default covariance arguments and their values ", 
            fill = TRUE)
        nlist <- as.character(names(x$args))
        NL <- length(nlist)
        for (k in 1:NL) {
            cat("   Argument:", nlist[k], " ")
            if (object.size(x$args[[k]]) <= 1024) {
                cat("has the value(s): ", fill = TRUE)
                print(x$args[[k]])
            }
            else {
                cat("too large to print value, size > 1K ...", 
                  fill = TRUE)
            }
        }
    }
    invisible(x)
}

summary.mKrig <- function(object, ...) {
    print.mKrig(object, ...)
}

predict.mKrig <- function(object, xnew = NULL, ynew = NULL, grid.list=NULL,
    derivative = 0, Z = NULL, drop.Z = FALSE, just.fixed = FALSE, 
    ...) {
    # the main reason to pass new args to the covariance is to increase
    # the temp space size for sparse multiplications
    # other optional arguments from mKrig are passed along in the
    # list object$args
    cov.args <- list(...)
    # predict at observation locations by default
    if( !is.null(grid.list)){
        xnew<- make.surface.grid(grid.list)
      }
    if (is.null(xnew)) {
        xnew <- object$x
    }
    if (is.null(Z)) {
        Z <- object$Tmatrix[, !object$ind.drift]
    }
    if (!is.null(ynew)) {
        coef.hold <- mKrig.coef(object, ynew)
        c.coef <- coef.hold$c
        d.coef <- coef.hold$d
    }
    else {
        c.coef <- object$c
        d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    #
    if (derivative == 0) {
        if (drop.Z | object$nZ == 0) {
            # just evaluate polynomial and not the Z covariate
            temp1 <- fields.mkpoly(xnew, m = object$m) %*% d.coef[object$ind.drift, 
                ]
        }
        else {
            temp1 <- cbind(fields.mkpoly(xnew, m = object$m), 
                Z) %*% d.coef
        }
    }
    else {
        if (!drop.Z & object$nZ > 0) {
            stop("derivative not supported with Z covariate included")
        }
        temp1 <- fields.derivative.poly(xnew, m = object$m, d.coef[object$ind.drift, 
            ])
    }
    if (just.fixed) {
        return(temp1)
    }
    # add nonparametric part. Covariance basis functions
    # times coefficients.
    # syntax is the name of the function and then a list with
    # all the arguments. This allows for different covariance functions
    # that have been passed as their name.
    if (derivative == 0) {
        # argument list are the parameters and other options from mKrig
        #  locations and coefficients,
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef), cov.args))
    }
    else {
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
            cov.args))
    }
    # add two parts together and coerce to vector
    return((temp1 + temp2))
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

# wrapper for Tps object 
"predict.Tps"<- function(object, ...){
  UseMethod("Krig")
 }

"predict.Krig" <- function(object, x = NULL, Z = NULL, 
    drop.Z = FALSE, just.fixed = FALSE, lambda = NA, df = NA, 
    model = NA, eval.correlation.model = TRUE, y = NULL, yM = NULL, 
    verbose = FALSE, ...) {
    #NOTE: most of this function is figuring out what to do!
    #
    # check that derivative is not called
    if (!is.null(list(...)$derivative)) {
        stop("For derivatives use predictDerivative")
    }
    # y is full data yM are the data collapsed to replicate means
    # if new data is not passed then copy from the object
    if (is.null(y) & is.null(yM)) {
        temp.c <- object$c
        temp.d <- object$d
    }
    # check for passed x but no Z -- this is an error
    # if there  are Z covariates in the model and drop.Z is FALSE
    ZinModel<- !is.null(object$Z)
    newX<- !is.null(x)
    missingZ<- is.null(Z)
    if( ZinModel&newX){
    if( missingZ  & !drop.Z) {
        stop("Need to specify drop.Z as TRUE or pass Z values")
    }
    }
    # default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    else {
        x <- as.matrix(x)
    }
    # default is to predict at data Z's
    if (is.null(Z)) {
        Z <- object$Z
    }
    else {
        Z <- as.matrix(Z)
    }
    if (verbose) {
        print(x)
        print(Z)
    }
    # transformations of x values used in Krig
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    x <- scale(x, xc, xs)
    # NOTE knots are already scaled in Krig object and are used
    # in transformed scale.
    #  i.e.   knots <- scale( object$knots, xc, xs)
    #
    # figure out if the coefficients for the surface needto be recomputed.
    find.coef <- (!is.null(y) | !is.null(yM) | !is.na(lambda) | 
        !is.na(df) | !is.na(model[1]))
    if (verbose) {
        cat("find.coef", find.coef, fill = TRUE)
    }
    #   convert effective degrees of freedom to equivalent lambda
    if (!is.na(df)) {
        lambda <- Krig.df.to.lambda(df, object$matrices$D)
    }
    if (!is.na(model)) {
        lambda <- model[1]
    }
    if (is.na(lambda)) 
        lambda <- object$lambda
    #
    # if the coefficients need to be recomputed  do it.
    if (find.coef) {
        if (verbose) {
            cat("new coefs found", fill = TRUE)
        }
        object3 <- Krig.coef(object, lambda = lambda, y = y, 
            yM = yM)
        temp.d <- object3$d
        temp.c <- object3$c
    }
    if (verbose) {
        cat(" d coefs", fill = TRUE)
        print(temp.d)
        cat("c coefs", fill = TRUE)
        print(temp.c)
    }
    
    # this is the fixed part of predictor
    #
    Tmatrix <- do.call(object$null.function.name, c(object$null.args, 
        list(x = x, Z = Z, drop.Z = drop.Z)))
    if (drop.Z) {
        temp <- Tmatrix %*% temp.d[object$ind.drift]
    }
    else {
        temp <- Tmatrix %*% temp.d
    }
    # add in spatial piece
    if (!just.fixed) {
        #
        # Now find sum of covariance functions times coefficients
        # Note that the multiplication of the cross covariance matrix
        # by the coefficients is done implicitly in the covariance function
        #
        # The covariance function is
        # evaluated by using its name, the do.call function, and any
        # additional arguments.
        #
        temp <- temp + do.call(object$cov.function.name, c(object$args, 
            list(x1 = x, x2 = object$knots, C = temp.c)))
    }
    #
    # transform back into raw scale if this is a correlation model.
    # if y's are in the scale of correlations
    # if so scale by sd and add back in mean
    correlation.model <- (object$correlation.model & eval.correlation.model)
    if (correlation.model) {
        if (!is.na(object$sd.obj[1])) {
            temp <- temp * predict(object$sd.obj, x)
        }
        if (!is.na(object$mean.obj[1])) {
            temp <- temp + predict(object$mean.obj, x)
        }
    }
    return(temp)
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


"predictSurface.Krig" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE,
       ZGrid=NULL, drop.Z= FALSE, just.fixed=FALSE,  ...) {
  
      if( is.null(ZGrid) & !drop.Z & (!is.null(object$Z)) ) {
      stop("Need to specify covariate (Z) values or set drop.Z==TRUE")
    }
# create a default grid if it is not passed    
    if (is.null(grid.list)) {
    # NOTE: 
    # without grid.list
    # default is 80X80 grid on first two variables
    # rest are set to median value of the x's
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    }
# do some checks on Zgrid and also reshape as a matrix
# rows index grid locations and columns  are the covariates
# (as Z in predict).
# if ZGrid is NULL just returns that back 
    Z<- unrollZGrid( grid.list, ZGrid) 
# here is the heavy lifting
    xg <- make.surface.grid(grid.list)
# NOTE: the predict function called will need to do some internal  the checks
# whether the evaluation of a large number of grid points (xg)  makes sense.
if( verbose){
print( dim( xg))
print( drop.Z)
print( dim( Z))
}
    out<-  predict(object, x=xg, Z=Z, drop.Z= drop.Z,   
                     just.fixed=just.fixed, ...)
# reshape as list with x, y and z components    
    out <-  as.surface( xg, out )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predict.surface" <- function(object, ...) {
    UseMethod("predict.surface")
}

predict.surface.default<- function(object,...){
   cat("predict.surface is now the function predictSurface")
 }

"predictSurface"<- function( object,...){
  UseMethod("predictSurface")
}

"predictSurface.default" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
    # NOTE: 
    # without grid.list
    # default is 80X80 grid on first two variables
    # rest are set to median value of x.
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
    # here is the heavy lifting
    xg <- make.surface.grid(grid.list)
# NOTE: the specific predict function called will need to do the checks
# whether the evaluation of a large number of grid points makes sense. 
    out <-  as.surface( xg, predict(object, xg,...) )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}

"predictSurface.mKrig" <- function( object, ...){
	NextMethod("predictSurface.Krig")
}

"predictSurface.fastTps" <- function(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...) {
# NOTE:  See  predictSurface.default for comments
    if (is.null(grid.list)) {
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    } 
# in the case of fastTps pass the grid list instead of the locations of grid points
#  (see xg in predictSurface.default)
    out <-  predict(object, grid.list=grid.list, xy=xy, ...)
    out <-  as.surface(grid.list, out )
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        xg<- make.surface.grid( grid.list)
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"quilt.plot" <- function(x, y, z, nx = 64, ny = 64, 
     grid = NULL, add.legend = TRUE, add = FALSE, nlevel=64, 
    col = tim.colors(nlevel), nrow = NULL, ncol = NULL, FUN=NULL,
    plot=TRUE, ...) {
    #
    # note that nrow and ncol refer to the resulting 'image format' for plotting.
    # here the x values are the rows and the y values are the columns
    # FUN = NULL means the weighted means are found for each grid cell
    if( !is.null(nrow)|!is.null(nrow)){
      nx<- nrow
      ny<- ncol
     }
    x <- as.matrix(x)
    if (ncol(x) == 2) {
        z <- y
    }
    if (ncol(x) == 1) {
        x <- cbind(x, y)
    }
    if (ncol(x) == 3) {
        z <- x[, 3]
        x <- x[, 1:2]
    }
    # at this point x should be a 2 column matrix of x-y locations
    #  z is a vector or one column matrix of the z values.
    #discretize data
    out.p <- as.image(z, x = x, nx = nx, ny = ny, na.rm = TRUE, 
        grid = grid, FUN=FUN)
    # besides the image information this list has the indices that 
    # map each z value to a grid box
    #    
    # plot it
    if( plot){
    if (add.legend) {
        image.plot(out.p, nlevel = nlevel, col = col, add = add, ...)
    }
    else {
        image(out.p, col = col, add = add, ...)
    }
    }
    invisible(out.p)
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"rdist" <- function(x1, x2=NULL ) {
    if (!is.matrix(x1)){
        x1 <- as.matrix(x1)
        }

    if (is.null(x2)){
    	 storage.mode(x1) <- "double"
 #     	 Rdist1C is  setup to only fill the upper triangle and diagonal 
 #       But the contents of the lower triangle is uncertain. 
 #        return( .Call("Rdist1C", x1, PACKAGE="fields") )
 # 
          return( .Call("RdistC", x1,x1, PACKAGE="fields") ) 
    }
    else{ 
    	if (!is.matrix(x2)){
        x2 <- as.matrix(x2)
        }
    	storage.mode(x1) <- "double"
    	storage.mode(x2) <- "double"
        return( .Call("RdistC", x1,x2, PACKAGE="fields") )         
    }                 
    
}
spatialProcess <- function(x, y, cov.function = "stationary.cov", 
	cov.args = list(Covariance = "Matern", smoothness = 1),
	 ngrid=10, theta.grid = NULL, ...) {
	MLEfit <- MLESpatialProcess(x, y, cov.function = cov.function, 
		cov.args = cov.args, ngrid=ngrid,theta.grid=theta.grid,
		 ...)
# now fit spatial model with MLE for theta (range parameter)
# reestimate the other parameters for simplicity to get the complete Krig
# object.		
	obj <- Krig(x, y, cov.function = cov.function, cov.args = cov.args, 
		theta = MLEfit$pars[1],  
		method = "REML", give.warnings=TRUE, 
		...)
	obj <- c(obj, MLEfit)
	obj$theta.MLE<- MLEfit$pars[1]
# replace call with this top level one
    obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","Krig")
 
	return(obj)
}
unrollZGrid<- function( grid.list, ZGrid){
  if( is.null(ZGrid)){
    return(ZGrid)
  }
  if( is.list( ZGrid) ){
     if( any(grid.list[[1]] != ZGrid[[1]]) |any(grid.list[[2]] != ZGrid[[2]]) ){
         stop("grid list does not match grid for covariates")
       }  
# wipe out the x and y components of ZGrid because grid.list will be used
  ZGrid<- ZGrid$z
  }
# check dimensions
    Zdim<- dim( ZGrid)
      nx<- length( grid.list[[1]])
      ny<- length( grid.list[[2]])
      if( (Zdim[1] != nx) | (Zdim[2] != ny) ){
         stop( "Dimension of ZGrid does not match dimensions of location grid list.")
      }
# reshape as a matrix where rows index locations.
# Note that this works whether Zdim[3] exists or not! 
      return( matrix( c(ZGrid),  nrow= Zdim[1]*Zdim[2] ))
 }