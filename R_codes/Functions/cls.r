#####################################################################
#Getting Started
#Step 1: Download R software from cran.r-project.org.
#Step 2: Load "lars" package from cran.r-project.org.
#Step 3: Paste the following two functions to R program.
# "lars" package is developed by Brad Efron and Trevor Hastie. The first three 
#functions are modified from the functions in "lars" in order to obtain positive estimators
#and use cross validation method to select the tuning parameter by
#randomly portioning the sampling locations. We acknowledge Brad Efron and
#Trevor Hastie for contributing this useful R package.
#####################################################################
library(lars)

"larspositive"<-
function (x, y, BPHI, is.CWLS=FALSE, is.CGLS=FALSE, type = c("lasso", "lar", "forward.stagewise",
    "stepwise"), trace = FALSE, normalize = TRUE, intercept = FALSE,
    Gram, eps = .Machine$double.eps, max.steps, use.Gram = TRUE)
{
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise",
        stepwise = "Forward Stepwise")
    if (trace)
        cat(paste(TYPE, "sequence\n"))
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]
    
    ###BPHI is related to the covariance matrix var(vec(S)).
    ###Default is positive Lasso approach and BPHI is not needed.
    ###BPHI=diag(var(vec(S))) and is.CWLS=TRUE for CWLS approach
    ###BPHI=var(vec(S)) and is.CGLS=TRUE for CGLS approach
    if(is.CWLS){
        R.m12.CWLS<-diag(diag(BPHI)^(-0.5))
        x<-R.m12.CWLS%*%x
        y<-R.m12.CWLS%*%y
    }
    if(is.CGLS){
        BPHI.eigen<-eigen(BPHI.fin)
        B.eigen<-BPHI.eigen$values
        R.m12.CGLS<-(BPHI.eigen$vectors)%*%diag(c(B.eigen^(-0.5)))%*%t(BPHI.eigen$vectors)
        x<-R.m12.CGLS%*%x
        y<-R.m12.CGLS%*%y
    }
    #############################################################

    if (intercept) {
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0, m)
        mu <- 0
        y <- drop(y)
    }
    if (normalize) {
        normx <- sqrt(drop(one %*% (x^2)))
        nosignal <- normx/sqrt(n) < eps
        if (any(nosignal)) {
            ignores <- im[nosignal]
            inactive <- im[-ignores]
            normx[nosignal] <- eps * sqrt(n)
            if (trace)
                cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
        }
        else ignores <- NULL
        names(normx) <- NULL
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, m)
        ignores <- NULL
    }
    if (use.Gram & missing(Gram)) {
        if (m > 500 && n < m)
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
        if (trace)
            cat("Computing X'X .....\n")
        Gram <- t(x) %*% x
    }
    Cvec <- drop(t(y) %*% x)
    ssy <- sum(y^2)
    residuals <- y
    if (missing(max.steps))
        max.steps <- 8 * min(m, n - intercept)
    beta <- matrix(0, max.steps + 1, m)
    lambda = double(max.steps)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while ((k < max.steps) & (length(active) < min(m - length(ignores),
        n - intercept))) {
        action <- NULL
        C <- Cvec[inactive]
        Cmax <- max(max(C),eps) ####### modification for positive estimators
        if (Cmax < eps * 100) {
            if (trace)
                cat("Max |corr| = 0; exiting...\n")
            break
        }
        k <- k + 1
        lambda[k] = Cmax
        if (!any(drops)) {
            new <- (C) >= Cmax - eps ####### modification for positive estimators
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                if (use.Gram) {
                  R <- updateR(Gram[inew, inew], R, drop(Gram[inew,
                    active]), Gram = TRUE, eps = eps)
                }
                else {
                  R <- updateR(x[, inew], R, x[, active], Gram = FALSE,
                    eps = eps)
                }
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace)
                    cat("LARS Step", k, ":\t Variable", inew,
                      "\tcollinear; dropped for good\n")
                }
                else {
                  if (first.in[inew] == 0)
                    first.in[inew] <- k
                  active <- c(active, inew)
                  Sign <- c(Sign, 1)  ####### modification for positive estimators
                  action <- c(action, inew)
                  if (trace)
                    cat("LARS Step", k, ":\t Variable", inew,
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, backsolvet(R, Sign))
        dropouts <- NULL
        if (type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if (!all(directions > 0)) {
                if (use.Gram) {
                  nnls.object <- nnls.lars(active, Sign, R, directions,
                    Gram[active, active], trace = trace, use.Gram = TRUE,
                    eps = eps)
                }
                else {
                  nnls.object <- nnls.lars(active, Sign, R, directions,
                    x[, active], trace = trace, use.Gram = FALSE,
                    eps = eps)
                }
                positive <- nnls.object$positive
                dropouts <- active[-positive]
                action <- c(action, -dropouts)
                active <- nnls.object$active
                Sign <- Sign[positive]
                Gi1 <- nnls.object$beta[positive] * Sign
                R <- nnls.object$R
                C <- Cvec[-c(active, ignores)]
            }
        }
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        if (!use.Gram)
            u <- drop(x[, active, drop = FALSE] %*% w)
        if ((length(active) >= min(n - intercept, m - length(ignores))) |
            type == "stepwise") {
            gamhat <- Cmax/A
        }
        else {
            if (use.Gram) {
                a <- drop(w %*% Gram[active, -c(active, ignores),
                  drop = FALSE])
            }
            else {
                a <- drop(u %*% x[, -c(active, ignores), drop = FALSE])
            }
            gam <- c((Cmax - C)/(A - a))####### modification for positive estimators
            gamhat <- min(gam[gam > eps & is.na(gam)!=1], Cmax/A)####### modification for positive estimated values
        }
        if (type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if (zmin < gamhat) {
                gamhat <- zmin
                drops <- z1 == zmin
            }
            else drops <- FALSE
        }
        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, active] <- beta[k + 1, active] + gamhat *
            w
        if (use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*%
                w
        }
        else {
            residuals <- residuals - gamhat * u
            Cvec <- drop(t(residuals) %*% x)
        }
        Gamrat <- c(Gamrat, gamhat/(Cmax/A))
        arc.length <- c(arc.length, gamhat)
        if (type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace)
                  cat("Lasso Step", k + 1, ":\t Variable", active[id],
                    "\tdropped\n")
                R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k + 1, dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn))
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]
        if (type == "stepwise")
            Sign = Sign * 0
    }
    beta <- beta[seq(k + 1), , drop = FALSE]
    lambda = lambda[seq(k)]
    dimnames(beta) <- list(paste(0:k), vn)
    if (trace)
        cat("Computing residuals, RSS etc .....\n")
    residuals <- y - x %*% t(beta)
    beta <- scale(beta, FALSE, normx)
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]
    actions = actions[seq(k)]
    netdf = sapply(actions, function(x) sum(sign(x)))
    df = cumsum(netdf)
    if (intercept)
        df = c(Intercept = 1, df + 1)
    else df = c(Null = 0, df)
    rss.big = rev(RSS)[1]
    df.big = n - rev(df)[1]
    if (rss.big < eps | df.big < eps)
        sigma2 = NaN
    else sigma2 = rss.big/df.big
    Cp <- RSS/sigma2 - n + 2 * df
    attr(Cp, "sigma2") = sigma2
    attr(Cp, "n") = n
    object <- list(call = call, type = TYPE, df = df, lambda = lambda,
        R2 = R2, RSS = RSS, Cp = Cp, actions = actions[seq(k)],
        entry = first.in, Gamrat = Gamrat, arc.length = arc.length,
        Gram = if (use.Gram) Gram else NULL, beta = beta, mu = mu,
        normx = normx, meanx = meanx)
    class(object) <- "lars"
    object
}


cv.lars.positive <-
function(x, y, BPHI,is.CWLS=FALSE,is.CGLS=FALSE,K = 10, fraction = seq(from = 0, to = 1, length = 100), trace = FALSE, plot.it = TRUE, se = TRUE)
{
  n.temp<-trunc(sqrt(2*length(y[,1])))
  all.folds <- cv.folds(n.temp, K)
  residmat <- matrix(0, length(fraction), K)

  for(i in seq(K)) {
  omit.t <- all.folds[[i]]
  omit<-c(1:length(y[,1]))
  loc.temp<-y[,1:2]
        for(j in 1:length(omit.t)){
	omit<-omit[loc.temp[,1]!=omit.t[j]&loc.temp[,2]!=omit.t[j]]
	loc.temp<-loc.temp[loc.temp[,1]!=omit.t[j]&loc.temp[,2]!=omit.t[j],]}
 	x.omit<-x[omit,]
 	y.omit<-y[omit,3]

		if(is.CWLS){
		BPHIO<-BPHI[omit,omit]
		BPHIO.m12<-diag(diag(BPHIO)^(-0.5))
		x.omit<-BPHIO.m12%*%x[omit,]
		y.omit<-BPHIO.m12%*%y[omit,3]
		}

		if(is.CGLS){
		BPHIO<-BPHI[omit,omit]
		BPHIO.eigen<-eigen(BPHIO)
		BPHIO.m12<-(BPHIO.eigen$vectors)%*%diag(c(BPHIO.eigen$values)^(-0.5))%*%t(BPHIO.eigen$vectors)
		x.omit<-BPHIO.m12%*%x[omit,]
		y.omit<-BPHIO.m12%*%y[omit,3]
		}
    fit <- larspositive(x.omit, y.omit, type="lasso", trace = trace)
    fit <- predict.larspositive(fit, x[-omit,  ,drop=FALSE], mode = "fraction", s = fraction
                   )$fit
    if(length(omit)==1)fit<-matrix(fit,nrow=1)
    residmat[, i] <- apply((y[-omit,3] - fit)^2, 2, mean)
    if(trace)
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(fraction = fraction, cv = cv, cv.error = cv.error)
}


"predict.larspositive"<-
function(object, newx, s, type = c("fit", "coefficients"), mode = c("step",
                                                               "fraction", "norm"), ...)
{
  mode <- match.arg(mode)
  type <- match.arg(type)
  if(missing(newx) & type == "fit") {
    warning("Type=fit with no newx argument; type switched to coefficients"
            )
    type <- "coefficients"
  }
  betas <- object$beta
  sbetas <- scale(betas, FALSE, 1/object$normx)#scaled for unit norm x
  kp <- dim(betas)
  k <- kp[1]
  p <- kp[2]
  steps <- seq(k)
  if(missing(s)) {
    s <- steps
    mode <- "step"
  }
  sbeta <- switch(mode,
                  step = {
                    if(any(s < 0) | any(s > k))
                      stop("Argument s out of range")
                    steps
                  }
                  ,
                  fraction = {
                    if(any(s > 1) | any(s < 0))
                      stop("Argument s out of range")
                    nbeta <- drop(abs(sbetas) %*% rep(1, p))
                    nbeta/nbeta[k]
                  }
                  ,
                  norm = {
                    nbeta <- drop(abs(sbetas) %*% rep(1, p))
                    if(any(s > nbeta[k]) | any(s < 0))
                      stop("Argument s out of range")
                    nbeta
                  }
                  )
  sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
  sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
  usbeta<-unique(sbeta)
  useq<-match(usbeta,sbeta)
  sbeta<-sbeta[useq]
  betas<-betas[useq,]
  coord <- approx(sbeta, seq(sbeta), sfrac)$y
  left <- floor(coord)
  right <- ceiling(coord)
  newbetas <- ((sbeta[right] - sfrac) * betas[left,  , drop = FALSE] + (sfrac -
                                                         sbeta[left]) * betas[right,  , drop = FALSE])/(sbeta[right] - sbeta[
                                                                                          left])
  newbetas[left == right,  ] <- betas[left[left == right],  ]
  robject <- switch(type,
                    coefficients = list(s = s, fraction = sfrac, mode = mode,
                      coefficients = drop(newbetas)),
                    fit = list(s = s, fraction = sfrac, mode = mode, fit = drop(
                                                                       scale(newx, object$meanx, FALSE) %*% t(newbetas)) + object$
                      mu))
  robject
}



 tensor<-function(data,data1){
	n.prod<-dim(data)[1]
	m.prod<-dim(data)[2]
	value.ten<-NULL
	for(i in 1:n.prod){
		value.i<-NULL
		for(j in 1:m.prod){
			ten.temp<-data[i,j]*data1
			value.i<-cbind(value.i,ten.temp)
		}
		value.ten<-rbind(value.ten,value.i)
	}
	return(value.ten)
}


Knn.func<-function(n){
	Knn<-matrix(0,n^2,n^2)
	for(i in 1:(n^2)){
		j.temp<-(i-floor((i-1)/n)*(n)-1)*n+1+floor((i-1)/n)
		Knn[j.temp,i]<-1
	}
	return(Knn)
}

n.ijkl.fast<-function(XX.obs){
n.obs<-dim(XX.obs)[2]
na.matrix<-NULL
for(i in 1:n.obs){
		XX.temp<-(1-is.na(XX.obs[,i]*XX.obs[,i:n.obs]))
		na.matrix<-cbind(na.matrix,XX.temp)

	}
	n1<-dim(na.matrix)[2]
	n2<-dim(na.matrix)[1]
	na.matrix.L<-numeric((n1^2))
	for(k in 1:n1){

		na.matrix.L[((k-1)*n1+1):(k*n1)]<-colSums(na.matrix[,k]*na.matrix)


		}
	na.matrix.L

}

