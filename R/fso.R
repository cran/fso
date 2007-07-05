fso <- function(...)
{
    UseMethod("fso")
}

mfso <- function(...)
{
    UseMethod("mfso")
}

core.fso <- function(x,sim)
{
    mu <- rep(0, length(x))
    a <- (x - min(x))/(max(x) - min(x))
    wghtb <- rep(sum(a), length(x))
    b <- (sim %*% a) - a
    b <- b/(wghtb - a)
    wghtc <- rep(sum(1 - a), length(x))
    c <- (sim %*% (1 - a)) - (1 - a)
    c <- c/(wghtc - (1 - a))
    mu <- as.vector((1. + (1. - c)^2 - (1. - b)^2)/2)
    list(mu = mu, data = x, r = cor(mu, x), var=deparse(substitute(x)))
}

fso.default <- function (x, dis, permute = FALSE, ...)
{
    if (missing(x)) {
        stop("You must provide a numeric vector or matrix")
    }
    if (class(dis) != "dist") {
        stop("You must provide and object of class 'dist'")
    }
    if (max(dis) > 1) {
        #cat("Rescaling to max=1.0 \n")
        dis <- dis/max(dis)
    }
    #cat("Converting dist object to similarity matrix \n")
    sim <- 1 - as.matrix(dis)
    if (is.vector(x)) {
        tmp <- core.fso(x, sim)
        tmp$var <- deparse(substitute(x))
        if (permute != FALSE && is.numeric(permute)) {
            cors <- rep(0, permute - 1)
            for (i in 1:(permute - 1)) {
                cors[i] <- core.fso(sample(x), sim)$r
            }
            tmp$p <- (sum(cors >= tmp$r) + 1)/(permute)
        }
        else {
            tmp$p <- cor.test(tmp$mu, tmp$data)$p.value
        }
        tmp$d <- cor(dist(tmp$mu),dis)
        class(tmp) <- "fso"
        return(tmp)
    }
    else if (is.matrix(x) || is.data.frame(x)) {
        if (is.matrix(x))
            x <- data.frame(x)
        r <- rep(0, ncol(x))
        mu <- matrix(0, nrow = nrow(x), ncol = ncol(x))
        p <- rep(0, ncol(x))
        d <- rep(0, ncol(x))
        for (i in 1:ncol(x)) {
            if (any(is.na(x[, i]))) {
                cat(paste("missing values in", names(x)[i], "\n"))
            }
            else if (is.factor(x[, i])) {
                cat(paste("variable ", names(x)[i], " is a factor\n"))
            }
            else {
                tmp <- core.fso(x[, i], sim)
                r[i] <- tmp$r
                mu[, i] <- tmp$mu
                if (permute != FALSE && is.numeric(permute)) {
                  cors <- rep(0, permute - 1)
                  for (j in 1:(permute - 1)) {
                    cors[j] <- core.fso(sample(x[, i]), sim)$r
                  }
                  p[i] <- (sum(cors >= tmp$r) + 1)/(permute)
                }
                else {
                  p[i] <- cor.test(tmp$mu, tmp$data)$p.value
                }
                d[i] <- cor(dist(tmp$mu),dis)
            }
        }
        out <- list(mu = mu, data = data.frame(x), r = r, p = p,
            d = d, var = names(x))
        class(out) <- "fso"
        return(out)
    }
}


fso.formula <- function (formula, dis, data, permute = FALSE, ...)
{
    if (!inherits(formula, "formula"))
        stop("Needs a model formula")
    if (missing(dis) || class(dis) != "dist") {
        stop("You must supply a (dis)similarity matrix\n")
    }
    if (missing(data)) {
        tmp <- data.frame(model.matrix(formula,data=parent.frame())[, -1])
        notmis <- as.numeric(row.names(tmp))
    }
    else {
        tmp <- data.frame(model.matrix(formula,data=data)[, -1])
        notmis <- match(row.names(tmp),row.names(data))
    }
    dis <- as.dist(as.matrix(dis)[notmis, notmis])
    fso.default(tmp, dis=dis, permute=permute)
}

mfso.default <- function(x,dis,permute=FALSE,scaling=1,lm=TRUE,notmis=NULL,...)
{
    relati <- function(z) {(z-min(z)) / (max(z)-min(z))}

    perm.fso <- function(dis,start,add,numitr,scaling) {
        rndcor <- rep(0,numitr)
        y <- as.dist(dis)
        for (i in 1:numitr) {
            if (is.null(start)) tmp <- data.frame(sample(add,replace=FALSE))
            else tmp <- data.frame(start,sample(add,replace=FALSE))
            tmp2 <- mfso(tmp,dis,scaling=scaling)
            a <- dist(tmp2$mu)
            rndcor[i] <- cor(y,a)
        }
        return(rndcor)
    }

    if (missing(x)){
        stop("You must provide a numeric matrix, or dataframe")
    }
    if (class(dis) != "dist") {
        stop("You must supply an object of class 'dist'")
    }
    if (max(dis) > 1.0) {
        cat("Rescaling to max=1.0 \n")
        dis <- dis/max(dis)
    }
    #cat("Converting dist object to similarity matrix \n")
    sim <- 1 - as.matrix(dis)
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop('Multi FSO only works on matrices or data frames')
    }
    if (is.matrix(x)) x <- data.frame(x)

    if (any(is.na(x))) {
        notmis <- apply(x,1,function(x){!any(is.na(x))})
        x <- x[notmis,]
        sim <- sim[notmis,notmis]
        dis <- as.dist(as.matrix(dis)[notmis,notmis])
    }

    r <- rep(0,ncol(x))
    mu <- matrix(0,nrow=nrow(x),ncol=ncol(x))
    gamma <- rep(1,ncol(x))
    p <- rep(0,ncol(x))
    for (i in 1:ncol(x)) {
        if (any(is.na(x[,i]))) {
            cat(paste('missing values in',names(x)[i],'\n'))
            is.na(r[i]) <- TRUE
        }
        else if (is.factor(x[,i])) {
            cat(paste('variable ',names(x)[i],' is a factor\n'))
            is.na(r[i]) <- TRUE
        }
        else {
            tmp <- core.fso(x[,i],sim)
            mu[,i] <- tmp$mu
        }
    }

    coor <- matrix(0,nrow=nrow(mu),ncol=ncol(mu))
    coor[,1] <- switch(scaling,
        mu[,1],
        relati(mu[,1]),
        relati(mu[,1]) * r[1] + (0.5 * (1-r[1])))
    r[1] <- cor(dist(coor[,1]),dis)
    if (permute) {
       rndcor <- perm.fso(dis=dis,start=NULL,add=x[,1],numitr=permute-1,scaling=scaling)
       p[1] <- (sum(rndcor>=r[1])+1)/(permute)
    }
    else {
       p[1] <- cor.test(dist(coor[,1]),dis)$p.value
    }
    if (ncol(mu) >= 2) {
        for (i in 2:ncol(mu)) {
            if (!is.na(r[i])) {
                coor[,i] <- switch(scaling,
                    mu[,i],
                    relati(mu[,i]),
                    relati(mu[,i]) * r[i] + (0.5 * (1-r[i])))
                if (lm) {
                    tmp <- lm(coor[,i]~coor[,1:(i-1)])
                    coor[,i] <- tmp$residuals + mean(coor[,i])
                    gamma[i] <- 1 - summary(tmp)$r.squared
                }
                r[i] <- cor(dist(coor[,1:i]),dis)
                if (permute) {
                    rndcor <- perm.fso(dis=dis,start=x[,1:(i-1)],
                                    add=x[,i],numitr=permute-1,scaling=scaling)
                    p[i] <- (sum(rndcor>=r[i])+1)/(permute)
                }
                else {
                    p[i] <- cor.test(dist(mu[,1:i]),dis)$p.value
                }
            }
        }
    }
    out <- list(mu=coor,data=data.frame(x),r=r,p=p,var=names(x),gamma=gamma,
                notmis=notmis)
    class(out) <- "mfso"
    attr(out,'permute') <- permute
    return(out)
}


mfso.formula <- function (formula, dis, data, permute = FALSE, lm = TRUE, scaling = 1, ...)
{
    if (!inherits(formula, "formula"))
        stop("Needs a model formula")
    if (missing(dis) || class(dis) != "dist") {
        stop("You must supply a (dis)similarity matrix\n")
    }
    if (missing(data)) {
        tmp <- data.frame(model.matrix(formula,data=parent.frame())[, -1])
        notmis <- as.numeric(row.names(tmp))
    }
    else {
        tmp <- data.frame(model.matrix(formula,data=data)[, -1])
        notmis <- match(row.names(tmp),row.names(data))
    }
    dis <- as.dist(as.matrix(dis)[notmis,notmis])
    mfso.default(tmp, dis=dis, permute=permute, lm=lm, scaling=scaling, notmis=notmis)
}


summary.fso <- function (object,...)
{
    if (class(object) != "fso") {
        stop("You must specify an object of class fso from fso()")
    }
    if (is.matrix(object$mu)) {
        num <- seq(1:length(object$r))
        out <- data.frame(num[rev(order(object$r, na.last = FALSE))],
               object$var[rev(order(object$r, na.last = FALSE))], 
               object$r[rev(order(object$r, na.last = FALSE))], 
               object$p[rev(order(object$r, na.last = FALSE))],
               object$d[rev(order(object$r, na.last = FALSE))])
        names(out) <- c("col", "variable", "r", "p","d")
        cat("\nFuzzy set statistics\n\n")
        print(format(out))
        cat("\nCorrelations of fuzzy sets\n\n")
        cormat <- data.frame(cor(object$mu))
        names(cormat) <- object$var
        row.names(cormat) <- object$var
        print(cormat)
    }
    else {
        tmp.mu <- summary(object$mu)
        tmp.dat <- summary(object$dat)
        cat(paste("fuzzy set ordination for", object$var, "\n\n"))
        cat(paste("r = ", formatC(object$r, width = 5), "\n"))
        cat(paste("p = ",formatC(object$p,width = 5), "\n\n"))
        cat(paste("            data       mu \n"))
        cat(paste("Min.   ", formatC(tmp.dat[[1]], width = 8),
            formatC(tmp.mu[[1]], width = 8), "\n"))
        cat(paste("1st Qu.", formatC(tmp.dat[[2]], width = 8),
            formatC(tmp.mu[[2]], width = 8), "\n"))
        cat(paste("Med.   ", formatC(tmp.dat[[3]], width = 8),
            formatC(tmp.mu[[3]], width = 8), "\n"))
        cat(paste("Mean   ", formatC(tmp.dat[[4]], width = 8),
            formatC(tmp.mu[[4]], width = 8), "\n"))
        cat(paste("3rd Qu.", formatC(tmp.dat[[5]], width = 8),
            formatC(tmp.mu[[5]], width = 8), "\n"))
        cat(paste("Max.   ", formatC(tmp.dat[[6]], width = 8),
            formatC(tmp.mu[[6]], width = 8), "\n"))
    }
}

summary.mfso <- function (object,...)
{
    if (class(object) != "mfso") {
        stop("You must specify an object of class mfso from mfso()")
    }
    num <- seq(1:length(object$r))
    incr <- rep(0, length(object$var))
    incr[1] <- object$r[1]
    if (length(incr) > 1) {
        for (i in 2:length(incr)) {
            incr[i] <- object$r[i] - object$r[i - 1]
        }
     }
     out <- data.frame(object$var, object$r, incr, object$p,
            object$gamma)
     names(out) <- c("variable", "cumulative_r", "increment",
            "p-value", "gamma")
     cat("\nStatistics for multidimensional fuzzy set ordination\n\n")
     print(format(out))
}

plot.fso <- function (x, which = "all", xlab = x$var, ylab = "mu(x)",
    title = "", r = TRUE, pch = 1, ...)
{
    if (class(x) != "fso") {
        stop("You must specify n object of class fso from fso()")
    }
    if (is.matrix(x$mu)) {
        if (which == "all") {
            for (i in 1:ncol(x$mu)) {
                if (is.na(x$r[i])) {
                  cat(paste("variable ", x$var[i], " has missing values \n"))
                }
                else {
                  plot(x$data[, i], x$mu[, i], xlab = x$var[i],
                    ylab = ylab, main = title)
                  if (r) {
                    ax <- min(x$data[, i])
                    ay <- max(x$mu[, i])
                    text(ax, ay, paste("r = ", format(x$r[i],
                      digits = 3)), pos = 4)
                  }
                }
                readline("\nHit Return for Next Plot\n")
            }
        }
        else if (is.numeric(which)) {
            for (i in which) {
                plot(x$data[, i], x$mu[, i], xlab = x$var[i],
                    ylab = ylab, main = title)
                if (r) {
                    ax <- min(x$data[, i])
                    ay <- max(x$mu[, i])
                    text(ax, ay, paste("r = ", format(x$r[i],
                      digits = 3)), pos = 4)
                }
                readline("\nHit Return for Next Plot\n")
            }
        }
        else {
            for (j in 1:ncol(x$mu)) {
                if (which == x$var[j])
                    plot(x$data[, j], x$mu[, j], xlab = x$var[j],
                      ylab = ylab, main = title)
                if (r) {
                    ax <- min(x$data[, j])
                    ay <- max(x$mu[, j])
                    text(ax, ay, paste("r = ", format(x$r[j],
                      digits = 3)), pos = 4)
                  }
            }
        }
    }
    else {
        plot(x$data, x$mu, xlab = x$var, ylab = ylab, main = title)
        if (r) {
            ax <- min(x$data)
            ay <- max(x$mu)
            text(ax, ay, paste("r = ", format(x$r, digits = 3)),
                pos = 4)
        }
    }
}

plot.mfso <- function (x,dis=NULL,pch=1,...)
{
    if (class(x) != 'mfso') {
        stop("You must pass an object of class 'mfso'")
    }
    for (i in 1:(ncol(x$mu)-1)) {
        for (j in (i+1):ncol(x$mu)) {
            plot(x$mu[,i],x$mu[,j],asp=1,
                xlim=range(apply(x$mu,2,range)),
                ylim=range(apply(x$mu,2,range)),
                xlab=names(x$data)[i],
                ylab=names(x$data)[j])
                readline("\nHit Return for Next Plot\n")
        }
    }
    if(!is.null(dis)) {
        tmp <- as.numeric(x$notmis)
        dis <- as.matrix(dis)[tmp,tmp]
        a <- dist(x$mu)
        y <- as.dist(dis)
        if (length(y) > 5000) pch <- "."
        plot(y,a,ylab="Ordination Distance",xlab="Matrix Dissimilarity",pch=pch)
        text(min(y),max(a),paste("r = ",format(cor(y,a),digits=3)),pos=4)
    }
}

plotid.fso <- function (x, which = "all", xlab = x$var, ylab = "mu(x)",
    title = "", r = TRUE, pch = 1, labels=NULL)
{
    if (class(x) != "fso") {
        stop("You must specify n object of class fso from fso()")
    }
    if (is.matrix(x$mu)) {
        if (which == "all") {
            for (i in 1:ncol(x$mu)) {
                if (is.na(x$r[i])) {
                  cat(paste("variable ", x$var[i], " has missing values \n"))
                }
                else {
                  plot(x$data[, i], x$mu[, i], xlab = x$var[i],
                    ylab = ylab, main = title)
                  if (r) {
                    ax <- min(x$data[, i])
                    ay <- max(x$mu[, i])
                    text(ax, ay, paste("r = ", format(x$r[i],
                      digits = 3)), pos = 4)
                  }
                }
                if (is.null(labels)) {
                    identify(x$data[, i], x$mu[, i])
                } 
                else {
                    identify(x$data[, i], x$mu[, i],labels)
                }
            }
        }
        else if (is.numeric(which)) {
            for (i in which) {
                plot(x$data[, i], x$mu[, i], xlab = x$var[i],
                    ylab = ylab, main = title)
                if (r) {
                    ax <- min(x$data[, i])
                    ay <- max(x$mu[, i])
                    text(ax, ay, paste("r = ", format(x$r[i],
                      digits = 3)), pos = 4)
                }
                if (is.null(labels)) {
                    identify(x$data[, i], x$mu[, i])
                }
                else {
                    identify(x$data[, i], x$mu[, i],labels)
                }
            }
        }
        else {
            for (j in 1:ncol(x$mu)) {
                if (which == x$var[j])
                    plot(x$data[, j], x$mu[, j], xlab = x$var[j],
                      ylab = ylab, main = title)
                if (r) {
                    ax <- min(x$data[, j])
                    ay <- max(x$mu[, j])
                    text(ax, ay, paste("r = ", format(x$r[j],
                      digits = 3)), pos = 4)
                  }
            }
        }
    }
    else {
        plot(x$data, x$mu, xlab = x$var, ylab = ylab, main = title)
        if (r) {
            ax <- min(x$data)
            ay <- max(x$mu)
            text(ax, ay, paste("r = ", format(x$r, digits = 3)),
                pos = 4)
        }
        if (is.null(labels)) {
            identify(x$data[, i], x$mu[, i])
        }
        else {
            identify(x$data[, i], x$mu[, i],labels)
        }
    }
}

plotid.mfso <- function (x,dis=NULL,labels=NULL,...) 
{
    if (class(x) != 'mfso') {
        stop("You must pass an object of class 'mfso'")
    }
    for (i in 1:(ncol(x$mu)-1)) {
        for (j in (i+1):ncol(x$mu)) {
            plot(x$mu[,i],x$mu[,j],asp=1,
                xlim=range(apply(x$mu,2,range)),
                ylim=range(apply(x$mu,2,range)),
                xlab=names(x$data)[i],
                ylab=names(x$data)[j])
            if (is.null(labels)) {
                identify(x$mu[,i],x$mu[,j])
            } 
            else {
                identify(x$mu[,i],x$mu[,j],labels)
            } 
        }
    }
    if(!is.null(dis)) {
        tmp <- as.numeric(row.names(x$data))
        dis <- as.matrix(dis)[tmp,tmp]
        a <- dist(x$mu)
        y <- as.dist(dis)
        if (length(y) > 5000) pch <- "."
        plot(y,a,ylab="Ordination Distance",xlab="Matrix Dissimilarity",pch=pch)
        text(min(y),max(a),paste("r = ",format(cor(y,a),digits=3)),pos=4)
    }
}

points.fso <- function (x, overlay, which = 1, col = 2, cex = 1, pch = 1, ...)
{
    if (class(x) != "fso") {
        stop("You must specify an object of class fso from fso()")
    }
    if (missing(overlay)) {
        stop("You must specify a overlay variable to highlight")
    }
    if (is.matrix(x$mu)) {
        points(x$data[, which][overlay], x$mu[, which][overlay],
            col = col, cex = cex)
    }
    else {
        points(x$data[overlay], x$mu[overlay], col = col,
            pch = pch, cex = cex)
    }
}

points.mfso <- function (x, overlay, col = 2, pch = 1, ...)
{
    if (class(x) != "mfso")
        stop("You must pass an object of class mfso")
    overlay <- overlay[as.numeric(row.names(x$data))]
    for (i in 1:(ncol(x$data)-1)) {
        for (j in (i+1):ncol(x$data)) {
            plot(x$mu[,i],x$mu[,j],asp=1,
                xlim=range(apply(x$mu,2,range)),
                ylim=range(apply(x$mu,2,range)),
                xlab=names(x$data)[i],
                ylab=names(x$data)[j])
            points(x$mu[overlay,i], x$mu[overlay, j],
                col = col, pch = pch)
            readline("\nHit Return for Next Plot\n")
        }
    }
}

boxplot.fso <- function (x, ...)
{
    if (class(x) != 'fso') {
        stop("You must pass an object of class 'fso' from fso()")
    }
    boxplot(data.frame(x$mu),names=x$var, ...)
}

boxplot.mfso <- function (x, ...)
{
    if (class(x) != 'mfso') {
        stop("You must pass an object of class 'mfso' from mfso()")
    }
    boxplot(data.frame(x$mu),names=x$var, ...)
}

chullord.fso <- function (x, overlay, which = 1, cols = c(2, 3, 4, 5,
    6, 7), ltys = c(1, 2, 3), ...)
{
    if (class(x) != "fso")
        stop("You must pass an object of class fso")
    if (inherits(overlay, c("partana", "pam", "slice")))
        overlay <- overlay$clustering
    else if (is.logical(overlay))
        overlay <- as.numeric(overlay)
    else if (is.factor(overlay))
        overlay <- as.numeric(overlay)
    pass <- 1
    layer <- 0
    lty <- ltys[pass]
    if (is.vector(x$data)) {
        for (i in 1:max(overlay, na.rm = TRUE)) {
            x <- x$dat[overlay == i & !is.na(overlay)]
            y <- x$mu[overlay == i & !is.na(overlay)]
            pts <- chull(x, y)
            layer <- layer + 1
            if (layer > length(cols)) {
                layer <- 1
                pass <- min(pass + 1, length(ltys))
            }
            col <- cols[layer]
            lty = ltys[pass]
            polygon(x[pts], y[pts], col = col, density = 0, lty = lty)
        }
    }
    else if (is.matrix(x$data) || is.data.frame(x$data)) {
        for (i in 1:max(overlay, na.rm = TRUE)) {
            x <- x$dat[, which][overlay == i & !is.na(overlay)]
            y <- x$mu[, which][overlay == i & !is.na(overlay)]
            pts <- chull(x, y)
            layer <- layer + 1
            if (layer > length(cols)) {
                layer <- 1
                pass <- min(pass + 1, length(ltys))
            }
            col <- cols[layer]
            lty = ltys[pass]
            polygon(x[pts], y[pts], col = col, density = 0, lty = lty)
        }
    }
}

chullord.mfso <- function (x, overlay, cols = c(2, 3, 4, 5,
    6, 7), ltys = c(1, 2, 3), ...)
{
    if (class(x) != "mfso")
        stop("You must pass an object of class mfso")
    if (inherits(overlay, c("partana", "pam", "slice")))
        overlay <- overlay$clustering
    else if (is.logical(overlay))
        overlay <- as.numeric(overlay)
    else if (is.factor(overlay))
        overlay <- as.numeric(overlay)
    overlay <- overlay[as.numeric(row.names(x$data))]
    for (i in 1:(ncol(x$data)-1)) {
        for (j in (i+1):ncol(x$data)) {
            pass <- 1
            layer <- 0
            lty <- ltys[pass]
            plot(x$mu[,i],x$mu[,j],asp=1,
                xlim=range(apply(x$mu,2,range)),
                ylim=range(apply(x$mu,2,range)),
                xlab=x$var[i],ylab=x$var[j],col=cols(overlay))
            for (k in 1:max(overlay, na.rm = TRUE)) {
                x <- x$mu[overlay==k,i]
                y <- x$mu[overlay==k,j]
                pts <- chull(x, y)
                layer <- layer + 1
                if (layer > length(cols)) {
                    layer <- 1
                    pass <- min(pass + 1, length(ltys))
                }
                col <- cols[layer]
                lty = ltys[pass]
                polygon(x[pts], y[pts], col = col, density = 0, lty = lty)
                points(x,y,col=col)
            }
        readline("\nHit Return for Next Plot\n")
        }
    }
}

hilight.fso <- function (x, overlay, which = 1, cols = c(2, 3, 4, 5,
    6, 7), symbol = c(1, 3, 5), origpch = 1, blank = "#FFFFFF",
    ...)
{
    if (class(x) != "fso")
        stop("You must pass an object of class fso")
    if (inherits(overlay, c("partana", "pam", "slice")))
        overlay <- overlay$clustering
    if (is.logical(overlay) || is.factor(overlay))
        overlay <- as.numeric(overlay)
    layer <- 0
    pass <- 1
    if (is.vector(x$data)) {
        points(x$data, x$mu,
            col = blank, pch = origpch)
        for (i in 1:max(overlay, na.rm = TRUE)) {
            layer <- layer + 1
            if (layer > length(cols)) {
                layer <- 1
                pass <- pass + 1
            }
            col <- cols[layer]
            pch <- symbol[pass]
            points(x$data[overlay==i], x$mu[overlay==i],
                col = col, pch = pch)
        }
    }
    else if (is.matrix(x$data) || is.data.frame(x$data)) {
        points(x$data[, which], x$mu[, which],
            col = blank, pch = origpch)
        for (i in 1:max(overlay, na.rm = TRUE)) {
            layer <- layer + 1
            if (layer > length(cols)) {
                layer <- 1
                pass <- pass + 1
            }
            col <- cols[layer]
            pch <- symbol[pass]
            points(x$data[overlay==i, which], x$mu[overlay==i, which],
                col = col, pch = pch)
        }
    }
}

hilight.mfso <- function (x, overlay, cols = c(2, 3, 4, 5,
    6, 7), symbol = c(1, 3, 5), origpch = 1, blank = "#FFFFFF", ...)
{
    if (class(x) != "mfso")
        stop("You must pass an object of class mfso")
    if (inherits(overlay, c("partana", "pam", "slice")))
        overlay <- overlay$clustering
    if (is.logical(overlay) || is.factor(overlay))
        overlay <- as.numeric(overlay)
    for (i in 1:(ncol(x$data)-1)) {
        for (j in (i+1):ncol(x$data)) {
        layer <- 0
        pass <- 1
            plot(x$mu[,i],x$mu[,j],asp=1,
                xlim=range(apply(x$mu,2,range)),
                ylim=range(apply(x$mu,2,range)),
                xlab=names(x$data)[i],
                ylab=names(x$data)[j])
            points(x$mu[,i], x$mu[, j],
                col = blank, pch = origpch)
                for (k in 1:max(overlay, na.rm = TRUE)) {
                    layer <- layer + 1
                    if (layer > length(cols)) {
                        layer <- 1
                        pass <- pass + 1
                    }
                    col <- cols[layer]
                    pch <- symbol[pass]
                    points(x$mu[overlay==k, i], x$mu[overlay==k, j],
                        col = col, pch = pch)
            }
        readline("\nHit Return for Next Plot\n")
        }
    }
}

step.mfso <- function (dis,start,add,numitr=100,scaling=1)
{
    perm.fso <- function(dis,start,add,numitr,scaling) {
        rndcor <- rep(0,numitr)
        for (i in 1:numitr) {
            if (is.null(start)) tmp <- data.frame(sample(add,replace=FALSE))
            else tmp <- data.frame(start,sample(add,replace=FALSE))
            tmp2 <- mfso(tmp,dis,scaling=scaling)
            a <- dist(tmp2$mu)
            y <- as.dist(dis)
            rndcor[i] <- cor(y,a)
        }
        return(rndcor)
    }

    if (is.null(start)) {
        missing <- seq(1,nrow(add))[apply(add,1,function(x){any(is.na(x))})]
        if (length(missing) > 0) {
            add <- data.frame(add[-missing,])
            dis <- as.dist(as.matrix(dis)[-missing,-missing])
        }
        basval <- 0

        res <- data.frame(names(add),rep(NA,ncol(add)),rep(NA,ncol(add)))
        names(res) <- c('variable','delta_cor','p_val')
        for (i in 1:ncol(add)) {
            tmp <- data.frame(add[,i])
            tmp <- mfso(tmp,dis,scaling=scaling)
            a <- dist(tmp$mu)
            y <- as.dist(dis)
            res[i,2] <- cor(a,y)
            rndcor <- perm.fso(dis,start,add[,i],numitr-1,scaling)
            res[i,3] <- (sum(rndcor>=res[i,2])+1)/(numitr)
        }
        cat(paste('Baseline = ',format(basval,3),'\n'))
        print(res)
    } else {
        full <- data.frame(start,add)
        missing <- seq(1,nrow(full))[apply(full,1,function(x){any(is.na(x))})]
        if (length(missing) > 0) {
            start <- data.frame(start[-missing,])
            add <- data.frame(add[-missing,])
            dis <- as.dist(as.matrix(dis)[-missing,-missing])
        }
        if (ncol(start) == 1) {
            basval <- cor(dist(fso(start[,1],dis)$mu),dis)
        } else {
            base <- mfso(start,dis,scaling=scaling)
            a <- dist(base$mu)
            y <- as.dist(dis)
            basval <- cor(a,y)
        }

        res <- data.frame(names(add),rep(NA,ncol(add)),rep(NA,ncol(add)))
        names(res) <- c('variable','delta_cor','p_val')
        for (i in 1:ncol(add)) {
            tmp <- data.frame(start,add[,i])
            tmp <- mfso(tmp,dis,scaling=scaling)
            if (is.vector(tmp$mu)) a <- dist(tmp$mu)
            else a <- dist(tmp$mu)
            y <- as.dist(dis)
            res[i,2] <- cor(a,y)
            rndcor <- perm.fso(dis,start,add[,i],numitr-1,scaling)
            res[i,3] <- (sum(rndcor>=res[i,2])+1)/(numitr)
            res[i,2] <- res[i,2] - basval
        }
        cat(paste('Baseline = ',format(basval,3),'\n'))
        print(res)
    }
}
