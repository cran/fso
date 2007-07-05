\name{mfso}
\alias{mfso}
\alias{mfso.default}
\alias{mfso.formula}
\alias{summary.mfso}
\title{Multidimensional Fuzzy Set Ordination}
\description{A multidimensional extension of fuzzy set ordination (FSO) that
constructs a multidimensional ordination by mapping samples from fuzzy
topological space to Euclidean space for
statistical analysis.  MFSO can be used in exploratory or testing modes.
}
\usage{
\method{mfso}{formula}(formula,dis,data,permute=FALSE,lm=TRUE,scaling=1,\dots)
\method{mfso}{default}(x,dis,permute=FALSE,scaling=1,lm=TRUE,notmis=NULL,\dots)
\method{summary}{mfso}(object,\dots)
}

\arguments{
    \item{formula}{Model formula, with no left hand side.  Right hand side gives
the independent variables to use in fitting the model}

    \item{dis}{a dist object of class \sQuote{dist} returned from
        \code{\link[mva]{dist}}, \code{\link[vegan]{vegdist}}, or
        \code{\link[labdsv]{dsvdis}}}

    \item{data}{a data frame containing the variables specified in the formula}

    \item{permute}{a switch to control how the probability of correlations
is calculated.  permute=FALSE (the default) uses a parametric Z distribution
approximation; permute=n permutes the independent variables (permute-1) times
and estimates the probability as (m+1)/(permute) where m is the number of
permuted correlations greater than or equal to the observed correlation.}

    \item{lm}{a switch to control scaling of axes after the first axis.  If
lm=TRUE (the default) each axis is constructed independently, and then
subjected to a Gram-Schmidt orthogonalization to all previous axes to preserve only the the variability
that is uncorrelated with all previous axes.  If lm=FALSE, the full extent of all axes is
preserved without correcting for correlation with previous axes.}

    \item{scaling}{a switch to control how the initial fuzzy set axes are
scaled: 1 = use raw $\mu$ membership values, 2 = relativize $\mu$ values [0,1],
3 = relativize $\mu$ values [0,1] and multiply by respective correlation
coefficient. 
}

    \item{x}{a quantitative matrix or dataframe.  One axis will be fit for each
column} 
    \item{notmis}{a vector passed from the formula version of mfso to control
               for missing values in the data}

    \item{object}{an object of class \sQuote{mfso}}
    \item{\dots}{generic arguments for future use}
}
\details{mfso performs individual fso calculations on each column of a
data frame or matrix, and then combines those fso axes into a higher dimensional
object.  The algorithm of fuzzy set ordination is described in the help
file for fso.  The key element in mfso is the Gram-Schmidt orthogonalization,
which ensures that
each axis is independent of all previous axes.  In practice, each axis is
regressed against all previous axes, and the residuals are retained as the data.}


\value{an object of class \sQuote{mfso} with components:
\item{mu}{a matrix of fuzzy set memberships of samples, analogous to the coordinates of
the samples along the axes, one column for each axis}
\item{data}{a dataframe containing the independent variables as columns}
\item{r}{a vector of correlation coefficients, one for each axis in order}
\item{p}{a vector of probabilities of observing correlations as high as
observed}
\item{var}{a vector of variables names used in fitting the model}
\item{gamma}{a vector of the fraction of variance for an axis that is
independent of all previous axes}
}
\note{MFSO is an extension of single dimensional fuzzy set ordination
designed to achieve low dimensional representations of a dissimilarity or
distance matrix as a function of environmental or experimental variables.
Although it is not technically a ``constrained ordination,'' in practice its use
is similar to \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.

If you set lm=FALSE, an mfso is
equivalent to an fso, but the plotting routines differ.  For an mfso, the plotting
routine plots each axis against all others in turn; for an fso the plotting routine
plots each axis against the environmental or experimental variable it is derived
from.}
}
\references{
Roberts, D.W.  2007.  Statistical analysis of multidimensional 
fuzzy set ordinations.  Ecology in press. 

\url{http://ecology.montana.msu.edu/R/labdsv/lab11/lab11.html}
}

\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\seealso{\code{\link[vegan]{cca}},\code{\link[vegan]{capscale}}}
\examples{
    require(labdsv)
    data(bryceveg) # returns a vegetation dataframe
    data(brycesite) # returns a dataframe of environmental variables
    dis.bc <- dsvdis(bryceveg,'bray/curtis')
                  # returns an object of class sQuote{dist}
    demo.mfso <- mfso(~elev+slope+av,dis.bc,data=brycesite) # creates the mfso
    plot(demo.mfso)
}
\keyword{multivariate}

