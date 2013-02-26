\name{fso}
\alias{fso}
\alias{fso.default}
\alias{fso.formula}
\alias{summary.fso}
\title{Fuzzy Set Ordination}
\description{
   Computes a fuzzy set for samples along a specified environmental or
   experimental gradient based on
   sample similarities and gradient values as weights.  The fuzzy set memberships
   represent the degree to which a sample is similar to one end of the
   gradient while not similar to the other.
}
\usage{
\method{fso}{formula}(formula,dis,data,permute=FALSE,\dots)
\method{fso}{default}(x,dis,permute=FALSE,\dots)
\method{summary}{fso}(object,\dots)
}
\arguments{
    \item{formula}{a formula in the form of ~x+y+z (no LHS)}
    \item{dis}{a dist object such as that returned by \code{dist},
        \code{\link[labdsv]{dsvdis}}, or \code{\link[vegan]{vegdist}}}
    \item{data}{a data frame that holds variables listed in the formula}
    \item{permute}{if FALSE, estimate probabilities from Z distribution for correlation;
          if numeric, estimate probabilities from permutation of input}
    \item{x}{a numerical vector, a matrix, or numeric dataframe}
    \item{object}{an object of class \sQuote{fso}}
    \item{\dots}{generic arguments for future use}
}
\details{The algorithm converts the input to a full symmetric similarity matrix
    and bounds [0,1] (if necessary).  It then calculates
    several fuzzy sets: 
    \deqn{\mu_a(i) = (x_i-min(x))/(max(x)-min(x))}{mu_a(i) = (x_i-min(x))/(max(x)-min(x))}
    \deqn{\mu_b(i) = 1 - \mu_a(i)}{mu_b(i) = 1 - mu_a(i)}
    \deqn{\mu_c(i) = \Bigl(\sum_{j\ne i} \mu_a(j) \times y_{ij}\Bigr) / \sum_{j\ne i} \mu_a(j)}{mu_c(i) = (sum_j mu_a(j) * y_{ij}) / sum_j mu_a(j)}
    \deqn{\mu_d(i) = \Bigl(\sum_{j\ne i} \mu_b(j) \times y_{ij}\Bigr) / \sum_{j\ne i} \mu_b(j)}{mu_d(i) = (sum_j mu_b(j) *  y_{ij}) / sum_j mu_b(j)}
    \deqn{\mu_e(i) = \Bigl(1 + \bigl(\mu_d(i))^2 - (\mu_c(i)\bigr)^2\Bigr)/2}{mu_e(i) = (1 + mu_d(i)^2 - mu_c(i)^2)/2}

    where \eqn{y_{i,j}} is the similarity of sample i to sample j.

    A separate fuzzy set ordination is calculated for each term in the
    formula.  If x is a matrix or dataframe a separate fuzzy set ordination is
    calculated for each column or field.

    If permute is numeric, the permutation is performed permute-1 times,
    and the probability is estimated as 
    \eqn{(correlations $>=$ observed + 1)/permute}{(correlations >= observed + 1)/permute.}
}

\value{
    An object of class \sQuote{fso} which has the following elements:
    \item{mu}{the fuzzy membership values for individual plots in the fuzzy
    set.  If x is a matrix or dataframe then mu is also a matrix of the 
    same dimension.}
    \item{data}{a copy of data vector or matrix y}
    \item{r}{the correlation between the original vector and the fuzzy
    set.  If x is a matrix or dataframe then r is a vector with length equal
    to the number of columns in the matrix or dataframe.}
    \item{p}{the probability of obtaining a correlation between the data
    and fuzzy set as large as observed}
    \item{d}{the correlation of pair-wise distances among each fuzzy set
             compared to the dissimilarity matrix from which the fso was constructed}
    \item{var}{the variable name(s) from matrix y}

}
\references{Roberts, D.W. 1986. Ordination on the basis of fuzzy set theory.
    Vegetatio 66:123-131.

    Roberts, D.W. 2007.  Statistical analysis of multidimensional fuzzy set 
    ordinations.  Ecology 89:1246-1260.

    Roberts, D.W.  2009.  Comparison of multidimensional fuzzy set ordination
    with CCA and DB-RDA.  Ecology.  90:2622-2634.

\url{http://ecology.msu.montana.edu/labdsv/R/labs/lab11/lab11.html}
} 
\note{Fuzzy set ordination is a method of multivariate analysis employed in 
vegetation analysis.

   fso can be run with the first argument either a dataframe or a formula
(with no left hand side).  The formula version has distinct advantages:

1) The data= argument allows the user to specify a data frame
          containing the variables of interest.  In this way variables
          need not be local.

2) The formula version handles categorical variables by converting
          them to dummy variables.  In the default version, all variables 
          must be quantitative or binary.

3) The formula version is somewhat more graceful about handling
          missing values in the data.
}
\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\examples{
    library(labdsv)
    data(bryceveg)
    data(brycesite)
    dis <- dsvdis(bryceveg,'bray/curtis')
    elev.fso <- fso(brycesite$elev,dis)
    elev.fso <- fso(~elev,dis,data=brycesite)
    plot(elev.fso)
    summary(elev.fso)
}
\keyword{nonparametric}

