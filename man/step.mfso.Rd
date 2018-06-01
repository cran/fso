\name{step.mfso}
\alias{step.mfso}
\title{Step-Wise Forward Variable Selection in a Multivariate Fuzzy Set Ordination}
\description{A simple routine to screen variables for addition to a multivariate
fuzzy set ordination (MFSO).  The routine operates by adding variables one at a time
to an existing MFSO (which can be NULL), and calculating the correlation coefficient
between the underlying dissimilarity matrix (object of class \sQuote{dist}) and the 
pair-wise distances in the MFSO ordination.
}
\usage{
step.mfso(dis,start,add,numitr=100,scaling=1)
}
\arguments{
	\item{dis}{a dissimilarity of distance object from \code{dist},
                   \code{\link{dsvdis}}, or \sQuote{vegdist}
                    or other \sQuote{dist} object}
	\item{start}{either NULL (to find the first variable to add) or a data.frame
                     of binary or quantitative variables to use in the base model}
        \item{add}{a data.frame of binary or quantitative variables to screen for 
                    addition to the model}
        \item{numitr}{the number of random permutations of a vector to use in establishing
                     the probability of observing as large an increase in correlation
                     as observed}
        \item{scaling}{the scaling parameter to pass along to \code{mfso}}
}

\details{

\sQuote{mfso} is intended as a tool for analysis of multiple competing hypotheses,
and the analyst is expected to have a priori models to compare.  Nonetheless,
\sQuote{mfso} can be used in a hypothesis generating variable screening mode by
maximizing the correlation between the underlying dissimilarity matrix and the
pair-wise distances in the \sQuote{mfso} ordination.

The step.mfso function is an inelegant approach to step-wise forward variable
selection in \code{mfso}.  It considers each variable offered in turn, calculates the
\code{mfso} resulting from adding that variable to the given \code{mfso}, permutes that
variable \sQuote{numitr} times, and determines a probability of observing as large
an increase in correlation as observed.  After testing all variables for inclusion, it
simply prints a table of the calculations, and the analyst has to rerun the routine
adding the selected variable to data.frame \sQuote{start} and deleting it from \sQuote{add}.

While it would be nice to automate the production of the step-wise \sQuote{mfso}, to date
I have only implemented this limited function.  In addition, model parsimony is ensured by
the permutation routine, rather than an AIC-based approach, and doesn't directly
penalize for degrees of freedom (number of variables).

}

\value{
Produces a table of the analysis but does not produce any objects
}
\references{Roberts, D.W.  2007.  Statistical analysis of multidimensional fuzzy
set ordinations.  Ecology in press

\url{http://ecology.msu.montana.edu/labdsv/R/labs/lab11/lab11.html}
}
\author{
	David W. Roberts
	\email{droberts@montana.edu}
}
\examples{
	\dontrun{require(labdsv) # make data available}
        \dontrun{data(bryceveg)  # get vegetation data}
        \dontrun{data(brycesite) # get environmental data}
        \dontrun{dis.bc <- dsvdis(bryceveg,'bray.curtis') # produce dist object}
        \dontrun{attach(brycesite)   # make variables easily available}
        \dontrun{step.mfso(dis.bc,start=NULL,add=data.frame(elev,slope,av))}
        \dontrun{step.mfso(dis.bc,start=data.frame(elev),add=data.frame(slope,av))}
}
\keyword{multivariate}

