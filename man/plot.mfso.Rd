\name{plot.mfso}
\alias{plot.mfso}
\alias{points.mfso}
\alias{plotid.mfso}
\alias{hilight.mfso}
\alias{chullord.mfso}
\alias{boxplot.mfso}

\title{Plotting Routines for Multidimensional Fuzzy Set Ordinations}
\description{A set of routines for plotting, identifying, or 
highlighting points in a multidimensional fuzzy set ordination (MFSO).}

\usage{
\method{plot}{mfso}(x, dis=NULL, pch=1, \dots)
\method{points}{mfso}(x, overlay, col=2, pch=1, ...)
\method{plotid}{mfso}(x, dis=NULL, labels=NULL, \dots)
\method{hilight}{mfso}(x, overlay, cols = c(2, 3, 4, 5,
    6, 7), symbol = c(1, 3, 5), origpch = 1, blank = "#FFFFFF", \dots)
\method{chullord}{mfso}(x, overlay, cols = c(2, 3, 4, 5,
    6, 7), ltys = c(1, 2, 3), \dots)
\method{boxplot}{mfso}(x, \dots)
}

\arguments{
    \item{x}{an object of class \sQuote{mfso}}
    \item{dis}{an object of class \sQuote{dist} from \code{dist},
        \code{\link[labdsv]{dsvdis}}, or 
        \code{\link[vegan]{vegdist}}}
    \item{overlay}{a logical vector of the same length as the number of 
        points in the plot}
    \item{labels}{a vector of labels to print next to the identified points}
    \item{symbol}{an integer or vector of integers to control which symbols
        are printed in which order on the plot by specifying values to
        \code{pch}}
    \item{origpch}{a switch to explain to hilight what symbol was used in the
        original plot to be blanked out}
    \item{blank}{a color specification to control blanking out of symbols to
        match the background color}
    \item{ltys}{an integer or vector of integers to control the line styles
        of convex hull polygons}
    \item{pch}{the symbol to plot}
    \item{col}{the color to use for plotted symbols}
    \item{cols}{an integer vector for color order}
    \item{\dots}{arguments to pass to function points}
}

\details{
Multidimensional fuzzy set ordinations (MFSO) are almost inherently graphical,
and routines to facilitate plotting and overlaying are essential to work
effectively with them.

A multidimensional fuzzy set ordination object (an object of class
\sQuote{mfso}) generally contains at least two axes, and may contain many more.
The \code{plot} routine plots all possible axis pairs in order.  If a
\sQuote{dist} object is passed with the \sQuote{dis=} argument, the final panel
is a plot of the dissimilarity or distance matrix values on the X axis and the
pair-wise ordination distances on the Y axis with the correlation coefficient in
the upper left hand corner.

The \sQuote{points} function can be used to highlight or identify specific
points in the plot.  The \sQuote{points} function requires a logical vector
(TRUE/FALSE) of the same length as the number of points in the plot.  The
default behavior is to color the points with a respective TRUE value red.  It is
possible to control the color (with col=), size (with cex=) and symbol (with
pch=) of the points.  

The \sQuote{plotid} function can be used to label or identify specific points
with the mouse.  Clicking the left mouse button adjacent to a point causes 
the point to be labeled, 
offset in the direction of the click relative to the point.  Clicking the right
mouse button exits the routine.  The default
(labels=NULL) is to label points with the row number in the data.frame (or 
position in the vector) for the point.  Alternatively, specifying a vector of
labels (labels=) prints the respective labels.  If the data were derived from a
data.frame, the row.names of the data.frame are often a good choice, but the
labels can also be used with a factor vector to identify the distribution of
values of a factor in the ordination (but see hilight as well).

The \sQuote{hilight} function identifies the factor values of points in the
ordination, using color and symbols to identify unique values (up to 18 values by
default).  The colors and symbols used can be specified by the \sQuote{col=} and
\sQuote{symbol=} arguments, which should both be integers or integer vectors.
The default of colors 2, 3, 4, 5, 6, 7 and symbols 1, 3, 5 shows well in most
cases, but on colored backgrounds you may need to adjust \sQuote{col=}.  If you
have a factor with more than 18 classes you will need to augment the
\sQuote{symbol=} vector with more values.

The \sQuote{chullord} function plots convex hulls around all points sharing the
same value for a factor variable, and colors all points of that value to match.  
The convention on colors follows
\code{hilight}. 

The \sQuote{boxplot} function plots boxplots of the $\mu$ membership values
in the MFSO.


}

\value{none}

\note{The plotting and highlighting routines for mfso are designed to match the
same routines for other ordinations in package \code{labdsv}.
}

\author{
    David W. Roberts
    \email{droberts@montana.edu}
    \url{http://ecology.msu.montana.edu/R/labdsv/lab11/lab11.html}
}
  
\seealso{
\code{\link[labdsv]{plot.pco}}, \code{\link[labdsv]{points.pco}},
\code{\link[labdsv]{hilight.pco}},\code{\link[labdsv]{chullord.pco}}
}

\examples{
    require(labdsv) # to obtain access to data sets and dissimilarity function
    data(bryceveg) # vegetation data
    data(brycesite) # environmental data
    dis.bc <- dsvdis(bryceveg,'bray/curtis') # produce \sQuote{dist} object
    demo.mfso <- mfso(~elev+slope+av,dis.bc,data=brycesite)
    plot(demo.mfso)
    \dontrun{hilight(demo.mfso,brycesite$quad) # requires interaction}
}

\keyword{hplot}
\keyword{aplot}
\keyword{iplot}
