\name{plot.fso}
\alias{plot.fso}
\alias{points.fso}
\alias{plotid.fso}
\alias{hilight.fso}
\alias{chullord.fso}
\alias{boxplot.fso}

\title{Plotting Routines for Fuzzy Set Ordinations}
\description{A set of routines for plotting, highlighting points,
or identifying the distribution of a third variable on an fso.} 

\usage{
\method{plot}{fso}(x, which="all", xlab = x$var, ylab="mu(x)",
    title="",r=TRUE,pch=1,\dots)
\method{points}{fso}(x, overlay, which="all", col=2, cex=1, pch=1, ...)
\method{plotid}{fso}(ord, which="all", xlab=ord$var, ylab="mu(x)", title="",
    r=TRUE, pch=1, labels=NULL, \dots)
\method{hilight}{fso}(ord, overlay, which=1, cols = c(2, 3, 4, 5,
    6, 7), symbol = c(1, 3, 5), \dots)
\method{chullord}{fso}(ord, overlay, which = 1, cols = c(2, 3, 4, 5,
    6, 7), ltys = c(1, 2, 3), \dots)
\method{boxplot}{fso}(x, \dots)
}

\arguments{
    \item{x}{an object of class \sQuote{fso}}
    \item{ord}{an object of class \sQuote{fso}}
    \item{which}{a switch to control which axis is plotted}
    \item{r}{a switch to control printing the correlation coefficient 
        in the plot}
    \item{fso}{an object of class \sQuote{fso} from \code{fso}}
    \item{overlay}{a logical vector of the same length as the number of points in the plot}
    \item{labels}{a vector of labels to print next to the identified points}
    \item{symbol}{an integer or vector of integers to control which symbols
        are printed in which order on the plot by specifying values to
        \code{pch}}
    \item{ltys}{an integer or vector of integers to control the line styles
        of convex hull polygons}
    \item{xlab}{text label for X axis}
    \item{ylab}{text label for Y axis}
    \item{title}{an overall title for the plot (equivalent to main)}
    \item{pch}{the symbol for plotting}
    \item{col}{the color for plotted symbols}
    \item{cex}{the character expansion factor (font size)}
    \item{cols}{an integer vector specifying color order}
    \item{\dots}{arguments to pass to the underlying plot function}
}

\details{
Fuzzy set ordinations (FSO) are almost inherently graphical, and routines to
facilitate plotting and overlaying are essential to work effectively with them.

A fuzzy set ordination object (an object of class \sQuote{fso}) may contain one
or more axes.  In the simplest case, for a single-axis fso, the plot routine
plots the underlying raw data on the X axis and the fuzzy set memberships on the Y
axis, including by default the correlation coefficient in the upper left corner.
For fsos containing multiple axes, the default (which="all") is to plot the raw
data on the X axis, the respective fuzzy set memberships on the Y axis, plotting all
axes in turn with a prompt to move to the next panel.  This is often effective.
It is also possible to plot a single panel out of the set of axes, specifying
the axis as an integer with, e.g., "which = 2."

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
default).  The colors and symbols used can be specified by the \sQuote{cols=} and
\sQuote{symbol=} arguments, which should both be integers or integer vectors.
The default of colors 2, 3, 4, 5, 6, 7 and symbols 1, 3, 5 shows well in most
cases, but on colored backgrounds you may need to adjust \sQuote{cols=}.  If you
have a factor with more than 18 classes you will need to augment the
\sQuote{symbol=} vector with more values.

The \sQuote{chullord} function plots convex hulls around all points sharing the
same value for a factor variable, and colors all points of that value to match.
The convention on colors follows \sQuote{hilight}. 

The \sQuote{boxplot} function plots boxplots of the \eqn{\mu} membership values
for the fuzzy sets in the fso.


}


\note{The plotting and highlighting routines for fso are designed to match the
same routines for other ordinations in package \code{labdsv}.
}

\author{
    David W. Roberts
    \email{droberts@montana.edu}
}
\references{
    \url{http://ecology.msu.montana.edu/labdsv/R/labs/lab11/lab11.html}
}

\examples{
    require(labdsv) # to obtain access to data sets and dissimilarity function
    data(bryceveg) # vegetation data
    data(brycesite) # environmental data
    dis.bc <- dsvdis(bryceveg,'bray/curtis') # produce \sQuote{dist} object
    demo.fso <- fso(~elev+slope+av,dis.bc,data=brycesite)
    \dontrun{plot(demo.fso)}
    \dontrun{hilight(demo.mfso,brycesite$quad)}
}

\keyword{hplot}
\keyword{aplot}
\keyword{iplot}
