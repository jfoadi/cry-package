\name{UnitCell-class}
\Rdversion{1.1}
\docType{class}
\alias{UnitCell-class}
\alias{changeCellParameters,UnitCell-method}
\alias{checkUnitCell,UnitCell-method}
\alias{computeCellVolume,UnitCell-method}
\alias{computeReciprocalUnitCell,UnitCell-method}
\alias{extractCellParameters,UnitCell-method}
\alias{print,UnitCell-method}

\title{Class \code{"UnitCell"}}
\description{
A \code{UnitCell} object includes all main features of a crystallographic unit cell. Routine crystallographic operations
that use unit cell parameters can be carried out using the methods of this class
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("UnitCell", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{a}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{c}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alpha}:}{Object of class \code{"Angle"} ~~ }
    \item{\code{beta}:}{Object of class \code{"Angle"} ~~ }
    \item{\code{gamma}:}{Object of class \code{"Angle"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{changeCellParameters}{\code{signature(object = "UnitCell")}: ... }
    \item{checkUnitCell}{\code{signature(object = "UnitCell")}: ... }
    \item{computeCellVolume}{\code{signature(object = "UnitCell")}: ... }
    \item{computeReciprocalUnitCell}{\code{signature(object = "UnitCell")}: ... }
    \item{extractCellParameters}{\code{signature(object = "UnitCell")}: ... }
    \item{print}{\code{signature(x = "UnitCell")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("UnitCell")
}
\keyword{classes}
