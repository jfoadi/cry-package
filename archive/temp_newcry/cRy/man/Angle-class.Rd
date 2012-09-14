\name{Angle-class}
\Rdversion{1.1}
\docType{class}
\alias{Angle-class}
\alias{changeAngle,Angle-method}
\alias{checkAngle,Angle-method}
\alias{degToRad,Angle-method}
\alias{print,Angle-method}
\alias{radToDeg,Angle-method}

\title{Class \code{"Angle"}}
\description{
Angles in degrees or radians. Their value is always forced withing the 0-360 or 0-2*pi range,
even when the input value lies outside this range
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Angle", ...)}. 

The user-friendly way to create an object of class
\code{Angle} is through the generic function \code{createAngle}
}
\section{Slots}{
  \describe{
    \item{\code{ang}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rad_flag}:}{Object of class \code{"logical"} ~~

      If rad_flag is TRUE, the slot value is intended to be the angle value in radians. Default is for rad_flag to
      be FALSE, i.e. for the angle value to be intended in degrees
    }
  }
}
\section{Methods}{
 \describe{
  \item{print}{\code{signature(x = "Angle")}:
   Prints value of angle and if the angle is in degrees or radians.}
  \item{checkAngle}{\code{signature(object = "Angle")}:
   Checks whether the object of class Angle has been correctly created (see ).}
  \item{changeAngle}{\code{signature(object = "Angle")}:
   Change angle's value and/or units (see ).}
  \item{degToRad}{\code{signature(object = "Angle")}:
   Transforms value of object from degrees to radians.}
  \item{radToDeg}{\code{signature(object = "Angle")}:
   Transforms value of object from radians to degrees.}
 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
James Foadi and David G. Waterman
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
showClass("Angle")
}
\keyword{classes}
