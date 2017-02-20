\name{d_efficiencysimple}
\alias{d_efficiencysimple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple d-efficiency calculation for a supplied model matrix
}
\description{
Calculates the d-efficiency of a supplied model matrix (det(M'M)^(1/ncol(ModelMatrix))) for a linear model
}
\usage{
d_efficiencysimple(CurrentMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CurrentMatrix}{
Model matrix with attribute factors coded using contr.sum coding. Matrix must be standardized before passing to function if standardization is desired.
}
  \item{returncov}{
TRUE/FALSE - Value specifying whether to return the variance-covariance matrix. Defaults to FALSE
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The d-efficiency is returned based on a linear model

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (CurrentMatrix, det_ref)
{
    det_calc <- det(t(CurrentMatrix) \%*\% CurrentMatrix)
    if (det_calc < 0) {
        det_calc <- 0
    }
    d_eff <- ((det_calc/det_ref)^(1/(ncol(CurrentMatrix))))
    returnvect <- c(d_eff, det_calc)
    names(returnvect) <- c("D efficiency", "Info Matrix Determinant")
    return(returnvect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
