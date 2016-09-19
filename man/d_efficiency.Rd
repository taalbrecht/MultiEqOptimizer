\name{d_efficiency}
\alias{d_efficiency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
D-efficiency calculation for a supplied model matrix
}
\description{
Calculates the d-efficiency of a supplied data frame of input variables based on a supplied formula and the ratio of the d-efficiency to the optimal d-efficiency. Model matrix is standardized by the function.
}
\usage{
d_efficiency(CurrentMatrix, det_ref, input_formula, Input_range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CurrentMatrix}{
Data frame with column for each basic variable contained in the input formula.
}
  \item{det_ref}{
Reference optimal d-efficiency for calculating the d-efficiency ratio.
}
  \item{input_formula}{
Formula to be used for model matrix creation
}
  \item{Input_range}{
Range of basic input variables in the input data frame. Column names must match input_formula term names. Format is a matrix or data frame with column in Input_range matching each column name in CurrentMatrix, minimum value in first row, and maximum value in second row.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector is returned containing:
  \item{d_eff}{
  Ratio of the d-efficiency of the supplied model to the optimal model, (det(M'M)/det_ref)^(1/k) where M is the supplied model matrix of the supplied model and  det_ref is the supplied optimal information matrix determinant
}

  \item{vcov}{
  Variance covariance matrix.
}
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
function (CurrentMatrix, det_ref, input_formula, Input_range)
{
    CurrentMatrix <- data.frame(CurrentMatrix)
    modelmat <- model.matrix(input_formula, data = CurrentMatrix)
    modelmat <- standardize_cols(modelmat, colnames(modelmat[,
        2:ncol(modelmat)]), Input_range = Input_range)
    det_calc <- det(t(modelmat) \%*\% modelmat)
    if (det_calc < 0) {
        det_calc <- 0
    }
    d_eff <- ((det_calc/det_ref)^(1/(ncol(modelmat))))
    returnvect <- c(d_eff, det_calc)
    names(returnvect) <- c("D efficiency", "Info Matrix Determinant")
    return(returnvect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
