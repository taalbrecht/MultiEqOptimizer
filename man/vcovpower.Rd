\name{vcovpower}
\alias{vcovpower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation for model coefficients
}
\description{
Calculates the power for model coefficients based on a supplied alpha error, beta error, covariance matrix, and desired resolution for each coefficient.
}
\usage{
vcovpower(vcovmat, detectdiff, test_alpha = 0.05, test_beta = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcovmat}{
Covariance matrix to use for power calculations. Alternatively, this can be a vector taken from the diagonal of a covariance matrix.
}
  \item{detectdiff}{
Vector of detectable difference to use for power calculations. This should be equal to the absolute difference in each coefficient that the user would like to detect.
}
  \item{test_alpha}{
Vector of alpha errors to use for each power calculation. Defaults to 0.05.
}
  \item{test_beta}{
Vector of beta errors to use for each power calculation. Defaults to 0.20.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Data frame with estimated minimum sample size required to estimate each parameter as well as the input values provided by the user for each calculation.
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
function (vcovmat, detectdiff, test_alpha = 0.05, test_beta = 0.2)
{
    if (length(test_alpha) == 1) {
        test_alpha <- rep(test_alpha, ncol(CurrentMatrix))
    }
    if (length(test_beta) == 1) {
        test_beta <- rep(test_beta, ncol(CurrentMatrix))
    }
    z_one_minus_alpha <- qnorm(1 - test_alpha)
    z_one_minus_beta <- qnorm(1 - test_beta)
    minsamplesize <- ((z_one_minus_beta + z_one_minus_alpha) *
        sqrt(diag(vcovmat))/abs(detectdiff))^2
    output <- data.frame(rownames(vcovmat), detectdiff, minsamplesize,
        test_alpha, test_beta)
    colnames(output) <- c("Effect Name", "Difference to Detect",
        "Minimum Sample Size", "Alpha", "Beta")
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
