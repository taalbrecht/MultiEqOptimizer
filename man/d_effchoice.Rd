\name{d_effchoice}
\alias{d_effchoice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
D-efficiency calculation for a supplied discrete choice model matrix
}
\description{
Calculates the d-efficiency of a supplied model matrix per the methods outlined in Zwerina et.al (1996) and the ratio of the d-efficiency to the optimal d-efficiency.
}
\usage{
d_effchoice(CurrentMatrix, altvect, paramestimates = NULL, returncov = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CurrentMatrix}{
Model matrix with attribute factors coded using contr.sum coding. Matrix must be standardized before passing to function if standardization is desired. Opt out choices, if included, should be coded as a row of all zeroes for every parameter.
}
  \item{altvect}{
Vector with integer corresponding to the choice set that each row of the model matrix belongs to.
}
  \item{paramestimates}{
Vector of estimates for each effect (column) of model matrix sized corresponding to standardized model matrix. If this is not supplied, all effect estimates will be set to 0.
}
  \item{returncov}{
TRUE/FALSE binary value indicating Whether to return the covariance matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{d_eff}{
  Ratio of the d-efficiency of the supplied model to the optimal model.
}

  \item{det_calc}{
  Determinanat of the supplied model information matrix, det(M'M)
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
function (CurrentMatrix, altvect, paramestimates = NULL, returncov = FALSE)
{
    altnames <- unique(altvect)
    if (is.null(paramestimates)) {
        paramestimates <- rep(0, ncol(CurrentMatrix))
    }
    info_mat = matrix(rep(0, ncol(CurrentMatrix) * ncol(CurrentMatrix)),
        ncol(CurrentMatrix), ncol(CurrentMatrix))
    exputilities = exp(CurrentMatrix \%*\% paramestimates)
    for (k_set in 1:length(altnames)) {
        alternatives = which(altvect == altnames[k_set])
        p_set = exputilities[alternatives]/sum(exputilities[alternatives])
        p_diag = diag(p_set)
        middle_term <- p_diag - p_set \%o\% p_set
        full_term <- t(CurrentMatrix[alternatives, ]) \%*\% middle_term \%*\%
            CurrentMatrix[alternatives, ]
        info_mat <- info_mat + full_term
    }
    sigma_beta <- tryCatch(solve(info_mat, diag(ncol(CurrentMatrix))),
        error = function(x) diag(x = Inf, nrow = 2, ncol = 2))
    if (returncov == TRUE) {
        output <- list(d_eff = det(sigma_beta)^(-1/ncol(CurrentMatrix)),
            vcov = sigma_beta)
    }
    else {
        output <- det(sigma_beta)^(-1/ncol(CurrentMatrix))
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
