\name{d_effchoiceupdate}
\alias{d_effchoiceupdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast multinomial logistic d-efficiency update function.
}
\description{
Fast d-efficiency update function for a single point change to a multinomial logistic model. Supply a subset of the existing design corresponding to a single question from the model matrix and the information matrix for the remaining design for all other questions.
}
\usage{
d_effchoiceupdate(CurrentMatrix, paramestimates, info_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CurrentMatrix}{
Subset of model matrix with attribute factors coded using contr.sum coding. Matrix must be standardized before passing to function if standardization is desired. Opt out choices, if included, should be coded as a row of all zeroes for every parameter. Subset should be 1 row per alternative and subset should only include all alternatives for one question only.
}
  \item{paramestimates}{
Vector of estimates for each effect (column) of model matrix sized corresponding to standardized model matrix. If this is not supplied, all effect estimates will be set to 0.
}
  \item{info_mat}{
The information matrix for all other questions in the model matrix that were not passed as the CurrentMatrix argument.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
D-efficiency of the supplied model.
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
