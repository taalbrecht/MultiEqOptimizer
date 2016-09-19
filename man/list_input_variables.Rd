\name{list_input_variables}
\alias{list_input_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List the basic input variables of a formula
}
\description{
Lists all unique fundamental input variables from formula passed to this function.
}
\usage{
list_input_variables(input_formula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_formula}{
A one or two-sided formula.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector of all unique basic variable names.
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
function (input_formula)
{
    library(nlme)
    Input_Func <- as.formula(paste(nlme::splitFormula(input_formula)))
    InputVarList <- all.vars(Input_Func)
    return(InputVarList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
