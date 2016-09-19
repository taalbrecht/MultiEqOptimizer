\name{algebraic_range}
\alias{algebraic_range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds the maximum and minimum values for all supplied terms in the supplied algebraic formula.
}
\description{
Uses the solver to attempt to find a minimum and maximum value for all supplied algebraic terms in the input formula based on the supplied basic variable ranges. If the algebraic terms are complex shapes or step functions, this may not find the absolute minima and maxima and may, instead, find a local minima or maxima.
}
\usage{
algebraic_range(base_var_range, algebraic_formula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base_var_range}{
Range of basic input variables that are contained in the algebraic formula. Column names must match basic variable names in the algebraic formula. Format is a matrix or data frame with column in Input_range matching each column name in StartingMat and column_names, minimum value in first row, and maximum value in second row. Values may also be supplied for some algebraic combinations. If supplied, these will be passed through to the output and will not be solved.
}
  \item{algebraic_formula}{
A one or two sided formula. Algebraic combinations on the input side of the equation will have ranges found.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with column names matching algebraic combinations on the input side of formula where row 1 = minimum values and row 2 = maximum values.
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
function (base_var_range, algebraic_formula)
{
    algebraic_input_terms <- colnames(attr(terms(algebraic_formula),
        "factors"))
    AlgebraicRange <- matrix(data = 0, nrow = 2, ncol = length(algebraic_input_terms))
    colnames(AlgebraicRange) <- algebraic_input_terms
    rownames(AlgebraicRange) <- c("minimum", "maximum")
    for (i in 1:length(algebraic_input_terms)) {
        transmult <- gsub(":", "*", algebraic_input_terms[i])
        optimexp <- parse(text = transmult)
        base_inputs <- all.vars(optimexp)
        optimfunc <- function(optimvect) {
            basevarmat <- matrix(nrow = 1, ncol = length(base_inputs),
                data = 0)
            basevarmat <- data.frame(basevarmat)
            colnames(basevarmat) <- base_inputs
            basevarmat[1, ] <- optimvect
            output <- eval(optimexp, basevarmat)
            return(output)
        }
        minvect <- as.vector(base_var_range[1, base_inputs])
        maxvect <- as.vector(base_var_range[2, base_inputs])
        AlgebraicRange["minimum", algebraic_input_terms[i]] <- optim(minvect,
            optimfunc, method = "L-BFGS-B", lower = minvect,
            upper = maxvect)$value
        AlgebraicRange["maximum", algebraic_input_terms[i]] <- optim(maxvect,
            optimfunc, method = "L-BFGS-B", lower = minvect,
            upper = maxvect, control = list(fnscale = -1))$value
    }
    return(AlgebraicRange)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
