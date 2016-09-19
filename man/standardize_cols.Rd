\name{standardize_cols}
\alias{standardize_cols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardize columns of supplied matrix
}
\description{
Standardizes the columns of a supplied matrix from -1 to 1 based on a supplied data frame of ranges. Also reverses standardization based on the same inputs.
}
\usage{
standardize_cols(StartingMat, column_names, Input_range, reverse_standard = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{StartingMat}{
Data frame or matrix with named columns to be standardized.
}
  \item{column_names}{
Vector of column names to be standardized by this operation. Columns that are not included in this vector will not be modified.
}
  \item{Input_range}{
Range of basic input variables in the input data frame. Column names must match StartingMat column names. Format is a matrix or data frame with column in Input_range matching each column name in StartingMat and column_names, minimum value in first row, and maximum value in second row.
}
  \item{reverse_standard}{
TRUE/FALSE binary value indicating Whether standardize (FALSE) or reverse a standardization (TRUE) of the supplied matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

Returns a matrix with the same column names as the supplied matrix where the columns in column_names have been standardized or reversed.

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
function (StartingMat, column_names, Input_range, reverse_standard = FALSE)
{
    StandardMat <- StartingMat
    for (i in 1:length(column_names)) {
        if (reverse_standard == FALSE) {
            if (missing(Input_range)) {
                currentmax <- max(StartingMat[, column_names[i]])
                currentmin <- min(StartingMat[, column_names[i]])
                currentavg <- (currentmax + currentmin)/2
                currentrange <- currentmax - currentmin
            }
            else {
                currentmax <- max(Input_range[, column_names[i]])
                currentmin <- min(Input_range[, column_names[i]])
                currentavg <- (currentmax + currentmin)/2
                currentrange <- currentmax - currentmin
            }
            if (currentrange == 0) {
                StandardMat[, column_names[i]] <- 0
            }
            else {
                StandardMat[, column_names[i]] <- 2 * (StartingMat[,
                  column_names[i]] - currentavg)/currentrange
            }
        }
        else {
            if (missing(Input_range)) {
            }
            else {
                currentmax <- max(Input_range[, column_names[i]])
                currentmin <- min(Input_range[, column_names[i]])
                currentavg <- (currentmax + currentmin)/2
                currentrange <- currentmax - currentmin
                StandardMat[, column_names[i]] <- (currentrange/2) *
                  (StartingMat[, column_names[i]]) + currentavg
            }
        }
    }
    return(StandardMat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
