\name{MultipleModelOptimize}
\alias{MultipleModelOptimize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimize a single experimental design for multiple discrete choice models and/or linear models
}
\description{
Optimizes a single experimental design for multiple discrete choice and/or linear models by joint optimization of the conditional logistic d-efficiency for all models. The algorithm attemtps to maximize the sum of the d-efficiency for each formula multiplied by the weight for each formula. If no weight vector is supplied, the algorithm will attempt to maximize lowest d-efficiency.
}
\usage{
MultipleModelOptimize(base_input_range, formulalist, questions, alts, blocks = NA, optout = FALSE, det_ref_list, mesh, tolerance, weight, candset = NA, priors = NA, searchstyle = "Fedorov")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base_input_range}{
Range of all base input variables used for all models in formulalist. Names must match those used in formulalist. Format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row. May also include ranges for algebraic combinations used in each model in formulalist. If algebraic combination range is not included, this function will attempt to estimate the maximum and minimum range for the algebraic combination based on the base input variable ranges provided.
}
  \item{formulalist}{
List of each model for joint d-efficiency optimization. Dependent variable does not need to be included.
Format = list(I(x/A)~ I(A/B^2), ~ A + B + A:B) etc
}
  \item{questions}{
Integer specifying the number of choice sets (questions) to be used in the experiment.
}
  \item{alts}{
Number of alternatives per choice set (question) not including an opt-out option.
}
  \item{blocks}{
Optional. An integer specifying the number of blocks to be used in the design. Will default to 1 block if not supplied.
}
  \item{optout}{
TRUE/FALSE - whether to include an opt-out in each question set (e.g. "Neither"). This will be added to number of alternates per question. Defaults to FALSE.
}
  \item{det_ref_list}{
List of reference maximum information matrix determinants for each formula in formulalist. Order of this list matches order of input formulas. Format = list(200, 216, ...)
}
  \item{mesh}{
This is required for Columnwise sampling and represents the number of steps each input factor is broken into for the optimization algorithm. Can be supplied in the following formats:
1. As a single value where it will be applied to all continuous variables which is interpreted as the number of steps to break each continuous variable into.
2. As a list with a single mesh size for each variable. It is then is interpreted as the number of steps to break each continuous variable into.
3. As a list of vectors with one vector for each variable. In this format the vector is interpreted as a list of all valid values for each variable.
}
  \item{tolerance}{
Numeric value indicating the minimum change in optimality criterion needed to terminate optimization function.
}
  \item{weight}{
Vector of weighting values to apply to the d-efficiency calculation for each formula in formulalist.
}
  \item{candset}{
Data frame or matrix of candidate points to search through for model creation. All basic variables in the formulalist must be contained in this matrix. This is required for Fedorov sampling but is not used for Columnwise sampling.
}
  \item{priors}{
List of known or estimated effect sizes (in model terms). Structure should be a list as follows: list(I(A^2) = 0.23, B = -0.5) where parameter names in the list match those in formulas. For attribute parameters, the list should be the effect size with reference to first factor level. Any effects not provided will be assumed to be equal to zero.
}
  \item{searchstyle}{
Character term defining which sampling style to use for optimization. The options are:

Fedorov:
Samples from a supplied matrix of available model points defined by candset.

Columnwise:
Samples across the range of each base variable using step sizes defined by mesh input.
}

  \item{startingdesign}{
Optional. Design to use as a starting point for optimization. Format is a data frame with one column for each base input variable and the proper number of rows corresponding to the number of model points. If not provided, a random starting design will be created.
}

  \item{eqtype}{
Optional. Vector of model types to assume for d-efficiency calculations. Defaults to choice model but can be set to linear calculation as well. Accepted values are "Choice" and "Linear". Vector length must equal the length of the list of formulas.
}

  \item{questionblockvars}{
Optional. Vector of named variables that should be the same for all alternatives within one question set. Currently only functional for Columnwise search strategy.
}

  \item{augment}{
Optional. TRUE/FALSE value that specifies whether the model should optimize all points in the design or only new points. To use this option a compatible design must be furnished to the startingdesign argument. In this case, the questions and alts arguments will be assumed to apply to the number of design points that should be ADDED to the design. Default is FALSE.
}

  \item{priorsnormalized}{
TRUE/FALSE value that specifies whether the supplied prior estimates for the coefficient values are scaled for a model effect range of -1 to 1. If FALSE, the priors supplied are with respect to the real-world values of the variable (e.g. a prior of 0.5 for a variable that can take values from 1 to 5 would provide an effect of 0.5 - 2.5). If TRUE, the priors are scaled for a model effect range of -1 to 1 (e.g. a prior of 0.5 for a variable will provide an effect of -0.5 to 0.5 for any supplied varlable range). Default is FALSE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list:

  \item{ModelMatrix}{
Matrix containing the single experiment optimized for all input formulas
}

  \item{Deff}{
List of d-efficiency values for each formula
}

  \item{DeffvsOptimal}{
List of d-efficiency values for each formula divided by the optimal d-efficiency for each formula (supplied as det_ref_list in the input of this function)
}

  \item{CovList}{
List of the covariance matrix for each formula
}

  \item{ObjectiveFunction}{
Final value of the objective function
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

##Define data needed to create design

#Create list of formulas
formulalist = list(~A+B, ~I(A/B)+B)

#Create a set of candidate points (for Fedorov optimization)
candset <- expand.grid(A = c(1, 2, 3), B = c(4, 5, 6))

#Define step sizes for each variable (for Columnwise optimization)
mesh <- list(20, 20)

#Define ranges of each basic variable in all of the formulas in formulalist
baseinputranges <- data.frame(A = c(1,3), B = c(4, 6))

#Define weights
weightvect <- c(0.5, 0.5)

#Define equation types
eqtypes <- c("Choice", "Linear")

#Define estimates of coefficient values for each formula (affects discrete choice optimization but not linear optimization)
priors <- list(c(-1, 1.5), c(0.75, 0.5))

#Define reference determinant values to use for scaling d-efficiency for each formula. Used for weighting purposes
deffrefs <- c(1,1)

##Create design using Fedorov optimization
fedorovdesign <- MultipleModelOptimize(base_input_range = baseinputranges,
formulalist = formulalist,
questions = 10,
alts = 2,
det_ref_list = deffrefs,
tolerance = 0.01,
candset = candset,
priors = priors,
weight = weightvect,
searchstyle = "Fedorov",
eqtype = eqtypes)

##Create design using Columnwise optimization
columnwisedesign <- MultipleModelOptimize(base_input_range = baseinputranges,
formulalist = formulalist,
questions = 10,
alts = 2,
det_ref_list = deffrefs,
tolerance = 0.01,
mesh = mesh,
priors = priors,
weight = weightvect,
searchstyle = "Columnwise",
eqtype = eqtypes)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
