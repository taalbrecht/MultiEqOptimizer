\name{DiscChoiceMultipleModel}
\alias{DiscChoiceMultipleModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimize a single experimental design for multiple discrete choice models
}
\description{
Optimizes a single experimental design for multiple discrete choice models by joint optimization of the conditional logistic d-efficiency for all models. The algorithm attemtps to maximize the sum of the d-efficiency for each formula multiplied by the weight for each formula. If no weight vector is supplied, the algorithm will attempt to maximize lowest d-efficiency.
}
\usage{
DiscChoiceMultipleModel(base_input_range, formulalist, questions, alts, blocks = NA, optout = FALSE, det_ref_list, mesh, tolerance, weight, candset = NA, priors = NA, searchstyle = "Federov")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base_input_range}{
Range of all base input variables used for all models in formulalist. Names must match those used in formulalist. Format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row. May also include ranges for algebraic combinations used in each model in formulalist. If algebraic combination range is not included, this function will attempt to estimate the maximum and minimum range for the algebraic combination based on the base input variable ranges provided.
}
  \item{formulalist}{
List of each model for joint d-efficiency optimization. Dependent variable does not need to be included.
Format = list(I(x/A)~ I(A/B^2), ~ A + B + A:B etc)
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
This is required for Gibbs sampling and represents the number of steps each input factor is broken into for the optimization algorithm. Can be supplied in the following formats:
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
Data frame or matrix of candidate points to search through for model creation. All basic variables in the formulalist must be contained in this matrix. This is required for Federov sampling but is not used for Gibbs sampling.
}
  \item{priors}{
List of known or estimated effect sizes (in model terms). Structure should be a list as follows: list(I(A^2) = 0.23, B = -0.5) where parameter names in the list match those in formulas. For attribute parameters, the list should be the effect size with reference to first factor level. Any effects not provided will be assumed to be equal to zero.
}
  \item{searchstyle}{
Character term defining which sampling style to use for optimization. The options are:

Federov:
Samples from a supplied matrix of available model points defined by candset.

Gibbs:
Samples across the range of each base variable using step sizes defined by mesh input.
}

  \item{startingdesign}{
Optional. Design to use as a starting point for optimization. Format is a data frame with one column for each base input variable and the proper number of rows corresponding to the number of model points. If not provided, a random starting design will be created.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (base_input_range, formulalist, questions, alts, blocks = NA,
    optout = FALSE, det_ref_list, mesh, tolerance, weight, candset = NA,
    priors = NA, searchstyle = "Federov")
{
    model_points <- questions * alts
    library(nlme)
    altvect <- c()
    for (i in 1:questions) {
        altvect <- c(altvect, rep(i, alts))
    }
    inputs_only <- lapply(formulalist, list_input_variables)
    input_list <- unique(unlist(inputs_only))
    list_input_ranges <- lapply(formulalist, algebraic_range,
        base_var_range = base_input_range)
    all_input_ranges <- data.frame(base_input_range, list_input_ranges)
    colnames(all_input_ranges) <- c(colnames(base_input_range),
        unlist(lapply(list_input_ranges, colnames)))
    all_input_ranges <- as.matrix(all_input_ranges)
    all_input_ranges <- all_input_ranges[, unique(colnames(all_input_ranges))]
    input_formulas <- sapply(formulalist, nlme::splitFormula)
    if (searchstyle == "Gibbs") {
        if (length(mesh) == 1) {
            mesh <- rep(mesh, times = length(input_list))
        }
        stepseq <- list()
        for (i in 1:length(input_list)) {
            step <- 2/(mesh[i] - 1)
            stepseq[[i]] <- seq(-1, 1, by = step)
        }
        ModelMatStand <- sapply(stepseq, sample, size = model_points,
            replace = TRUE)
        colnames(ModelMatStand) <- input_list
        d_efficiency_vect <- Vectorize(d_efficiency, c("input_formula",
            "det_ref"))
        ModelMatReal <- standardize_cols(ModelMatStand, input_list,
            all_input_ranges, reverse_standard = TRUE)
        eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal,
            det_ref = det_ref_list, input_formula = input_formulas,
            Input_range = all_input_ranges)
        if (missing(weight)) {
            obj_current <- min(eff_vect[1, ])
        }
        else {
            obj_current <- sum(eff_vect[1, ] * weight)
        }
        objective_change <- tolerance + 1
        while (objective_change > tolerance) {
            obj_prev <- obj_current
            for (i in 1:nrow(ModelMatStand)) {
                for (j in 1:ncol(ModelMatStand)) {
                  for (s in stepseq[[j]]) {
                    oldvalue <- ModelMatStand[i, j]
                    ModelMatStand[i, j] <- s
                    ModelMatReal <- standardize_cols(ModelMatStand,
                      input_list, all_input_ranges, reverse_standard = TRUE)
                    eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal,
                      det_ref = det_ref_list, input_formula = input_formulas,
                      Input_range = all_input_ranges)
                    if (missing(weight)) {
                      obj_temp <- min(eff_vect[1, ])
                    }
                    else {
                      obj_temp <- sum(eff_vect[1, ] * weight)
                    }
                    if (obj_temp <= obj_current) {
                      ModelMatStand[i, j] <- oldvalue
                    }
                    else {
                      obj_current <- obj_temp
                    }
                  }
                }
            }
            objective_change <- obj_current - obj_prev
        }
        ModelMatReal <- standardize_cols(ModelMatStand, input_list,
            all_input_ranges, reverse_standard = TRUE)
        eff_vect_final <- d_efficiency_vect(CurrentMatrix = ModelMatReal,
            det_ref = det_ref_list, input_formula = input_formulas,
            Input_range = all_input_ranges)
    }
    if (searchstyle == "Federov") {
        candset <- candset[, input_list]
        factcols <- rep(FALSE, ncol(candset))
        for (i in 1:ncol(candset)) {
            factcols[i] <- is.factor(candset[, i])
        }
        for (i in 1:length(factcols)) {
            if (factcols[i] == TRUE) {
                contrasts(candset[, i]) <- contr.sum(nlevels(candset[,
                  i]))
            }
        }
        candexpand <- lapply(input_formulas, model.matrix, data = candset)
        if (sum(factcols == FALSE) > 0) {
            candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x,
                column_names = colnames(x)[colnames(x) \%in\% colnames(inrange)],
                Input_range = inrange))
        }
        rownums <- sample(nrow(candset), size = model_points,
            replace = TRUE)
        d_efficiency_vect <- Vectorize(d_effchoice, c("CurrentMatrix",
            "paramestimates"))
        eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand,
            function(x, rowind = rownums) x[rowind, 2:ncol(x)]),
            paramestimates = priors, altvect = altvect)
        if (missing(weight)) {
            obj_current <- min(eff_vect)
        }
        else {
            obj_current <- sum(eff_vect * weight)
        }
        objective_change <- tolerance + 1
        while (objective_change > tolerance) {
            obj_prev <- obj_current
            for (i in 1:length(rownums)) {
                for (j in 1:nrow(candset)) {
                  oldvalue <- rownums[i]
                  rownums[i] <- j
                  eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand,
                    function(x, rowind = rownums) x[rowind, 2:ncol(x)]),
                    paramestimates = priors, altvect = altvect)
                  if (missing(weight)) {
                    obj_temp <- min(eff_vect)
                  }
                  else {
                    obj_temp <- sum(eff_vect * weight)
                  }
                  if (obj_temp <= obj_current) {
                    rownums[i] <- oldvalue
                  }
                  else {
                    obj_current <- obj_temp
                  }
                }
            }
            objective_change <- obj_current - obj_prev
        }
        ModelMatReal <- candset[rownums, ]
        eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand,
            function(x, rowind = rownums) x[rowind, 2:ncol(x)]),
            paramestimates = priors, altvect = altvect, returncov = TRUE)
    }
    outputfinal <- list(data.frame(Question = altvect, ModelMatReal),
        eff_vect_final, obj_current)
    names(outputfinal) <- c("ModelMatrix", "D-Efficiency", "ObjectiveFunction")
    return(outputfinal)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
