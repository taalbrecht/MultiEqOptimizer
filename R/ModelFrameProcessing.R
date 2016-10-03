############################Convert to sum contrasts while preserving factor order####################

#Converts a contrast to contrast sum while preserving factor order and applying appropriate column names

#Inputs
#factorin - vector of factor values
#speclevels - vector of factor level names in desired order

#Output
#vector as a factor with contrast sum encoding that respects the order in speclevels

customcontrsum <- function(factorin, speclevels){

      #Correct factor levels to proper order
  factorin <- factor(factorin, levels = speclevels)

      #Apply contrast sum encoding
      contrasts(factorin) <- contr.sum(nlevels(factorin))

      #Apply contrast column names
      colnames(contrasts(factorin)) <- rownames(contrasts(factorin))[1:ncol(contrasts(factorin))]

      return(factorin)

}


#################################Standardize Columns in Matrix Function##########################

#Calculates: Standardizes (from {-1 to 1}) or reverses standardization of list of column names in matrix

#Output: matrix(StandardMat), where StandardMat matches StartingMat but each column name in column_names has been standardized from -1 to 1 based on min and max values in StartingMat or Input_range if provided

#Inputs, required:
#StartingMat - a raw input matrix of continuous variables with column names matching variable names in input_formula. May also include outputs or variables to not be analyzed
#column_names - a list that contains all of the column names to be standardized. format = list("Response", "X1", ...)

#Inputs, optional:
#Input_range - matrix or data frame listing the range of columns used for standardization. Names must match those used in column_names
#			format is matrix or data frame with column name for each column to be standardized and minimum value in first row and maximum value in second row.
#reverse_standard - logical(TRUE, FALSE) indicating whether standardization should be reversed using Input_range. Default is FALSE

standardize_cols <- function(StartingMat, column_names, Input_range, reverse_standard = FALSE){

  #Initialize StandardMat by setting equal to StartingMat
  StandardMat <- StartingMat

  #For loop to standardize each variable column used in the formula from -1 to 1
  for(i in 1:length(column_names)){

    if(reverse_standard == FALSE){
      #If Input_range was not provided, pull min and max values from StartingMat. Otherwise, use min and max from Input_range
      if(missing(Input_range)){
        currentmax <- max(StartingMat[,column_names[i]])
        currentmin <- min(StartingMat[,column_names[i]])
        currentavg <- (currentmax + currentmin)/2
        currentrange <- currentmax-currentmin
      }
      else{currentmax <- max(Input_range[,column_names[i]])
      currentmin <- min(Input_range[,column_names[i]])
      currentavg <- (currentmax + currentmin)/2
      currentrange <- currentmax-currentmin
      }
      if(currentrange == 0){
        StandardMat[,column_names[i]] <- 0
      }
      else{
        StandardMat[,column_names[i]] <- 2*(StartingMat[,column_names[i]]-currentavg)/currentrange
      }
    }
    else{
      #Reverse standardization. If Input_range was not provided, return same matrix as StartingMat
      if(missing(Input_range)){}
      else{currentmax <- max(Input_range[,column_names[i]])
      currentmin <- min(Input_range[,column_names[i]])
      currentavg <- (currentmax + currentmin)/2
      currentrange <- currentmax-currentmin

      StandardMat[,column_names[i]] <- (currentrange/2)*(StartingMat[,column_names[i]])+currentavg
      }
    }
  }
  #Return Standardized Matrix
  return(StandardMat)}

##################################################################################################################

###############################Pull input variables from formula function#######################################

#Calculates: list of all unique fundamental input variables from formula passed to function

#Output: list(InputVarList), where:
#InputVarList - list of all input variables in function, format = list("x", "y", "z") where input formula may include combinations of x, y, and z

#Inputs:
#input_formula - a formula (one or two sided) that contains all of the input variables to be modeled. format = x~ A + B + A:B + I(A^2/B), etc

list_input_variables <- function(input_formula){

  #library(nlme)


  #Pull inputs only from input_formula and count
  Input_Func <- as.formula(paste(nlme::splitFormula(input_formula)))
  InputVarList <- all.vars(Input_Func)

  # #Pull inputs only from input_formula and count - work on this
  # InputVarList <- attributes(terms(input_formula))$term.labels

  #Return list of input variables
  return(InputVarList)
}

#########################################################################################################################

#################################Function to return min and max range for algebraic combinations of base variables###################

#Calculates: Max and min for algebraic combination of base variables when range of base variables is known

#Output: matrix(AlgebraicRange), where Algebraic_Range is a matrix with column names matching algebraic combos on input side of formula. Row 1 = min and row 2 = max.

#Inputs, required:
#base_var_range - matrix or data frame listing the range of base input variables used for algebraic model. Names must match those used in algebraic_formula
#			format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.
#algebraic_formula - a formula (one or two sided) that contains all of the algebraic combos to be evaluated. format = x~ A + B + A:B + I(A^2/B)+ I(ln(A)), etc

algebraic_range <- function(base_var_range, algebraic_formula){

#   #Old code that does not work with factors and continuous variables
#   #Pull names of each algebraic term on input side of formula
#   algebraic_input_terms <- colnames(attr(terms(algebraic_formula),"factors"))

  #Pull names of each algebraic term and factor combination on input side of formula
  algebraic_input_terms <- colnames(model.matrix(algebraic_formula, base_var_range), "factors")[-1]

  #Create matrix to store min and max values for each algebraic term
  AlgebraicRange <- matrix(data = 0, nrow = 2, ncol = length(algebraic_input_terms))
  colnames(AlgebraicRange) <- algebraic_input_terms
  rownames(AlgebraicRange) <- c("minimum", "maximum")

  #Loop to be executed for each algebraic term found

  for (i in 1:length(algebraic_input_terms)){

    #Replace : multiplication operators with * operators to ensure expression will work correctly
    transmult <- gsub(":","*",algebraic_input_terms[i])

    #Translate algebraic term to expression
    optimexp <- parse(text = transmult)

    #Pull base variables from algebraic term
    base_inputs <- all.vars(optimexp)

    #Create function using current algebraic term for use with optimization routine
    optimfunc <- function(optimvect){

      basevarmat <- matrix(nrow = 1, ncol = length(base_inputs), data = 0)
      basevarmat <- data.frame(basevarmat)
      colnames(basevarmat) <- base_inputs
      basevarmat[1,] <- optimvect
      output <- eval(optimexp,basevarmat)
      return(output)
    }

    #Create expanded grid using all possible combinations of base input variables to use as starting guesses

    trymat <- expand.grid(base_var_range[base_inputs])

    #Define vector of minimum and maximum values for base variables used in this algebraic term
    minvect <- as.vector(base_var_range[1,base_inputs])
    maxvect <- as.vector(base_var_range[2,base_inputs])

    #Find minimum and maximum values by trying all starting combinations in expanded grid
    AlgebraicRange["minimum",algebraic_input_terms[i]] <- min(unlist(lapply(apply(trymat, MARGIN = 1, optim, fn = optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect), "[", "value")))
    AlgebraicRange["maximum",algebraic_input_terms[i]] <- max(unlist(lapply(apply(trymat, MARGIN = 1, optim, fn = optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect, control = list(fnscale = -1)), "[", "value")))

#     #Old, less accurate code that doesn't find min and max for factor combinations
#     #Find minimum value for algebraic term
#     AlgebraicRange["minimum",algebraic_input_terms[i]] <- optim(minvect, optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect)$value
#     AlgebraicRange["maximum",algebraic_input_terms[i]] <- optim(maxvect, optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect, control = list(fnscale = -1))$value
  }

  #Return Min and Max for Algebraic Functions
  return(AlgebraicRange)}

#################################################################################################################################
