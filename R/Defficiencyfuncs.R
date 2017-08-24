########################D-efficiency calc function######################################
#Calculates determinant of supplied model matrix and scales to d-efficiency

#Output: vector(d_eff,vcov) where:
#d_eff - d-efficiency of model matrix assuming linear model
#vcov - the variance-covariance matrix of the model

#Inputs:
#CurrentMatrix - model matrix with appropriate factor numerical coding as created by the model.matrix function. Data.frame formats not accepted

d_efficiencysimple <- function(CurrentMatrix, returncov = FALSE){

  #Calculate information matrix
  infomat <- t(CurrentMatrix)%*%CurrentMatrix
  #Calculate determinant of information matrix
  det_calc <- det(infomat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Calculate d-efficiency
  d_eff <- ((det_calc)^(1/(ncol(CurrentMatrix))))/nrow(CurrentMatrix)

  if(returncov == TRUE){

    #Calculate variance-covariance matrix
    vocvmat <- tryCatch(solve(infomat), error = function(x) diag(x = Inf, nrow = 2, ncol = 2))
    output <- list(d_eff = d_eff, vcov = vocvmat)

  }else{output <- d_eff}

  #Return objective function and determinants for both current models
  return(output)}

###########################################################################################

########################D-efficiency calc function######################################
#Calculates determinant of supplied model matrix and scales to d-efficiency

#Output: d-efficiency of model matrix assuming linear model

#Inputs:
#CurrentMatrix - model matrix with appropriate factor numerical coding as created by the model.matrix function. Data.frame formats not accepted

d_efflinearupdate <- function(CurrentMatrix){

  #Calculate information matrix
  infomat <- t(CurrentMatrix)%*%CurrentMatrix

  #Calculate determinant of information matrix
  det_calc <- det(infomat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Return d-efficiency
  return(((det_calc)^(1/(ncol(CurrentMatrix))))/nrow(CurrentMatrix))}

###########################################################################################


########################Optimization Calc Function, D-Optimal######################################
#Calculates determinant ratio based on supplied det_ref. Will convert a base data.frame/matrix to a model matrix and standardize based on supplied input_formula and Input_range.

#Output: vector(d_eff,det_calc) where:
#d_eff - d-efficiency of
#det_calc - determinanat of info matrix created using CurrentMatrix and input_formula

#Inputs:
#CurrentMatrix - either a matrix or data.frame of base variables only (ex: Col for A, B, C) that will be converted to a model matrix using input_formula
#det_ref - reference determinant for calculation of d_efficiency. Should be max attainable determinant of info matrix
#input_formula - one-sided formula (inputs only) for model to calculate d efficiency from. Format = ~ A + B + A:B + I(A/B^2), etc
####required if CurrentMatrix is a data frame of base variables that needs to be converted to a model matrix
#Input_range - matrix or data frame listing the range of columns used for standardization. Column names must match input_formula term names
####format is matrix or data frame with column in Input_range matching each column name in CurrentMatrix and minimum value in first row and maximum value in second row.


d_efficiency <- function(CurrentMatrix, det_ref, input_formula, Input_range){

  #Convert input data to data.frame in case it was passed as matrix or array
  CurrentMatrix <- data.frame(CurrentMatrix)

  #Get model matrix
  modelmat <- model.matrix(input_formula,data = CurrentMatrix)

  #Standardize model matrix
  modelmat <- standardize_cols(modelmat, colnames(modelmat[,2:ncol(modelmat)]), Input_range = Input_range)

  #Calculate determinant of information matrix
  det_calc <- det(t(modelmat)%*%modelmat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Calculate ratio of current determinant to optimal determinanat for additive and mechanistic model
  d_eff <- ((det_calc/det_ref)^(1/(ncol(modelmat))))

  #Construct return vector with named elements
  returnvect <- c(d_eff, det_calc)
  names(returnvect) <- c("D efficiency", "Info Matrix Determinant")

  #Return objective function and determinants for both current models
  return(returnvect)}

###########################################################################################


########################D-efficiency calc function for discrete choice experiments######################################
#Calculates d-error of supplied model matrix. Will also calculate probability centered d-error if vector of parameter estimates is supplied

#Output: list containing:
#d_error - d-error of supplied design with respect to parameter estimates
#covmat - covariance matrix

#Inputs:
#CurrentMatrix - model matrix as created by model.matrix function with all continuous parameters centered and standardized
#### and appropriate factor numerical coding (generated using contr.sum for factors) as created by the model.matrix function.
####Data.frame formats not accepted
####Opt out choices, if included, should be coded as a row of all zeroes for every parameter
#altvect - vector with integer corresponding to each row that indicates which choice set it is a member of
#paramestimates - estimates for each effect (column) of model matrix sized corresponding to standardized model matrix.
####If not supplied, parameter estimates will be assumed equal to zero for all parameters
#returninfomat - whether to return the covariance matrix as well. Default is FALSE

d_effchoice <- function(CurrentMatrix, altvect, paramestimates = NULL, returninfomat = FALSE){

  #get all unique alternate names
  altnames <- unique(altvect)

  #check supplied parameters and set = 0 if not supplied
  if(is.null(paramestimates)){
    paramestimates <- rep(0, ncol(CurrentMatrix))
  }

  info_mat=matrix(rep(0,ncol(CurrentMatrix)*ncol(CurrentMatrix)), ncol(CurrentMatrix), ncol(CurrentMatrix))
  # compute exp(design matrix times initial parameter values)
  exputilities=exp(CurrentMatrix%*%paramestimates)

  #Initialize vector that will be sum of product of each set probability (variance)
  p_var <- 0

  #Initialize vector that holds calculated probability of each alternative being selected
  probvect <- rep(0,nrow(CurrentMatrix))

  #Loop over each choice set
  for (k_set in 1:length(altnames)) {
    # select row numbers corresponding to current loop alternatives in the choice set
    alternatives= which(altvect == altnames[k_set])
    # obtain vector of choice shares within the choice set
    p_set=exputilities[alternatives]/sum(exputilities[alternatives])
    # also put these probabilities on the diagonal of a matrix that only contains zeros
    p_diag=diag(p_set)
    # compute middle term P-pp'
    middle_term<-p_diag-p_set%o%p_set
    # pre- and postmultiply with the Xs from the design matrix for the alternatives in this choice set
    full_term<-t(CurrentMatrix[alternatives,])%*%middle_term%*%CurrentMatrix[alternatives,]
    # Add contribution of this choice set to the information matrix
    info_mat<-info_mat+full_term

    #Calculate product of all probabilities in set and add to p_var
    p_var <- p_var + prod(p_set)

    #Enter all calculated probabilities for this set into the probability vector
    probvect[alternatives] <- p_set

  }
  #get the inverse of the information matrix (i.e., gets the variance-covariance matrix)
  #Use "try" wrapper to prevent unsolvable matrices from crashing. Return 2x2 diagonal infinite matrix on crash
  #sigma_beta<- tryCatch(solve(info_mat,diag(ncol(CurrentMatrix))), error = function(x) diag(x = Inf, nrow = ncol(CurrentMatrix), ncol = ncol(CurrentMatrix)))


  #Calculate determinant
  det_calc <- det(info_mat)
  #If determinant is negative (should not be possible but sometimes happens), return zero to prevent an error
  if(det_calc < 0){det_calc <- 0}

  if(returninfomat == TRUE){

    output <- list(d_eff = det_calc^(1/ncol(CurrentMatrix)), info_mat = info_mat, p_var = p_var, probvect = probvect)

  }else{output <- det_calc^(1/ncol(CurrentMatrix))}

  #Return objective function and determinants for both current models
  #return(list(d_eff = det(sigma_beta)^(-1/ncol(CurrentMatrix)), vcov = sigma_beta))}
  return(output)}


###########################################################################################


########################D-efficiency update function for single question for multinomial logit model#######################################Updates d-efficiency of supplied sub-model matrix and information matrix.
#Fast d-efficiency update function for use during search algorithms to reduce search time. Requires supplied information matrix for all questions not included
#Output: numeric value for d-efficiency

#Inputs:
#CurrentMatrix - subsection of model matrix for one question only. Should include all alternatives for that questionas created by the model.matrix function with all continuous parameters centered and standardized
#### and appropriate factor numerical coding (generated using contr.sum for factors) as created by the model.matrix function.
####Data.frame formats not accepted
####Opt out choices, if included, should be coded as a row of all zeroes for every parameter
#paramestimates - standardized effect estimates for each effect (column) of model matrix corresponding to standardized model matrix.
#info_mat - the information matrix for all other questions in the model matrix

d_effchoiceupdate <- function(CurrentMatrix, paramestimates, info_mat){

  # compute exp(design matrix times initial parameter values)
  exputilities <- c(exp(CurrentMatrix%*%paramestimates))

  # obtain vector of choice shares within the choice set
  p_set <- exputilities/sum(exputilities)

  #Calculate product of all probabilities in set for utility balance use
  p_var <- prod(p_set)

  # calculate information matrix of this choice set and add it to the info_matrix
  info_mat<-info_mat + t(CurrentMatrix)%*%(diag(p_set)-p_set%o%p_set)%*%CurrentMatrix

  #Calculate determinant
  det_calc <- det(info_mat)

  #If determinant is negative (should not be possible but sometimes happens due to numerical processing), return zero to prevent an error
  if(det_calc < 0){det_calc <- 0}

  #Return Determinant
  return(list(d_eff = det_calc^(1/ncol(CurrentMatrix)), p_var = p_var))}


###########################################################################################
