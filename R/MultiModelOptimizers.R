#########################CEX Function to optimize for all models in list of formulas###############################

#Calculates: Design that optimizes for multiple models with constraint to efficiency target for each model

#Output: list(ModelMat,det_mech,det_add) where:
#ModelMat, where ModelMat has been optimized for mechanistic model
#obj - objective function
#det_mech - determinant of info matrix of current mechanistic model
#det_add - determinanat of info matrix of current additive model

#Inputs:
#base_input_range - range of all base input variables used for all models in formulalist. Names must match those used in formulalist.
#			format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.
#			May also include ranges for algebraic combinations used in each model in formulalist. If algebraic combination range is not included, this function will
#			attemtp to estimate the maximum and minimum range for the algebraic combination based on the base input variable ranges provided.
#formulalist - list of each model to assess design against, output does not need to be included format = list(I(x/A)~ I(A/B^2), ~ A + B + A:B etc)
#model_points - integer specifying the number of unique model points to use
#blocks - number of experimental blocks (survey versions). Not currently implemented
#mesh - number of steps each factor broken into for optimization algorithm. Can be supplied as a single value where it will be applied to all variables or as a vector with a mesh size for each variable
#tolerance - minimum change in optimality criterion needed to terminate optimization function. Must be greater than 0
#det_ref_list - list of reference max determinants for each formula in formulalist. Order of this list matches order of input formulas. Format = list(200, 216, etc...)
#searchstyle - option for searching experimental space. "Federov" samples from a supplied matrix of available designs (candset). "Gibbs" samples across the range of each base variable by step sizes defined by mesh.
#candset - set of candidate points to search through for model creation. This is required if "Federov" searchstyle is used but is not used otherwise

#Inputs, optional:
#weight - vector of weighting values to apply to each formula in formulalist. Order of this vector matches order of input formulas. Format = c(0.5,0.25, etc)
#		if weight is not supplied, function will attempt to maximize lowest d_efficiency using det_ref_list

CEX_MultipleModel <- function(base_input_range, formulalist, model_points, blocks = NA, det_ref_list, mesh, tolerance, weight, candset = NA, searchstyle = "Federov"){

  #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
  # step <- 2/(mesh-1)
  # stepseq <- seq(-1,1, by = step)

  #library(nlme)

  #Pull inputs only from all formulas and place into a list
  inputs_only <- lapply(formulalist, list_input_variables)

  #Create list of all base input variables used in formulalist
  input_list <- unique(unlist(inputs_only))

  #Find min and max range for all base vars and algebraic combinations in formulalist (WARNING: If surface is very complex, THIS MAY FAIL)
  list_input_ranges <- lapply(formulalist, algebraic_range, base_var_range = base_input_range)
  all_input_ranges <- data.frame(base_input_range, list_input_ranges)
  colnames(all_input_ranges) <- c(colnames(base_input_range), unlist(lapply(list_input_ranges,colnames)))
  all_input_ranges <- as.matrix(all_input_ranges) #Convert from list to matrix
  all_input_ranges <- all_input_ranges[,unique(colnames(all_input_ranges))] # Remove redundant columns from matrix

  #Create list of one-sided expressions of input only side for all formulas in formulalist
  input_formulas <- sapply(formulalist, nlme::splitFormula)

  #Optimize via Gibbs search if that option is selected

  if(searchstyle == "Gibbs"){
    #If mesh supplied as a single value, apply it to all base input values
    if(length(mesh) == 1){

      mesh <- rep(mesh, times = length(input_list))

    }


    #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
    stepseq <- list()

    for(i in 1:length (input_list)){

      step <- 2/(mesh[i]-1)

      stepseq[[i]] <- seq(-1,1, by = step)

    }

    # #Create random model matrix using base vars (standardized {-1 to 1}) to start optimization function and assign input names to columns
    # ModelMatStand <- matrix(data = sample(stepseq, length(input_list)*model_points, replace = TRUE), nrow = model_points, ncol = length(input_list))
    # colnames(ModelMatStand) <- input_list

    #Create random model matrix using base vars (standardized {-1 to 1}) to start optimization function and assign input names to columns
    ModelMatStand <- sapply(stepseq, sample, size = model_points, replace = TRUE)
    colnames(ModelMatStand) <- input_list

    #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
    d_efficiency_vect <- Vectorize(d_efficiency, c("input_formula", "det_ref"))

    #Calculate initial D-efficiencies for randomly seeded model
    ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE) #Convert to real numbers before passing to d_efficiency
    eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)

    #Calculate starting objective function using vector of weights
    if(missing(weight))
    {obj_current <- min(eff_vect[1,])}
    else{obj_current <- sum(eff_vect[1,]*weight)}

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one parameter at a time
      for (i in 1:nrow(ModelMatStand)){
        for(j in 1:ncol(ModelMatStand)){

          #Step through -1 to 1 by step value
          for (s in stepseq[[j]]){

            #Replace current value in ModelMatStand with new temporary point
            oldvalue <- ModelMatStand[i,j]
            ModelMatStand[i,j] <- s

            #Calculate D-efficiencies for ModelMatStand with new temporary point
            ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE) #Convert to real numbers before passing to d_efficiency
            eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)

            #Calculate new objective function value using vector of weights and temporary point
            if(missing(weight))
            {obj_temp <- min(eff_vect[1,])}
            else{obj_temp <- sum(eff_vect[1,]*weight)}

            #If objective of new point is greater than old point, use new point in matrix. Otherwise, put old point back into matrix
            if(obj_temp <= obj_current)

            {ModelMatStand[i,j] <- oldvalue}

            else {obj_current <- obj_temp}
          }
        }
      }
      #Check objective function after replacement to see if improv
      objective_change <- obj_current - obj_prev
    }

    #Convert final matrix to real values
    ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE)

    #Calculate final d-efficiency
    eff_vect_final <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)

  }

  if(searchstyle == "Federov"){

    #Reduce candidate set to only include base variables in the supplied formulas to eliminate extraneous data
    candset <- candset[,input_list]

    #Create model matrix of candidate set for each input formula
    candexpand <- lapply(input_formulas, model.matrix, data = candset)

    #Standardize all candidate sets within each variable range
    candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x[,2:ncol(x)]), Input_range = inrange))

    #Convert determinant list to vector
    det_refs <- unlist(det_ref_list)

    #Select random row indices from supplied candidate set of model points. Reduce matrix
    rownums <- sample(nrow(candset), size = model_points, replace = TRUE)
    #ModelMatReal <- candset[sample(nrow(candset),size = model_points, replace = TRUE),]

    #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
    d_efficiency_vect <- Vectorize(d_efficiencysimple, c("CurrentMatrix", "det_ref"))

    #Calculate initial D-efficiencies for randomly seeded model
    eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)

    #Calculate starting objective function using vector of weights
    if(missing(weight))
    {obj_current <- min(eff_vect[1,])}
    else{obj_current <- sum(eff_vect[1,]*weight)}

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one row at a time
      for (i in 1:length(rownums)){

        #Try every possible point in candidate set of points for each row
        for(j in 1:nrow(candset)){

          #Replace current value in ModelMatReal with new temporary point
          oldvalue <- rownums[i]
          rownums[i] <- j

          #Calculate D-efficiencies for ModelMatReal with new temporary point
          eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)

          #Calculate new objective function value using vector of weights and temporary point
          if(missing(weight))
          {obj_temp <- min(eff_vect[1,])}
          else{obj_temp <- sum(eff_vect[1,]*weight)}

          #If objective of new point is greater than old point, use new point in matrix. Otherwise, put old point back into matrix
          if(obj_temp <= obj_current)

          {rownums[i] <- oldvalue}

          else {obj_current <- obj_temp}
        }
      }
      #Check objective function after replacement to see if improv
      objective_change <- obj_current - obj_prev
    }

    #Convert final matrix to real values
    ModelMatReal <- candset[rownums,]

    #Calculate final d-efficiency
    eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)

  }

  colnames(eff_vect_final) <- input_formulas #Name efficiency vector formulas

  #Name final output
  outputfinal <- list(ModelMatReal, eff_vect_final, obj_current)
  names(outputfinal) <- c("ModelMatrix","D-Efficiency","ObjectiveFunction")

  #Return model matrix, objective function value, and 1/D-optimality for add, mech model

  return(outputfinal)
}

##########################################################################################################################################

#########################CEX Function to optimize for all models in list of formulas###############################

#Calculates: Design that optimizes for multiple models with constraint to efficiency target for each model

#Output: list(ModelMat,det_mech,det_add) where:
#ModelMat, where ModelMat has been optimized for mechanistic model
#obj - objective function
#det_mech - determinant of info matrix of current mechanistic model
#det_add - determinanat of info matrix of current additive model

#Inputs:
#base_input_range - range of all base input variables used for all models in formulalist. Names must match those used in formulalist.
#			format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.
#			May also include ranges for algebraic combinations used in each model in formulalist. If algebraic combination range is not included, this function will
#			attemtp to estimate the maximum and minimum range for the algebraic combination based on the base input variable ranges provided.
#formulalist - list of each model to assess design against, output does not need to be included format = list(I(x/A)~ I(A/B^2), ~ A + B + A:B etc)
#questions - integer specifying the number of survey/study questions to use
#alts - number of alternates presented per choice question set not including opt out (e.g. "Neither"). Must be integer  >= 2
#blocks - number of experimental blocks (survey versions). Not currently implemented
#optout - whether to include an opt-out in each question set (e.g. "Neither"). Will be added to number of alternates per question.
#mesh - number of steps each factor broken into for optimization algorithm. Can be supplied as a single value where it will be applied to all variables or as a list with a mesh size for each numeric variable and a list of valid factor levels for each categorical variable
#tolerance - minimum change in optimality criterion needed to terminate optimization function. Must be greater than 0
#det_ref_list - list of reference max determinants for each formula in formulalist. Order of this list matches order of input formulas. Format = list(200, 216, etc...)
#searchstyle - option for searching experimental space. "Federov" samples from a supplied matrix of available designs (candset). "Gibbs" samples across the range of each base variable by step sizes defined by mesh.
#candset - set of candidate points to search through for model creation. This is required if "Federov" searchstyle is used but is not used otherwise

#Inputs, optional:
#weight - vector of weighting values to apply to each formula in formulalist. Order of this vector matches order of input formulas. Format = c(0.5,0.25, etc)
#		if weight is not supplied, function will attempt to maximize lowest d_efficiency using det_ref_list
#priors - list of known or estimated effect sizes (in model terms). Structure should be a list as follows: list(I(A^2) = 0.23, B = -0.5) where parameter names in list match those in formulas. For attribute parameters, list should be effect size with reference to first factor level. Any effects not provided will be assumed to be equal to zero.
#startingdesign - design used as startingpoint for optimization. If not provided, a random design will be selected as the starting point. Will cause an error if this is not properly formatted.
##best practice is to use a design created by this function as a starting desgin for another round of optimization

DiscChoiceMultipleModel <- function(base_input_range, formulalist, questions, alts, blocks = NA, optout = FALSE, det_ref_list, mesh, tolerance, weight, candset = NA, priors = NA, searchstyle = "Federov", startingdesign = NULL){

  #Calculate number of model points to use
  model_points <- questions*alts

  #library(nlme)

  #Generate vector that describes sequence of question vs model points
  altvect <- c()

  for(i in 1:questions){
    altvect <- c(altvect, rep(i, alts))
  }

  #Convert determinant list to vector
  det_ref_list <- unlist(det_ref_list)

  #Pull inputs only from all formulas and place into a list
  inputs_only <- lapply(formulalist, list_input_variables)

  #Create list of all base input variables used in formulalist
  input_list <- unique(unlist(inputs_only))

  #Create list of one-sided expressions of input only side for all formulas in formulalist
  input_formulas <- sapply(formulalist, nlme::splitFormula)

  #Optimize via Gibbs search if that option is selected
  if(searchstyle == "Gibbs"){

    #If mesh supplied as a single value, apply it to all base input values
    if(length(mesh) == 1){

      mesh <- rep(mesh, times = length(input_list))

    }


    #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
    stepseq <- list()

    for(i in 1:length (input_list)){

      #For numeric values, create steps based on mesh value from minimum supplied value to maximum supplied value
      if(is.numeric(base_input_range[[input_list[[i]]]]) == TRUE){
        step <- (max(base_input_range[[input_list[[i]]]]) - min(base_input_range[[input_list[[i]]]]))/(mesh[[i]]-1)

        stepseq[[i]] <- seq(min(base_input_range[[input_list[[i]]]]),max(base_input_range[[input_list[[i]]]]), by = step)

      }

      #For attribute values already classified as factors, either use supplied factor levels or use all levels present in the factor in base_input_range
      if((is.factor(base_input_range[[input_list[[i]]]]) == TRUE) | (is.character(base_input_range[[input_list[[i]]]]) == TRUE)){

        #If mesh values were provided for this factor, use them. Otherwise return all levels of base_input_range as possibilities.
        if((is.factor(mesh[[i]]) == TRUE) | (is.character(mesh[[i]]) == TRUE)){

          stepseq[[i]] <- as.factor(mesh[[i]])

        }else{

          stepseq[[i]] <- levels(as.factor(base_input_range[[input_list[[i]]]]))

        }

      }

    }

    #Apply contrast sum encoding to base_input_range
    for(i in 1:ncol(base_input_range)){

      if(is.numeric(base_input_range[[i]]) == FALSE){

        #Check factor levels in base_input_range against supplied mesh value

        base_input_range[,i] <- customcontrsum(factorin = base_input_range[,i], speclevels = stepseq[[i]])
      }

    }


    #Translate all factor or character input columns in base input range to -1 to 1 range for algebraic range finding
    base_input_range_numeric <- base_input_range

    for(i in 1:ncol(base_input_range)){


      #Expand factor into contrast sum coding for range of -1 to 1 for all variables
      if(is.numeric(base_input_range_numeric[[i]]) == FALSE){

        #Expand coefficients based on number of factor levels
        tempframe <- model.matrix(as.formula(paste0("~", colnames(base_input_range)[i])), base_input_range)

        #Set min to -1 and max to 1 for all factors
        tempframe[1,] <- -1
        tempframe[2,] <- 1

        #Add expanded factor names to base_input_range_numeric. Exclude intercept column
        removecol <- ncol(base_input_range_numeric) + 1
        base_input_range_numeric <- cbind.data.frame(base_input_range_numeric, tempframe)[,-removecol]

      }

    }

    #Find min and max range for all base vars and algebraic combinations in formulalist (WARNING: If surface is very complex, THIS MAY FAIL)
    list_input_ranges <- lapply(formulalist, algebraic_range, base_var_range = base_input_range_numeric)
    all_input_ranges <- data.frame(list_input_ranges)
    colnames(all_input_ranges) <- c(unlist(lapply(list_input_ranges,colnames)))
    all_input_ranges <- as.matrix(all_input_ranges) #Convert from list to matrix
    all_input_ranges <- all_input_ranges[,unique(colnames(all_input_ranges))] # Remove redundant columns from matrix

    if(is.null(startingdesign) == TRUE){
    #Create random model matrix using base vars to start optimization function and assign input names to columns
    ModelMatReal <- data.frame(lapply(stepseq, sample, size = model_points, replace = TRUE))
    colnames(ModelMatReal) <- input_list

    }else{

      #Use provided design as starting point
      ModelMatReal <- startingdesign[,input_list]
    }

    #Change coding for all factor columns to contr.sum for proper balancing for d-error
    for(i in 1:ncol(ModelMatReal)){

      if(is.factor(ModelMatReal[,i]) == TRUE){

        ModelMatReal[,i] <- customcontrsum(factorin = ModelMatReal[,i], speclevels = stepseq[[i]])

      }

    }

    #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
    d_efficiency_vect <- Vectorize(d_effchoice, c("CurrentMatrix", "paramestimates"))

    #Expand model matrix for each formula
    candexpand <- lapply(formulalist, model.matrix, data = ModelMatReal)

    #Standardize numeric columns
    candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))

    #Calculate initial D-error for randomly seeded model
    eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x) x[,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

    #Calculate starting objective function using vector of weights
    if(missing(weight))
    {obj_current <- min(eff_vect)}else{
      obj_current <- sum(eff_vect*weight)}

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one parameter at a time
      for (i in 1:nrow(ModelMatReal)){
        for(j in 1:ncol(ModelMatReal)){

          #Step through -1 to 1 by step value
          for (s in stepseq[[j]]){

            #Replace current value in ModelMatReal with new temporary point
            oldvalue <- ModelMatReal[i,j]
            ModelMatReal[i,j] <- s

            #Expand model matrix for each formula
            candexpand <- lapply(formulalist, model.matrix, data = ModelMatReal)

            #Standardize numeric columns
            candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))

            #Calculate D-efficiency
            eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x) x[,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list


            #Calculate new objective function value using vector of weights and temporary point
            if(missing(weight))
            {obj_temp <- min(eff_vect)}else{
              obj_temp <- sum(eff_vect*weight)}

            #If objective of new point is greater than or equal to old point, use new point in matrix. Otherwise, put old point back into matrix
            if(obj_temp < obj_current)

            {ModelMatReal[i,j] <- oldvalue}else{
              obj_current <- obj_temp}
          }
        }
      }
      #Check objective function after replacement to see if improv
      objective_change <- obj_current - obj_prev
    }


    #Expand final model matrix for each formula
    candexpand <- lapply(formulalist, model.matrix, data = ModelMatReal)

    #Standardize numeric columns
    candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))

    #Calculate final D-efficiencies
    eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x) x[,2:ncol(x)]), paramestimates = priors, altvect = altvect, returncov = TRUE)
  }

  if(searchstyle == "Federov"){

    #Reduce candidate set to only include base variables in the supplied formulas to eliminate extraneous data
    candset <- candset[,input_list]

    #Identify all columns that are coded as factors
    factcols <- rep(FALSE, ncol(candset))

    for(i in 1:ncol(candset)){

      factcols[i] <- is.factor(candset[,i])

    }


    #Change frame contrasts on all factors to contr.sum for proper balancing for d-error
    for (i in 1:length(factcols)){

      if(factcols[i] == TRUE){

        candset[,i] <- customcontrsum(factorin = candset[,i], speclevels = levels(candset[,i]))
        #contrasts(candset[,i]) <- contr.sum(nlevels(candset[,i]))

      }
    }

    #Create model matrix of candidate set for each input formula
    candexpand <- lapply(input_formulas, model.matrix, data = candset)

    #Get all input ranges from candidate set
    list_input_ranges <- lapply(candexpand, function(x) apply(x, MARGIN = 2, function(y) c(min(y), max(y)))[,2:ncol(x)])
    all_input_ranges <- data.frame(list_input_ranges)
    colnames(all_input_ranges) <- c(unlist(lapply(list_input_ranges,colnames)))
    all_input_ranges <- as.matrix(all_input_ranges) #Convert from list to matrix
    all_input_ranges <- all_input_ranges[,unique(colnames(all_input_ranges))] # Remove redundant columns from matrix

    #Standardize all numeric or integer variables

    if(sum(factcols == FALSE) > 0){
      candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))
    }

    if(is.null(startingdesign) == TRUE){

      #Select random row indices from supplied candidate set of model points. Reduce matrix
      rownums <- sample(nrow(candset), size = model_points, replace = TRUE)

    }else{

      #Find corresponding rows based on starting design
      rownums <- apply(startingdesign[,input_list], MARGIN = 1, function(x){

        rowmatch <- sapply(input_list, function(y) candset[,y] == x[[y]])
        rowmatch <- which(apply(rowmatch, MARGIN = 1, prod) == 1)
        rowmatch <- rowmatch[1]

return(rowmatch)
})

    }

#     #Select random row indices from supplied candidate set of model points. Reduce matrix
#     rownums <- sample(nrow(candset), size = model_points, replace = TRUE)

    #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
    d_efficiency_vect <- Vectorize(d_effchoice, c("CurrentMatrix", "paramestimates"))

    #Calculate initial D-error for randomly seeded model
    eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

    #Calculate starting objective function using vector of weights
    if(missing(weight))
    {obj_current <- min(eff_vect)}
    else{obj_current <- sum(eff_vect*weight)}

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one row at a time
      for (i in 1:length(rownums)){

        #Try every possible point in candidate set of points for each row
        for(j in 1:nrow(candset)){

          #Replace current value in ModelMatReal with new temporary point
          oldvalue <- rownums[i]
          rownums[i] <- j

          #Calculate D-efficiencies for ModelMatReal with new temporary point
          eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

          #Calculate new objective function value using vector of weights and temporary point
          if(missing(weight))
          {obj_temp <- min(eff_vect)}
          else{obj_temp <- sum(eff_vect*weight)}

          #If objective of new point is greater than or equal to old point, use new point in matrix. Otherwise, put old point back into matrix
          if(obj_temp < obj_current)

          {rownums[i] <- oldvalue}

          else {obj_current <- obj_temp}
        }
      }
      #Check objective function after replacement to see if improv
      objective_change <- obj_current - obj_prev
    }

    #Convert final matrix to real values
    ModelMatReal <- candset[rownums,]

    #Calculate final d-efficiency
    eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect, returncov = TRUE)

  }

  #colnames(eff_vect_final) <- input_formulas #Name efficiency vector formulas

  #Name final output
  outputfinal <- list("ModelMatrix" = data.frame(Question = altvect, ModelMatReal),
                      "Deff" = unlist(eff_vect_final[1,]),
                      "DeffvsOptimal" = unlist(eff_vect_final[1,])/det_ref_list,
                      "CovList" = eff_vect_final[2,],
                      "ObjectiveFunction" = obj_current)
  #names(outputfinal) <- c("ModelMatrix","D-Efficiency","ObjectiveFunction")

  #Return model matrix, objective function value, and 1/D-optimality for add, mech model

  return(outputfinal)
}

##########################################################################################################################################
