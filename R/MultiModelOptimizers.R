############################DEPRECATED FUNCTION. ALL FUNCTIONS AND MORE CAN NOW BE EXECUTED USING MultipleModelOptimize#############################################

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
#searchstyle - option for searching experimental space. "Fedorov" samples from a supplied matrix of available designs (candset). "Columnwise" samples across the range of each base variable by step sizes defined by mesh.
#candset - set of candidate points to search through for model creation. This is required if "Fedorov" searchstyle is used but is not used otherwise

#Inputs, optional:
#weight - vector of weighting values to apply to each formula in formulalist. Order of this vector matches order of input formulas. Format = c(0.5,0.25, etc)
#		if weight is not supplied, function will attempt to maximize lowest d_efficiency using det_ref_list

# CEX_MultipleModel <- function(base_input_range, formulalist, model_points, blocks = NA, det_ref_list, mesh, tolerance, weight, candset = NA, searchstyle = "Fedorov"){
#
#   #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
#   # step <- 2/(mesh-1)
#   # stepseq <- seq(-1,1, by = step)
#
#   #library(nlme)
#
#   #Pull inputs only from all formulas and place into a list
#   inputs_only <- lapply(formulalist, list_input_variables)
#
#   #Create list of all base input variables used in formulalist
#   input_list <- unique(unlist(inputs_only))
#
#   #Find min and max range for all base vars and algebraic combinations in formulalist (WARNING: If surface is very complex, THIS MAY FAIL)
#   list_input_ranges <- lapply(formulalist, algebraic_range, base_var_range = base_input_range)
#   all_input_ranges <- data.frame(base_input_range, list_input_ranges)
#   colnames(all_input_ranges) <- c(colnames(base_input_range), unlist(lapply(list_input_ranges,colnames)))
#   all_input_ranges <- as.matrix(all_input_ranges) #Convert from list to matrix
#   all_input_ranges <- all_input_ranges[,unique(colnames(all_input_ranges))] # Remove redundant columns from matrix
#
#   #Create list of one-sided expressions of input only side for all formulas in formulalist
#   input_formulas <- sapply(formulalist, nlme::splitFormula)
#
#   #Optimize via Columnwise search if that option is selected
#
#   if(searchstyle == "Columnwise"){
#     #If mesh supplied as a single value, apply it to all base input values
#     if(length(mesh) == 1){
#
#       mesh <- rep(mesh, times = length(input_list))
#
#     }
#
#
#     #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
#     stepseq <- list()
#
#     for(i in 1:length (input_list)){
#
#       step <- 2/(mesh[i]-1)
#
#       stepseq[[i]] <- seq(-1,1, by = step)
#
#     }
#
#     # #Create random model matrix using base vars (standardized {-1 to 1}) to start optimization function and assign input names to columns
#     # ModelMatStand <- matrix(data = sample(stepseq, length(input_list)*model_points, replace = TRUE), nrow = model_points, ncol = length(input_list))
#     # colnames(ModelMatStand) <- input_list
#
#     #Create random model matrix using base vars (standardized {-1 to 1}) to start optimization function and assign input names to columns
#     ModelMatStand <- sapply(stepseq, sample, size = model_points, replace = TRUE)
#     colnames(ModelMatStand) <- input_list
#
#     #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
#     d_efficiency_vect <- Vectorize(d_efficiency, c("input_formula", "det_ref"))
#
#     #Calculate initial D-efficiencies for randomly seeded model
#     ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE) #Convert to real numbers before passing to d_efficiency
#     eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)
#
#     #Calculate starting objective function using vector of weights
#     if(missing(weight))
#     {obj_current <- min(eff_vect[1,])}
#     else{obj_current <- sum(eff_vect[1,]*weight)}
#
#     #Initialize objective change value to start while loop
#     objective_change <- tolerance + 1
#
#     #Outer loop to run until tolerance achieved
#     while (objective_change > tolerance){
#
#       #Store objective function value before executing replacement loop
#       obj_prev <- obj_current
#
#       #Loop to make single pass through design, one parameter at a time
#       for (i in 1:nrow(ModelMatStand)){
#         for(j in 1:ncol(ModelMatStand)){
#
#           #Step through -1 to 1 by step value
#           for (s in stepseq[[j]]){
#
#             #Replace current value in ModelMatStand with new temporary point
#             oldvalue <- ModelMatStand[i,j]
#             ModelMatStand[i,j] <- s
#
#             #Calculate D-efficiencies for ModelMatStand with new temporary point
#             ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE) #Convert to real numbers before passing to d_efficiency
#             eff_vect <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)
#
#             #Calculate new objective function value using vector of weights and temporary point
#             if(missing(weight))
#             {obj_temp <- min(eff_vect[1,])}
#             else{obj_temp <- sum(eff_vect[1,]*weight)}
#
#             #If objective of new point is greater than old point, use new point in matrix. Otherwise, put old point back into matrix
#             if(obj_temp <= obj_current)
#
#             {ModelMatStand[i,j] <- oldvalue}
#
#             else {obj_current <- obj_temp}
#           }
#         }
#       }
#       #Check objective function after replacement to see if improv
#       objective_change <- obj_current - obj_prev
#     }
#
#     #Convert final matrix to real values
#     ModelMatReal <- standardize_cols(ModelMatStand, input_list, all_input_ranges, reverse_standard = TRUE)
#
#     #Calculate final d-efficiency
#     eff_vect_final <- d_efficiency_vect(CurrentMatrix = ModelMatReal, det_ref = det_ref_list, input_formula = input_formulas, Input_range = all_input_ranges)
#
#   }
#
#   if(searchstyle == "Fedorov"){
#
#     #Reduce candidate set to only include base variables in the supplied formulas to eliminate extraneous data
#     candset <- candset[,input_list]
#
#     #Create model matrix of candidate set for each input formula
#     candexpand <- lapply(input_formulas, model.matrix, data = candset)
#
#     #Standardize all candidate sets within each variable range
#     candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x[,2:ncol(x)]), Input_range = inrange))
#
#     #Convert determinant list to vector
#     det_refs <- unlist(det_ref_list)
#
#     #Select random row indices from supplied candidate set of model points. Reduce matrix
#     rownums <- sample(nrow(candset), size = model_points, replace = TRUE)
#     #ModelMatReal <- candset[sample(nrow(candset),size = model_points, replace = TRUE),]
#
#     #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
#     d_efficiency_vect <- Vectorize(d_efficiencysimple, c("CurrentMatrix", "det_ref"))
#
#     #Calculate initial D-efficiencies for randomly seeded model
#     eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)
#
#     #Calculate starting objective function using vector of weights
#     if(missing(weight))
#     {obj_current <- min(eff_vect[1,])}
#     else{obj_current <- sum(eff_vect[1,]*weight)}
#
#     #Initialize objective change value to start while loop
#     objective_change <- tolerance + 1
#
#     #Outer loop to run until tolerance achieved
#     while (objective_change > tolerance){
#
#       #Store objective function value before executing replacement loop
#       obj_prev <- obj_current
#
#       #Loop to make single pass through design, one row at a time
#       for (i in 1:length(rownums)){
#
#         #Try every possible point in candidate set of points for each row
#         for(j in 1:nrow(candset)){
#
#           #Replace current value in ModelMatReal with new temporary point
#           oldvalue <- rownums[i]
#           rownums[i] <- j
#
#           #Calculate D-efficiencies for ModelMatReal with new temporary point
#           eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)
#
#           #Calculate new objective function value using vector of weights and temporary point
#           if(missing(weight))
#           {obj_temp <- min(eff_vect[1,])}
#           else{obj_temp <- sum(eff_vect[1,]*weight)}
#
#           #If objective of new point is greater than old point, use new point in matrix. Otherwise, put old point back into matrix
#           if(obj_temp <= obj_current)
#
#           {rownums[i] <- oldvalue}
#
#           else {obj_current <- obj_temp}
#         }
#       }
#       #Check objective function after replacement to see if improv
#       objective_change <- obj_current - obj_prev
#     }
#
#     #Convert final matrix to real values
#     ModelMatReal <- candset[rownums,]
#
#     #Calculate final d-efficiency
#     eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,]), det_ref = det_ref_list)
#
#   }
#
#   colnames(eff_vect_final) <- input_formulas #Name efficiency vector formulas
#
#   #Name final output
#   outputfinal <- list(ModelMatReal, eff_vect_final, obj_current)
#   names(outputfinal) <- c("ModelMatrix","D-Efficiency","ObjectiveFunction")
#
#   #Return model matrix, objective function value, and 1/D-optimality for add, mech model
#
#   return(outputfinal)
# }

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
#searchstyle - option for searching experimental space. "Fedorov" samples from a supplied matrix of available designs (candset). "Columnwise" samples across the range of each base variable by step sizes defined by mesh.
#candset - set of candidate points to search through for model creation. This is required if "Fedorov" searchstyle is used but is not used otherwise

#Inputs, optional:
#weight - vector of weighting values to apply to each formula in formulalist. Order of this vector matches order of input formulas. Format = c(0.5,0.25, etc)
#		if weight is not supplied, function will attempt to maximize lowest d_efficiency using det_ref_list
#priors - list of known or estimated effect sizes (in model terms). Structure should be a list as follows: list(I(A^2) = 0.23, B = -0.5) where parameter names in list match those in formulas. For attribute parameters, list should be effect size with reference to first factor level. Any effects not provided will be assumed to be equal to zero.
#startingdesign - design used as startingpoint for optimization. If not provided, a random design will be selected as the starting point. Will cause an error if this is not properly formatted.
##best practice is to use a design created by this function as a starting desgin for another round of optimization
#Optimality - character defining which optimality criterion to use. More to be added soon Choices are:
##D - d-optimality
##D-Util - blend of d-optimality and utility balance in choice sets (for choice models only)
#OptBlend - vector with scaling factor to apply to each efficiency type before adding together for total optimality when there are multiple (D-Util only currently)

MultipleModelOptimize <- function(base_input_range, formulalist, questions, alts, blocks = NA, optout = FALSE, det_ref_list, mesh, tolerance, weight, candset = NA, priors = NA, searchstyle = "Fedorov", startingdesign = NULL, eqtype = NULL, questionblockvars = NULL, augment = FALSE, priorsnormalized = FALSE, Optimality = "D", OptBlend = c(0.5,0.5)){

  #Calculate number of model points to use
  model_points <- questions*alts

  #library(nlme)

  #Generate vector that describes sequence of question vs model points
  altvect <- c()

  for(i in 1:questions){
    altvect <- c(altvect, rep(i, alts))
  }

  #Define starting row for new designs as row 1 for both the columnwise and fedorov optimizer below
  optimizerowstart <- 1


  if(augment == TRUE){
    #For designs that will be augmented, add the alternative vector generated above to the existing vector from staringdesign
    #Format expected for this is a column labeled "Question"

    altvect <- c(startingdesign$Question, altvect + max(startingdesign$Question))

    #Define starting row as row after the starting design
    optimizerowstart <- nrow(startingdesign)+1

  }

  #Convert determinant list to vector
  det_ref_list <- unlist(det_ref_list)

  #Create list of largest possible variance for all choice equations (assumes basic multinomial logit)
  p_var_ref <- 0

  #Loop through all choice sets and calculate maximum possible variance of each set
  for(i in unique(altvect)){

    #Calculate maximum variance of set and add to total variance
    p_var_ref <- p_var_ref + (1/length(which(altvect == i)))^length(which(altvect == i))

  }

  #Pull inputs only from all formulas and place into a list
  inputs_only <- lapply(formulalist, list_input_variables)

  #Create list of all base input variables used in formulalist
  input_list <- unique(unlist(inputs_only))

  #Create list of one-sided expressions of input only side for all formulas in formulalist
  input_formulas <- sapply(formulalist, nlme::splitFormula)

  #If eqtype is not populated, assume all equations are choice equations
  if(is.null(eqtype) == TRUE){

    eqtype <- rep("Choice", length(formulalist))

  }

  #Create vector of choice equation numbers
  choiceeqs <- which(eqtype == "Choice")

  #Create vector of linear equation numbers
  lineareqs <- which(eqtype == "Linear")

  #Optimize via Columnwise search if that option is selected
  if(searchstyle == "Columnwise"){

    #If mesh supplied as a single value, apply it to all base input values
    if(length(mesh) == 1){

      mesh <- rep(mesh, times = length(input_list))

    }


    #Determine step number based on mesh input and construct sequence from -1 to 1 by step size
    stepseq <- list()

    for(i in 1:length (input_list)){

      #For numeric values, create steps based on mesh value from minimum supplied value to maximum supplied value
      if(is.numeric(base_input_range[[input_list[[i]]]]) == TRUE){

        #If the length of the mesh entry is equal to 1, it refers to the number of steps to break the search into.
        if(length(mesh[[i]]) == 1){
          step <- (max(base_input_range[[input_list[[i]]]]) - min(base_input_range[[input_list[[i]]]]))/(mesh[[i]]-1)

          stepseq[[i]] <- seq(min(base_input_range[[input_list[[i]]]]),max(base_input_range[[input_list[[i]]]]), by = step)
        }

        #If the length of the mesh entry is greater than 1, it is interpreted as a vector of valid values to use in the search
        if(length(mesh[[i]]) > 1){

          stepseq[[i]] <- mesh[[i]]

        }

      }

      #For attribute values already classified as factors, either use supplied factor levels or use all levels present in the factor in base_input_range
      if((is.factor(base_input_range[[input_list[[i]]]]) == TRUE) | (is.character(base_input_range[[input_list[[i]]]]) == TRUE)){

        #If mesh values were provided for this factor, use them. Otherwise return all levels of base_input_range as possibilities.
        if((is.factor(mesh[[i]]) == TRUE) | (is.character(mesh[[i]]) == TRUE)){

          stepseq[[i]] <- as.factor(mesh[[i]])

        }else{

          stepseq[[i]] <- levels(as.factor(base_input_range[[input_list[[i]]]]))

        }

        #Add factor levels present in startingdesign that are not in base_input_range or mesh
        if(is.null(startingdesign) == FALSE){

          levels(stepseq[[i]]) <- c(levels(stepseq[[i]]),
                                    levels(as.factor(startingdesign[[input_list[[i]]]]))[which(!levels(as.factor(startingdesign[[input_list[[i]]]])) %in% stepseq[[i]])])

        }

      }

    }


    for(i in 1:ncol(base_input_range)){

      #Apply contrast sum encoding to base_input_range
      if(is.numeric(base_input_range[[i]]) == FALSE){

        #Check factor levels in base_input_range against supplied mesh value

        #base_input_range[,i] <- customcontrsum(factorin = base_input_range[,i], speclevels = stepseq[[i]])

        #New code to handle levels in startingdesign. Note that vector concatenation at end is necessary to preserve order
        base_input_range[,i] <- customcontrsum(factorin = base_input_range[,i], speclevels = c(paste(stepseq[[i]]), levels(stepseq[[i]])[which(!levels(stepseq[[i]])%in%stepseq[[i]])]))
      }else{

        ####Finish Here#######Do I actually want to redefine the range based on the previous design? Should we consider leaving base_input_range as is to reflect the new design space?
        #Expand numeric base input range to include range in startingdesign
        base_input_range[1,i] <- min(c(base_input_range[1,i], startingdesign[[colnames(base_input_range)[i]]]))
        base_input_range[2,i] <- max(c(base_input_range[2,i], startingdesign[[colnames(base_input_range)[i]]]))
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

      if(augment == TRUE){

        ModelMatReal <- data.frame(lapply(stepseq, sample, size = model_points, replace = TRUE))
        colnames(ModelMatReal) <- input_list

        ModelMatReal <- rbind.data.frame(startingdesign[,input_list], ModelMatReal)

      }else{
        #Use provided design as starting point
        ModelMatReal <- startingdesign[,input_list]
      }
    }

    #Change coding for all factor columns to contr.sum for proper balancing for d-error
    for(i in 1:ncol(ModelMatReal)){

      if((is.factor(ModelMatReal[,i]) == TRUE)|(is.character(ModelMatReal[,i]) == TRUE)){

        ModelMatReal[,i] <- customcontrsum(factorin = ModelMatReal[,i], speclevels = levels(base_input_range[[colnames(ModelMatReal)[i]]]))

      }

    }

    #Vectorize full D-efficiency function to be able to call later using lists of model matrix and parameter estimates
    d_efficiency_vect <- Vectorize(d_effchoice, c("CurrentMatrix", "paramestimates"))

    #Vectorize d-efficiency update function for multinomial models to be able to call later using lists of model matrix, parameter estimates, and existing information matrix
    d_efficiency_update_vect <- Vectorize(d_effchoiceupdate, c("CurrentMatrix", "paramestimates", "info_mat"))

    #Expand model matrix for each formula
    candexpand <- lapply(formulalist, model.matrix, data = ModelMatReal)

    #Standardize numeric columns
    candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))

    #Rescale priors to model matrix if they were not supplied in scaled format

    if(priorsnormalized == FALSE){

      for(i in 1:length(candexpand)){

        for(j in 2:ncol(candexpand[[i]])){

          #Recalculate prior for -1 to 1 normalization of related variable
          priors[[i]][j-1] <- priors[[i]][j-1]*0.5*(max(all_input_ranges[1:2,colnames(candexpand[[i]])[j]])-min(all_input_ranges[1:2,colnames(candexpand[[i]])[j]]))


        }

      }

    }

    #Initialize efficiency vector
    eff_vect <- rep(0, length(formulalist))

    ##############################New Code#########################################

    #Ensure that blocking variables by question for choice designs are set to the same value across the question choice sets
    if(is.null(questionblockvars) == FALSE){
      for(i in 1:length(questionblockvars)){
        for(j in 1:questions){

          ModelMatReal[[questionblockvars]][altvect == j] <- ModelMatReal[[questionblockvars]][altvect == j][1]

        }
      }
    }

    #######################

    #Calculate initial D-efficiency for randomly seeded model for Choice equations
    if(length(choiceeqs) > 0){
      eff_vect[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)/det_ref_list[choiceeqs]
    }

    #Calculate d-efficiency for linear equations
    if(length(lineareqs) > 0){

      eff_vect[lineareqs] <- sapply(candexpand[lineareqs], d_efficiencysimple)/det_ref_list[lineareqs]

    }

    #Calculate starting objective function using vector of weights
    if(missing(weight)){obj_current <- min(eff_vect)}else{
      obj_current <- sum(eff_vect*weight)}

    #Initialize the current question for choice experiments
    currentquestion <- 0

    #Initialize info_mat list for current information matrix of non-updated question
    info_mat <- as.list(rep(0, length(formulalist)))

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one parameter at a time
      for (i in optimizerowstart:nrow(ModelMatReal)){

        #Calculate d-efficiency for choice models
        if(length(choiceeqs) > 0){

          #Check current question number. If it has changed, recalculate the logistic multinomial model information matrix for all questions not being updated in this loop
          if(currentquestion != altvect[i]){

            #Update current question to match current row of design being updated
            currentquestion <- altvect[i]

            #Calculate information matrix for all of design except current question being updated
            info_mat[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[(altvect != currentquestion),2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect[altvect != currentquestion], returninfomat = TRUE)[2,]

          }

        }

        #Calculate d-efficiency for linear models
        if(length(lineareqs) > 0){

          eff_vect[lineareqs] <- sapply(candexpand[lineareqs], d_efficiencysimple)/det_ref_list[lineareqs]

        }



        for(j in 1:ncol(ModelMatReal)){

          #Step through -1 to 1 by step value
          for (s in stepseq[[j]]){

            #Replace current value in ModelMatReal with new temporary point
            oldvalue <- ModelMatReal[i,j]

            #Check variables against blocking variable and make sure this value stays the same across this question if it is a blocking variable
            if((colnames(ModelMatReal)[j] %in% questionblockvars) == FALSE){
              #Replace only the current point if this is not a blocking variable
              ModelMatReal[i,j] <- s
            }else{
              #Replace all sets in this question if this is a blocking variable
              ModelMatReal[[questionblockvars]][altvect == altvect[i]] <- s
            }

            #Expand model matrix for each formula
            candexpand <- lapply(formulalist, model.matrix, data = ModelMatReal)

            #Standardize numeric columns
            candexpand <- lapply(candexpand, function(x, inrange = all_input_ranges) standardize_cols(StartingMat = x, column_names = colnames(x)[colnames(x) %in% colnames(inrange)], Input_range = inrange))

            ##Old method that only considered one form of d-efficiency
            #             #Calculate D-efficiency
            #             eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x) x[,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

            #Calculate d-efficiency for choice models
            if(length(choiceeqs) > 0){

              #Update d-efficiency with new point (must send entire alternative set for this question)
              eff_vect[choiceeqs] <- d_efficiency_update_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[(altvect == currentquestion),2:ncol(x)]), paramestimates = priors[choiceeqs], info_mat = info_mat[choiceeqs])/det_ref_list[choiceeqs]

              #eff_vect[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)/det_ref_list[choiceeqs]
            }

            #Calculate d-efficiency for linear models
            if(length(lineareqs) > 0){

              eff_vect[lineareqs] <- sapply(candexpand[lineareqs], d_efficiencysimple)/det_ref_list[lineareqs]

            }

            #Calculate new objective function value using vector of weights and temporary point
            if(missing(weight)){obj_temp <- min(eff_vect)}else{
              obj_temp <- sum(eff_vect*weight)}

            #If objective of new point is greater than or equal to old point, use new point in matrix. Otherwise, put old point back into matrix
            if(obj_temp <= obj_current){

              if((colnames(ModelMatReal)[j] %in% questionblockvars) == FALSE){
                #Replace only the current point if this is not a blocking variable
                ModelMatReal[i,j] <- oldvalue
              }else{
                #Replace all sets in this question if this is a blocking variable
                ModelMatReal[[questionblockvars]][altvect == altvect[i]] <- oldvalue
              }
            }else{
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

    #Initialize final d-efficiency, optimality, covariance, and additional diagnostic lists
    eff_vect_final <- as.list(eff_vect)

    opt_vect_final <- as.list(eff_vect)

    cov_vect_final <- as.list(rep(0, length(formulalist)))


    additionaldiags <- list()
    #Store probabilities where applicable
    additionaldiags$probvect <- as.list(rep(NA, length(formulalist)))
    #Store probability variance where applicable
    additionaldiags$p_var <- as.list(rep(NA, length(formulalist)))
    #Store probability variance optimality (ratio of variance/max variance)
    additionaldiags$p_opt <- as.list(rep(NA, length(formulalist)))


    #Calculate final d-efficiency for choice models
    if(length(choiceeqs) > 0){

      #Calculate final d-efficiency

      #if(Optimality == "D"){
      eff_vect_final[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)
      #}

#       if(Optimality == "D-Util"){
#
#         temp_eff <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)
#
#         #Blend d-optimality with utility balance (ratio of variance/max variance) based on OptBlend
#         eff_vect_final[choiceeqs] <- OptBlend[1]*unlist(temp_eff[1,])/det_ref_list[choiceeqs]+OptBlend[2]*unlist(temp_eff[3,])/p_var_ref
#
#       }

      #Calculate final covlist

      cov_vect_final[choiceeqs] <- lapply(d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x) x[,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)[2,], function(x) tryCatch(solve(x), error = function(x) diag(x = Inf, ncol = 2, nrow = 2)))
    }

    #Calculate final d-efficiency for linear models
    if(length(lineareqs) > 0){

      #Calculate final d-efficiency
      eff_vect_final[lineareqs] <- sapply(candexpand[lineareqs], d_efficiencysimple)
      #Calculate final covlist
      cov_vect_final[lineareqs] <- sapply(candexpand[lineareqs], d_efficiencysimple, returncov = TRUE)[2,]

    }

    ##########################Placeholder section prior to implementing D-Util calcs for columnwise search
    opt_vect_final <- unlist(eff_vect_final)/det_ref_list
    #######################################################################

  }

  if(searchstyle == "Fedorov"){

    #Reduce candidate set to only include base variables in the supplied formulas to eliminate extraneous data
    candset <- candset[,input_list]

    if(is.null(startingdesign) == FALSE){

      #If there is a starting design provided, combine this design and the candidate set into one candidate set
      candsetnew <- rbind.data.frame(startingdesign[,input_list], candset)

      #Loop through all factor columns and ensure they are recoded as factors with new levels added as necessary
      for(i in 1:ncol(candsetnew)){

        if((is.factor(candsetnew[,i]) == TRUE)|(is.character(candsetnew[,i]) == TRUE)){

          candsetnew[,i] <- customcontrsum(factorin = candsetnew[,i], speclevels = c(levels(candset[[colnames(candsetnew)[i]]]),
                                                                                     levels(as.factor(candsetnew[,i]))[which(!levels(as.factor(candsetnew[,i])) %in% levels(candset[[colnames(candsetnew)[i]]]))]))

        }

      }

      #Use this as candset moving forward
      candset <- candsetnew
    }

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

    #Rescale priors to model matrix if they were not supplied in scaled format
    if(priorsnormalized == FALSE){

      for(i in 1:length(candexpand)){

        for(j in 2:ncol(candexpand[[i]])){

          #Recalculate prior for -1 to 1 normalization of related variable
          priors[[i]][j-1] <- priors[[i]][j-1]*0.5*(max(all_input_ranges[1:2,colnames(candexpand[[i]])[j]])-min(all_input_ranges[1:2,colnames(candexpand[[i]])[j]]))


        }

      }

    }

    if(is.null(startingdesign) == TRUE){

      #Select random row indices from supplied candidate set of model points. Reduce matrix
      rownums <- sample(c(optimizerowstart:nrow(candset)), size = model_points, replace = TRUE)

    }else{

      #If a starting design was provided, use the row numbers of candset as these are equal to startingdesign
      rownums <- c(1:nrow(startingdesign))

      if(augment == TRUE){
        #Randomly sample additional points from the supplied candidate set of model points
        rownums <- c(1:nrow(startingdesign), sample(c(optimizerowstart:nrow(candset)), size = model_points, replace = TRUE))


        #       ##Matching code is no longer necessary as existing points are added in their entirety to the first rows of the candset
        #       for(i in 1:nrow(startingdesign)){
        #
        #         rowmatch <- sapply(input_list, function(y) candset[,y] == startingdesign[i,input_list][[y]])
        #         rowmatch <- which(apply(rowmatch, MARGIN = 1, prod) == 1)
        #         rownums[i] <- as.integer(rowmatch[1])
        #
        #       }

      }
    }

    #Vectorize D-efficiency function to be able to call later using lists of formulas and reference determinants
    d_efficiency_vect <- Vectorize(d_effchoice, c("CurrentMatrix", "paramestimates"))

    #Vectorize d-efficiency update function for multinomial models to be able to call later using lists of model matrix, parameter estimates, and existing information matrix
    d_efficiency_update_vect <- Vectorize(d_effchoiceupdate, c("CurrentMatrix", "paramestimates", "info_mat"))

    #Initialize efficiency vector
    eff_vect <- rep(0, length(formulalist))
    # #Original single equation type code
    #     #Calculate initial D-error for randomly seeded model
    #     eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

    #Calculate initial d-efficiency for both model types
    if(length(choiceeqs) > 0){

      #browser()
      #Calculate choice d-efficiency
      if(Optimality == "D"){
      eff_vect[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)/det_ref_list[choiceeqs]
      }
      if(Optimality == "D-Util"){

        temp_eff <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)

        #Blend d-optimality with utility balance (ratio of variance/max variance) based on OptBlend
        eff_vect[choiceeqs] <- OptBlend[1]*unlist(temp_eff[1,])/det_ref_list[choiceeqs]+OptBlend[2]*unlist(temp_eff[3,])/p_var_ref

      }
    }
    if(length(lineareqs) > 0){

      #Calculate linear d-efficiency
      eff_vect[lineareqs] <- sapply(lapply(candexpand[lineareqs], function(x, rowind = rownums) x[rowind,]), d_efficiencysimple)/det_ref_list[lineareqs]
    }

    #Calculate starting objective function using vector of weights
    if(missing(weight))
    {obj_current <- min(eff_vect)}
    else{obj_current <- sum(eff_vect*weight)}

    #Initialize the current question for choice experiments
    currentquestion <- 0

    #Initialize info_mat list for current information matrix of non-updated question
    info_mat <- as.list(rep(0, length(formulalist)))

    #Initialize p_var vector for current variance of non-updated question
    p_var_static <- rep(p_var_ref,length(formulalist))

    #Initialize objective change value to start while loop
    objective_change <- tolerance + 1

    #Outer loop to run until tolerance achieved
    while (objective_change > tolerance){

      #Store objective function value before executing replacement loop
      obj_prev <- obj_current

      #Loop to make single pass through design, one row at a time
      for (i in optimizerowstart:length(rownums)){

        #Check current question number. If it has changed, recalculate the logistic multinomial model information matrix for all questions not being updated in this loop
        if(currentquestion != altvect[i]){

          #Update current question to match current row of design being updated
          currentquestion <- altvect[i]
          #browser()

          #Calculate information matrix for all of design except current question being updated
          temp_eff <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind[(altvect != currentquestion)],2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect[altvect != currentquestion], returninfomat = TRUE)
          info_mat[choiceeqs] <- temp_eff[2,]

          #Calculate p_var for all of design except current question being updated for D-Util optimization
          p_var_static[choiceeqs] <- unlist(temp_eff[3,])


        }

        #Try every possible point in candidate set of points supplied by user for each row
        ##Note that this starts at optimizerowstart to account for the rows that were added to the start of the candidate set when a startingdesign is supplied
        for(j in optimizerowstart:nrow(candset)){

          #Replace current value in ModelMatReal with new temporary point
          oldvalue <- rownums[i]
          rownums[i] <- j

          #           #Old single equation type code
          #           #Calculate D-efficiencies for ModelMatReal with new temporary point
          #           eff_vect <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect)/det_ref_list

          #Calculate D-efficiencies for ModelMatReal with new temporary point
          if(length(choiceeqs) > 0){

            #Update d-efficiency with new point (must send entire alternative set for this question)
            #browser()
            temp_eff <- d_efficiency_update_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind[(altvect == currentquestion)],2:ncol(x)]), paramestimates = priors[choiceeqs], info_mat = info_mat[choiceeqs])
            #D-efficiency
            if(Optimality == "D"){
              eff_vect[choiceeqs] <- unlist(temp_eff[1,])/det_ref_list[choiceeqs]
            }

            #D-efficiency and utility balance
            if(Optimality == "D-Util"){

              eff_vect[choiceeqs] <- OptBlend[1]*unlist(temp_eff[1,])/det_ref_list[choiceeqs]+OptBlend[2]*(unlist(temp_eff[2,])+p_var_static[choiceeqs])/p_var_ref

            }
            #eff_vect[choiceeqs] <- unlist(d_efficiency_update_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind[(altvect == currentquestion)],2:ncol(x)]), paramestimates = priors[choiceeqs], info_mat = info_mat[choiceeqs])[1,])/det_ref_list[choiceeqs]

            ##Calculate choice d-efficiency
            #eff_vect[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)/det_ref_list[choiceeqs]
          }

          if(length(lineareqs) > 0){

            #Calculate linear d-efficiency
            eff_vect[lineareqs] <- sapply(lapply(candexpand[lineareqs], function(x, rowind = rownums) x[rowind,]), d_efficiencysimple)/det_ref_list[lineareqs]
          }

          #Calculate new objective function value using vector of weights and temporary point
          if(missing(weight))
          {obj_temp <- min(eff_vect)}
          else{obj_temp <- sum(eff_vect*weight)}

          #If objective of new point is greater than old point, use new point in matrix. Otherwise, put old point back into matrix
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

    #     #Old code for single equation type
    #     #Calculate final d-efficiency
    #     eff_vect_final <- d_efficiency_vect(CurrentMatrix = lapply(candexpand, function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors, altvect = altvect, returncov = TRUE)

    #Initialize final d-efficiency, optimality, covariance, and additional diagnostic lists
    eff_vect_final <- as.list(eff_vect)

    opt_vect_final <- as.list(eff_vect)

    cov_vect_final <- as.list(rep(0, length(formulalist)))

    additionaldiags <- list()
    #Store probabilities where applicable
    additionaldiags$probvect <- as.list(rep(NA, length(formulalist)))
    #Store probability variance where applicable
    additionaldiags$p_var <- as.list(rep(NA, length(formulalist)))
    #Store probability variance optimality (ratio of variance/max variance)
    additionaldiags$p_opt <- as.list(rep(NA, length(formulalist)))

    #Calculate final d-efficiency for choice models
    if(length(choiceeqs) > 0){

      #Calculate final model

      temp_eff <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)

      #Calculate final d-efficiency for D-efficiency model

      if(Optimality == "D"){

        #Extract d-efficiency
        eff_vect_final[choiceeqs] <- unlist(temp_eff[1,])

        #Scale to d-optimality
        opt_vect_final[choiceeqs] <- unlist(eff_vect_final[choiceeqs])/det_ref_list[choiceeqs]
      #eff_vect_final[choiceeqs] <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect)
      }

      if(Optimality == "D-Util"){

        temp_eff <- d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)

        #Pull d-efficiency into list
        eff_vect_final[choiceeqs] <- unlist(temp_eff[1,])

        #Blend d-optimality with utility balance (ratio of variance/max variance) based on OptBlend
        opt_vect_final[choiceeqs] <- OptBlend[1]*unlist(temp_eff[1,])/det_ref_list[choiceeqs]+OptBlend[2]*unlist(temp_eff[3,])/p_var_ref

      }

      #Calculate final covlist
      cov_vect_final[choiceeqs] <- lapply(d_efficiency_vect(CurrentMatrix = lapply(candexpand[choiceeqs], function(x, rowind = rownums) x[rowind,2:ncol(x)]), paramestimates = priors[choiceeqs], altvect = altvect, returninfomat = TRUE)[2,], function(x) tryCatch(solve(x), error = function(x) diag(x = Inf, ncol = 2, nrow = 2)))

      #Extract probability vectors
      additionaldiags$probvect[choiceeqs] <- temp_eff[4,]
      additionaldiags$p_var[choiceeqs] <- unlist(temp_eff[3,])
      additionaldiags$p_opt[choiceeqs] <- unlist(temp_eff[3,])/p_var_ref
    }

    #Calculate final d-efficiency for linear models
    if(length(lineareqs) > 0){

      #Calculate final d-efficiency
      eff_vect_final[lineareqs] <- sapply(lapply(candexpand[lineareqs], function(x, rowind = rownums) x[rowind,]), d_efficiencysimple)
      #Translate to final d-optimality
      opt_vect_final[lineareqs] <- unlist(eff_vect_final[lineareqs])/det_ref_list[lineareqs]
      #Calculate final covlist
      cov_vect_final[lineareqs] <- sapply(lapply(candexpand[lineareqs], function(x, rowind = rownums) x[rowind,]), d_efficiencysimple, returncov = TRUE)[2,]

    }
  }

  #colnames(eff_vect_final) <- input_formulas #Name efficiency vector formulas

  #Name final output
  outputfinal <- list("ModelMatrix" = data.frame(Question = altvect, ModelMatReal),
                      "Deff" = unlist(eff_vect_final),
                      "DeffvsOptimal" = opt_vect_final,
                      "CovList" = cov_vect_final,
                      "ObjectiveFunction" = obj_current,
                      "AdditionalDiagnostics" = additionaldiags)
  #names(outputfinal) <- c("ModelMatrix","D-Efficiency","ObjectiveFunction")

  #Return model matrix, objective function value, and 1/D-optimality for add, mech model
 #browser()
  return(outputfinal)
}

##########################################################################################################################################
