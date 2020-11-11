######################################################################
#' Function to estimate the expected value of partial perfect information
#' @param param_interest the main parameter of interest
#' @param value_param_interest value of the parameter of interest
#' @param names_params_needed names of needed parameters  from param
#' matrix returned
#' @param names_params_model names of parameters in the model
#' @param params_passed parameters passed while running the markov model
#' @param param_file all parameters required to run the model,provided with
#' name of the parameter, distribution, and parameters that define the
#' probability distribution
#' @param colnames_paramdistr col names where the parameter distribution
#' is defined
#' @return current_param_list list of parameters
#' @keywords internal
#' @details
#' this function gets all the parameters except the parameter of interest
#' if the parameter is fixed, just read from file or a distribution,
#' then gets it from a distribution, or if it to be calculated, just
#' give it back as it is. If the parameter is not found, it will be taken
#' same that as that of in the model parameter.
get_parameter_list <- function(param_interest, value_param_interest,
                               names_params_needed, names_params_model,
                               params_passed,
                               param_file, colnames_paramdistr) {
  list_checks <- list(param_interest, value_param_interest,
                  names_params_needed, names_params_model,
                  params_passed, param_file, colnames_paramdistr)

  results <- sapply(list_checks, packDAMipd::check_null_na)
  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")

  param_distr_data <- utils::read.csv(param_file, row.names = NULL)
  param_colno <- IPDFileCheck::get_columnno_fornames(param_distr_data,
                                                     "parameter")
  if (length(param_colno) > 1) {
    stop("Parameter column should be only one")
  }

  param_colname <- colnames(param_distr_data)[param_colno]
  distr_col_no <- IPDFileCheck::get_columnno_fornames(param_distr_data,
                                                     "distribution")
  if (length(distr_col_no) > 1) {
    stop("Distribution column should be only one")
  }

  current_param_needed <- c()
  current_param_names <- c()
  for (j in 2:length(names_params_needed)) {
    # fixing the param_interest , so that we don't need to sample
    if (names_params_needed[j] != param_interest) {
      # sample others - but identify those require calculations
      #if the selected parameters is needed to run the model
      if (names_params_needed[j] %in% names_params_model) {
        #if the params_passed for the current parameter considered
        # is a numerical value, that has been either read from a file or
        # directly assigned. It might not be calculated. so those can be
        # found from the parameter file
        res <- suppressWarnings(as.numeric(params_passed[names_params_needed[j]]))
        if (!is.na(res)) {
          # check if the parameter is defined using a distribution or actually
          # fixed (for death rates)
          row_cor_param <- param_distr_data[param_distr_data[param_colname] ==
                                              names_params_needed[j], ]
          if (nrow(row_cor_param) > 1) {
            stop("One row should correspond to the paramter in parameter file")
          }
          if (nrow(row_cor_param) == 0) {
            curr_parm_value <- params_passed[names_params_needed[j]]
          }else{
            distr_col_val_param <- row_cor_param[[distr_col_no]]
            # if the entry corresponding to distribution is null, na or empty
            # keep this param fixed
            if (packDAMipd::check_null_na(distr_col_val_param) < 0 |
                distr_col_val_param == "") {
              curr_parm_value <-
                packDAMipd::get_parameter_read(names_params_needed[j],
                                               param_file)
            }else{
              curr_parm_value <-
                packDAMipd::get_parameter_def_distribution(names_params_needed[j],
                                                           param_file,
                                                          colnames_paramdistr)
            }
          }
        }else{
          curr_parm_value <- params_passed[names_params_needed[j]]
        }
        curr_parm <- names_params_needed[j]
      }

    }else{
      curr_parm_value <- value_param_interest
      curr_parm <- param_interest
    }
    current_param_needed <- append(current_param_needed, curr_parm_value)
    current_param_names <- append(current_param_names, curr_parm)
  }
  current_param_list <- current_param_needed
  names(current_param_list) <- current_param_names
  return(current_param_list)
}
######################################################################
#' Function to get NMB for a set of parameter
#' This is the inner loop of the two stage Monte Carlo process
#' @param param_interest the main parameter of interest
#' @param value_param_interest value for the parameter of interest
#' @param param_file all parameters required to run the model,provided with
#' name of parameter, distribution and parameters that define the probability
#' distribution
#' @param colnames_paramdistr column names where parameter distribution and
#' values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold threshold values of WTP
#' @param comparator optional parameter if need ICER as the output
#' @return icer_nmb ICER and NMBs for the set of parameters and given
#' strategies
#' @keywords internal
#' @details
#' For the given set of parameter find the nmb from the packDAMipd function
#' calculate_icer_nmb
find_nmb_setofparams <- function(param_interest, value_param_interest,
                                 param_file, colnames_paramdistr,
                                 list_markov, threshold, comparator = NULL) {
  list_checks <- list(param_interest, value_param_interest,
                  threshold, param_file, colnames_paramdistr)
  results <- sapply(list_checks, packDAMipd::check_null_na)
  if (sum(results) != 0)
    stop("Some of the parameters may by NULL or NA - please check")

  if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }
  new_list_markov <- list()
  for (i in 1:nrow(list_markov)) {
    this <- list_markov[i, ]
    # names of strategies extracted - control, intervention etc
    # parameters needed get it from param_matrix
    names_params_needed <- colnames(list_markov[i, ]$param_matrix)
    names_params_model <- names(list_markov[i, ]$list_param_values)
    params_passed <- list_markov[i, ]$list_param_values
    current_param_list <- get_parameter_list(param_interest,
                                             value_param_interest,
                                             names_params_needed,
                                             names_params_model,
                                             params_passed, param_file,
                                             colnames_paramdistr)

    this_markov  <- packDAMipd::markov_model(current_strategy = this$strategy,
                                             cycles = this$cycles,
                                             initial_state = this$initial_state,
                                             discount = this$discount,
                                             parameter_values = current_param_list,
                                             half_cycle_correction =
                                               this$half_cycle_correction,
                                             state_cost_only_prevalent =
                                               this$state_cost_only_prevalent,
                                             state_util_only_prevalent =
                                               this$state_util_only_prevalent,
                                             method = this$method,
                                             startup_cost = this$startup_cost,
                                             startup_util = this$startup_util)

    new_list_markov[[length(new_list_markov) + 1]] <- this_markov
  }
  new_list_markov <- packDAMipd::combine_markov(new_list_markov)
  icer_nmb <- packDAMipd::calculate_icer_nmb(new_list_markov, threshold,
                                            comparator)
  return(icer_nmb)
}
######################################################################
#' Function to estimate the expected value of partial perfect information
#' This is the outer loop of the two stage Monte Carlo process
#' @param param_interest the main parameter of interest
#' @param param_file all parameters required to run the model,provided with
#' name of parameter, distribution and parameters that define the probability
#' distribution
#' @param colnames_paramdistr column names where parameter distribution
#' and values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold threshold values of WTP
#' @param outer_iterations number of iterations for outer loop
#' @param inner_iterations number of iterations for inner loop
#' @param comparator optional parameter if need ICER as the output
#' @return result result after evppi calculation, including three forms
#' of evppi
#' @source https://www.sciencedirect.com/science/article/pii/S1098301510605888
#' @export
estimate_evppi <- function(param_interest, param_file, colnames_paramdistr,
                        list_markov, threshold, outer_iterations,
                        inner_iterations,
                        comparator = NULL) {

  list_checks <- list(param_interest, param_file, colnames_paramdistr,
                  threshold, outer_iterations, inner_iterations)

  results <- sapply(list_checks, packDAMipd::check_null_na)

  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")

    if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }

  cols_needed <- nrow(list_markov)
  all_mean_thetac_nmb_thetai <- matrix(0, nrow = outer_iterations,
                                      ncol = cols_needed)
  all_mean_thetac_maxT_nmbs_thetai <- matrix(0, nrow = outer_iterations,
                                            ncol = 1)
  all_maxT_mean_thetac_maxT_nmbs_thetai <- matrix(0, nrow = outer_iterations,
                                                 ncol = 1)

  #step 9
  for (outer_loop in 1:outer_iterations) {
    #Step 1
    value_param_interest <- packDAMipd::get_parameter_def_distribution(
                            param_interest, param_file, colnames_paramdistr)

    all_nmbs_thetai <- matrix(0, nrow = inner_iterations, ncol = cols_needed)
    all_max_nmb_iteration <- list()
    #step  5
    for (inner_loop in 1:inner_iterations) {
      # step 2,3
      this_result <- find_nmb_setofparams(param_interest, value_param_interest,
                                         param_file, colnames_paramdistr,
                                      list_markov, threshold, comparator)
      nmbs_thetai <- (as.numeric(this_result$NMB))
      # store them
      all_nmbs_thetai[inner_loop, 1:cols_needed] <- nmbs_thetai
      # step 4
      maxT_nmbs_thetai <- max(as.numeric(nmbs_thetai))
      # store
      all_max_nmb_iteration[length(all_max_nmb_iteration) + 1] <-
        maxT_nmbs_thetai
    }
    # step 6
    mean_thetac_nmb_thetai <- colSums(all_nmbs_thetai) / inner_iterations
    # step 7
    mean_thetac_maxT_nmbs_thetai  <- mean(unlist(all_max_nmb_iteration))
    # step 8
    maxT_mean_thetac_maxT_nmbs_thetai <- max(mean_thetac_nmb_thetai)
    #store them
    all_mean_thetac_nmb_thetai[outer_loop, 1:cols_needed] <-
      mean_thetac_nmb_thetai
    all_mean_thetac_maxT_nmbs_thetai[outer_loop] <-
      mean_thetac_maxT_nmbs_thetai
    all_maxT_mean_thetac_maxT_nmbs_thetai[outer_loop] <-
      maxT_mean_thetac_maxT_nmbs_thetai
  }
  #step 10
  mean_thetai_mean_thetac_nmb_thetai <-
    colSums(all_mean_thetac_nmb_thetai) / outer_iterations
  mean_thetai_mean_thetac_maxT_nmb_thetai <-
    mean(all_mean_thetac_maxT_nmbs_thetai)
  mean_thetai_maxT_mean_thetac_nmbs_thetai <-
    mean(all_maxT_mean_thetac_maxT_nmbs_thetai)

  evppi <- mean_thetai_mean_thetac_maxT_nmb_thetai -
    mean_thetai_maxT_mean_thetac_nmbs_thetai
  result <- structure(list(
    max_mean_thetai_mean_thetac_nmb =
      max(mean_thetai_mean_thetac_nmb_thetai),
    mean_thetai_mean_thetac_nmb = mean_thetai_mean_thetac_nmb_thetai,
    mean_thetai_mean_thetac_maxT_nmb = mean_thetai_mean_thetac_maxT_nmb_thetai,
    mean_thetai_maxT_mean_thetac_nmbs =
      mean_thetai_maxT_mean_thetac_nmbs_thetai,
    Chilcott_evppi = evppi
  ))
  return(result)

}

######################################################################
#' Function to estimate the expected value of partial perfect information
#' This is the outer loop of the two stage Monte Carlo process
#' @param param_file all parameters required to run the model,provided
#' with name of parameter, distribution and parameters that define the
#' probability distribution
#' @param colnames_paramdistr column names where parameter distribution and
#' values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold threshold values of WTP
#' @param iterations number of iterations needed
#' @param comparator optional parameter if need ICER as the output
#' @return result result after evpi calculation
#' @source https://www.sciencedirect.com/science/article/pii/S1098301510605888
#' @export
estimate_evpi <- function(param_file, colnames_paramdistr,
                           list_markov, threshold, iterations,
                          comparator = NULL) {

  list_checks <- list(param_file, colnames_paramdistr,
                  threshold, iterations)

  results <- sapply(list_checks, packDAMipd::check_null_na)
  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")
  if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }
  names_strategy <- list()
  for (i in 1:nrow(list_markov)) {
    names_strategy <- append(names_strategy, list_markov[[i]]$name_strategy)
  }
  names_strategy <- unlist(names_strategy)


  cols_needed <- nrow(list_markov)
  all_nmbs_thetai <- matrix(0, nrow = iterations, ncol = cols_needed)
  all_max_nmb_iteration <- list()
  for (loop in 1:iterations) {
      this_result <- find_nmb_allparams(param_file, colnames_paramdistr,
                                         list_markov, threshold, comparator)
      nmbs_thetai <- (as.numeric(this_result$NMB))
      # store them
      all_nmbs_thetai[loop, 1:cols_needed] <- nmbs_thetai
      # step 4
      maxT_nmbs_thetai <- max(as.numeric(nmbs_thetai))
      # store
      all_max_nmb_iteration[length(all_max_nmb_iteration) + 1] <-
        maxT_nmbs_thetai

    }
    mean_nmb_thetai <- colSums(all_nmbs_thetai) / iterations
    mean_maxT_nmbs_thetai  <- mean(unlist(all_max_nmb_iteration))
    maxT_mean_nmbs_thetai <- max(mean_nmb_thetai)
    col_with_max_mean_nmbs <- which(mean_nmb_thetai == maxT_mean_nmbs_thetai)
    trt_with_max_mean_nmbs <- names_strategy[col_with_max_mean_nmbs]
    evpi <- mean_maxT_nmbs_thetai - maxT_mean_nmbs_thetai
    result <- structure(list(
      mean_maxT_nmbs = mean_maxT_nmbs_thetai,
      maxT_mean_nmbs = maxT_mean_nmbs_thetai,
      treatment_max_mean_nmbs = trt_with_max_mean_nmbs,
      evpi = evpi
    ))
    return(result)
}

######################################################################
#' Function to get NMB for a set of parameter
#' This is the inner loop of the two stage Monte Carlo process
#' @param param_file all parameters required to run the model,provided with
#' name of parameter, distribution and parameters that define the probability
#' distribution
#' @param colnames_paramdistr column names where parameter distribution and
#' values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold threshold values of WTP
#' @param comparator optional parameter if need ICER as the output
#' @return icer_nmb ICER and NMB using all sampled parameters
#' @keywords internal
#' @details
#' For the given set of parameter find the nmb from the packDAMipd function
#' calculate_icer_nmb
find_nmb_allparams <- function(param_file, colnames_paramdistr,
                               list_markov, threshold, comparator = NULL) {

  list_checks <- list(param_file, colnames_paramdistr,
                  threshold)
  results <- sapply(list_checks, packDAMipd::check_null_na)
  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")

  if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }
  new_list_markov <- list()
  for (i in 1:nrow(list_markov)) {
    this <- list_markov[i, ]
    # names of strategies extracted - control, intervention etc
    # parameters needed get it from param_matrix
    names_params_needed <- colnames(list_markov[i, ]$param_matrix)
    names_params_model <- names(list_markov[i, ]$list_param_values)
    params_passed <- list_markov[i, ]$list_param_values
    current_param_list <- get_all_parameter_list(names_params_needed,
                                                 names_params_model,
                                             params_passed, param_file,
                                             colnames_paramdistr)
    this_markov  <- packDAMipd::markov_model(current_strategy = this$strategy,
                                             cycles = this$cycles,
                                             initial_state = this$initial_state,
                                             discount = this$discount,
                                             parameter_values = current_param_list,
                                             half_cycle_correction =
                                               this$half_cycle_correction,
                                             state_cost_only_prevalent =
                                               this$state_cost_only_prevalent,
                                             state_util_only_prevalent =
                                               this$state_util_only_prevalent,
                                             method = this$method,
                                             startup_cost = this$startup_cost,
                                             startup_util = this$startup_util)
    new_list_markov[[length(new_list_markov) + 1]] <- this_markov
  }
  new_list_markov <- packDAMipd::combine_markov(new_list_markov)
  icer_nmb <- packDAMipd::calculate_icer_nmb(new_list_markov, threshold)
  return(icer_nmb)
}
######################################################################
#' Function to estimate the expected value of partial perfect information
#' @param names_params_needed names of needed parameters  from param matrix
#' returned
#' @param names_params_model names of parameters in the model
#' @param params_passed parameters passed while running the markov model
#' @param param_file all parameters required to run the model,provided
#' with name of parameter, distribution and parameters that define the
#' probability distribution
#' @param colnames_paramdistr col names where the parameter distribution is
#' defined
#' @return current_parameter_list set of all parameters
#' @keywords internal
#' @details
#' this function gets all the parameters except the parameter of interest
#' if they parameter is fixed, just read from file or a distribution,
#' then gets it from a distribution, or if it to be calculated, just
#' give it back as it is
get_all_parameter_list <- function(names_params_needed, names_params_model,
                               params_passed,
                               param_file, colnames_paramdistr) {
  list_checks <- list(names_params_needed, names_params_model,
                  params_passed, param_file, colnames_paramdistr)
  results <- sapply(list_checks, packDAMipd::check_null_na)
  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")


  param_distr_data <- utils::read.csv(param_file, row.names = NULL)
  param_colno <- IPDFileCheck::get_columnno_fornames(param_distr_data,
                                                     "parameter")
  if (length(param_colno) > 1) {
    stop("Parameter column should be only one")
  }
  param_colname <- colnames(param_distr_data)[param_colno]
  distr_col_no <- IPDFileCheck::get_columnno_fornames(param_distr_data,
                                                     "distribution")
  if (length(distr_col_no) > 1) {
    stop("Distribution column should be only one")
  }
  current_param_needed <- c()
  current_param_names <- c()
  for (j in 2:length(names_params_needed)) {
      # sample all - but identify those require calculations
      #if the selected parameters is needed to run the model
    if (names_params_needed[j] %in% names_params_model) {
        #if the params_passed for the current parameter considered
        # is a numerical value, that has been either read from a file
        # or directly assigned. It might not be calculated. so those can be
        # found from the parameter file
        res <- suppressWarnings(as.numeric(params_passed[names_params_needed[j]]))
        if (!is.na(res)) {
          # check if the parameter is defined using a distribution or actually
          # fixed (for death rates)
          row_cor_param <- param_distr_data[param_distr_data[param_colname] ==
                                             names_params_needed[j], ]
          if (nrow(row_cor_param) > 1) {
            stop("One row should correspond to the paramter in parameter file")
          }
          if (nrow(row_cor_param) == 0) {
            curr_parm_value <- params_passed[names_params_needed[j]]
          }else{
            distr_col_val_param <- row_cor_param[[distr_col_no]]
            # if the entry corresponding to distribution is null, na or
            # empty keep this param fixed
            if (packDAMipd::check_null_na(distr_col_val_param) < 0 |
                distr_col_val_param == "") {
              curr_parm_value <-
                packDAMipd::get_parameter_read(names_params_needed[j],
                                               param_file)
            }else{
              curr_parm_value <-
            packDAMipd::get_parameter_def_distribution(names_params_needed[j],
                                                           param_file,
                                                        colnames_paramdistr)
            }
          }
        }else{
          curr_parm_value <- params_passed[names_params_needed[j]]
        }
        curr_parm <- names_params_needed[j]
    }
    current_param_needed <- append(current_param_needed, curr_parm_value)
    current_param_names <- append(current_param_names, curr_parm)
  }
  current_param_list <- current_param_needed
  names(current_param_list) <- current_param_names
  return(current_param_list)
}
######################################################################
#' Function to estimate the expected value of partial perfect information
#' This is the outer loop of the two stage Monte Carlo process
#' @param param_interest the main parameter of interest
#' @param param_file all parameters required to run the model,provided with
#' name of the parameter, distribution and parameters that define the
#' probability distribution
#' @param colnames_paramdistr column names where parameter distribution
#' and values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold threshold values of WTP
#' @param outer_iterations number of iterations for outer loop
#' @param inner_iterations number of iterations for inner loop
#' @param comparator optional parameter if need ICER as the output
#' @return result evpi and evppi for a single param
#' @source https://www.sciencedirect.com/science/article/pii/S1098301510605888
#' @export
estimate_evpi_evppi_single_param <- function(param_interest, param_file,
                                            colnames_paramdistr, list_markov,
                                            threshold, outer_iterations,
                                            inner_iterations,
                                            comparator = NULL) {

  list_checks <- list(param_interest, param_file, colnames_paramdistr,
                  threshold, outer_iterations, inner_iterations)
  results <- sapply(list_checks, packDAMipd::check_null_na)

  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")

  if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }

  results_evpi <- estimate_evpi(param_file, colnames_paramdistr,
                                list_markov, threshold, outer_iterations,
                               comparator)
  results_evppi <- estimate_evppi(param_interest, param_file,
                                 colnames_paramdistr, list_markov,
                                 threshold, outer_iterations,
                                 inner_iterations, comparator)

  Chilcott_evppi <- results_evppi$Chilcott_evppi
  Brennan_evppi <- results_evppi$mean_thetai_maxT_mean_thetac_nmbs -
                  results_evpi$maxT_mean_nmbs
  Claxton_evppi <- results_evppi$mean_thetai_maxT_mean_thetac_nmbs -
                  results_evppi$max_mean_thetai_mean_thetac_nmb

  result <- structure(list(
    evpi = results_evpi,
    evppi_2level = results_evppi,
    Chilcott_evppi = Chilcott_evppi,
    Brennan_evppi = Brennan_evppi,
    Claxton_evppi = Claxton_evppi
  ))
  return(result)
}

######################################################################
#' Function to estimate the expected value of partial perfect information
#' This is the outer loop of the two stage Monte Carlo process
#' @param parameter_of_interest parameter of interest
#' @param param_file all parameters required to run the model,provided
#' with name of the parameter, distribution and parameters that define
#' the probability distribution
#' @param colnames_paramdistr column names where parameter distribution and
#' values are defined
#' @param list_markov list of markov models to estimate the NMB
#' @param threshold_values threshold values of WTP
#' @param outer_iterations number of iterations for outer loop
#' @param inner_iterations number of iterations for inner loop
#' @param comparator optional parameter if need ICER as the output
#' @return result result for evpi and evppi for different values of threshold
#' @source https://www.sciencedirect.com/science/article/pii/S1098301510605888
#' @export
estimate_evpi_evppi_diff_threshold <- function(parameter_of_interest,
                                              param_file,
                                              colnames_paramdistr,
                                              list_markov, threshold_values,
                                              outer_iterations = NULL,
                                              inner_iterations= NULL,
                                              comparator = NULL) {
  list_checks <- list(parameter_of_interest, param_file, colnames_paramdistr,
                  threshold_values)
  results <- sapply(list_checks, packDAMipd::check_null_na)

  if (sum(results) != 0)
    stop("Some of the parameters may by null or na - please check")

  if (is.null(list_markov)) {
    stop("Markov models can not be NULL")
  }
  if (is.null(outer_iterations))
    outer_iterations <- 500
  if (is.null(inner_iterations))
    inner_iterations <- 1000

  all_evpis <- as.data.frame(0, nrow = length(threshold_values), ncol = 4)
  num_trt_strategies <- nrow(list_markov)

  all_evppis <- as.data.frame(0, nrow = length(threshold_values),
                             ncol = num_trt_strategies + 3)
  named_evppis <- as.data.frame(0, nrow = length(threshold_values), ncol = 3)

  number_threshold <- length(threshold_values)
  for (jj in 1:number_threshold) {
    this_threshold <- threshold_values[jj]
    res_a_param <- estimate_evpi_evppi_single_param(parameter_of_interest,
                                                   param_file,
                                                   colnames_paramdistr,
                                                   list_markov,
                                                   this_threshold,
                                                   outer_iterations,
                                                   inner_iterations,
                                                   comparator)
    all_evpis[jj, 1] <- as.numeric(res_a_param$evpi[1])
    all_evpis[jj, 2] <- as.numeric(res_a_param$evpi[2])
    all_evpis[jj, 3] <- (res_a_param$evpi[3])
    all_evpis[jj, 4] <- as.numeric(res_a_param$evpi[4])
    all_evppis[jj, 1] <- as.numeric(res_a_param$evppi[1])
    names_strategy <- list()
    for (i in 1:num_trt_strategies) {
      all_evppis[jj, i + 1] <- as.numeric(res_a_param$evppi[[2]][i])
      names_strategy <- append(names_strategy,
                               list_markov[i, ]$strategy$name_strategy)
    }
    names_strategy <- unlist(names_strategy)
    all_evppis[jj, num_trt_strategies + 2] <-
      as.numeric(res_a_param$evppi[3])
    all_evppis[jj, num_trt_strategies + 3] <-
      as.numeric(res_a_param$evppi[4])

    named_evppis[jj, 1] <- as.numeric(res_a_param$Chilcott_evppi)
    named_evppis[jj, 2] <- as.numeric(res_a_param$Brennan_evppi)
    named_evppis[jj, 3] <- as.numeric(res_a_param$Claxton_evppi)

  }
  all_evpis[["thresholds"]] <- threshold_values
  all_evppis[["thresholds"]] <- threshold_values
  named_evppis[["thresholds"]] <- threshold_values
  colnames(all_evpis) <- c("Mean max NMB", "Max Mean NMB",
                           "Treatment with Max Mean NMB",
                           "EVPI", "Thresholds")
  temp_col_names <- list()
  for (i in 1:num_trt_strategies) {
    temp_col_names <- append(temp_col_names,
                             paste("Mean Thetai Mean Thetac NMB ",
                                   names_strategy[i]))
  }
  temp_col_names <- unlist(temp_col_names)
  colnames(all_evppis) <- c("Max Mean Thetai Mean Thetac NMB",
                            temp_col_names,
                           "Mean Thetai Mean Thetac Max NMB",
                           "Mean Thetai Max Mean Thetac NMB",
                           "Thresholds")
  colnames(named_evppis) <- c("Chilcott", "Brennan", "Claxton", "Thresholds")
  result <- list(Parameter = parameter_of_interest,
                EVPIs = all_evpis,
                EVPPIs = all_evppis,
                ThreeEVPPIs = named_evppis)

  return(result)
}
######################################################################
#' Function to estimate the expected value of partial perfect information
#' This is the outer loop of the two stage Monte Carlo process
#' @param result_evpi_evppi result from estimation of evpi and evppi
#' @return plot the plot
#' @export
plot_evpi_threshold <- function(result_evpi_evppi) {
  if (is.null(result_evpi_evppi)) {
    stop("result from EVPI /EVPPI calculation can not be NULL")
  }
  required_names <- c("Parameter", "EVPIs", "EVPPIs", "ThreeEVPPIs")
  current_names <- names(result_evpi_evppi)
  if (sum(current_names %in% required_names) != length(required_names)) {
    stop("Error - Given result do not contain all the required infomration")
  }
  EVPIs <- result_evpi_evppi$EVPIs
  EVPPIs <- result_evpi_evppi$ThreeEVPPIs

  name_file_plot <- paste0("EVPI for threshold values for ",
                           result_evpi_evppi$Parameter, ".pdf", sep = "")
  grDevices::pdf(name_file_plot)
  p <- ggplot2::ggplot(EVPIs, ggplot2::aes_(x = ~Thresholds,
                                                        y = ~EVPI)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Threshold") +
    ggplot2::labs(y = "EVPI") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  graphics::plot(p) # plot result
  grDevices::dev.off()

  name_file_plot <- paste0("EVPPI for threshold values for ",
                           result_evpi_evppi$Parameter, ".pdf", sep = "")
  grDevices::pdf(name_file_plot)

  EVPPIs_melted <- reshape2::melt(EVPPIs, id.var = "Thresholds")
  p <- ggplot2::ggplot(EVPPIs_melted, ggplot2::aes_(x = ~Thresholds,
                                            y = ~value, color = ~variable)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Threshold") +
    ggplot2::labs(y = "EVPPI") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  graphics::plot(p) # plot result
  grDevices::dev.off()
}
