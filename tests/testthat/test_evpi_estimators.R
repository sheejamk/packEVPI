param_file = system.file("extdata", "table_param.csv", package = "packEVPI")

well <- packDAMipd::health_state("well", cost = "cost_well_co", utility = 1)
disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_co", utility = "utility_dis_co")
dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
tm <- packDAMipd::populate_transition_matrix(3, tmat, c("tp_well_well_co", "tp_well_dis_co",
                                            "tp_well_dead", "tp_dis_dis_co",
                                            "tp_dis_dead", "tp_dead_dead"),
                                 colnames(tmat))
health_states <- packDAMipd::combine_state(well, disabled, dead)
this.strategy <- packDAMipd::strategy(tm, health_states, "control")
param_list = packDAMipd::define_parameters(
  tp_well_dis_co = packDAMipd::get_parameter_read("tp_well_dis_co", param_file),
  tp_well_dis_in =  packDAMipd::get_parameter_read("tp_well_dis_in", param_file),
  tp_well_dead =  packDAMipd::get_parameter_read("tp_well_dead", param_file),
  tp_dis_dead =  packDAMipd::get_parameter_read("tp_dis_dead", param_file),
  tp_dead_dead =  1,
  cost_well_co =  packDAMipd::get_parameter_read("cost_well_co", param_file),
  cost_well_in =  packDAMipd::get_parameter_read("cost_well_in", param_file),
  cost_dis_co =  packDAMipd::get_parameter_read("cost_dis_co", param_file),
  cost_dis_in =  packDAMipd::get_parameter_read("cost_dis_in", param_file),
  utility_dis_co =  packDAMipd::get_parameter_read("utility_dis_co", param_file),
  utility_dis_in =  packDAMipd::get_parameter_read("utility_dis_in", param_file),
  tp_well_well_co = "1-(tp_well_dis_co + tp_well_dead)",
  tp_well_well_in = "1-(tp_well_dis_in + tp_well_dead)",
  tp_dis_dis_co = "1-( tp_dis_dead)",
  tp_dis_dis_in = "1-( tp_dis_dead)"
)

this_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0), discount = c(0, 0),
                            method = "half cycle correction", param_list)

well <- packDAMipd::health_state("well", cost = "cost_well_in", utility = 1)
disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_in", utility = "utility_dis_in")
dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
tm <- packDAMipd::populate_transition_matrix(3, tmat, c("tp_well_well_in", "tp_well_dis_in",
                                            "tp_well_dead", "tp_dis_dis_in",
                                            "tp_dis_dead", "tp_dead_dead"),
                                 colnames(tmat))
health_states <- packDAMipd::combine_state(well, disabled, dead)
this.strategy <- packDAMipd::strategy(tm, health_states, "intervention")
sec_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0), discount = c(0, 0),
                           method = "half cycle correction",param_list)
list_markov <- packDAMipd::combine_markov(list(this_markov, sec_markov))


param_interest = "tp_well_dis_co"
colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")

parameter_of_interest =  "tp_well_dis_in"
threshold_values = c(5000,10000,15000,20000)

res = estimate_evpi_evppi_diff_threshold(parameter_of_interest, param_file,
                                         colnames_paramdistr,
                                         list_markov,threshold_values,
                                         outer_iterations = 2,
                                         inner_iterations = 3)

###############################################################################
context("testing getting parameter list")
test_that("testing getting parameter list", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  param_interest = "tp_well_dis_co"
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  value_param_interest = 0.2
  names_params_needed = colnames(list_markov[1,]$param_matrix)
  names_params_model = names(list_markov[1,]$list_param_values)
  params_passed = list_markov[1,]$list_param_values

  expect_error(get_parameter_list(NULL,value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))

  expect_error(get_parameter_list(param_interest,NULL,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  NA, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))

  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  names_params_needed, NA,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  names_params_needed, names_params_model,
                                  NA,
                                  param_file,colnames_paramdistr))
  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  NA,colnames_paramdistr))
  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,NA))
  param_file = system.file("extdata", "table_param_name_error.csv", package = "packEVPI")
  expect_error(get_parameter_list(param_interest,value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  param_file = system.file("extdata", "table_param_name_col_twice.csv", package = "packEVPI")
  expect_error(get_parameter_list(param_interest,value_param_interest,
                     names_params_needed, names_params_model,
                     params_passed,
                     param_file,colnames_paramdistr))
  param_file = system.file("extdata", "table_param_name_twice.csv", package = "packEVPI")
  param_interest = "tp_well_dis_co"
  expect_error(get_parameter_list(param_interest,value_param_interest,
                     names_params_needed, names_params_model,
                     params_passed,
                     param_file,colnames_paramdistr))
})
# ###############################################################################
context("testing finding nmb for a set of parameters")
test_that("testing finding nmb for a set of parameters", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  param_interest = "tp_well_dis_co"
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  value_param_interest = 0.2
  expect_error(find_nmb_setofparams(NULL,value_param_interest,
                                   param_file,colnames_paramdistr,
                                   list_markov,threshold))
  expect_error(find_nmb_setofparams(param_interest,NULL,
               param_file,colnames_paramdistr,
               list_markov,threshold))
  expect_error(find_nmb_setofparams(param_interest,value_param_interest,
                                    NULL,colnames_paramdistr,
                                    list_markov,threshold))
  expect_error(find_nmb_setofparams(param_interest,value_param_interest,
                                    param_file,NULL,
                                    list_markov,threshold))
  expect_error( find_nmb_setofparams(param_interest,value_param_interest,
                                     param_file,colnames_paramdistr,
                                     NULL,threshold))

  expect_error(find_nmb_setofparams(param_interest,value_param_interest,
                                   param_file,colnames_paramdistr,
                                   list_markov,NULL))

})
###############################################################################
context("testing estimating evppi")
test_that("testing estimating evppi", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  param_interest = "tp_well_dis_co"
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  value_param_interest = 0.2
  expect_error(estimate_evppi(param_interest, param_file,
                                    colnames_paramdistr,
                                    list_markov,threshold, outer_iterations,
                                    inner_iterations))
  expect_error(estimate_evppi(NULL, param_file,
                 colnames_paramdistr,
                 list_markov,threshold, outer_iterations,
                 inner_iterations))
  expect_error(estimate_evppi(param_interest, NULL,
                              colnames_paramdistr,
                              list_markov,threshold, outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              NULL,
                              list_markov,threshold, outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              NULL,threshold, outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              list_markov,NULL, outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                 colnames_paramdistr,
                 list_markov,threshold, NULL,
                 inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              list_markov,threshold, outer_iterations,
                              NULL))

  expect_error(estimate_evppi("cost_IT", param_file,
                              colnames_paramdistr,
                              list_markov,threshold, outer_iterations,
                              inner_iterations))

})
###############################################################################
context("testing estimating evpi")
test_that("testing estimating evpi", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  threshold  =  20000
  iterations = 1
  expect_error( estimate_evpi(NULL,
                colnames_paramdistr,
                list_markov,threshold, iterations,
  ))
  expect_error( estimate_evpi(param_file,NULL,
                              list_markov,threshold, iterations,
  ))
  expect_error( estimate_evpi(param_file,colnames_paramdistr,
                              NULL,threshold, iterations,
  ))
  expect_error( estimate_evpi(param_file,colnames_paramdistr,
                              list_markov,NULL, iterations,
  ))
  expect_error( estimate_evpi(param_file,colnames_paramdistr,
                              list_markov,threshold, NULL,
  ))
})
###############################################################################
context("testing finding nmb for all parameters")
test_that("testing finding nmb for all parameters", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  expect_error(find_nmb_allparams(NULL,colnames_paramdistr,
                                    list_markov,threshold))
  expect_error(find_nmb_allparams(param_file,NULL,
                                    list_markov,threshold))
  expect_error(find_nmb_allparams(param_file,colnames_paramdistr,
                                     NULL,threshold))
  expect_error(find_nmb_allparams(param_file,colnames_paramdistr,
                                    list_markov,NULL))
})
###############################################################################
context("testing getting all parameter list")
test_that("testing getting allparameter list", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  names_params_needed = colnames(list_markov[1,]$param_matrix)
  names_params_model = names(list_markov[1,]$list_param_values)
  params_passed = list_markov[1,]$list_param_values

  expect_error(get_all_parameter_list(NA, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))

  expect_error(get_all_parameter_list(names_params_needed, NA,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  NULL,
                                  param_file,colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  params_passed,
                                  NULL,colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,NULL))
  param_file = system.file("extdata", "table_param_name_error.csv", package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  param_file = system.file("extdata", "table_param_name_col_twice.csv", package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))
  param_file = system.file("extdata", "table_param_name_twice.csv", package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                  params_passed,
                                  param_file,colnames_paramdistr))
})
###############################################################################

context("testing evpi and evppi for single parameter")
test_that("testing evpi and evppi for single parameter", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  param_interest = "tp_well_dis_co"
  outer_iterations = 1
  inner_iterations = 2
  expect_error(estimate_evpi_evppi_single_param(NULL, param_file,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))

  expect_error(estimate_evpi_evppi_single_param(param_interest, NULL,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))

  expect_error(estimate_evpi_evppi_single_param(param_interest, param_file,
                                                NULL,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))
  expect_error(estimate_evpi_evppi_single_param(param_interest, param_file,
                                                colnames_paramdistr,
                                                NULL,threshold,
                                                outer_iterations, inner_iterations))
  expect_error(estimate_evpi_evppi_single_param(param_interest, param_file,
                                                colnames_paramdistr,
                                                list_markov,NULL,
                                                outer_iterations, inner_iterations))

  expect_error(estimate_evpi_evppi_single_param(param_interest, param_file,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                NULL, inner_iterations))

  expect_error(estimate_evpi_evppi_single_param(param_interest, param_file,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                outer_iterations, NULL))

})

###############################################################################

context("testing evpi and evppi for single parameter but for mulitple values of threshold")
test_that("testing evpi and evppi for single parameter but for mulitple values of threshold", {
  param_file = system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  = c("Param1_name", "Param1_value","Param2_name", "Param2_value")
  param_interest = "tp_well_dis_co"
  outer_iterations = 1
  inner_iterations = 2

  expect_error(estimate_evpi_evppi_diff_threshold(NULL, param_file,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))

  expect_error(estimate_evpi_evppi_diff_threshold(param_interest, NULL,
                                                colnames_paramdistr,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))

  expect_error(estimate_evpi_evppi_diff_threshold(param_interest, param_file,
                                                NULL,
                                                list_markov,threshold,
                                                outer_iterations, inner_iterations))
  expect_error(estimate_evpi_evppi_diff_threshold(param_interest, param_file,
                                                colnames_paramdistr,
                                                NULL,threshold,
                                                outer_iterations, inner_iterations))
  expect_error(estimate_evpi_evppi_diff_threshold(param_interest, param_file,
                                                colnames_paramdistr,
                                                list_markov,NULL,
                                                outer_iterations, inner_iterations))
})

