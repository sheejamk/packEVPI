
param_file <- system.file("extdata", "table_param.csv", package = "packEVPI")

well <- packDAMipd::health_state("well", cost = "cost_well_co", utility = 1)
disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_co",
                                     utility = "utility_dis_co")
dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
tm <- packDAMipd::populate_transition_matrix(3, tmat,
                                             c("tp_well_well_co",
                                              "tp_well_dis_co",
                                              "tp_well_dead", "tp_dis_dis_co",
                                            "tp_dis_dead", "tp_dead_dead"),
                                             colnames(tmat))
health_states <- packDAMipd::combine_state(well, disabled, dead)
this.strategy <- packDAMipd::strategy(tm, health_states, "control")
param_list <- packDAMipd::define_parameters(
  tp_well_dis_co = packDAMipd::get_parameter_read("tp_well_dis_co",
                                                  param_file),
  tp_well_dis_in =  packDAMipd::get_parameter_read("tp_well_dis_in",
                                                   param_file),
  tp_well_dead =  packDAMipd::get_parameter_read("tp_well_dead", param_file),
  tp_dis_dead =  packDAMipd::get_parameter_read("tp_dis_dead", param_file),
  tp_dead_dead =  1,
  cost_well_co =  packDAMipd::get_parameter_read("cost_well_co", param_file),
  cost_well_in =  packDAMipd::get_parameter_read("cost_well_in", param_file),
  cost_dis_co =  packDAMipd::get_parameter_read("cost_dis_co", param_file),
  cost_dis_in =  packDAMipd::get_parameter_read("cost_dis_in", param_file),
  utility_dis_co =  packDAMipd::get_parameter_read("utility_dis_co",
                                                   param_file),
  utility_dis_in =  packDAMipd::get_parameter_read("utility_dis_in",
                                                   param_file),
  tp_well_well_co = "1-(tp_well_dis_co + tp_well_dead)",
  tp_well_well_in = "1-(tp_well_dis_in + tp_well_dead)",
  tp_dis_dis_co = "1-( tp_dis_dead)",
  tp_dis_dis_in = "1-( tp_dis_dead)"
)

this_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0),
                                        discount = c(0, 0),
                                        method = "half cycle correction",
                                        param_list)

well <- packDAMipd::health_state("well", cost = "cost_well_in", utility = 1)
disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_in",
                                     utility = "utility_dis_in")
dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
tm <- packDAMipd::populate_transition_matrix(3, tmat,
                                             c("tp_well_well_in",
                                              "tp_well_dis_in",
                                            "tp_well_dead", "tp_dis_dis_in",
                                              "tp_dis_dead", "tp_dead_dead"),
                                             colnames(tmat))
health_states <- packDAMipd::combine_state(well, disabled, dead)
this.strategy <- packDAMipd::strategy(tm, health_states, "intervention")
sec_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0),
                                       discount = c(0, 0),
                                method = "half cycle correction", param_list)
list_markov <- packDAMipd::combine_markov(list(this_markov, sec_markov))


A <- packDAMipd::health_state("A", cost = 100, utility = 1)
B <- packDAMipd::health_state("B", cost = 200, utility = 0.8)
C <- packDAMipd::health_state("C", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("A", "B", "C")
tm <- packDAMipd::populate_transition_matrix(3, tmat, c(0.5, 0.4, 0.1,
                                                        0.6, 0.4, 1),
                                             colnames(tmat))
health_states <- packDAMipd::combine_state(A, B, C)
control_strategy <- packDAMipd::strategy(tm, health_states, "control")
control_markov <- packDAMipd::markov_model(control_strategy, 10,
                                           c(1000, 0, 0),
                                           discount = c(0, 0))

A <- packDAMipd::health_state("A", cost = 200, utility = 1)
B <- packDAMipd::health_state("B", cost = 400, utility = 1)
C <- packDAMipd::health_state("C", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("A", "B", "C")
tm <- packDAMipd::populate_transition_matrix(3, tmat, c(0.5, 0.4, 0.1,
                                                        0.6, 0.4, 1),
                                             colnames(tmat))
health_states <- packDAMipd::combine_state(A, B, C)
inter_strategy <- packDAMipd::strategy(tm, health_states, "inter")
inter_markov <- packDAMipd::markov_model(inter_strategy, 10, c(1000, 0, 0),
                                         discount = c(0, 0))
list_markov2 <- packDAMipd::combine_markov(control_markov, inter_markov)

A <- packDAMipd::health_state("A", cost = "cost_a", utility = 1)
B <- packDAMipd::health_state("B", cost = 400, utility = 1)
C <- packDAMipd::health_state("C", cost = 0, utility = 0)
tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
colnames(tmat) <- rownames(tmat) <- c("A", "B", "C")
tm <- packDAMipd::populate_transition_matrix(3, tmat, c(0.5, 0.4, 0.1,
                                                        0.6, 0.4, 1),
                                             colnames(tmat))
param2_list <- packDAMipd::define_parameters(cost_a = 100)
health_states <- packDAMipd::combine_state(A, B, C)
inter2_strategy <- packDAMipd::strategy(tm, health_states, "inter")
inter2_markov <- packDAMipd::markov_model(inter2_strategy, 24, c(1000, 0, 0),
                                          discount = c(0, 0), param2_list)

list_markov3 <-  packDAMipd::combine_markov(this_markov, inter2_markov)
###############################################################################
context("testing getting parameter list")
test_that("testing getting parameter list", {
  param_file <- system.file("extdata", "table_param.csv",
                            package = "packEVPI")
  param_interest <- "tp_well_dis_co"
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  value_param_interest <- 0.2
  names_params_needed <- colnames(list_markov[1, ]$param_matrix)
  names_params_model <- names(list_markov[1, ]$list_param_values)
  params_passed <- list_markov[1, ]$list_param_values

  parameters <- get_parameter_list(param_interest, value_param_interest,
                                   names_params_needed, names_params_model,
                                   params_passed,
                                   param_file, colnames_paramdistr)
  expect_equal(parameters$tp_well_dis_co, 0.2)
  expect_error(get_parameter_list(NULL, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))

  expect_error(get_parameter_list(param_interest, NULL,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  NA, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))

  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, NA,
                                  params_passed,
                                  param_file, colnames_paramdistr))
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  NA,
                                  param_file, colnames_paramdistr))
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  NA, colnames_paramdistr))
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, NA))
  param_file <- system.file("extdata", "table_param_name_error.csv",
                            package = "packEVPI")
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_name_col_twice.csv",
                            package = "packEVPI")
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_name_twice.csv",
                            package = "packEVPI")
  param_interest <- "tp_well_dis_co"
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_distribution_twice.csv",
                            package = "packEVPI")
  param_interest <- "tp_well_dis_co"
  expect_error(get_parameter_list(param_interest, value_param_interest,
                                  names_params_needed, names_params_model,
                                  params_passed,
                                  param_file, colnames_paramdistr))

})
#############################################################################
context("testing finding nmb for a set of parameters")
test_that("testing finding nmb for a set of parameters", {
  param_file <- system.file("extdata", "table_param.csv",
                            package = "packEVPI")
  param_interest <- "tp_well_dis_co"
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  value_param_interest <- 0.2
  threshold <- 2000
  res <- find_nmb_setofparams(param_interest, value_param_interest,
                              param_file, colnames_paramdistr,
                              list_markov, threshold)

  expect_equal(res$Strategy[1], "control")

  expect_error(find_nmb_setofparams(NULL, value_param_interest,
                                    param_file, colnames_paramdistr,
                                    list_markov, threshold))
  expect_error(find_nmb_setofparams(param_interest, NULL,
                                    param_file, colnames_paramdistr,
                                    list_markov, threshold))
  expect_error(find_nmb_setofparams(param_interest, value_param_interest,
                                    NULL, colnames_paramdistr,
                                    list_markov, threshold))
  expect_error(find_nmb_setofparams(param_interest, value_param_interest,
                                    param_file, NULL,
                                    list_markov, threshold))
  expect_error(find_nmb_setofparams(param_interest, value_param_interest,
                                    param_file, colnames_paramdistr,
                                    NULL, threshold))

  expect_error(find_nmb_setofparams(param_interest, value_param_interest,
                                    param_file, colnames_paramdistr,
                                    list_markov, NULL))
  # only cycle as the parameter - error
  expect_error(find_nmb_setofparams("cost", value_param_interest,
                                    param_file, colnames_paramdistr,
                                    list_markov2, threshold))
  find_nmb_setofparams("cost_a", 200,
                       param_file, colnames_paramdistr,
                       list_markov3, threshold)
})
###############################################################################
context("testing estimating evppi")
test_that("testing estimating evppi", {
  param_file <- system.file("extdata", "table_param.csv",
                            package = "packEVPI")
  param_interest <- "tp_well_dis_co"
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  value_param_interest <- 0.2
  inner_iterations <- 5
  outer_iterations <- 2
  threshold <- 2000
  res <- estimate_evppi(param_interest, param_file,
                        colnames_paramdistr,
                        list_markov, threshold,
                        outer_iterations,
                        inner_iterations)
  expect_error(estimate_evppi(NULL, param_file,
                              colnames_paramdistr,
                              list_markov, threshold, outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, NULL,
                              colnames_paramdistr,
                              list_markov, threshold,
                              outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              NULL,
                              list_markov, threshold,
                              outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              NULL, threshold,
                              outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              list_markov, NULL,
                              outer_iterations,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              list_markov, threshold, NULL,
                              inner_iterations))
  expect_error(estimate_evppi(param_interest, param_file,
                              colnames_paramdistr,
                              list_markov, threshold,
                              outer_iterations,
                              NULL))

  expect_error(estimate_evppi("cost_IT", param_file,
                              colnames_paramdistr,
                              list_markov, threshold,
                              outer_iterations,
                              inner_iterations))
})
###############################################################################
context("testing estimating evpi")
test_that("testing estimating evpi", {
  param_file <- system.file("extdata", "table_param.csv",
                            package = "packEVPI")
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  threshold  <-  20000
  iterations <- 1
  estimate_evpi(param_file,
                colnames_paramdistr,
                list_markov, threshold, iterations)
  expect_error(estimate_evpi(NULL,
                             colnames_paramdistr,
                             list_markov, threshold, iterations))
  expect_error(estimate_evpi(param_file, NULL,
                             list_markov, threshold, iterations))
  expect_error(estimate_evpi(param_file, colnames_paramdistr,
                             NULL, threshold, iterations))
  expect_error(estimate_evpi(param_file, colnames_paramdistr,
                             list_markov, NULL, iterations))
  expect_error(estimate_evpi(param_file, colnames_paramdistr,
                             list_markov, threshold, NULL))
})
###############################################################################
context("testing finding nmb for all parameters")
test_that("testing finding nmb for all parameters", {
  param_file <- system.file("extdata", "table_param.csv",
                            package = "packEVPI")
  threshold <- 2000
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  nmbs <- find_nmb_allparams(param_file, colnames_paramdistr,
                             list_markov, threshold)
  expect_error(find_nmb_allparams(NULL, colnames_paramdistr,
                                  list_markov, threshold))
  expect_error(find_nmb_allparams(param_file, NULL,
                                  list_markov, threshold))
  expect_error(find_nmb_allparams(param_file, colnames_paramdistr,
                                  NULL, threshold))
  expect_error(find_nmb_allparams(param_file, colnames_paramdistr,
                                  list_markov, NULL))
})
###############################################################################
context("testing getting all parameter list")
test_that("testing getting allparameter list", {
  param_file <- system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  names_params_needed <- colnames(list_markov[1, ]$param_matrix)
  names_params_model <- names(list_markov[1, ]$list_param_values)
  params_passed <- list_markov[1, ]$list_param_values
  parameters <- get_all_parameter_list(names_params_needed, names_params_model,
                                       params_passed,
                                       param_file, colnames_paramdistr)

  expect_error(get_all_parameter_list(NA, names_params_model,
                                      params_passed,
                                      param_file, colnames_paramdistr))

  expect_error(get_all_parameter_list(names_params_needed, NA,
                                      params_passed,
                                      param_file, colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      NULL,
                                      param_file, colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      NULL, colnames_paramdistr))
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      param_file, NULL))
  param_file <- system.file("extdata", "table_param_name_error.csv",
                            package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_name_col_twice.csv",
                            package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_name_twice.csv",
                            package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      param_file, colnames_paramdistr))
  param_file <- system.file("extdata", "table_param_distribution_twice.csv",
                            package = "packEVPI")
  expect_error(get_all_parameter_list(names_params_needed, names_params_model,
                                      params_passed,
                                      param_file, colnames_paramdistr))
})
###############################################################################

context("testing evpi and evppi for single parameter")
test_that("testing evpi and evppi for single parameter", {
  param_file <- system.file("extdata", "table_param.csv", package = "packEVPI")
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  param_interest <- "tp_well_dis_co"
  outer_iterations <- 1
  inner_iterations <- 2
  threshold <- 20000
  res <- estimate_evpi_evppi(param_interest, param_file,
                             colnames_paramdistr,
                             list_markov, threshold,
                             outer_iterations,
                             inner_iterations)
  expect_error(estimate_evpi_evppi(NULL, param_file,
                                   colnames_paramdistr,
                                   list_markov, threshold,
                                   outer_iterations,
                                   inner_iterations))

  expect_error(estimate_evpi_evppi(param_interest, NULL,
                                   colnames_paramdistr,
                                   list_markov, threshold,
                                   outer_iterations,
                                   inner_iterations))

  expect_error(estimate_evpi_evppi(param_interest, param_file,
                                   NULL,
                                   list_markov, threshold,
                                   outer_iterations,
                                   inner_iterations))
  expect_error(estimate_evpi_evppi(param_interest, param_file,
                                   colnames_paramdistr,
                                   NULL, threshold,
                                   outer_iterations,
                                   inner_iterations))
  expect_error(estimate_evpi_evppi(param_interest, param_file,
                                   colnames_paramdistr,
                                   list_markov, NULL,
                                   outer_iterations,
                                   inner_iterations))

  expect_error(estimate_evpi_evppi(param_interest, param_file,
                                   colnames_paramdistr,
                                   list_markov, threshold,
                                   NULL, inner_iterations))

  expect_error(estimate_evpi_evppi(param_interest, param_file,
                                   colnames_paramdistr,
                                   list_markov, threshold,
                                   outer_iterations, NULL))

})

###############################################################################

context("testing evpi and evppi for single parameter but for mulitple
        values of threshold")
test_that("testing evpi and evppi for single parameter but for mulitple
          values of threshold", {
            param_file <- system.file("extdata", "table_param.csv", package = "packEVPI")
            colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                                      "Param2_value")
            param_interest <- "tp_well_dis_co"
            outer_iterations <- 1
            inner_iterations <- 2
            threshold <- 20000
            res <- estimate_evpi_evppi_diff_threshold(param_interest,
                                                      param_file,
                                                      colnames_paramdistr,
                                                      list_markov, threshold,
                                                      outer_iterations,
                                                      inner_iterations)

            expect_error(estimate_evpi_evppi_diff_threshold(NULL, param_file,
                                                            colnames_paramdistr,
                                                            list_markov,
                                                            threshold,
                                                            outer_iterations,
                                                            inner_iterations))

            expect_error(estimate_evpi_evppi_diff_threshold(param_interest,
                                                            NULL,
                                                            colnames_paramdistr,
                                                            list_markov,
                                                            threshold,
                                                            outer_iterations,
                                                            inner_iterations))

            expect_error(estimate_evpi_evppi_diff_threshold(param_interest,
                                                            param_file,
                                                            NULL,
                                                            list_markov,
                                                            threshold,
                                                            outer_iterations,
                                                            inner_iterations))
            expect_error(estimate_evpi_evppi_diff_threshold(param_interest,
                                                            param_file,
                                                          colnames_paramdistr,
                                                            NULL, threshold,
                                                            outer_iterations,
                                                            inner_iterations))
            expect_error(estimate_evpi_evppi_diff_threshold(param_interest,
                                                            param_file,
                                                          colnames_paramdistr,
                                                            list_markov, NULL,
                                                            outer_iterations,
                                                            inner_iterations))
          })
###############################################################################

context("plotting evpi threshold")
test_that("plotting evpi threshold", {
  param_file <- system.file("extdata", "table_param.csv", package = "packEVPI")

  well <- packDAMipd::health_state("well", cost = "cost_well_co", utility = 1)
  disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_co",
                                       utility = "utility_dis_co")
  dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
  tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
  colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
  tm <- packDAMipd::populate_transition_matrix(3, tmat,
                                               c("tp_well_well_co",
                                                 "tp_well_dis_co",
                                                 "tp_well_dead",
                                                 "tp_dis_dis_co",
                                                 "tp_dis_dead",
                                                 "tp_dead_dead"),
                                               colnames(tmat))
  health_states <- packDAMipd::combine_state(well, disabled, dead)
  this.strategy <- packDAMipd::strategy(tm, health_states, "control")
  param_list <- packDAMipd::define_parameters(
    tp_well_dis_co = packDAMipd::get_parameter_read("tp_well_dis_co",
                                                    param_file),
    tp_well_dis_in =  packDAMipd::get_parameter_read("tp_well_dis_in",
                                                     param_file),
    tp_well_dead =  packDAMipd::get_parameter_read("tp_well_dead", param_file),
    tp_dis_dead =  packDAMipd::get_parameter_read("tp_dis_dead", param_file),
    tp_dead_dead =  1,
    cost_well_co =  packDAMipd::get_parameter_read("cost_well_co", param_file),
    cost_well_in =  packDAMipd::get_parameter_read("cost_well_in", param_file),
    cost_dis_co =  packDAMipd::get_parameter_read("cost_dis_co", param_file),
    cost_dis_in =  packDAMipd::get_parameter_read("cost_dis_in", param_file),
    utility_dis_co =  packDAMipd::get_parameter_read("utility_dis_co",
                                                     param_file),
    utility_dis_in =  packDAMipd::get_parameter_read("utility_dis_in",
                                                     param_file),
    tp_well_well_co = "1-(tp_well_dis_co + tp_well_dead)",
    tp_well_well_in = "1-(tp_well_dis_in + tp_well_dead)",
    tp_dis_dis_co = "1-( tp_dis_dead)",
    tp_dis_dis_in = "1-( tp_dis_dead)"
  )

  this_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0),
                                          discount = c(0, 0),
                                          method = "half cycle correction",
                                          param_list)

  well <- packDAMipd::health_state("well", cost = "cost_well_in", utility = 1)
  disabled <- packDAMipd::health_state("disabled", cost = "cost_dis_in",
                                       utility = "utility_dis_in")
  dead <- packDAMipd::health_state("dead", cost = 0, utility = 0)
  tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
  colnames(tmat) <- rownames(tmat) <- c("well", "disabled", "dead")
  tm <- packDAMipd::populate_transition_matrix(3, tmat, c("tp_well_well_in",
                                                          "tp_well_dis_in",
                                                          "tp_well_dead",
                                                          "tp_dis_dis_in",
                                                          "tp_dis_dead",
                                                          "tp_dead_dead"),
                                               colnames(tmat))
  health_states <- packDAMipd::combine_state(well, disabled, dead)
  this.strategy <- packDAMipd::strategy(tm, health_states, "intervention")
  sec_markov <- packDAMipd::markov_model(this.strategy, 24, c(1000, 0, 0),
                                         discount = c(0, 0),
                                         method = "half cycle correction",
                                         param_list)
  list_markov <- packDAMipd::combine_markov(list(this_markov, sec_markov))

  param_interest <- "tp_well_dis_co"
  colnames_paramdistr  <- c("Param1_name", "Param1_value", "Param2_name",
                            "Param2_value")
  threshold_values <- c(5000, 10000, 150000, 20000)
  res <- estimate_evpi_evppi_diff_threshold(param_interest, param_file,
                                            colnames_paramdistr,
                                            list_markov, threshold_values,
                                            outer_iterations = 3,
                                            inner_iterations = 5)
  plot_evpi_threshold(res)
  expect_error(plot_evpi_threshold(NULL))
  error_names <- c("one", "two", "three", "four")
  error_res <- res
  names(error_res) <- error_names
  expect_error(plot_evpi_threshold(error_res))
})

###############################################################################
context("testing evpi and evppi for subset of parameters but for mulitple
        values of threshold")
test_that("testing evpi and evppi for subset of parameters  but for mulitple
          values of threshold", {

            H <- packDAMipd::health_state("H", cost = "c_H ", utility = "u_H")
            S <- packDAMipd::health_state("S", cost = "c_S", utility = "u_S")
            D <- packDAMipd::health_state("D", cost = "c_D ", utility = "u_D")

            tmat <- rbind(c(1, 2, 3), c(NA, 4, 5), c(NA, NA, 6))
            colnames(tmat) <- rownames(tmat) <- c("H", "S", "D")

            tm <- packDAMipd::populate_transition_matrix(3, tmat,
                                                         c("p_HH", "p_HS",
                                                          "p_HD", "p_SS",
                                                          "p_SD", "p_DD"),
                                                         colnames(tmat))

            param_file <- system.file("extdata", "table_param_2.csv",
                                      package = "packEVPI")
            param_list <-
              packDAMipd::define_parameters(
                p_HD = packDAMipd::get_parameter_read("p_HD", param_file),
                p_HS  = packDAMipd::get_parameter_def_distribution("p_HS",
                                                                   param_file,
              c("Param1_name", "Param1_value", "Param2_name", "Param2_value")),
                hr_S   = packDAMipd::get_parameter_read("hr_S", param_file),
                p_SD  = "1 - exp(log(1 - p_HD) * hr_S)",
                p_HH = "1 - (p_HS + p_HD)",
                p_SS = "1 - (p_SD)",
                p_DD = 1,
                c_H   = packDAMipd::get_parameter_read("c_H", param_file),
                c_S  = packDAMipd::get_parameter_read("c_S", param_file),
                c_D   = packDAMipd::get_parameter_read("c_D", param_file),
                c_Trt = packDAMipd::get_parameter_read("c_Trt", param_file),
                u_H   = packDAMipd::get_parameter_read("u_H", param_file),
                u_S  = packDAMipd::get_parameter_read("u_S", param_file),
                u_D   = packDAMipd::get_parameter_read("u_D", param_file),
                u_Trt = packDAMipd::get_parameter_read("u_Trt", param_file),
                du_HS = packDAMipd::get_parameter_read("du_HS", param_file),
                ic_HS = packDAMipd::get_parameter_read("ic_HS", param_file),
                ic_D   = packDAMipd::get_parameter_read("ic_D", param_file))
            health_states <- packDAMipd::combine_state(H, S, D)
            uc_strategy <- packDAMipd::strategy(tm, health_states, "Usual care")
            uc_markov <- packDAMipd::markov_model(uc_strategy, 85, c(1, 0, 0),
                                                  c(0, 0),
                                                  param_list, TRUE,
                                          method = "half cycle correction")
            packDAMipd::plot_model(uc_markov)

            H <- packDAMipd::health_state("H", cost = "c_H ", utility = "u_H")
            S <- packDAMipd::health_state("S1", cost = "c_S + c_Trt ",
                                          utility = "u_Trt")
            D <- packDAMipd::health_state("D", cost = "c_D ", utility = "u_D")
            health_states <- packDAMipd::combine_state(H, S, D)
            trt_strategy <- packDAMipd::strategy(tm, health_states,
                                                 "New treatment")
            trt_markov <-
              packDAMipd::markov_model(current_strategy = trt_strategy,
                                                   cycles = 85,
                                                   initial_state = c(1, 0, 0),
                                                   discount = c(0, 0),
                                                   parameter_values =
                                                     param_list, TRUE, method =
                                                     "half cycle correction")

            list_markov <- packDAMipd::combine_markov(uc_markov, trt_markov)


            colnames_paramdistr  <- c("Param1_name", "Param1_value",
                                      "Param2_name",
                                      "Param2_value")

            parameters_of_interest <-  c("p_HS", "c_Trt")
            threshold_values <- c(10000)
            res <- `estimate_evpi_evppi_diff_threshold`(parameters_of_interest,
                                                        param_file,
                                                        colnames_paramdistr,
                                                        list_markov,
                                                        threshold_values,
                                                        outer_iterations = 1,
                                                        inner_iterations = 1)
            expect_equal(res$Parameter, c("p_HS", "c_Trt"))

          })
