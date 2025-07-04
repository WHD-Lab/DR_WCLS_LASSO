#' DR_WCLS_LASSO
#'
#' A function for basic LASSO DR-WCLS
#'
#' @param data raw data without pseudo-outcome, ptSt
#' @param fold # of folds hope to split when generating pesudo outcome
#' @param ID the name of column where participants' ID are stored
#' @param time the name of column where time in study are stored
#' @param Ht a vector that contains column names of control variables
#' @param St a vector that contains column names of moderator variables; St should be a subset of Ht
#' @param At column names of treatment (At)
#' @param prob column names of \eqn{p_t(A_t = 1|H_t)}, the experiment design treatment probability
#' @param outcome column names of outcome variable
#' @param method_pesu the machines learning method used when generate estimates of the nuisance parameters,
#' and those values will be used to calculate the pseudo outcome. The available machine learning algorithms include
#' cross validation LASSO ("CVLASSO"), Random Forest ("RandomForest"), and gradient boosting ("GradientBoosting").
#' @param lam the value of penalty term of randomized LASSO. If it is not provided, the default value will be used
#' @param noise_scale Scale of Gaussian noise added to objective. Default is \eqn{\sqrt{\frac{(1 - \rho)}{\rho}*number\ of\ rows}*sd(y)} where \eqn{\rho} is the split rate.
#' The random noises, \eqn{\omega}, are iid drawn from normal distribution with standard deviation noise_scale
#' @param splitrat this value is corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
#' This value will be used only when user doesn't provide the noise_scale.
#' @param virtualenv_path Python virtual environment path
#' @param beta the true coefficient value (if simulation is conducted)
#' @param level the CI significant level
#' @param core_num the number of cores will be used for parallel calculation
#' @param CI_algorithm when calculate CI can choice using for loop, built-in parallel function, and doParallel function; "lapply", "parallel", "doParallel"
#' @param max_iterate: the max iteration number when searching CI
#' @param max_tol the max error tolerance when calculate pivot value. The default value is \eqn{10^{-3}} i.e. the target pivot value is 5% for lower bound when calculating
#' 90% confidence interval, and the provided results is within 4.999% to 5.001%.
#' @param varSelect_program the user can decide using which program to do variable selection. If it is "Python", a valid virtual environment path must be provided, i.e.
#' virtualenv_path can't be NULL. If it is "R", no virtual environment path is required.
#'
#' @return A table with the selected variables is returned. The returned table contains GEE estimate,
#' the post selection true value (if simulation is conducted; Otherwise, NA is provided.)
#' the p value,
#' the confidence interval,
#' the true corresponding pivot value for the lower bound and the upper bound.
#' @examples
#'
#' UI_return = DR_WCLS_LASSO(data = data, fold = 5, ID = "id", time = "decision_point",
#' Ht = Ht, St = St, At = "action", prob = "prob", outcome = "outcome", method_pesu = "CVLASSO",
#' lam = NULL, noise_scale = NULL, splitrat = 0.8,
#' virtualenv_path = "path to selective-inference folder/env3",
#' beta =  c(-0.2, 0.8, 0.3, 0.7, 0.3, rep(0, 21)), level = 0.9, core_num = 3,
#' CI_algorithm = "parallel", varSelect_program = "R")
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import devtools
#' @import zoo
#'
#' @export

DR_WCLS_LASSO = function(data, fold, ID, time, Ht, St, At, prob, outcome, method_pesu,
                                 lam = NULL, noise_scale = NULL, splitrat = 0.8, virtualenv_path,
                                 beta = NULL, level = 0.9, core_num = NULL, CI_algorithm = "lapply",
                         max_iterate = 10^{6}, max_tol = 10^{-3}, varSelect_program = "Python"){
  # data: raw data without pesudo-outcome, ptSt.
  # fold: # of folds hope to split when generating pesudo outcome
  # ID: the name of column where participants' ID are stored
  # time: the name of column where time in study are stored
  # Ht: a vector that contains column names of control variables
  # St: a vector that contains column names of moderator variables; St should be a subset of Ht
  # At: column names of treatment (At)
  # outcome: column names of outcome variable
  # method_pesu: the machines learning method used when generate estimates of the nuisance parameters, and those values will be used to calculate the
  #             pesudo outcome
  # lam: the value of penalty term of randomized LASSO. If it is not provided, the default value will be used
  # noise_scale: Scale of Gaussian noise added to objective. Default is sqrt((1 - splitrat)/splitrat*NT)*sd(y).
  #             The random noises, omega, are iid drawn from normal distribution with standard deviation noise_scale
  # splitrat: the corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
  #           This value will be used only when user doesn't provide the noise_scale.
  # virtualenv_path: Python virtual environment path
  # beta: the true coefficient value (if simulation is conducted)
  # level: the CI significant level
  # core_num: the number of cores will be used for parallel calculation
  # CI_algorithm: when calculate CI can choice using for loop, built-in parallel function, and doParallel function; "lapply", "parallel", "doParallel"
  # max_iterate: the max iteration number when searching CI
  # max_tol: the max error tolerance when calculate pivot value
  # varSelect_program: the user can decide using which program to do variable selection. If it is "Python", a valid virtual environment path must be provided, i.e.
  # virtualenv_path can't be NULL. If it is "R", no virtual environment path is required.

  require(devtools)
  if(method_pesu == "CVLASSO") {
    ps = pseudo_outcome_generator_CVlasso(fold, ID, data, Ht, St, At, prob=prob, outcome, core_num)
  }

  if(method_pesu == "RandomForest") {
    ps = pseudo_outcome_generator_rf_v2(fold, ID, data, Ht, St, At, prob=prob, outcome, core_num)
  }

  if(method_pesu == "GradientBoosting") {
    ps = pseudo_outcome_generator_gbm(fold, ID, data, Ht, St, At, prob=prob, outcome, core_num)
  }

  my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))

  if(is.null(lam) & is.null(noise_scale)) {
    if(varSelect_program == "Python") {select = variable_selection_PY_penal_int(ps, ID, my_formula, splitrat=splitrat, virtualenv_path= virtualenv_path, beta = beta)}

    if(varSelect_program == "R") {select = FISTA_backtracking(ps, ID, my_formula, splitrat=splitrat, beta = beta)}
  }

  if(!is.null(lam) & is.null(noise_scale)) {
    if(varSelect_program == "Python") {select = variable_selection_PY_penal_int(ps, ID, my_formula, lam = lam, splitrat = splitrat,
                                                                                virtualenv_path= virtualenv_path, beta = beta)}

    if(varSelect_program == "R") {select = FISTA_backtracking(ps, ID, my_formula, lam = lam, splitrat = splitrat, beta = beta)}

  }

  if(is.null(lam) & !is.null(noise_scale)) {
    if(varSelect_program == "Python") {select = variable_selection_PY_penal_int(ps, ID, my_formula, noise_scale = noise_scale, splitrat = splitrat,
                                                                                virtualenv_path= virtualenv_path, beta = beta)}
    if(varSelect_program == "R") {
      select = FISTA_backtracking(ps, ID, my_formula, noise_scale = noise_scale, splitrat = splitrat, beta = beta)
    }

  }

  if(!is.null(lam) & !is.null(noise_scale)) {
    if(varSelect_program == "Python") {select = variable_selection_PY_penal_int(ps, ID, my_formula, lam = lam, noise_scale = noise_scale,
                                             splitrat = splitrat, virtualenv_path= virtualenv_path, beta = beta)}

    if(varSelect_program == "R") {select = FISTA_backtracking(ps, ID, my_formula, lam = lam, noise_scale = noise_scale,
                                                              splitrat = splitrat, beta = beta)}

  }

  # print selection
  print(paste("select predictors:", select$E))
  print(length(setdiff(select$E[select$E != "(Intercept)"], St)) == 0)
  print(setdiff(select$E[select$E != "(Intercept)"], St))

  if(length(setdiff(select$E[select$E != "(Intercept)"], St)) == 0) {
    AsyNormbeta_shared = joint_dist_selectAll(E = select$E, pes_outcome = "yDR", data = ps, id = ID, time = time)
    PQR_shared = PQR_shared_selectALL(AsyNormbeta_shared, select)

    CI_per_select_var = function(ej) {
      AsyNormbeta_ej = joint_dist_Penal_Int_ej(AsyNormbeta_shared, ej)
      PQR_ej = PQR_ej_selectALL(PQR_shared, AsyNormbeta_ej, AsyNormbeta_shared)
      condition = conditional_dist_selectALL(PQR_shared, PQR_ej, AsyNormbeta_shared, AsyNormbeta_ej, select)
      CI = pivot_split_update_sellectALL(PQR_shared, PQR_ej, condition, AsyNormbeta_shared, AsyNormbeta_ej, select, level = level, pes_outcome = "yDR", data = ps,
                              id = ID, time = time)

      return(CI)
    }
  } else {
    AsyNormbeta_shared = joint_dist_Penal_Int_shared(E = select$E, NE = select$NE, pes_outcome = "yDR", data = ps, id = ID, time = time)
    PQR_shared = PQR_Pint_shared(AsyNormbeta_shared, select)

    CI_per_select_var = function(ej) {
      AsyNormbeta_ej = joint_dist_Penal_Int_ej(AsyNormbeta_shared, ej)
      PQR_ej = PQR_Pint_ej(PQR_shared, AsyNormbeta_ej, AsyNormbeta_shared)
      condition = conditional_dist(PQR_shared, PQR_ej, AsyNormbeta_shared, AsyNormbeta_ej, select)
      CI = pivot_split_update(PQR_shared, PQR_ej, condition, AsyNormbeta_shared, AsyNormbeta_ej, select, level = level, pes_outcome = "yDR", data = ps,
                              id = ID, time = time)

      return(CI)
    }
  }


  #the list that stores the ej vectors for selected variables
  ejs = list()
  for(i in 1:length(select$E)) {
    vec_ej = rep(0,length(select$E))
    vec_ej[i] = 1

    ejs[[i]] = vec_ej
  }

  if(CI_algorithm == "parallel") {
    require(parallel)

    if(!is.null(core_num)) {cl = makeCluster(core_num)} else {cl = makeCluster(detectCores())}

    # List all objects in the global environment
    all_functions <- ls(envir = globalenv())

    # Filter the list to keep only the function names (not variables)
    function_names <- all_functions[sapply(all_functions, function(x) is.function(get(x)))]

    clusterExport(cl, varlist = function_names, envir = environment())
    # do parallel calculation
    results = parLapply(cl, ejs, CI_per_select_var)

    stopCluster(cl)
  }

  if(CI_algorithm == "lapply") {results <- lapply(ejs, CI_per_select_var)}

  if(CI_algorithm == "doParallel") {
    require(doParallel)
    require(foreach)

    if(!is.null(core_num)) {cl = makeCluster(core_num)} else {cl = makeCluster(detectCores())}
    registerDoParallel(cl)

    results = foreach(i = 1:length(select$E)) %dopar% {
      ej = ejs[[i]]
      CI_per_select_var(ej)
    }

    stopCluster(cl)

  }


  final_results <- do.call(rbind, lapply(results, function(x) {
    as.data.frame(as.list(x), stringsAsFactors = FALSE)
  }))

  num_cols <- setdiff(names(final_results), c("E", "message_low", "message_up"))
  final_results[num_cols] <- lapply(final_results[num_cols], as.numeric)

  # now makes the returned outcomes look like a table
  final_results = data.frame(
    E = character(),
    GEE_est = numeric(),
    post_beta = numeric(),
    pvalue = numeric(),
    lowCI = numeric(),
    upperCI = numeric(),
    prop_low = numeric(),
    prop_up = numeric(),
    message_low_CI = character(),
    message_up_CI = character()
  )

  for(i in 1:length(results)) {final_results[i,] = unlist(results[[i]])}

  return(final_results)

  # Output:
  # A table
  # E: the selected variables for which CI is calculated
  # GEE_est: the GEE estimate for this predictor
  # post_beta: the post selection true value for this predictor if simulation is conducted. Otherwise, NA is provided.
  # pvalue: the p value
  # lowCI: the lower bound of the confidence interval
  # upperCI: the upper bound of the confidence interval
  # prop_low: the true corresponding pivot value for the lower bound. For example, if it's 90% CI, this
  #          value will close to 0.05
  # prop_up: the true corresponding pivot value for the upper bound. For example, if it's 90% CI, this
  #          value will close to 0.95
}




