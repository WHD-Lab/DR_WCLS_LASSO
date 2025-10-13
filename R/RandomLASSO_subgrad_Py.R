# This version allow the program penalize the intercept


variable_selection_PY_penal_int = function(data,ID, moderator_formula, lam = NULL, noise_scale = NULL,
                                 splitrat = 0.8, virtualenv_path, ridge_term = 0, beta = NULL) {
  # data: the output of pesudo_outcomecal function
  # ID: the name of column where participants' ID are stored
  # moderator_formula: determines the formula for the f(St)T*beta function
  # lam: the value of lambda used to compute beta.
  # noise_scale: Scale of Gaussian noise added to objective. Default is 0.5*sd(y) times the sqrt of the mean of the trace of X^TX.
  #             where omega is drawn from IID normals with standard deviation noise_scale
  # splitrat: the corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
  #           This value will be used only when user doesn't provide the noise_scale.
  # virtualenv_path: Python virtual environment path
  # ridge_term

  ptSt = data[,"ptSt"]
  n = dplyr::n_distinct(data[,ID])
  wssqrt = c(sqrt(ptSt*(1-ptSt)*2/n^(1/2)))

  X = data.matrix(modelr::model_matrix(data, moderator_formula))
  Xw = X * wssqrt # this gives the same result as diag(wssqrt) %*% X
  ystring <- toString(moderator_formula[[2]])
  val_names <- all.vars(moderator_formula)
  y <- data.matrix(dplyr::select(stats::na.omit(data[,val_names]),tidyselect::all_of(ystring)))
  yw = c(wssqrt*y)

  #if(is.null(lam)) {lam = sqrt(2*log(dim(Xw)[2]))*sd(yw)*splitrat*sqrt(dim(Xw)[1])}
  mod = lm(yw ~ Xw)
  dispersion = sum(mod[["residuals"]]^2)/mod[["df.residual"]] # this value is slightly different from the python
  # because in python they use (Moore-Penrose) pseudo-inverse of a matrix

  if(is.null(noise_scale)) {
    noise_scale = sqrt(dispersion * (1 - splitrat)/splitrat * dim(Xw)[1])
  }
  if(is.null(lam)) {lam = sqrt(2*log(dim(Xw)[2]))*sd(yw)*splitrat*sqrt(dim(Xw)[1])}

  # load python virtualenv
  require(reticulate)
  use_virtualenv(virtualenv_path)
  #use_condaenv(condaenv = 'env3', conda = "/opt/anaconda3/bin/conda", required = TRUE)
  # load required modules and functions
  np = import("numpy")
  selectinf = import("selectinf")
  lassopy = selectinf$randomized$lasso$lasso
  selected_targets = selectinf$base$selected_targets
  const = lassopy$gaussian
  exact_grid_inference = selectinf$randomized$exact_reference$exact_grid_inference

  # convert data to Python array
  X1 = array_reshape(Xw, c(dim(Xw)[1], dim(Xw)[2]))
  Y1 = np_array(yw, dtype = "float64")
  # penalize the intercept
  wt = np_array(c(0,rep(lam,dim(Xw)[2]-1)))


  # plug in function
  conv = const(X1, Y1, feature_weights = wt, randomizer_scale = noise_scale, ridge_term = ridge_term)

  # get sign solution
  signs = conv$fit()
  nonzero = (signs!=0)
  # beta_lambdaE and subgradient
  perturb = py_get_attr(conv, "_initial_omega")
  soln = c(conv$observed_soln)
  subgrad = c(conv$observed_subgrad/lam)
  perturb = c(perturb$T)

  E = colnames(X)[nonzero]
  NE = colnames(X)[!nonzero] # NE may include intercept
  Z = subgrad[!nonzero]

  # calculate the post selection beta
  if(!is.null(beta) & sum(nonzero) > 1) {
    postbeta = np$dot(np$linalg$pinv(Xw[, nonzero]), np$dot(Xw, beta))
  } else if(!is.null(beta) & sum(nonzero) == 1) {
    postbeta = solve(t(Xw[, nonzero]) %*% Xw[, nonzero]) %*% t(Xw[, nonzero]) %*% Xw %*% beta
  } else {postbeta = rep(NA, sum(nonzero))}
  # below calculation matched with the above
  # But be careful if you change the weights, have to double check two calculation
  #wt = ptSt * (1-ptSt)
  #postbeta2 = solve(t(X[,nonzero] * c(wt / sqrt(n))) %*% X[,nonzero]) %*% t(X[,nonzero] * c(wt / sqrt(n))) %*% X %*% beta

  # get all Python code return
  #disperson = np$linalg$norm(Y1 - np$dot(X1[,nonzero], np$dot(np$linalg$pinv(X1[,nonzero]), Y1)))^2 / (n - sum(nonzero))
  #conv$setup_inference(dispersion=dispersion)
  #target_spec = selected_targets(conv$loglike,
  #                               conv$observed_soln,
  #                               dispersion=dispersion)
  #result = conv$inference(target_spec,
  #                        method='exact',
  #                        level=0.9)
  #G = exact_grid_inference(query_spec = conv$specification, target_spec = target_spec)
  #grids = as.matrix(G$stat_grid)

  return(list(formula = moderator_formula, E = E, NE = NE, n = n,
              perturb = perturb/(-2), lam = lam/(-2), Z = Z, OMEGA = (noise_scale)^2/4,
              sign_soln = signs[nonzero], soln = soln, postbeta = postbeta, nonzero = nonzero
              #Python_Output = list(conv = conv,
              #                     target_spec = target_spec,
              #                     result = result,
              #                     dispersion = dispersion,
              #                     dispersion_point = np$dot(np$linalg$pinv(X1[,nonzero]), Y1),
              #                     HEE = t(X1[,nonzero]) %*% X1[,nonzero]/c(n),
              #                     grids = grids)
              ))

  # Output:
  # E: Selected variables that remain in the model after penalization.
  # NE: Non-selected variables that were excluded from the model.
  # n: Number of unique subjects in the dataset.
  # perturb: Random noise added for each variables.
  # lam: Regularization parameter (lambda) used in the penalization.
  # Z: Subgradient of unselected variables.
  # OMEGA: (noise_scale)^2/4. Divide 4 because in theoretical calculation -1/2 is hidden inside penalty and random noise
  # sign_soln: Signs of the estimated coefficients for selected variables.
  # soln: Estimated coefficients of selected variables from the penalized regression.
  # postbeta: True beta projected on space of post-selection predictors.
  # nonzero: Boolean vector indicating which variables were selected.
}

