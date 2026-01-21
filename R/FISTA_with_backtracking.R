# FISTA with backtracking

# Lipschitz constant for l1 regularization depens on the maximum eigenvalue of
# inner product of design matrix.
# Sometimes this value is hard to obtain
# so we can use below algorithm to estimate the Lipschitz constant


#' FISTA_backtracking
#'
#' This fuction applies the FISTA algorithm with backtracking to solve randomized LASSO problem. Reference paper is
#' "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems" by Amir Beck and Marc Teboulle.
#'
#' @param data the output of pseudo_outcomecal function
#' @param ID the name of column where participants' ID are stored
#' @param moderator_formula determines the formula for the f(St)T*beta function
#' @param lam the value of penalty term of randomized LASSO. If it is not provided, the default value will be used.
#' Default is
#' \eqn{\sqrt{2n* logp} \rho sd(y)} where \eqn{\rho}
#'   is the split rate and \eqn{n} is the number of rows.
#' @param noise_scale Standard deviation of Gaussian noise added to objective.
#' Default is \eqn{\sqrt{\frac{1-\rho}{\rho}\, n}\,\mathrm{sd}(y)} where \eqn{\rho}
#'   is the split rate and \eqn{n} is the number of rows.
#' The random noises, \eqn{\omega}, are iid drawn from normal distribution with standard deviation noise_scale
#' @param splitrat this value is corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
#' This value will be used only when user doesn't provide the `noise_scale` or `lam`.
#' @param max_ite the maximum iteration for searching predictors
#' @param tol when the improvement of residual sum of square less than this number, stop searching
#' @param beta the true coefficient value (if simulation is conducted)
#'
#' @return
#' \describe{
#' \item{formula}{the raw model before selection}
#' \item{E}{Selected variables}
#' \item{NE}{Not selected variables}
#' \item{n}{the number of participants in the study}
#' \item{soln}{Estimated coefficients of selected variables from the penalized randomized regression.}
#' \item{Z}{Subgradient of unselected variables.}
#' \item{OMEGA}{The scaled variance for the added random noise. This value is scaled by 4.}
#' \item{lam}{Scale orginal regularization parameter by -2 for later inference purpose.}
#' \item{ori_lam}{Regularization parameter (lambda) used in the penalization.}
#' \item{perturb}{Random noise added for each variables. This value is scaled by -2 for later inference purpose.}
#' \item{nonzero}{Boolean vector indicating which variables were selected.}
#' \item{postbeta}{True beta projected on space spinned by post-selection predictors}
#' \item{sign_soln}{Signs of the estimated coefficients for selected variables}
#' }
#'
#' @examples
#'
#' sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
#'   beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
#'   beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
#'   theta1 = 0.8)
#'
#' Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
#' St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))
#'
#' pseudo_outcome_CVlasso = pseudo_outcome_generator_CVlasso(fold = 5,ID = "id", sim_data, Ht,
#'   St, "action", "outcome",core_num = 5)
#'
#' my_formula = as.formula(paste("yDR ~ ", paste(St, collapse = " + ")))
#'
#' FISTA_backtracking(data = sim_data, moderator_formula = my_formula, lam = NULL, noise_scale = NULL,
#'   splitrat = 0.7, beta = matrix(c(-1, 1.7, 1.5, -1.3, -1, rep(0,21)), ncol = 1))
#'
#' @export
#' 

FISTA_backtracking = function(data,ID, moderator_formula, lam = NULL, noise_scale = NULL,
                              splitrat = 0.8, max_ite = 10^(5), tol = 10^(-4), beta = NULL){

  # data: the output of pesudo_outcomecal function
  # ID: the name of column where participants' ID are stored
  # moderator_formula: determines the formula for the f(St)T*beta function
  # lam: the value of lambda used to compute beta.
  # noise_scale: Scale of Gaussian noise added to objective. Default is 0.5*sd(y) times the sqrt of the mean of the trace of X^TX.
  #             where omega is drawn from IID normals with standard deviation noise_scale
  # splitrat: the corresponding to the data splitting rate. Details can read "Exact Selective Inference with Randomization" page 15 equation (10).
  #           This value will be used only when user doesn't provide the noise_scale.
  # max_ite: the maximum iteration for searching predictors
  # tol: when the improvement of residual sum of square less than this number, stop searching
  # beta: the true coefficient value (if simulation is conducted)

  ptSt = data[,"ptSt"]
  n = dplyr::n_distinct(data[,ID])
  wssqrt = c(sqrt(ptSt*(1-ptSt)/sqrt(n)))
  wt = ptSt * (1 - ptSt)

  X = data.matrix(modelr::model_matrix(data, moderator_formula))
  Xw = X * wssqrt
  ystring <- toString(moderator_formula[[2]])
  val_names <- all.vars(moderator_formula)
  y <- data.matrix(dplyr::select(stats::na.omit(data[,val_names]),tidyselect::all_of(ystring)))
  yw = c(wssqrt*y)

  # calcualate penalty term and noise scale
  mod = lm(yw ~ Xw-1)
  dispersion = sum(mod[["residuals"]]^2)/mod[["df.residual"]]

  if(is.null(noise_scale)) {
    noise_scale = sqrt(dispersion * (1 - splitrat)/splitrat * dim(Xw)[1])
  }

  if(is.null(lam)) {
    lam = sqrt(2*log(dim(Xw)[2]))*sd(yw)*splitrat* sqrt(dim(Xw)[1])
  }

  lam_vec = c(0, rep(lam, dim(Xw)[2] - 1))

  # simulate random variables
  perturb = matrix(rnorm(n = dim(Xw)[2], mean = 0, sd = noise_scale), ncol = 1)

  # intial y_k, t_k
  k = 1
  L_kmins1 = 3
  eta = 2
  x_kmins1 = matrix(mod[["coefficients"]], ncol = 1)
  y_k = matrix(mod[["coefficients"]], ncol = 1)
  t_k = 1
  resSumSqDiff = Inf

  while(k < max_ite & resSumSqDiff > tol) {
    L_k = i_finder(L_kmins1, y_k, design_matrix = X, outcome = y, wt, n, perturb, lam_vec, eta)

    p_L_returns = p_L(y_k, L_k, design_matrix = X, outcome = y, wt, n, perturb, lam_vec)
    x_k = p_L_returns[["x_k"]]
    t_kplus1 = (1 + sqrt(1 + 4 * t_k^2))/2
    y_kplus1 = x_k + (t_k - 1)/(t_kplus1) * (x_k - x_kmins1)

    y_k = y_kplus1
    t_k = t_kplus1
    L_kmins1 = L_k

    resSumSq_kmins1 = sum(c(y - X %*% x_kmins1)^2 * ptSt * (1 - ptSt))/sqrt(n) +
      lam*sum(abs(x_kmins1)) - c(t(perturb) %*% x_kmins1)
    resSumSq_k = sum(c(y - X %*% x_k)^2 * ptSt * (1 - ptSt))/sqrt(n) +
      lam*sum(abs(x_k)) - c(t(perturb) %*% x_k)
    resSumSqDiff = abs(resSumSq_kmins1 - resSumSq_k)

    x_kmins1 = x_k

    k = k + 1
  }

  sign_soln = sign(x_k)
  nonzero = (sign_soln != 0)

  if(!is.null(beta)) {
    postbeta = solve(t(Xw[, nonzero]) %*% Xw[, nonzero]) %*% t(Xw[, nonzero]) %*% Xw %*% beta
  } else {postbeta = rep(NA, sum(nonzero))}

  return(list(formula = moderator_formula,
              E = colnames(X)[nonzero], NE = colnames(X)[!nonzero],
              n = n, soln = c(x_k),
              Z = p_L_returns[["Z"]][!nonzero],
              OMEGA = (noise_scale^2)/4,
              lam = lam/(-2), ori_lam = lam,
              perturb = c(perturb/(-2)),
              nonzero = c(nonzero), postbeta = c(postbeta),
              sign_soln = sign_soln[nonzero]))
}

p_L = function(y_k, L, design_matrix, outcome, wt, n, perturb, lam_vec){
  # y_k: the value will be used to obtain x_k
  # L: estimater for the Lipschitz constant
  # design_matrix: in matrix form
  # outcome: in matrix form
  # wt: the weights. It is ptSt * (1 - ptSt)
  # n: the number of participants
  # perturb: the random value added to the loss function

  # first order derivative
  der_designmatrix = c(-2/sqrt(n)) *  t(design_matrix) %*%
    matrix((outcome - design_matrix %*% y_k)* wt, ncol = 1) - perturb

  x_k = pmax(abs(y_k - 1/L * der_designmatrix) - lam_vec/L, 0) * sign(y_k - der_designmatrix/L)

  # subgradients
  Z = L/lam_vec * (y_k - der_designmatrix/L - x_k)
  Z[1] = sign(x_k)[1]


  return(list(x_k = x_k,
              y_k = y_k,
               n = n,
              wt = wt,
              design_matrix = design_matrix,
              outcome = outcome,
              lam = lam_vec,
              L = L, perturb = perturb, Z = Z))
}

F_loss_function = function(p_L_results) {
  # p_L_results: the results of function p_L

  x_k = p_L_results[["x_k"]]
  n = p_L_results[["n"]]
  wt = p_L_results[["wt"]]
  design_matrix = p_L_results[["design_matrix"]]
  outcome = p_L_results[["outcome"]]
  lam = p_L_results[["lam"]]
  perturb = p_L_results[["perturb"]]

  loss_F = sum((outcome - design_matrix %*% x_k)^2 * wt)/sqrt(n) + sum(lam * abs(x_k)) - c(t(perturb) %*% x_k)

  return(loss_F = loss_F)
}

Q_loss_function = function(p_L_results) {
  # this is a quadratic approximation of the target loss function

  # p_L_results: the results of function p_L
  y_k = p_L_results[["y_k"]]
  n = p_L_results[["n"]]
  wt = p_L_results[["wt"]]
  design_matrix = p_L_results[["design_matrix"]]
  outcome = p_L_results[["outcome"]]
  perturb = p_L_results[["perturb"]]
  x_k = p_L_results[["x_k"]]
  L = p_L_results[["L"]]
  lam = p_L_results[["lam"]]

  f_y_k = sum((outcome - design_matrix %*% y_k)^2 * wt)/sqrt(n) - t(perturb) %*% y_k
  der_designmatrix_y_k = c(-2/sqrt(n)) *  t(design_matrix) %*%
    matrix((outcome - design_matrix %*% y_k)* wt, ncol = 1) - perturb

  QL_x_y = f_y_k + t(x_k - y_k) %*% der_designmatrix_y_k + L/2 * t(x_k - y_k) %*% (x_k - y_k) +  sum(lam * abs(x_k))

  return(QL_x_y = QL_x_y)
}

i_finder = function(L_kmins1, y_k, design_matrix, outcome, wt, n, perturb, lam_vec, eta){
  # y_k: the value will be used to obtain x_k
  # L_kmins1: previous estimater for the Lipschitz constant
  # design_matrix: in matrix form
  # outcome: in matrix form
  # wt: the weights. It is ptSt * (1 - ptSt)
  # n: the number of participants
  # perturb: the random value added to the loss function

  i_k = 0
  L_bar = eta^(i_k) * L_kmins1
  p_L_return = p_L(y_k, L_bar, design_matrix, outcome, wt, n, perturb, lam_vec)
  loss_F = F_loss_function(p_L_return)
  QL_x_y = Q_loss_function(p_L_return)

  while(loss_F > QL_x_y){
    i_k = i_k + 1
    L_bar = eta^(i_k) * L_kmins1
    p_L_return = p_L(y_k, L_bar, design_matrix, outcome, wt, n, perturb, lam_vec)
    loss_F = F_loss_function(p_L_return)
    QL_x_y = Q_loss_function(p_L_return)
  }
  
  return(L_bar)
}


