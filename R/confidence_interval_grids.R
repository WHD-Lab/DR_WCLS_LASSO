# calculate TP^{[I_,I+]} or S^[a,b] values

truncate_normal_weights = function(PQR_shared, PQR_ej, joint_distcal_shared, joint_distcal_ej, select_E,
                                  condition, grid_values) {

  lowerbound = condition[["pivot_related"]][["bounds"]][["lowerbond"]]
  upperbound = condition[["pivot_related"]][["bounds"]][["upperbond"]]
  eta = condition[["eta"]]
  etaLAMBDAeta = condition[["pivot_related"]][["etaLAMBDA"]]
  HE = condition[["HE"]]
  omega = condition[["omega"]]
  P = PQR_ej[["P"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  lam = select_E[["lam"]]
  ZNE = matrix(select_E[["Z"]], ncol = 1)
  R = PQR_shared[["R"]]
  n = joint_distcal_shared[["n"]]
  sigmajsq = condition[["pivot_related"]][["sigmasq_final"]]
  zetaj = condition[["pivot_related"]][["zetaj"]]
  lambdaj = condition[["pivot_related"]][["lambdaj"]]
  betaEjhat = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]

  NEnum = pnum - Enum

  # calculate weights
  Zero = matrix(0,nrow = Enum, ncol = pnum - Enum)
  lammatrix = diag(rep(lam, pnum - Enum))

  # log transformed truncated cdf value
  log_trunc_cdf = sapply(grid_values, FUN = function(betaEj) {
    mu = solve(t(HE) %*% solve(omega) %*% HE) %*% (t(HE) %*% solve(omega)) %*%
      (P %*% rbind(betaEj, GammaEjPerp) + rbind(Zero, lammatrix) %*% ZNE + R) /c(sqrt(n))
    etamu = t(eta) %*% mu

    logtrunc = log(pnorm(q= upperbound, mean = etamu, sd = sqrt(etaLAMBDAeta)) -
                     pnorm(q= lowerbound, mean = etamu, sd = sqrt(etaLAMBDAeta)))
    return(logtrunc)
  })

  logW = log_trunc_cdf - (grid_values - c(betaEjhat))^2/(2 * c(sigmajsq))
  # stablize calculate all values mins the max value
  weights = unlist(logW - max(logW))

  return(list(grid_values = grid_values,
              weights = weights))
}


# define discrete pdf function of the pivot
pdf_discrete = function(theta, grids, weights) {
  M = max(weights + grids * theta) - 5
  numerator = exp(weights + grids*theta - M)
  denominator = sum(numerator)
  pdfs = numerator/denominator

  return(list(grid_values = grids,
              pdf_value = pdfs))
}

# define ecdf function of the pivot here
cdf_discrete = function(theta, upper, grids, weights) {

  pdf_function_return = pdf_discrete(theta, grids, weights)

  pdfs = pdf_function_return[["pdf_value"]]

  return(sum(pdfs[grids <= c(upper)]))
}

# mimic the python code to decide the naive searching range
E_Var_null = function(theta, grids, weights) {
  pdf_function_return = pdf_discrete(theta, grids, weights)

  density_value = pdf_function_return[["pdf_value"]]
  E = sum(grids * density_value)
  Var = sum((grids - E)^2 * density_value)

  return(list(E = E, Var = Var))
}


##############################################################
# function that actually calculates the upper and lower bound

CI_grids = function(logW, betaEjhat, var_name, condition, tol = 10^(-3), level){
  # logW: the results of function truncate_normal_wights
  # betaEjhat: the point estimate from data
  # var_name: the name of variable, for which the inference is conducted
  # condition: the return of function conditional_dist
  # tol: the maximum accept difference between target quntile and estimated one
  # level: the nominate level of CI

  betaEjhat = c(betaEjhat)
  grid_values = logW[["grid_values"]]
  weights = logW[["weights"]]
  lambdaj = condition[["pivot_related"]][["lambdaj"]]
  zetaj = condition[["pivot_related"]][["zetaj"]]
  sigmajsq = condition[["pivot_related"]][["sigmasq_final"]]

  # initial range
  E_Var_theta0 = E_Var_null(theta = 0, grid_values, weights)
  E_null = E_Var_theta0[["E"]]
  sigma_null = sqrt(E_Var_theta0[["Var"]])
  lb = E_null - 20 * sigma_null
  ub = E_null + 20 * sigma_null

  lowertarget = (1 - level)/2
  uppertarget = level + (1 - level)/2

  thetalb = find_root(grids = grid_values, weights = weights,
                      lb = lb, ub = ub, target = uppertarget, betaEjhat)

  thetaub = find_root(grids = grid_values, weights = weights,
                      lb = lb, ub = ub, target = lowertarget, betaEjhat)

  # convert back to the betaEj scale instead of theta
  adjusted_value = (betaEjhat - zetaj)/lambdaj
  betaEjlb = thetalb$theta_values/lambdaj * sigmajsq + adjusted_value
  betaEjub = thetaub$theta_values/lambdaj * sigmajsq + adjusted_value

  # calculate p value
  theta_null = (zetaj - betaEjhat)/sigmajsq
  cdf_theta_null = cdf_discrete(theta = c(theta_null), upper = betaEjhat, grid_values, weights)
  pvalue = 2 * min(cdf_theta_null, 1 - cdf_theta_null)

  return(list(var_name = var_name,
              point_est = betaEjhat,
              lb_CI = betaEjlb,
              ub_CI = betaEjub,
              lb_quntile_CI = thetaub$pivot_values,
              ub_quntile_CI = thetalb$pivot_values,
              pvalue = pvalue))

}


find_root = function(grids, weights, lb, ub, target, betaEjhat, tol = 10^(-3)) {
  # target is the value we hope the CI to achieve like 0.05 or 0.95 for 90% CI

  # check pivot value at lb and ub
  pivot_lb = cdf_discrete(theta = lb, upper = betaEjhat, grids, weights)
  pivot_ub = cdf_discrete(theta = ub, upper = betaEjhat, grids, weights)


  # because pivot is monotonic decreasing function of theta
  # so when theta is smaller the pivot value is closer to 1
  if(pivot_lb < target) {
    # we need to adjust lower bound
    a = lb - (ub - lb)
    pivot_a = cdf_discrete(theta = a, upper = betaEjhat, grids, weights)
    while(pivot_a < target) {
      a = a - (ub - a)
      pivot_a = cdf_discrete(theta = a, upper = betaEjhat, grids, weights)
    }
  } else {
    a = lb
  }

  if(pivot_ub > target) {
    b = ub + (ub - lb)
    pivot_b = cdf_discrete(theta = b, upper = betaEjhat, grids, weights)

    while(pivot_b > target) {
      b = b + (b - lb)
      pivot_b = cdf_discrete(theta = b, upper = betaEjhat, grids, weights)
    }
  } else {
    b = ub
  }

  c = (a + b)/2
  pivot_c =  cdf_discrete(theta = c, upper = betaEjhat, grids, weights)

  while(abs(pivot_c - target) > tol) {
    if(pivot_c > target) {
      a = c
    } else {
      b = c
    }

    c = (a + b)/2
    pivot_c =  cdf_discrete(theta = c, upper = betaEjhat, grids, weights)
  }

  return(list(theta_values = c,
              pivot_values = pivot_c))
}








