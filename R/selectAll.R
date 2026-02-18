# this file is used when NE = NULL

joint_dist_selectAll = function(E, pes_outcome, data, id, time) {
  # E: vector of selected predictors (don't include intercept)
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # id: column name, a vector which identifies individuals
  # time: column name, a vector that records the decision points for each individual

  require(dplyr)
  if("(Intercept)" %in% E) {ftStE = t(cbind(1,as.matrix(data[,E[E!="(Intercept)"]])))} else {ftStE = t(as.matrix(data[,E]))}

  id = data[,id]
  n = n_distinct(id)

  # get betaE point estiamtes
  wt = data$ptSt * (1-data$ptSt)
  idf = as.factor(id)
  time = data[,time]

  if("(Intercept)" %in% E & length(E) > 1) {
    formula = as.formula(paste(pes_outcome, "~", paste(E[which(E != "(Intercept)")], collapse = "+")))
  } else if("(Intercept)" %in% E & length(E) == 1) {
    formula = as.formula(paste(pes_outcome, "~1"))
  } else {
    formula = as.formula(paste(pes_outcome, "~-1+", paste(E, collapse = "+")))
  }

  betaEM = geepack::geeglm(formula, data = data, weights = wt/sqrt(n), corstr = "independence", id = idf,
                           waves = time)

  res = betaEM[["residuals"]] # this is unweighted residuals

  #sigmaTT
  HEE = H(ftStE, ftStNE, wt, "EE", n)
  KEE = K(ftStE, ftStNE, wt, "EE", betaEM, n)
  sigmaTT = solve(HEE) %*% KEE %*% solve(HEE)


  # point estimate betaE
  betaE = matrix(betaEM[["coefficients"]], ncol = 1)


  return(list(wt = wt,
              betaEM = betaEM,
              H = list(HEE = HEE),
              K = list(KEE = KEE),
              Sigma = list(sigmaTT = sigmaTT),
              n = n,
              E = E
  ))

  # Output:
  # wt: ptSt * (1 - ptSt).
  # betaEM: fitted model for selected predictors.
  # H: HEE
  # K: KEE
  # Sigma: sigmaTT
  # n: # of unique subjects in the dataset.
  # E: selected predictors.
}

#################################################################################
# can still use the original "joint_dist_Penal_Int_ej" function
#################################################################################

PQR_shared_selectALL = function(joint_distcal_shared, select_E) {
  # joint_distcal_shared: the result of function joint_dist
  # select_E: the result of function variable_selection_PY_penal_int

  HEE = joint_distcal_shared[["H"]][["HEE"]]
  sigmaTT = joint_distcal_shared[["Sigma"]][["sigmaTT"]]
  E = joint_distcal_shared[["E"]]
  n = joint_distcal_shared[["n"]]

  Enum = length(E)
  pnum = Enum

  ###### Below code need to edit after you build randomized lasso ######
  lam = select_E[["lam"]]
  se = select_E[["sign_soln"]]

  # matrix P shared part
  P_share = HEE * sqrt(n)

  # matrix Q
  Q = HEE * c(-sqrt(n))

  # matrix R
  R = matrix(se*c(0, rep(lam, Enum-1)), ncol = 1)

  return(list(P_share = P_share,
              Q = Q,
              R = R,
              Enum = Enum,# number of selected parameters
              pnum = pnum, # number of all potential parameters
              se = se))

  # Output:
  # P_share: the part of P matrix that will not change with ej value
  # Q: Q matrix
  # R: R matrix
  # Enum: number of selected parameters
  # pnum: number of all potential parameters
  # se: Signs of the estimated coefficients for selected variables.

}

PQR_ej_selectALL = function(PQR_Pint_shared, joint_dist_Penal_Int_ej, joint_distcal_shared) {

  betaEjperp = joint_dist_Penal_Int_ej[["betaEjperp_cal"]][["betaEjperp"]]
  ej = joint_dist_Penal_Int_ej[["ej"]]
  sigmaTT = joint_distcal_shared[["Sigma"]][["sigmaTT"]]
  P_shared = PQR_Pint_shared[["P_share"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  n = joint_distcal_shared[["n"]]

  # Gamma Ej perp
  GammaEjPerp = betaEjperp

  # matrix P depends on ej part
  p1 = HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  P = cbind(p1 * sqrt(n), P_shared)


  return(list(GammaEjPerp = GammaEjPerp,
              p1 = p1* sqrt(n),
              P = P))
  # Output:
  # GammaEjPerp: the GammaEjPerp matrix
  # p1: the 1st block of P matrix. This value will be used to calculate eta
  # P: the P matrix
}

#####################################################################################
conditional_dist_selectALL = function(PQR_shared, PQR_ej, joint_distcal_shared, joint_distcal_ej, select_E) {
  # PQR_shared: output of function PQR_Pint_shared
  # PQR_ej: output of function PQR_Pint_ej
  # joint_distcal_shared: output of function joint_dist_Penal_Int_shared
  # joint_distcal_ej: output of function joint_dist_Penal_Int_ej
  # select_E: output of function variable_selection_PY_penal_int

  OMEGA = select_E[["OMEGA"]] # need to change it when get randomized lasso down
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  P = PQR_ej[["P"]]
  Q = PQR_shared[["Q"]]
  R = PQR_shared[["R"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  betaEj = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]
  sigmasq1 = joint_distcal_ej[["betaEj_cal"]][["sigmasq1"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  lam = select_E[["lam"]] # need to change it when get randomized lasso down
  n = joint_distcal_shared[["n"]]
  p1 = PQR_ej[["p1"]]
  hat_betaE_lambda = select_E[["soln"]]
  se = PQR_shared[["se"]]

  # create omega
  if(length(OMEGA) == 1) {
    omega = diag(rep(OMEGA, pnum))
  } else {
    omega = OMEGA
  }

  # conditional distribution of hat{beta}^{lambda}_E
  LAMBDA = solve(t(Q) %*% solve(omega) %*% Q)

  # calculate values related to pivot
  eta = -t(Q) %*% solve(omega) %*% p1
  Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)

  hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]

  Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)

  mu_final = betaEj
  sigmasq_final = sigmasq1 / n

  # calculate bounds for truncated normal distribution
  bounds = support(se, eta, LAMBDA, Aeta, Qn)
  #etamu = t(eta) %*% mu
  etaLAMBDA = t(eta) %*% LAMBDA %*% eta

  return(list(omega = omega,
              #mu = mu,
              LAMBDA = LAMBDA,
              eta = eta,
              Aeta = Aeta,
              Qn = Qn,
              se = se,
              n = n,hat_betaE_lambda = hat_betaE_lambda,
              pivot_related = list(mu_final = mu_final,
                                   sigmasq_final = sigmasq_final,
                                   bounds = bounds,
                                   #etamu = etamu,
                                   etaLAMBDA = etaLAMBDA,
                                   zetaj = 0,
                                   lambdaj = 1)))
  # Output:
  # omega: Matrix version of OMEGA. It the variance of added random noised.
  # mu: the mean of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # LAMBDA: the variance of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # eta: the specifically designed vector that will be used to reduce integration dimension when calculate pivot
  # Aeta: perpendicular to betaEj_lambda
  # Qn: value used to obtain support
  # se: Signs of the estimated coefficients for selected variables.
  # n: # of unique subjects in the dataset.
}

truncate_normal_weights_selectALL = function(PQR_shared, PQR_ej, joint_distcal_shared, joint_distcal_ej, select_E,
                                   condition, grid_values) {

  lowerbound = condition[["pivot_related"]][["bounds"]][["lowerbond"]]
  upperbound = condition[["pivot_related"]][["bounds"]][["upperbond"]]
  eta = condition[["eta"]]
  etaLAMBDAeta = condition[["pivot_related"]][["etaLAMBDA"]]
  Q = PQR_shared[["Q"]]
  omega = condition[["omega"]]
  P = PQR_ej[["P"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  lam = select_E[["lam"]]
  R = PQR_shared[["R"]]
  n = joint_distcal_shared[["n"]]
  sigmajsq = condition[["pivot_related"]][["sigmasq_final"]]
  zetaj = condition[["pivot_related"]][["zetaj"]]
  lambdaj = condition[["pivot_related"]][["lambdaj"]]
  betaEjhat = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]

  NEnum = pnum - Enum

  # calculate weights
  # log transformed truncated cdf value
  log_trunc_cdf = sapply(grid_values, FUN = function(betaEj) {
    mu = -solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(betaEj, GammaEjPerp))

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





















