# This function calculate all conditional distribution for betaEj, GammaPerp which will be used
# to built later pivot

conditional_dist = function(PQR_shared, PQR_ej, joint_distcal_shared, joint_distcal_ej, select_E) {
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
  HNEE = joint_distcal_shared[["H"]][["HNEE"]]
  lam = select_E[["lam"]] # need to change it when get randomized lasso down
  n = joint_distcal_shared[["n"]]
  p1 = PQR_ej[["p1"]]
  p4 = PQR_ej[["p4"]]
  ZNE = matrix(select_E[["Z"]], ncol = 1)
  hat_betaE_lambda = select_E[["soln"]]
  se = select_E[["sign_soln"]]

  NEnum = pnum - Enum

  # create omega
  if(length(OMEGA) == 1) {
    omega = diag(rep(OMEGA, pnum))
  } else {
    omega = OMEGA
  }

  # conditional distribution of hat{beta}^{lambda}_E, Z_{-E}
  delta = -solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(betaEj, GammaEjPerp))
  theta = solve(t(Q) %*% solve(omega) %*% Q)

  # conditional distribution of Z_{-E}
  delta2 = matrix(delta[(Enum+1):pnum,], ncol = 1)
  theta22 = theta[(Enum+1):pnum, (Enum+1):pnum]

  # conditional distribution of hat{beta}^{lambda}_E
  HE = rbind(HEE, HNEE)
  Zero = matrix(0,nrow = Enum, ncol = pnum - Enum)

  if(pnum - Enum > 1) {
    lammatrix = diag(rep(lam, pnum - Enum))
  } else {
    lammatrix = lam
  }

  mu = solve(t(HE) %*% solve(omega) %*% HE) %*% (t(HE) %*% solve(omega)) %*%
    (P %*% rbind(betaEj, GammaEjPerp) + rbind(Zero, lammatrix) %*% ZNE + R) /c(sqrt(n))
  LAMBDA = solve(c(n)*t(HE) %*% solve(omega) %*% HE)

  # calculate values related to pivot
  eta = t(HE) %*% solve(omega) %*% rbind(p1, p4)
  Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)

  hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]

  Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)

  zeros = matrix(0, nrow = NEnum, ncol = Enum)
  I = diag(rep(1, NEnum))
  shared = cbind(zeros, I)
  zero1p = matrix(0, nrow = 1, ncol = pnum)
  Ip = diag(rep(1, pnum))
  A = -shared %*% solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% P %*% matrix(c(1,rep(0,pnum)), ncol = 1)
  C = -shared %*% solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(zero1p, Ip) %*% GammaEjPerp)

  mu_final = (n*betaEj + c(sigmasq1) * (t(ZNE - C) %*% solve(theta22) %*% A))/(n + c(sigmasq1) * (t(A) %*% solve(theta22) %*% A))
  sigmasq_final = sigmasq1 / (n + c(sigmasq1) * (t(A) %*% solve(theta22) %*% A))

  # calculate bounds for truncated normal distribution
  bounds = support(se, eta, LAMBDA, Aeta, Qn)
  etamu = t(eta) %*% mu
  etaLAMBDA = t(eta) %*% LAMBDA %*% eta

  return(list(omega = omega,
              delta = delta,
              theta = theta,
              delta2 = delta2,
              theta22 = theta22,
              HE = HE,
              subgradient_unselected = ZNE,
              mu = mu,
              LAMBDA = LAMBDA,
              eta = eta,
              Aeta = Aeta,
              Qn = Qn,
              se = PQR_shared[["se"]],
              n = n,hat_betaE_lambda = hat_betaE_lambda,
              pivot_related = list(mu_final = mu_final,
                                   sigmasq_final = sigmasq_final,
                                   bounds = bounds,
                                   etamu = etamu,
                                   etaLAMBDA = etaLAMBDA,
                                   zetaj = sigmasq_final*(t(ZNE - C) %*% solve(theta22) %*% A),
                                   lambdaj = n * sigmasq_final/sigmasq1)
              ))
  # Output:
  # omega: Matrix version of OMEGA. It the variance of added random noised.
  # delta: the mean of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta: the variance of conditional distribution (hat{beta}^{lambda}_E, Z_{-E}) | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # delta2: the mean of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # theta22: the variance of conditional distribution Z_{-E} | (hat{beta}_Ej, hat{Gamma}EjPerp).
  # HE: the combined vector of HEE and HNEE.
  # subgradient_unselected: the subgradient vector of unselected variables
  # mu: the mean of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # LAMBDA: the variance of conditional distribution hat{beta}^{lambda}_E | (hat{beta}_Ej, hat{Gamma}EjPerp, Z_{-E}).
  # eta: the specifically designed vector that will be used to reduce integration dimension when calculate pivot
  # Aeta: perpendicular to betaEj_lambda
  # Qn: value used to obtain support
  # se: Signs of the estimated coefficients for selected variables.
  # n: # of unique subjects in the dataset.
}


# below function calculate support for truncated normal distribution
support = function(se, eta, LAMBDA, Aeta, Qn) {
  
  diagSe = diag(se)

  lowerbond = -Inf
  upperbond = Inf

  for(i in 1:length(se)) {
    # warning: I manually set i starts from 2 to skip intercept
    # Because we lack se and point estimate for intercept from randomized lasso
    # to avoid manipulated value impact later result I skip it

    frac = (t(diagSe[,i]) %*% Aeta)/(-t(diagSe[,i]) %*% Qn)

    if(t(diagSe[,i]) %*% LAMBDA %*% eta > 0) {
      lowerbond = max(frac, lowerbond)
    }

    if(t(diagSe[,i]) %*% LAMBDA %*% eta < 0) {
      upperbond = min(frac, upperbond)
    }
  }
  # what if exactly equal to 0

  return(list(upperbond = upperbond,
              lowerbond = lowerbond))
}
