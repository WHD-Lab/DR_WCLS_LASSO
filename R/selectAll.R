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
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  lam = select_E[["lam"]] # need to change it when get randomized lasso down
  n = joint_distcal_shared[["n"]]
  p1 = PQR_ej[["p1"]]
  hat_betaE_lambda = select_E[["soln"]]

  # create omega
  if(length(OMEGA) == 1) {
    omega = diag(rep(OMEGA, pnum))
  } else {
    omega = OMEGA
  }

  # conditional distribution of hat{beta}^{lambda}_E
  mu = -solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q) %*% solve(omega) %*% (R + P %*% rbind(betaEj, GammaEjPerp))
  LAMBDA = solve(t(Q) %*% solve(omega) %*% Q)

  # calculate values related to pivot
  eta = -t(Q) %*% solve(omega) %*% p1
  Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)

  hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
  # if(length(hat_betaE_lambda) < Enum) {hat_betaE_lambda = matrix(c(0.001,hat_betaE_lambda), ncol = 1)}
  # WARNING: hat{beta}^lambda_E doesn't provide estimate for intercept. Here I use 0.001 for default
  # try to use fitted intercept from gee
  # if(length(hat_betaE_lambda) < Enum) {hat_betaE_lambda = matrix(c(joint_distcal[["betaEM"]][["coefficients"]][["(Intercept)"]]
  #                                                                 ,hat_betaE_lambda), ncol = 1)}

  Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)

  return(list(omega = omega,
              mu = mu,
              LAMBDA = LAMBDA,
              eta = eta,
              Aeta = Aeta,
              Qn = Qn,
              se = PQR_shared[["se"]],
              n = n,hat_betaE_lambda = hat_betaE_lambda))
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
####################################################################################

# I recalculate all previous values under this hypothesis value beta_E,j = b.

# step 1. Get an estimate for the theoretical beta value under
# the hypothesis

Asynorm_beta_under_hy_selectALL = function(E, pes_outcome, data, ej, id, time, null_value,
                                 HEE) {
  # E: vector of selected predictors
  # pes_outcome: column name for pesudo-outcome
  # data: the output of pesudo_outcomecal function
  # ej: jth standard basis vector. used to get betaE,j, don't need in matrix form
  # id: a vector which identifies individuals
  # time: a vector that records the decision points for each individual
  # null_value: the value of betaE,j under the null hypothesis
  # HEE: can recycle from previous calculation

  require(dplyr)

  id = data[,id]
  n = n_distinct(id)

  if(("(Intercept)" %in% E) & (length(E) > 1)) {
    ftStE = t(cbind(1,as.matrix(data[,E[E!="(Intercept)"]])))
  } else if (("(Intercept)" %in% E) & (length(E) == 1)) {
    ftStE = matrix(1, nrow = 1, ncol = dim(data)[1])
  } else {
    ftStE = t(as.matrix(data[,E]))
  }


  # modify the pesudo-outcome to remove the impact of betaEj
  data$yDR_modify = data[,pes_outcome] - t(ftStE) %*% matrix(ej, ncol = 1) * c(null_value)

  wt = data$ptSt * (1-data$ptSt)
  idf = as.factor(id)
  time = data[,time]

  if(length(E) > 1) {
    # modify E to remove betaE,j from selection
    E_modify = E[which(ej != 1)]
    other_terms = setdiff(E_modify, "(Intercept)")

    if (length(other_terms) == 0) {
      formula = as.formula("yDR_modify ~ 1")
    } else if("(Intercept)" %in% E_modify) {
      formula = as.formula(paste("yDR_modify ~", paste(E_modify[which(E_modify != "(Intercept)")], collapse = "+")))
    } else {
      formula = as.formula(paste("yDR_modify~-1+", paste(E_modify, collapse = "+")))
    }

    # get modified betaE point estiamtes
    betaEM_modify = geepack::geeglm(formula, data = data, weights = wt/sqrt(n), corstr = "independence", id = idf,
                                    waves = time)
  } else {
    betaEM_modify = list(residuals = data$yDR_modify,
                         id = id)
  }

  #sigmaTT
  KEE = K(ftStE, ftStNE, wt, "EE", betaEM_modify, n)
  sigmaTT = solve(HEE) %*% KEE %*% solve(HEE)


  # point estimate betaE
  if(length(E) > 1) {
    matrix_pre = diag(c(1-ej))
    matrix_pre = matrix_pre[,-which(ej == 1)]
    betaE_rest = matrix(betaEM_modify[["coefficients"]], ncol = 1)
  } else{
    matrix_pre = matrix(0, ncol = 1)
    betaE_rest = matrix(0, ncol = 1)
  }
  ej = matrix(ej, ncol = 1)
  betaE = ej*c(null_value) + matrix_pre %*% betaE_rest

  # betaEj
  betaEj = null_value
  sigmasq1 = t(ej) %*% sigmaTT %*% ej

  # betaEjperp
  betaEjperp = betaE - sigmaTT %*% ej %*% solve(t(ej) %*% sigmaTT %*% ej) * c(betaEj)
  sigma3 = sigmaTT - sigmaTT %*% ej %*% t(ej) %*% sigmaTT * c(solve(t(ej) %*% sigmaTT %*% ej))

  # return those updated value
  return(list(betaEM_modify = betaEM_modify,
              H = list(HEE = HEE),
              K_modify = list(KEE = KEE),
              Sigma_modify = list(sigmaTT = sigmaTT),
              betaEj_cal_modify = list(betaEj = betaEj, sigmasq1 = sigmasq1),
              betaEjperp_cal_modify = list(betaEjperp = betaEjperp, sigma3 = sigma3),
              n = n,
              ej = ej,
              Enum = length(E)
  ))
}

# now update PQR calculation
# only P matrix will change with the new betaE,j value
# due to SigmaTT and SigmaST is changed
P_modify_selectALL = function(Asynorm_beta_under_hy){
  # Asynorm_beta_under_hy: the result from function Asynorm_beta_under_hy

  HEE = Asynorm_beta_under_hy[["H"]][["HEE"]]
  sigmaTT = Asynorm_beta_under_hy[["Sigma_modify"]][["sigmaTT"]]
  n = Asynorm_beta_under_hy[["n"]]
  ej = Asynorm_beta_under_hy[["ej"]]
  Enum = Asynorm_beta_under_hy[["Enum"]]

  # matrix P
  p1 = HEE %*% sigmaTT %*% ej/c(t(ej) %*% sigmaTT %*% ej)
  p2 = HEE


  P = cbind(p1, p2)* sqrt(n)

  return(list(P = P,
              p1 = p1 * sqrt(n),
              p2 = p2 * sqrt(n),
              n = n))
}

###########
# can still use the same support update
pivot_split_update_sellectALL = function(PQR_shared, PQR_ej, cond_dist, joint_distcal_shared, joint_distcal_ej, select_E, level = 0.9, pes_outcome, data,
                              id, time, max_iterate = 10^{6}, max_tol = 10^{-3}) {
  # PQR_shared: result of function PQR_Pint_shared
  # PQR_ej: result of function PQR_Pint_ej
  # cond_dist: result of function conditional_dist
  # joint_distcal_shared: result of function joint_dist_Penal_Int_shared
  # joint_distcal_ej: result of function joint_dist_Penal_Int_ej
  # select_E: result of function variable_selection_PY_penal_int
  # level: the significant level of confidence interval
  # pes_outcome: the column names of pseudo outcome
  # data: the dataset with pseudo outcome
  # id: column names of participant id
  # time: column names of decidion time points
  # max_iterate: the max iteration number when searching CI
  # max_tol: the max error tolerance when calculate pivot value


  E = select_E[["E"]]
  ej = joint_distcal_ej[["ej"]]
  HEE = joint_distcal_shared[["H"]][["HEE"]]
  omega = cond_dist[["omega"]]
  n = cond_dist[["n"]]
  GammaEjPerp = PQR_ej[["GammaEjPerp"]]
  lam = select_E[["lam"]]
  R = PQR_shared[["R"]]
  LAMBDA = cond_dist[["LAMBDA"]]
  pnum = PQR_shared[["pnum"]]
  Enum = PQR_shared[["Enum"]]
  se = cond_dist[["se"]]
  hat_betaE_lambda = select_E[["soln"]]
  betaEjhat = joint_distcal_ej[["betaEj_cal"]][["betaEj"]]
  sigmasq1 = joint_distcal_ej[["betaEj_cal"]][["sigmasq1"]]
  Q = PQR_shared[["Q"]]


  lowp = (1-level)/2
  upp = 1 - (1-level)/2


  zeroEPNE = matrix(0, Enum, pnum-Enum)
  if(pnum - Enum > 1) {
    lamPNE = diag(rep(lam, pnum - Enum))
  } else {
    lamPNE = lam
  }


  # recreate self-define function to construct pivot
  my_function = function(b, sigmasq1, n, eta,p1,p2,null_value, Var_etabeta,lower, upper) {

    mu_etabeta = sapply(b, function(b) {-t(eta) %*% solve(t(Q) %*% solve(omega) %*% Q) %*% t(Q)%*% solve(omega) %*%
        (b * p1 + p2 %*% GammaEjPerp + R)})

    out = dnorm(b, mean = null_value, sqrt(sigmasq1/n)) *
      (pnorm(upper,mu_etabeta, sqrt(Var_etabeta)) - pnorm(lower,mu_etabeta, sqrt(Var_etabeta)))

    return(out)
  }

  pivot_prop_value = function(null_value) {
    # update joint beta distribution
    update_beta_dist = Asynorm_beta_under_hy_selectALL(E, pes_outcome, data, ej, id, time, null_value = null_value,
                                             HEE)
    sigmasq1 = update_beta_dist[["betaEj_cal_modify"]][["sigmasq1"]]
    # update P matrix
    update_P = P_modify_selectALL(update_beta_dist)
    p1 = update_P[["p1"]]
    p2 = update_P[["p2"]]
    P = update_P[["P"]]

    eta = -t(Q) %*% solve(omega) %*% p1
    Var_etabeta  = t(eta) %*% LAMBDA %*% eta
    Qn = LAMBDA %*% eta/ c(t(eta) %*% LAMBDA %*% eta)
    hat_betaE_lambda = hat_betaE_lambda[which(hat_betaE_lambda != 0)]
    Aeta = hat_betaE_lambda - LAMBDA %*% eta %*% t(eta) %*% hat_betaE_lambda/c(t(eta) %*% LAMBDA %*% eta)

    # calculate upper and lower bound
    temp_support = support_update(se = se, eta = eta, LAMBDA = LAMBDA, Aeta = Aeta, Qn = Qn)
    lower = temp_support[["lowerbond"]]
    upper = temp_support[["upperbond"]]

    require(calculus)
    numerator = integral(my_function, bounds = list(b = c(betaEjhat - 15*sqrt(sigmasq1/n), betaEjhat)),
                         params = list(sigmasq1 = sigmasq1, n = n, eta = eta,
                                       p1 = p1, p2 = p2, null_value = null_value, Var_etabeta = Var_etabeta,lower = lower, upper = upper))

    denominator_p2 = integral(my_function, bounds = list(b = c(betaEjhat, betaEjhat + 15*sqrt(sigmasq1/n))),
                              params = list(sigmasq1 = sigmasq1, n = n, eta = eta,
                                            p1 = p1, p2 = p2, null_value = null_value, Var_etabeta = Var_etabeta,lower = lower, upper = upper))


    prop_null_value = numerator$value/(numerator$value + denominator_p2$value)

    return(prop_null_value)
  }

  # function that search betaE,j value using half split
  betaEj_search = function(range0, max_tol, max_iterate) {
    # check the pivot value under the GEE point estimate
    prop_initial = pivot_prop_value(betaEjhat)

    if(is.na(prop_initial)) {
      print("the pivot is undefined even when the null value equals to the point estimate, the denominator of pivot is 0." )
      return(list(low_bound = NA,
                  up_bound = NA,
                  prop_up = NA,
                  prop_low = NA,
                  message_low_CI = NA,
                  message_up_CI = NA))
    }

    require(zoo)

    # check whether the trend is increasing first then decreasing
    # or monotone decreasing
    check_trend_downward <- function(x, k = 1) {
      if (length(x) < k + 2) return(FALSE)

      # Smooth with moving average
      smoothed <- rollmean(x, k = k, fill = NA)
      smoothed <- na.omit(smoothed)

      if (length(smoothed) < 3) return(FALSE)

      # Check monotone decreasing
      if (all(diff(smoothed) <= 0)) return(TRUE)

      # Check increasing then decreasing
      peak_idx <- which.max(smoothed)
      if (peak_idx == 1 || peak_idx == length(smoothed)) return(FALSE)

      left <- smoothed[1:peak_idx]
      right <- smoothed[peak_idx:length(smoothed)]

      increasing_before <- mean(diff(left)) > 0
      decreasing_after <- mean(diff(right)) < 0

      return(increasing_before && decreasing_after)
    }

    # check whether the trend is decreasing then increasing or monotone increasing
    check_trend_upward <- function(x, k = 1) {
      if (length(x) < k + 2) return(FALSE)

      # Smooth with moving average
      smoothed <- rollmean(x, k = k, fill = NA)
      smoothed <- na.omit(smoothed)

      if (length(smoothed) < 3) return(FALSE)

      # Check monotone increasing
      if (all(diff(smoothed) >= 0)) return(TRUE)

      # Check decreasing then increasing
      valley_idx <- which.min(smoothed)
      if (valley_idx == 1 || valley_idx == length(smoothed)) return(FALSE)

      left <- smoothed[1:valley_idx]
      right <- smoothed[valley_idx:length(smoothed)]

      decreasing_before <- mean(diff(left)) < 0
      increasing_after <- mean(diff(right)) > 0

      return(decreasing_before && increasing_after)
    }

    var = E[which(ej != 0)]


    #######################
    # find the lower bound#
    #######################
    range = range0

    if(prop_initial > upp & prop_initial < 1) {
      bottom = betaEjhat
      top = betaEjhat + range

      # check the pivot value under top
      top_prop = pivot_prop_value(top)
      while(top_prop > upp & top_prop < 1){
        range = range*1.25
        top = betaEjhat + range
        top_prop = pivot_prop_value(top)
      }

      #print(list(var = var, top = top, top_prop = top_prop, "the upper range is updating; Lower CI"))

    } else {
      bottom = betaEjhat - range
      top = betaEjhat
      bottom_prop = pivot_prop_value(bottom)

      # store bottom and its prop value for tend analysis
      bottom_value = c()
      bottom_prop_value = c()

      while(bottom_prop < upp & !is.na(bottom_prop) & bottom_prop != 0) {
        range = range*1.25
        bottom = betaEjhat - range
        bottom_prop = pivot_prop_value(bottom)

        bottom_value = c(bottom_value, bottom)
        bottom_prop_value = c(bottom_prop_value, bottom_prop)

        #print(list(var = var, bottom = bottom, bottom_prop = bottom_prop, "the lower range is updating; Lower CI"))
      }

      trends = check_trend_downward(bottom_prop_value)

      if(trends == T) {
        bottom = bottom_value[which.max(bottom_prop_value)]
        bottom_prop = bottom_prop_value[which.max((bottom_prop_value))]
        range = betaEjhat - bottom
      }

      #print(list(var = var, betaEjhat = betaEjhat, bottom = bottom, bottom_prop = bottom_prop, top = top, "first adjustment for bottom; Lower CI"))

      if(trends == T & bottom_prop < (upp - (1 - upp)/2)) {
        # store bottom and its prop value for tend analysis
        bottom_value = c()
        bottom_prop_value = c()
        j = 0

        while(bottom_prop < (2*upp - 1) & j < 20) {
          range = range/1.15
          bottom = betaEjhat - range
          bottom_prop = pivot_prop_value(bottom)

          bottom_value = c(bottom_value, bottom)
          bottom_prop_value = c(bottom_prop_value, bottom_prop)

          #print(list(var = var, bottom = bottom, bottom_prop = bottom_prop, "the lower range is updating; Lower CI"))

          j = j+ 1
        }

        # might see increasing then decreasing trend or monotone decreasing trend
        # because we come back to reasonable range
        # when get close to betaEj_hat the pivot value will decrease
        trends = check_trend_downward(bottom_prop_value)

        if(trends == T) {
          bottom = bottom_value[which.max(bottom_prop_value)]
          bottom_prop = bottom_prop_value[which.max((bottom_prop_value))]
        }

      } else {
        while(bottom_prop > (upp + (1 - upp)/2) | is.na(bottom_prop) | bottom_prop == 0) {
          range = range/1.15
          bottom = betaEjhat - range
          bottom_prop = pivot_prop_value(bottom)


          #print(list(var = var, bottom = bottom, bottom_prop = bottom_prop, "the lower range is updating; Lower CI"))
        }
      }




      #print(list(var = var, betaEjhat = betaEjhat, bottom = bottom, bottom_prop = bottom_prop, top = top, "second adjustment for bottom; Lower CI"))

    }


    i = 0
    prop = 0
    message_low_CI = NULL
    while((abs(upp - prop) > max_tol) & (i < max_iterate)) {
      if (abs(top - bottom) < 1e-4) {
        message_low_CI = "stops lower bound calculation due to the little change for each loop"
        break
      }

      if(i == 0) {
        cal = (top + bottom)/2
      } else {
        if(prop > upp) {
          bottom = cal
          cal = (bottom + top)/2
        } else {
          top = cal
          cal = (bottom + top)/2
        }
      }

      prop = pivot_prop_value(cal)
      i = i + 1
    }

    # the lower bound value
    low_bound = cal
    prop_up = prop
    if(i == max_iterate) {message_low_CI = "stops due to iteration limitation"
    } else if(!is.null(message_low_CI)) {message_low_CI = message_low_CI
    } else {message_low_CI = "Normal"}



    #######################
    # find the upper bound#
    #######################
    ## reset top and bottom
    prop = 0
    range = range0
    if(prop_initial < lowp) {
      bottom = betaEjhat - range
      top = betaEjhat

      # check the pivot value under bottom
      bottom_prop = pivot_prop_value(bottom)
      while(bottom_prop < lowp) {
        range = range*1.25
        bottom = betaEjhat - range
        bottom_prop = pivot_prop_value(bottom)

        print(list(var = var, bottom = bottom, bottom_prop = bottom_prop, "the lower range is updating; Upper CI"))

      }


    } else {
      bottom = betaEjhat
      top = betaEjhat + range

      # check the pivot value under top
      top_prop = pivot_prop_value(top)

      # store bottom and its prop value for tend analysis
      top_value = c()
      top_prop_value = c()

      while(top_prop > lowp & top_prop < 1 & !is.na(top_prop)) {
        range = range * 1.25
        top = betaEjhat + range
        top_prop = pivot_prop_value(top)

        top_value = c(top_value, top)
        top_prop_value = c(top_prop_value, top_prop)

        print(list(var = var, top = top, top_prop = top_prop, "the top range is updating; Upper CI"))

      }


      trends = check_trend_upward(top_prop_value)
      #print(list(top_prop_value = top_prop_value, trends = trends))

      if(trends == T) {
        top = top_value[which.min(top_prop_value)]
        top_prop = top_prop_value[which.min(top_prop_value)]
        range = top - betaEjhat
      }

      print(list(var = var, betaEjhat = betaEjhat, bottom = bottom, top = top, top_prop = top_prop, "first adjustment for top; Upper CI"))

      if(trends == T & top_prop > 1.5*lowp) {
        # store bottom and its prop value for tend analysis
        top_value = c()
        top_prop_value = c()
        j = 0

        while(top_prop > 1.5*lowp & j < 20) {
          range = range/1.15
          top = betaEjhat + range
          top_prop = pivot_prop_value(top)

          top_value = c(top_value, top)
          top_prop_value = c(top_prop_value, top_prop)

          print(list(var = var, top = top, top_prop = top_prop, "the top range is updating; Upper CI"))

          j = j+ 1
        }
        # might see decreasing then increasing trend or monotone increasing trend
        # because we come back to reasonable range
        # when get close to betaEj_hat the pivot value will increase
        trends = check_trend_upward(top_prop_value)

        if(trends == T) {
          top = top_value[which.min(top_prop_value)]
          top_prop = top_prop_value[which.min(top_prop_value)]
        }

      } else {
        while(is.na(top_prop) | top_prop < lowp/2 | top_prop == 1) {
          range = range/1.15
          top = betaEjhat + range
          top_prop = pivot_prop_value(top)

          print(list(var = var, top = top, top_prop = top_prop, "the top range is updating; Upper CI"))
        }
      }



      print(list(var = var, betaEjhat = betaEjhat, bottom = bottom, top = top, top_prop = top_prop, "second adjustment for top; Upper CI"))


    }

    i = 0
    message_up_CI = NULL
    while((abs(lowp - prop) > max_tol) & (i < max_iterate)) {
      if (abs(top - bottom) < 1e-4) {
        message_up_CI = "stops upper bound calculation due to the little change for each loop"
        break
      }

      if(i == 0) {
        cal = (top + bottom)/2
      } else {
        if(prop > lowp) {
          bottom = cal
          cal = (bottom + top)/2
        } else {
          top = cal
          cal = (bottom + top)/2
        }
      }
      prop = pivot_prop_value(cal)
      i = i + 1
    }

    # the upper bound value
    up_bound = cal
    prop_low = prop
    if(i == max_iterate) {message_up_CI = "stops due to iteration limitation"
    } else if (!is.null(message_up_CI)) {message_up_CI = message_up_CI
    } else {message_up_CI = "Normal"}


    return(list(low_bound = low_bound,
                up_bound = up_bound,
                prop_up = prop_up,
                prop_low = prop_low,
                message_low_CI = message_low_CI,
                message_up_CI = message_up_CI))
  }

  CI = betaEj_search(sqrt(sigmasq1/n), max_tol, max_iterate)
  # CI = betaEj_search(sqrt(sigmasq1/n)*5, 10^{-2}, 10^{2})
  low = CI$low_bound
  up = CI$up_bound

  ############################
  # now calculate the p value#
  ############################

  prop0 = pivot_prop_value(0)
  if(prop0 > 1|is.na(prop0)) {pvalue = NA} else {pvalue = min(prop0, 1 - prop0)*2}

  return(list(E = E[which(ej == 1)],
              GEE_est = betaEjhat,
              post_beta = select_E[["postbeta"]][which(ej == 1)],
              pvalue = pvalue,
              lowCI = low,
              upperCI = up,
              prop_low = CI$prop_low,
              prop_up = CI$prop_up,
              message_low_CI = CI$message_low_CI,
              message_up_CI = CI$message_up_CI
  ))

  # Output:
  # E: the selected variables for which CI is calculated
  # GEE_est: the GEE estimate for this predictor
  # post_beta: the post selection true value for this predictor if simulation is conducted
  # pvalue: the p value
  # lowCI: the lower bound of the confidence interval
  # upperCI: the upper bound of the confidence interval
  # prop_low: the true corresponding pivot value for the lower bound. For example, if it's 90% CI, this
  #          value will close to 0.05
  # prop_up: the true corresponding pivot value for the upper bound. For example, if it's 90% CI, this
  #          value will close to 0.95
}




























