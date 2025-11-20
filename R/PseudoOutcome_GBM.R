#' pseudo_outcome_generator_gbm
#'
#' Cross-fitting is applied when generate the target pseudo-outcome.
#' This function uses Gradient Boosting to conduct model training on folds then estimates \eqn{E[Y_{t+1}|H_t, A_t], E[A_t|H_t], E[A_t|S_t]} for the reserved fold. Then
#' the function returns a dataset with column named "yDR" for DR-WCLS pseudo outcome.
#'
#' @param fold number of folds to split when do corss-fitting
#' @param ID the name of column where participants' ID are stored
#' @param data dataset name
#' @param Ht a vector that contains column names of control variables
#' @param St a vector that contains column names of moderator variables; St should be a subset of Ht
#' @param At column names of treatment (At)
#' @param prob column names of \eqn{p_t(A_t = 1|H_t)}, the experiment design treatment probability
#' @param outcome column names of outcome variable
#' @param core_num number of cores will be used for calculation
#'
#' @return This function returns a dataset with pseudo outcome. It learns appropriate working models with Gradient Boosting and
#' generates pseudo outcome using the DR-WCLS.
#'
#' @examples
#'
#'  sim_data = generate_dataset(N = 1000, T = 40, P = 50, sigma_residual = 1.5, sigma_randint = 1.5, main_rand = 3, rho = 0.7,
#'   beta_logit = c(-1, 1.6 * rep(1/50, 50)), model = ~ state1 + state2 + state3 + state4,
#'   beta = matrix(c(-1, 1.7, 1.5, -1.3, -1),ncol = 1),
#'   theta1 = 0.8)
#'
#' Ht = unlist(lapply(1:50, FUN = function(X) paste0("state",X)))
#' St = unlist(lapply(1:25, FUN = function(X) paste0("state",X)))
#'
#' pseudo_outcome_generator_gbm(fold = 5,ID = "id", sim_data, Ht, St, "action", "outcome",core_num = 5)
#'
#' @import xgboost
#' @import parallel
#' @import dplyr
#' @import modelr
#'
#' @export




pseudo_outcome_generator_gbm = function(fold, ID, data, Ht, St, At, prob, outcome, core_num = NULL) {
  fold_ind = split_data(data[[ID]], fold = fold)
  MRT_gbm = ps_gradient_boosting(fold_indices = fold_ind, fold = fold, ID = ID,
                                 data = data, Ht = Ht, St = St, At = At, outcome = outcome, core_num)
  pseudo = pseudo_outcome_cal_gbm(MRT_gbm, At, prob, outcome)
  return(pseudo)
}


ps_gradient_boosting = function(fold_indices, fold, ID, data, Ht, St, At, outcome, core_num = NULL) {
  data_withpred = data.frame()

  expectation_cal = function(i) {
    require(xgboost)
    require(parallel)
    require(dplyr)

    reserve = data[data[[ID]] %in% fold_indices[[i]], ]
    train_fold = data[!(data[[ID]] %in% fold_indices[[i]]), ]

    fomula = as.formula(paste("~-1+", paste(paste0(Ht,":",At), collapse = "+")))
    X_test_gt =  data.matrix(modelr::model_matrix(reserve, fomula))
    colnames(X_test_gt) <- gsub(":", ".", colnames(X_test_gt))
    X_test_gt = cbind(reserve, X_test_gt)

    X_train_gt =  data.matrix(modelr::model_matrix(train_fold, fomula))
    colnames(X_train_gt) <- gsub(":", ".", colnames(X_train_gt))
    names_val = colnames(X_train_gt)
    X_train_gt = cbind(train_fold, X_train_gt)


    train_matrix = xgb.DMatrix(data = as.matrix(X_train_gt[, c(Ht, At, names_val)]), label = X_train_gt[[outcome]])
    test_matrix = xgb.DMatrix(data = as.matrix(X_test_gt[, c(Ht, At, names_val)]))

    params_reg = list(objective = "reg:squarederror", booster = "gbtree", eta = 0.1, max_depth = 6, nrounds = 500)
    params_class = list(objective = "binary:logistic", booster = "gbtree", eta = 0.1, max_depth = 6, nrounds = 500)

    ## gt(Ht,At)
    gt_gbm = xgb.train(params = params_reg, data = train_matrix, nrounds = params_reg$nrounds)
    reserve$gt_pred_gbm = predict(gt_gbm, test_matrix)

    ## gt(Ht, At = 1) and gt(Ht, At = 0)
    X_testfixAt1 = reserve[, c(Ht, At)]; X_testfixAt1[[At]] = 1
    X_test_gtfixAt1 =  data.matrix(modelr::model_matrix(X_testfixAt1, fomula))
    colnames(X_test_gtfixAt1) <- gsub(":", ".", colnames(X_test_gtfixAt1))
    X_testfixAt1 = cbind(X_testfixAt1, X_test_gtfixAt1)


    X_testfixAt0 = reserve[, c(Ht, At)]; X_testfixAt0[[At]] = 0
    X_test_gtfixAt0 =  data.matrix(modelr::model_matrix(X_testfixAt0, fomula))
    colnames(X_test_gtfixAt0) <- gsub(":", ".", colnames(X_test_gtfixAt0))
    X_testfixAt0 = cbind(X_testfixAt0, X_test_gtfixAt0)


    reserve$gtAt1_pred_gbm = predict(gt_gbm, xgb.DMatrix(as.matrix(X_testfixAt1)))
    reserve$gtAt0_pred_gbm = predict(gt_gbm, xgb.DMatrix(as.matrix(X_testfixAt0)))

    ## E[At|Ht] and E[At|St]
    ptHt_gbm = xgb.train(params = params_class, data = xgb.DMatrix(as.matrix(train_fold[, Ht]), label = train_fold[[At]]), nrounds = params_class$nrounds)
    ptSt_gbm = xgb.train(params = params_class, data = xgb.DMatrix(as.matrix(train_fold[, St]), label = train_fold[[At]]), nrounds = params_class$nrounds)

    reserve$ptHt_pred_gbm = predict(ptHt_gbm, xgb.DMatrix(as.matrix(reserve[, Ht])))
    reserve$ptSt_pred_gbm = predict(ptSt_gbm, xgb.DMatrix(as.matrix(reserve[, St])))


    return(reserve)
  }

  folds_list = 1:fold
  if (is.null(core_num)) { cl = parallel::makeCluster(detectCores()) } else { cl = parallel::makeCluster(core_num) }
  clusterExport(cl, varlist = ls(envir= environment()), envir = environment())
  results = parLapply(cl, folds_list, expectation_cal)
  stopCluster(cl)
  data_withpred = dplyr::bind_rows(results)

  colnames(data_withpred)[ncol(data_withpred)+1 - 5:1] = c("gt_pred_gbm","gtAt1_pred_gbm","gtAt0_pred_gbm", "ptHt", "ptSt")
  data_withpred$ptHtobs = data_withpred[,At]*data_withpred$ptHt + (1-data_withpred[,At])*(1-data_withpred$ptHt)
  data_withpred$ptStobs = data_withpred[,At]*data_withpred$ptSt + (1-data_withpred[,At])*(1-data_withpred$ptSt)

  return(data_withpred)
}

pseudo_outcome_cal_gbm = function(data_withpred, At, prob, outcome) {
  At = data_withpred[,At]
  ptSt = data_withpred$ptSt
  y = data_withpred[,outcome]
  gt = data_withpred$gt_pred_gbm
  gtAt1 = data_withpred$gtAt1_pred_gbm
  gtAt0 = data_withpred$gtAt0_pred_gbm

  wt = data_withpred$ptStobs / (data_withpred[,prob] * At + (1 - data_withpred[,prob]) * (1 - At))
  data_withpred$yDR = wt * (At - ptSt) * (y - gt) / (ptSt * (1 - ptSt)) + (gtAt1 - gtAt0)

  return(data_withpred)
}
