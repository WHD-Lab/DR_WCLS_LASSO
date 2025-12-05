run_simple_lasso_test_snigdha_fixed <- function(
        venv,
        n = 200L,
        p = 20L,
        s = 5L,
        sigma = 1,
        rho = 0.3,
        signal = 3,
        lambda_frac = 1.0
) {
    
    use_virtualenv(venv, required = TRUE)
    
    np        <- import("numpy", convert = FALSE)
    lasso_mod <- import("selectinf.randomized.lasso", convert = FALSE)$lasso
    
    # ---- MANUAL SIMULATION ----
    # covariance matrix
    Sigma <- rho ^ abs(outer(1:p, 1:p, "-"))  # Toeplitz autocorrelation
    
    # Cholesky factor
    L <- np$linalg$cholesky(np$array(Sigma))
    
    # design matrix X (n x p)
    Z <- np$random$randn(n, p)
    X_np <- Z$dot(L)
    
    # true beta
    beta <- np$zeros(p)
    beta[0:s] <- signal
    
    # noise and Y
    eps <- sigma * np$random$randn(n)
    Y_np <- X_np$dot(beta) + eps
    
    # ---- weight ----
    w <- lambda_frac * sigma * sqrt(2 * log(p))
    
    conv <- lasso_mod$gaussian(
        X_np, Y_np,
        feature_weights  = w * np$ones(p),
        randomizer_scale = sigma,
        ridge_term       = 0.0
    )
    
    signs   <- conv$fit()
    nonzero <- signs != 0L
    
    list(
        selected      = np$sum(nonzero),
        selected_idx  = which(py_to_r(nonzero)),
        true_support  = 1:s,
        beta_true     = beta
    )
}
