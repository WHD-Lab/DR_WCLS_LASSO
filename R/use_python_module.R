#' venv_config
#' 
#' @description
#' This function helps to build a Python virtual environment in R
#' 
#' 
#' @export

venv_config = function(salt = "") {
  install_python("3.9:latest")
  venv_info = get_venv_info(salt)
  repo_dir <- clone_selective_inference()
  
  virtualenv_create(venv_info$hash, version = venv_info$deps[1])
  virtualenv_install(venv_info$hash, venv_info$deps[-1])
  virtualenv_install(venv_info$hash, "numpy==1.22.4")
  
  virtualenv_install(venv_info$hash,
                     c("git+https://github.com/regreg/regreg.git"),
                     pip_options = "--no-build-isolation")
  
  virtualenv_install(venv_info$hash, c("scikit-learn", "nose"))
  
  virtualenv_install(
    venv_info$hash,
    repo_dir,
    options = c("-e"),
    pip_options = "--no-build-isolation"
  )
  
  # unlisted deps??? If any other unlisted deps, list here
  virtualenv_install(venv_info$hash, c("pandas", "mpmath"))
  
  use_virtualenv(venv_info$hash, required = TRUE)
  
  return(venv_info)
}

hash_vec <- function(..., length = 8) {
    args <- list(...)
    # Convert all vectors to character and flatten
    flat <- unlist(lapply(args, as.character), use.names = FALSE)
    substr(digest::digest(flat, algo = "md5"), 1, length)
}

get_venv_info = function(salt = "") {
    py_version = "3.9"
    # pkg_deps = c(
    #          "numpy==1.22.4",
    #          "cython>=0.18",
    #          "scipy>=0.16.0",
    #          "joblib>=0.15.1",
    #          "scikit-learn"
    #         )
    pkg_deps = c("cython>=0.18",
                 "numpy>=1.7.1",
                 "scipy>=0.16.0",
                 "joblib>=0.15.1")
    
    deps = c(py_version, pkg_deps)
    
    return(
        list(
            deps = deps,
            hash = hash_vec(deps, salt)
        )
    )
    
}

clone_selective_inference <- function(branch = "refactor_names") {
    tmp <- tempfile("selective-inference-")
    
    # Step 1 — git clone specific branch
    status <- system2(
        "git",
        c(
            "clone",
            "--branch", branch,
            "--single-branch",
            "https://github.com/snigdhagit/selective-inference.git",
            tmp
        ),
        stdout = TRUE,
        stderr = TRUE
    )
    
    if (!is.null(attr(status, "status"))) {
        stop("git clone failed:\n", paste(status, collapse = "\n"))
    }
    
    # Step 2 — init submodules
    status2 <- system2(
        "git",
        c("-C", tmp, "submodule", "update", "--init"),
        stdout = TRUE,
        stderr = TRUE
    )
    
    if (!is.null(attr(status2, "status"))) {
        stop("git submodule update failed:\n", paste(status2, collapse = "\n"))
    }
    
    tmp
}
