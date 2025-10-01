#' Conditional Mixed Process Models
#'
#' Estimate conditional mixed process (CMP) models using maximum likelihood.
#' CMP models allow for systems of equations with different model types.
#'
#' @param formula A list of formulas or a single formula object specifying the model equations
#' @param data A data frame containing the variables
#' @param indicators A character vector or list specifying the model type for each equation.
#'   Available types: "continuous", "probit", "oprobit", "left", "right", "interval"
#' @param weights Optional weights vector
#' @param subset Optional subset specification
#' @param na.action How to handle missing values
#' @param method Estimation method ("ML" for maximum likelihood)
#' @param start Optional starting values
#' @param control Control parameters for optimization
#' @param ... Additional arguments
#'
#' @details
#' The CMP model allows estimation of systems of equations where each equation
#' can follow different model types:
#' \itemize{
#'   \item "continuous": Continuous outcomes (normal/tobit-style)
#'   \item "probit": Binary probit model
#'   \item "oprobit": Ordered probit model
#'   \item "left": Left-censored (Tobit) model with user-definable bounds
#'   \item "right": Right-censored model with user-definable bounds
#'   \item "interval": Interval regression with general lower/upper bounds
#' }
#'
#' @return An object of class "cmp" containing:
#' \itemize{
#'   \item coefficients: Estimated coefficients
#'   \item vcov: Variance-covariance matrix for the full parameter vector
#'   \item loglik: Log-likelihood value at the optimum
#'   \item fitted.values: Linear predictors for each equation
#'   \item residuals: Residuals for continuous-type equations
#'   \item call: Original function call
#'   \item formula: Model formulas
#'   \item data: Data used in estimation
#'   \item indicators: Model indicators
#'   \item nobs: Number of observations contributing to the likelihood
#'   \item npar: Number of estimated parameters
#' }
#'
#' @export
cmp <- function(formula, data, indicators, weights = NULL, subset = NULL,
                na.action = na.omit, method = "ML", start = NULL,
                control = list(), ...) {

  cl <- match.call()

  inputs <- validate_cmp_inputs(formula, indicators, data)
  formula_list <- inputs$formula
  data_full <- inputs$data
  indicators <- inputs$indicators

  weights_expr <- substitute(weights)
  weights_vec <- NULL
  if (!missing(weights) && !is.null(weights_expr)) {
    weights_vec <- eval(weights_expr, data_full, parent.frame())
    if (length(weights_vec) != nrow(data_full)) {
      stop("weights must have length equal to number of rows in data")
    }
    weights_vec <- as.numeric(weights_vec)
  }

  subset_idx <- NULL
  if (!is.null(subset)) {
    subset_expr <- substitute(subset)
    subset_idx <- eval(subset_expr, data_full, parent.frame())
    if (!is.logical(subset_idx)) {
      subset_idx <- as.logical(subset_idx)
    }
    subset_idx[is.na(subset_idx)] <- FALSE
    data_use <- data_full[subset_idx, , drop = FALSE]
    if (!is.null(weights_vec)) {
      weights_vec <- weights_vec[subset_idx]
    }
  } else {
    data_use <- data_full
    if (!is.null(weights_vec)) {
      weights_vec <- weights_vec
    }
  }

  prep <- prepare_cmp_models(formula_list, data_use, indicators, na.action)
  models <- prep$models
  n_obs <- prep$n_obs

  if (!is.null(weights_vec) && length(weights_vec) != n_obs) {
    stop("weights vector must align with the estimation sample")
  }

  param_info <- build_param_info(models)

  if (param_info$total_par == 0L) {
    stop("No parameters to estimate; check model specification")
  }

  if (is.null(start)) {
    start <- cmp_start_values(models, param_info, method = "ols")
  } else {
    if (length(start) != param_info$total_par) {
      stop(
        "start vector must have length ", param_info$total_par,
        " (received ", length(start), ")"
      )
    }
  }

  default_control <- list(fnscale = -1, maxit = 1000, reltol = 1e-8)
  control <- utils::modifyList(default_control, control, keep.null = TRUE)

  opt <- optim(
    par = start,
    fn = cmp_loglik,
    gr = cmp_gradient,
    models = models,
    param_info = param_info,
    weights = weights_vec,
    method = "BFGS",
    control = control,
    hessian = TRUE
  )

  names(opt$par) <- param_info$names

  total_coef <- param_info$total_coef
  total_sigma <- param_info$total_sigma
  total_cut <- param_info$total_cut

  coefficients <- opt$par[seq_len(total_coef)]
  names(coefficients) <- param_info$names[seq_len(total_coef)]

  sigma_par <- if (total_sigma) opt$par[total_coef + seq_len(total_sigma)] else numeric()
  corr_par <- if (param_info$corr$count) opt$par[(total_coef + total_sigma + total_cut) + seq_len(param_info$corr$count)] else numeric()
  cut_par <- if (total_cut) opt$par[(total_coef + total_sigma) + seq_len(total_cut)] else numeric()

  sigma_vec <- rep(1, param_info$n_eq)
  if (total_sigma) {
    sigma_vec[param_info$sigma$eq] <- exp(sigma_par)
  }

  rho_vec <- if (length(corr_par)) tanh(corr_par) else numeric()

  Sigma <- diag(sigma_vec^2)
  if (length(rho_vec)) {
    for (k in seq_len(param_info$corr$count)) {
      idx <- param_info$corr$pairs[k, ]
      i <- idx[1]
      j <- idx[2]
      val <- rho_vec[k] * sigma_vec[i] * sigma_vec[j]
      Sigma[i, j] <- Sigma[j, i] <- val
    }
  }

  cut_list <- vector("list", param_info$n_eq)
  if (total_cut) {
    for (i in seq_len(param_info$n_eq)) {
      count <- param_info$cut$counts[i]
      if (count) {
        offset <- param_info$cut$offsets[i]
        alpha <- cut_par[offset + seq_len(count)]
        cut_list[[i]] <- build_cutpoints(alpha)
      }
    }
  }

  vcov_matrix <- tryCatch({
    solve(-opt$hessian)
  }, error = function(e) {
    warning("Hessian is not invertible. Using pseudo-inverse.")
    MASS::ginv(-opt$hessian)
  })
  dimnames(vcov_matrix) <- list(param_info$names, param_info$names)

  fitted_values <- cmp_fitted(coefficients, models, param_info)
  residuals <- cmp_residuals(models, fitted_values)

  mask_mat <- vapply(models, function(m) as.integer(m$mask), integer(n_obs))
  active_obs <- rowSums(mask_mat) > 0
  n_active <- sum(active_obs)

  result <- list(
    coefficients = coefficients,
    vcov = vcov_matrix,
    loglik = opt$value,
    fitted.values = fitted_values,
    residuals = residuals,
    call = cl,
    formula = formula_list,
    data = data_use,
    indicators = vapply(models, function(m) m$indicator, character(1)),
    models = models,
    param_info = param_info,
    cutpoints = cut_list,
    Sigma = Sigma,
    sigma = sigma_vec,
    sigma_par = sigma_par,
    corr = rho_vec,
    corr_par = corr_par,
    nobs = n_active,
    npar = length(opt$par),
    convergence = opt$convergence,
    weights = weights_vec,
    optim = opt,
    mask = mask_mat,
    control = control
  )

  class(result) <- "cmp"
  result
}

#' Create starting values for CMP estimation
#'
#' @param models List of model specifications
#' @param param_info Parameter bookkeeping object
#' @param method Method used for starting values
#' @keywords internal
cmp_start_values <- function(models, param_info, method = "ols") {
  n_eq <- param_info$n_eq
  start_coef <- numeric(param_info$total_coef)
  start_sigma <- if (param_info$total_sigma) numeric(param_info$total_sigma) else numeric()
  start_cut <- if (param_info$total_cut) numeric(param_info$total_cut) else numeric()
  sigma_positions <- if (param_info$total_sigma) setNames(seq_len(param_info$total_sigma), param_info$sigma$eq) else integer()

  tol <- sqrt(.Machine$double.eps)

  for (i in seq_len(n_eq)) {
    model <- models[[i]]
    mask <- model$mask
    X <- model$X[mask, , drop = FALSE]
    offsets <- param_info$coef$offsets
    coef_idx <- offsets[i] + seq_len(model$n_coef)

    if (!any(mask)) {
      start_coef[coef_idx] <- 0
      next
    }

    if (model$indicator %in% c("continuous", "left", "right", "interval")) {
      y <- as.numeric(model$y[mask])
      if (model$indicator == "left") {
        limit <- model$lower[mask]
        uncensored <- y > limit + tol
        if (sum(uncensored) >= ncol(X)) {
          fit <- tryCatch(stats::lm.fit(X[uncensored, , drop = FALSE], y[uncensored]), error = function(e) NULL)
        } else {
          fit <- tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
        }
      } else if (model$indicator == "right") {
        limit <- model$upper[mask]
        uncensored <- y < limit - tol
        if (sum(uncensored) >= ncol(X)) {
          fit <- tryCatch(stats::lm.fit(X[uncensored, , drop = FALSE], y[uncensored]), error = function(e) NULL)
        } else {
          fit <- tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
        }
      } else if (model$indicator == "interval") {
        lower <- model$lower[mask]
        upper <- model$upper[mask]
        midpoint <- ifelse(is.finite(lower) & is.finite(upper), 0.5 * (lower + upper),
                           ifelse(is.finite(lower), lower, upper))
        fit <- tryCatch(stats::lm.fit(X, midpoint), error = function(e) NULL)
      } else {
        fit <- tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
      }

      if (!is.null(fit)) {
        start_coef[coef_idx] <- fit$coefficients
        if (model$has_sigma) {
          sigma_pos <- sigma_positions[as.character(i)]
          if (!is.na(sigma_pos)) {
            resid <- y - X %*% fit$coefficients
            sd_val <- stats::sd(resid, na.rm = TRUE)
            if (!is.finite(sd_val) || sd_val <= 0) {
              sd_val <- 1
            }
            start_sigma[sigma_pos] <- log(sd_val)
          }
        }
      } else {
        start_coef[coef_idx] <- 0
        if (model$has_sigma) {
          sigma_pos <- sigma_positions[as.character(i)]
          if (!is.na(sigma_pos)) {
            start_sigma[sigma_pos] <- 0
          }
        }
      }

    } else if (model$indicator == "probit") {
      mf_valid <- model$mf[mask, , drop = FALSE]
      formula <- stats::formula(model$terms)
      fit <- tryCatch(stats::glm(formula, data = mf_valid, family = stats::binomial(link = "probit")), error = function(e) NULL)
      if (!is.null(fit)) {
        beta_hat <- stats::coef(fit)
        beta_vec <- setNames(numeric(length(coef_idx)), colnames(model$X))
        beta_vec[names(beta_hat)] <- beta_hat
        start_coef[coef_idx] <- beta_vec
      } else {
        start_coef[coef_idx] <- 0
      }

    } else if (model$indicator == "oprobit") {
      mf_valid <- model$mf[mask, , drop = FALSE]
      formula <- stats::formula(model$terms)
      fit <- tryCatch(MASS::polr(formula, data = mf_valid, method = "probit", Hess = FALSE), error = function(e) NULL)
      if (!is.null(fit)) {
        beta_hat <- stats::coef(fit)
        beta_vec <- setNames(numeric(length(coef_idx)), colnames(model$X))
        beta_vec[names(beta_hat)] <- beta_hat
        start_coef[coef_idx] <- beta_vec
        zeta <- fit$zeta
        count <- length(zeta)
        offset <- param_info$cut$offsets[i]
        if (count) {
          alpha <- numeric(count)
          alpha[1L] <- zeta[1L]
          if (count > 1L) {
            diffs <- diff(zeta)
            diffs[diffs <= 0] <- tol
            alpha[-1L] <- log(diffs)
          }
          start_cut[offset + seq_len(count)] <- alpha
        }
      } else {
        start_coef[coef_idx] <- 0
      }
    } else {
      start_coef[coef_idx] <- 0
    }
  }

  c(start_coef, start_sigma, start_cut, rep(0, param_info$corr$count))
}
