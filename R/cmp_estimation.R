#' Log-likelihood function for CMP models
#'
#' @param par Parameter vector
#' @param models List of model specifications
#' @param param_info Parameter bookkeeping object
#' @param weights Optional weights
#' @keywords internal
cmp_loglik <- function(par, models, param_info, weights = NULL) {

  n_eq <- param_info$n_eq
  total_coef <- param_info$total_coef
  total_sigma <- param_info$total_sigma
  total_cut <- param_info$total_cut

  coef_par <- par[seq_len(total_coef)]
  cursor <- total_coef

  sigma_par <- if (total_sigma) par[cursor + seq_len(total_sigma)] else numeric()
  cursor <- cursor + total_sigma

  cut_par <- if (total_cut) par[cursor + seq_len(total_cut)] else numeric()
  cursor <- cursor + total_cut

  corr_par <- if (param_info$corr$count) par[cursor + seq_len(param_info$corr$count)] else numeric()

  sigma_vec <- rep(1, n_eq)
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

  chol_check <- safe_chol_logdet(Sigma)
  if (!chol_check$ok) {
    return(-Inf)
  }

  n_obs <- length(models[[1]]$mask)
  if (is.null(weights)) {
    weights <- rep(1, n_obs)
  } else {
    if (length(weights) != n_obs) {
      stop("weights vector must align with data used in estimation")
    }
    weights <- as.numeric(weights)
    weights[is.na(weights)] <- 0
  }

  coef_list <- vector("list", n_eq)
  linear_preds <- vector("list", n_eq)
  for (i in seq_len(n_eq)) {
    start <- param_info$coef$offsets[i] + 1L
    end <- param_info$coef$offsets[i + 1L]
    beta <- coef_par[start:end]
    coef_list[[i]] <- beta

    mask <- models[[i]]$mask
    lp <- rep(NA_real_, n_obs)
    if (any(mask)) {
      lp[mask] <- as.numeric(models[[i]]$X[mask, , drop = FALSE] %*% beta)
    }
    linear_preds[[i]] <- lp
  }

  cut_list <- vector("list", n_eq)
  if (total_cut) {
    for (i in seq_len(n_eq)) {
      count <- param_info$cut$counts[i]
      if (!count) {
        next
      }
      offset <- param_info$cut$offsets[i]
      alpha <- cut_par[offset + seq_len(count)]
      cut_list[[i]] <- build_cutpoints(alpha)
    }
  }

  ll_total <- 0
  active_any <- FALSE
  for (obs in seq_len(n_obs)) {
    active_eq <- which(vapply(models, function(m) m$mask[obs], logical(1)))
    if (!length(active_eq)) {
      next
    }
    active_any <- TRUE
    ll_obs <- cmp_obs_loglik(
      obs = obs,
      active_eq = active_eq,
      models = models,
      linear_preds = linear_preds,
      Sigma = Sigma,
      sigma_vec = sigma_vec,
      cut_list = cut_list
    )
    if (!is.finite(ll_obs)) {
      return(-Inf)
    }
    w <- weights[obs]
    if (!is.finite(w) || w == 0) {
      next
    }
    ll_total <- ll_total + w * ll_obs
  }

  if (!active_any) {
    return(-Inf)
  }

  if (!is.finite(ll_total)) {
    return(-Inf)
  }

  ll_total
}

#' Observation-level contribution to the CMP log-likelihood
#'
#' @param obs Observation index
#' @param active_eq Active equation indices
#' @param models Model list
#' @param linear_preds Linear predictor list
#' @param Sigma Covariance matrix for latent errors
#' @param sigma_vec Vector of error standard deviations
#' @param cut_list Cut point list for ordered probit equations
#' @keywords internal
cmp_obs_loglik <- function(obs, active_eq, models, linear_preds, Sigma, sigma_vec, cut_list) {
  tol <- sqrt(.Machine$double.eps)

  m <- length(active_eq)
  lower <- rep(-Inf, m)
  upper <- rep(Inf, m)
  value_flags <- logical(m)

  for (k in seq_len(m)) {
    idx <- active_eq[k]
    model <- models[[idx]]
    mu <- linear_preds[[idx]][obs]
    if (is.na(mu)) {
      return(-Inf)
    }
    sigma_i <- sigma_vec[idx]
    indicator <- model$indicator

    if (indicator == "continuous") {
      y_val <- as.numeric(model$y[obs])
      resid <- (y_val - mu) / sigma_i
      lower[k] <- upper[k] <- resid
      value_flags[k] <- TRUE

    } else if (indicator == "probit") {
      y_val <- as.numeric(model$y[obs])
      if (y_val == 1) {
        lower[k] <- -mu
        upper[k] <- Inf
      } else {
        lower[k] <- -Inf
        upper[k] <- -mu
      }

    } else if (indicator == "left") {
      y_val <- as.numeric(model$y[obs])
      limit <- model$lower[obs]
      if (!is.finite(limit)) {
        limit <- -Inf
      }
      if (y_val <= limit + tol) {
        lower[k] <- -Inf
        upper[k] <- (limit - mu) / sigma_i
      } else {
        resid <- (y_val - mu) / sigma_i
        lower[k] <- upper[k] <- resid
        value_flags[k] <- TRUE
      }

    } else if (indicator == "right") {
      y_val <- as.numeric(model$y[obs])
      limit <- model$upper[obs]
      if (!is.finite(limit)) {
        limit <- Inf
      }
      if (y_val >= limit - tol) {
        lower[k] <- (limit - mu) / sigma_i
        upper[k] <- Inf
      } else {
        resid <- (y_val - mu) / sigma_i
        lower[k] <- upper[k] <- resid
        value_flags[k] <- TRUE
      }

    } else if (indicator == "interval") {
      lower_y <- model$lower[obs]
      upper_y <- model$upper[obs]
      if (is.finite(lower_y)) {
        lower[k] <- (lower_y - mu) / sigma_i
      } else {
        lower[k] <- -Inf
      }
      if (is.finite(upper_y)) {
        upper[k] <- (upper_y - mu) / sigma_i
      } else {
        upper[k] <- Inf
      }
      if (is.finite(lower_y) && is.finite(upper_y) && abs(upper_y - lower_y) <= tol) {
        resid <- (lower_y - mu) / sigma_i
        lower[k] <- upper[k] <- resid
        value_flags[k] <- TRUE
      }

    } else if (indicator == "oprobit") {
      cat_idx <- as.integer(model$y[obs])
      cuts <- cut_list[[idx]]
      lower_cut <- if (cat_idx == 1L) -Inf else cuts[cat_idx - 1L]
      upper_cut <- if (cat_idx > length(cuts)) Inf else cuts[cat_idx]
      lower[k] <- lower_cut - mu
      upper[k] <- upper_cut - mu

    } else {
      stop("Unhandled indicator type: ", indicator)
    }
  }

  Sigma_sub <- Sigma[active_eq, active_eq, drop = FALSE]
  value_idx <- which(abs(upper - lower) <= tol)
  interval_idx <- setdiff(seq_len(m), value_idx)
  sigma_active <- sigma_vec[active_eq]
  log_jacobian <- if (any(value_flags)) sum(log(sigma_active[value_flags])) else 0

  ll <- 0

  if (length(value_idx)) {
    e_vals <- lower[value_idx]
    Sigma_vv <- Sigma_sub[value_idx, value_idx, drop = FALSE]
    density_v <- mvtnorm::dmvnorm(
      x = e_vals,
      mean = rep(0, length(value_idx)),
      sigma = Sigma_vv,
      log = TRUE
    )
    ll <- density_v - log_jacobian
  }

  if (length(interval_idx)) {
    if (length(value_idx)) {
      Sigma_vv <- Sigma_sub[value_idx, value_idx, drop = FALSE]
      Sigma_vi <- Sigma_sub[value_idx, interval_idx, drop = FALSE]
      Sigma_iv <- t(Sigma_vi)
      Sigma_ii <- Sigma_sub[interval_idx, interval_idx, drop = FALSE]
      Sigma_vv_inv <- tryCatch(solve(Sigma_vv), error = function(e) NULL)
      if (is.null(Sigma_vv_inv)) {
        return(-Inf)
      }
      mu_cond <- drop(Sigma_iv %*% Sigma_vv_inv %*% lower[value_idx])
      cov_cond <- Sigma_ii - Sigma_iv %*% Sigma_vv_inv %*% Sigma_vi
      if (!safe_chol_logdet(cov_cond)$ok) {
        return(-Inf)
      }
      prob <- mvtnorm::pmvnorm(
        lower = lower[interval_idx],
        upper = upper[interval_idx],
        mean = mu_cond,
        sigma = cov_cond
      )
    } else {
      prob <- mvtnorm::pmvnorm(
        lower = lower,
        upper = upper,
        mean = rep(0, m),
        sigma = Sigma_sub
      )
    }
    prob_val <- as.numeric(prob)
    if (!is.finite(prob_val) || prob_val <= 0) {
      return(-Inf)
    }
    ll <- ll + log(prob_val)
  }

  ll
}

#' Gradient function for CMP models
#'
#' @param par Parameter vector
#' @param models List of model specifications
#' @param param_info Parameter bookkeeping object
#' @param weights Optional weights
#' @keywords internal
cmp_gradient <- function(par, models, param_info, weights = NULL) {
  grad_func <- function(x) cmp_loglik(x, models = models, param_info = param_info, weights = weights)
  numDeriv::grad(grad_func, par)
}

#' Calculate fitted values on the linear predictor scale
#'
#' @param coefficients Coefficient vector
#' @param models List of model specifications
#' @param param_info Parameter bookkeeping object
#' @keywords internal
cmp_fitted <- function(coefficients, models, param_info) {
  n_eq <- param_info$n_eq
  n_obs <- length(models[[1]]$mask)
  fitted <- matrix(NA_real_, nrow = n_obs, ncol = n_eq)
  colnames(fitted) <- vapply(models, function(m) m$eq_name, character(1))

  for (i in seq_len(n_eq)) {
    start <- param_info$coef$offsets[i] + 1L
    end <- param_info$coef$offsets[i + 1L]
    beta <- coefficients[start:end]
    mask <- models[[i]]$mask
    if (any(mask)) {
      fitted[mask, i] <- as.numeric(models[[i]]$X[mask, , drop = FALSE] %*% beta)
    }
  }
  fitted
}

#' Calculate residuals for CMP models
#'
#' @param models List of model specifications
#' @param fitted_values Fitted matrix on linear predictor scale
#' @keywords internal
cmp_residuals <- function(models, fitted_values) {
  n_obs <- nrow(fitted_values)
  n_eq <- length(models)
  residuals <- matrix(NA_real_, nrow = n_obs, ncol = n_eq)
  colnames(residuals) <- vapply(models, function(m) m$eq_name, character(1))
  tol <- sqrt(.Machine$double.eps)

  for (i in seq_len(n_eq)) {
    model <- models[[i]]
    indicator <- model$indicator
    mask <- model$mask

    if (indicator == "continuous") {
      y <- as.numeric(model$y)
      residuals[mask, i] <- y[mask] - fitted_values[mask, i]
    } else if (indicator == "left") {
      y <- as.numeric(model$y)
      lower <- model$lower
      uncensored <- mask & (y > lower + tol)
      residuals[uncensored, i] <- y[uncensored] - fitted_values[uncensored, i]
    } else if (indicator == "right") {
      y <- as.numeric(model$y)
      upper <- model$upper
      uncensored <- mask & (y < upper - tol)
      residuals[uncensored, i] <- y[uncensored] - fitted_values[uncensored, i]
    }
  }

  residuals
}
