#' Print method for CMP objects
#'
#' @param x A cmp object
#' @param ... Additional arguments
#' @export
print.cmp <- function(x, ...) {
  cat("Conditional Mixed Process Model

")
  cat("Call:
")
  print(x$call)
  cat("
")

  cat("Model Types:
")
  for (i in seq_along(x$models)) {
    cat(sprintf("%s: %s
", x$models[[i]]$eq_name, x$models[[i]]$indicator))
  }
  cat("
")

  cat("Coefficients:
")
  print(x$coefficients)
  cat("
")

  if (length(x$sigma_par)) {
    cat("Log-standard deviations:
")
    sigma_names <- paste0("lnsigma_", vapply(x$models[x$param_info$sigma$eq], function(m) m$eq_name, character(1)))
    sigma_show <- setNames(x$sigma_par, sigma_names)
    print(sigma_show)
    cat("
")
  }

  if (length(x$corr_par)) {
    cat("Fisher-z correlations:
")
    corr_names <- apply(x$param_info$corr$pairs, 1, function(idx) {
      paste0("atanhrho_", x$models[[idx[1]]]$eq_name, "_", x$models[[idx[2]]]$eq_name)
    })
    corr_show <- setNames(x$corr_par, corr_names)
    print(corr_show)
    cat("
")
  }

  cat("Log-likelihood:", x$loglik, "
")
  cat("Observations:", x$nobs, "
")
  cat("Parameters:", x$npar, "
")

  if (x$convergence != 0) {
    cat("Warning: Optimization did not converge (code:", x$convergence, ")
")
  }

  invisible(x)
}

#' Summary method for CMP objects
#'
#' @param object A cmp object
#' @param ... Additional arguments
#' @export
summary.cmp <- function(object, ...) {
  coef_len <- length(object$coefficients)
  se <- sqrt(diag(object$vcov)[seq_len(coef_len)])
  names(se) <- names(object$coefficients)

  t_stat <- object$coefficients / se
  p_values <- 2 * stats::pnorm(-abs(t_stat))

  coef_table <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = se,
    `z value` = t_stat,
    `Pr(>|z|)` = p_values,
    check.names = FALSE
  )

  coef_table$Signif <- cut(
    p_values,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " " ),
    right = FALSE
  )

  aic_val <- -2 * object$loglik + 2 * object$npar
  bic_val <- -2 * object$loglik + log(object$nobs) * object$npar

  sigma_eq <- object$param_info$sigma$eq
  sigma_vals <- if (length(sigma_eq)) object$sigma[sigma_eq] else numeric()
  sigma_names <- if (length(sigma_eq)) {
    paste0("sigma_", vapply(object$models[sigma_eq], function(m) m$eq_name, character(1)))
  } else character()
  corr_names <- if (length(object$corr_par)) {
    apply(object$param_info$corr$pairs, 1, function(idx) {
      paste0("rho_", object$models[[idx[1]]]$eq_name, "_", object$models[[idx[2]]]$eq_name)
    })
  } else character()

  result <- list(
    call = object$call,
    coefficients = coef_table,
    loglik = object$loglik,
    aic = aic_val,
    bic = bic_val,
    nobs = object$nobs,
    npar = object$npar,
    indicators = vapply(object$models, function(m) m$indicator, character(1)),
    convergence = object$convergence,
    sigma = setNames(sigma_vals, sigma_names),
    corr = setNames(object$corr, corr_names)
  )

  class(result) <- "summary.cmp"
  result
}

#' Print method for summary.cmp objects
#'
#' @param x A summary.cmp object
#' @param ... Additional arguments
#' @export
print.summary.cmp <- function(x, ...) {
  cat("Conditional Mixed Process Model

")
  cat("Call:
")
  print(x$call)
  cat("
")

  cat("Coefficients:
")
  print(x$coefficients)
  cat("---
")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

")

  if (length(x$sigma)) {
    cat("Error Standard Deviations:
")
    print(x$sigma)
    cat("
")
  }

  if (length(x$corr)) {
    cat("Error Correlations:
")
    print(x$corr)
    cat("
")
  }

  cat("Log-likelihood:", x$loglik, "
")
  cat("AIC:", x$aic, "  BIC:", x$bic, "
")
  cat("Observations:", x$nobs, "
")

  if (x$convergence != 0) {
    cat("Warning: Optimization did not converge (code:", x$convergence, ")
")
  }

  invisible(x)
}

#' Predict method for CMP objects
#'
#' @param object A cmp object
#' @param newdata Optional new data for predictions
#' @param type Type of prediction ("link", "response", "prob")
#' @param ... Additional arguments
#' @export
predict.cmp <- function(object, newdata = NULL, type = "link", ...) {
  type <- match.arg(type, c("link", "response", "prob"))

  link_vals <- if (is.null(newdata)) {
    object$fitted.values
  } else {
    cmp_predict_newdata(object, newdata)
  }

  if (type == "link") {
    return(link_vals)
  }

  indicators <- vapply(object$models, function(m) m$indicator, character(1))

  if (type == "response") {
    response <- link_vals
    for (i in seq_along(indicators)) {
      if (indicators[i] == "probit") {
        response[, i] <- stats::pnorm(link_vals[, i])
      }
    }
    return(response)
  }

  # type == "prob"
  probs <- vector("list", length(indicators))
  names(probs) <- vapply(object$models, function(m) m$eq_name, character(1))

  for (i in seq_along(indicators)) {
    indicator <- indicators[i]
    if (indicator == "probit") {
      p1 <- stats::pnorm(link_vals[, i])
      probs[[i]] <- cbind(`Pr(0)` = 1 - p1, `Pr(1)` = p1)
    } else if (indicator == "oprobit") {
      cuts <- object$cutpoints[[i]]
      n_cat <- length(cuts) + 1L
      prob_mat <- matrix(NA_real_, nrow = nrow(link_vals), ncol = n_cat)
      colnames(prob_mat) <- levels(object$models[[i]]$y)
      for (j in seq_len(nrow(link_vals))) {
        mu <- link_vals[j, i]
        lower <- c(-Inf, cuts)
        upper <- c(cuts, Inf)
        probs_j <- stats::pnorm(upper - mu) - stats::pnorm(lower - mu)
        probs_j <- pmax(probs_j, 0)
        total <- sum(probs_j)
        if (total > 0) {
          prob_mat[j, ] <- probs_j / total
        } else {
          prob_mat[j, ] <- NA_real_
        }
      }
      probs[[i]] <- prob_mat
    } else {
      probs[[i]] <- matrix(NA_real_, nrow = nrow(link_vals), ncol = 0)
    }
  }

  probs
}

#' Variance-covariance matrix for CMP objects
#'
#' @param object A cmp object
#' @param ... Additional arguments
#' @export
vcov.cmp <- function(object, ...) {
  coef_len <- length(object$coefficients)
  vcov_coef <- object$vcov[seq_len(coef_len), seq_len(coef_len), drop = FALSE]
  rownames(vcov_coef) <- colnames(vcov_coef) <- names(object$coefficients)
  vcov_coef
}

#' Extract coefficients from CMP objects
#'
#' @param object A cmp object
#' @param ... Additional arguments
#' @export
coef.cmp <- function(object, ...) {
  object$coefficients
}

#' Log-likelihood for CMP objects
#'
#' @param object A cmp object
#' @param ... Additional arguments
#' @export
logLik.cmp <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$npar
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Make predictions for new data
#'
#' @param object A cmp object
#' @param newdata New data frame
#' @keywords internal
cmp_predict_newdata <- function(object, newdata) {
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame")
  }

  n_eq <- length(object$models)
  n_obs <- nrow(newdata)
  preds <- matrix(NA_real_, nrow = n_obs, ncol = n_eq)
  colnames(preds) <- vapply(object$models, function(m) m$eq_name, character(1))

  for (i in seq_len(n_eq)) {
    model <- object$models[[i]]
    terms <- model$terms
    mf <- model.frame(terms, data = newdata, na.action = stats::na.pass)
    X_new <- stats::model.matrix(terms, mf, na.action = stats::na.pass)
    template_cols <- colnames(model$X)
    missing_cols <- setdiff(template_cols, colnames(X_new))
    if (length(missing_cols)) {
      add_mat <- matrix(0, nrow = n_obs, ncol = length(missing_cols))
      colnames(add_mat) <- missing_cols
      X_new <- cbind(X_new, add_mat)
    }
    X_new <- X_new[, template_cols, drop = FALSE]
    start <- object$param_info$coef$offsets[i] + 1L
    end <- object$param_info$coef$offsets[i + 1L]
    beta <- object$coefficients[start:end]
    preds[, i] <- as.numeric(X_new %*% beta)
  }

  preds
}
