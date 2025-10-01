#' Utility functions for CMP package
#'
#' @name cmp-utils
#' @keywords internal
NULL

#' Validate CMP model specification
#'
#' @param formula List of formulas or single formula
#' @param indicators Vector or list of indicators
#' @param data Data frame
#' @keywords internal
validate_cmp_inputs <- function(formula, indicators, data) {
  if (!is.list(formula) && !inherits(formula, "formula")) {
    stop("formula must be a formula or list of formulas")
  }

  if (inherits(formula, "formula")) {
    formula <- list(formula)
  }

  if (is.null(indicators)) {
    stop("indicators argument is required")
  }

  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }

  if (length(formula) != length(indicators)) {
    stop(
      "Number of formulas (", length(formula),
      ") must match number of indicators (", length(indicators), ")"
    )
  }

  list(formula = formula, indicators = indicators, data = data)
}

#' Parse indicator definition
#'
#' @param indicator Indicator definition
#' @param eq_index Equation index
#' @param data Data frame
#' @keywords internal
parse_indicator_definition <- function(indicator, eq_index, data) {
  make_name <- function(idx) paste0("eq", idx)

  if (is.list(indicator)) {
    if (is.null(indicator$type)) {
      stop("Indicator definition for equation ", eq_index, " must include a 'type' entry")
    }
    type <- tolower(indicator$type)
    name <- indicator$name %||% make_name(eq_index)
    lower <- indicator$lower %||% indicator$limit %||% indicator$bounds
    upper <- indicator$upper %||% indicator$limit
  } else {
    type <- tolower(as.character(indicator))
    name <- make_name(eq_index)
    lower <- NULL
    upper <- NULL
  }

  supported <- c("continuous", "probit", "oprobit", "left", "right", "interval")
  if (!type %in% supported) {
    stop(
      "Unsupported indicator type '", type, "' for equation ", eq_index,
      ". Supported types: ", paste(supported, collapse = ", ")
    )
  }

  has_sigma <- type %in% c("continuous", "left", "right", "interval")
  n_cutpoints <- if (type == "oprobit") NA_integer_ else 0L

  lower_vec <- upper_vec <- NULL
  if (type == "left") {
    lower_vec <- resolve_bound(lower, data, default = 0)
  } else if (type == "right") {
    upper_vec <- resolve_bound(upper, data, default = 0)
  }

  list(
    type = type,
    name = name,
    has_sigma = has_sigma,
    lower = lower_vec,
    upper = upper_vec,
    n_cutpoints = n_cutpoints
  )
}

#' Resolve constant/vector bound definition
#'
#' @param spec Specification (NULL, numeric, or character)
#' @param data Data frame
#' @param default Default value if spec is NULL
#' @keywords internal
resolve_bound <- function(spec, data, default = 0) {
  if (is.null(spec)) {
    return(rep(default, nrow(data)))
  }
  if (length(spec) == 1L && is.numeric(spec)) {
    return(rep(as.numeric(spec), nrow(data)))
  }
  if (is.character(spec) && length(spec) == 1L) {
    if (!spec %in% names(data)) {
      stop("Bound variable '", spec, "' not found in data")
    }
    return(as.numeric(data[[spec]]))
  }
  if (length(spec) == nrow(data)) {
    return(as.numeric(spec))
  }
  stop("Unable to resolve bound specification; provide scalar, vector, or column name")
}

#' Prepare per-equation model objects
#'
#' @param formula_list List of formulas
#' @param data Data frame
#' @param indicators Indicator definitions
#' @param na.action NA handling function
#' @keywords internal
prepare_cmp_models <- function(formula_list, data, indicators, na.action = na.omit) {
  n_eq <- length(formula_list)
  models <- vector("list", n_eq)
  n_obs <- nrow(data)

  for (i in seq_len(n_eq)) {
    info <- parse_indicator_definition(indicators[[i]], i, data)
    mf <- model.frame(formula_list[[i]], data = data, na.action = stats::na.pass)
    trms <- stats::terms(mf)
    X <- stats::model.matrix(trms, mf, na.action = stats::na.pass)

    response <- model.response(mf)
    mask <- rep(TRUE, n_obs)

    if (info$type == "continuous" || info$type %in% c("left", "right")) {
      response <- as.numeric(response)
    }

    if (info$type == "probit") {
      if (is.factor(response)) {
        if (nlevels(response) != 2L) {
          stop("Binary probit outcomes must have exactly two levels")
        }
        response <- as.numeric(response == levels(response)[2L])
      } else if (is.logical(response)) {
        response <- as.numeric(response)
      } else {
        response <- as.numeric(response)
      }
      uniq <- unique(stats::na.omit(response))
      if (!all(uniq %in% c(0, 1))) {
        stop("Binary probit outcomes must take values 0/1")
      }
    }

    if (info$type == "oprobit") {
      if (!is.factor(response)) {
        response <- stats::as.ordered(response)
      } else if (!is.ordered(response)) {
        response <- stats::as.ordered(response)
      }
      if (nlevels(response) < 3L) {
        stop("Ordered probit equations require at least three outcome categories")
      }
      info$n_cutpoints <- nlevels(response) - 1L
    }

    if (!is.null(response)) {
      mf[[1L]] <- response
    }

    if (info$type == "interval") {
      if (!(is.matrix(response) || is.data.frame(response)) || ncol(response) != 2L) {
        stop("Interval outcomes must be supplied as cbind(lower, upper)")
      }
      lower <- as.numeric(response[, 1L])
      upper <- as.numeric(response[, 2L])
      mask <- !(is.na(lower) & is.na(upper))
    } else {
      lower <- upper <- NULL
      mask <- !is.na(response)
    }

    if (anyNA(X)) {
      mask <- mask & stats::complete.cases(X)
    }

    if (info$type == "left") {
      lower <- info$lower
      if (length(lower) != n_obs) {
        stop("Left-censoring bound must have length equal to number of observations")
      }
    }

    if (info$type == "right") {
      upper <- info$upper
      if (length(upper) != n_obs) {
        stop("Right-censoring bound must have length equal to number of observations")
      }
    }

    models[[i]] <- list(
      eq_index = i,
      eq_name = info$name,
      indicator = info$type,
      X = X,
      terms = trms,
      mf = mf,
      y = response,
      lower = if (!is.null(lower)) as.numeric(lower) else rep(-Inf, n_obs),
      upper = if (!is.null(upper)) as.numeric(upper) else rep(Inf, n_obs),
      mask = mask,
      has_sigma = info$has_sigma,
      n_cutpoints = info$n_cutpoints %||% 0L,
      n_coef = ncol(X)
    )
  }

  list(models = models, n_obs = n_obs)
}

#' Build parameter bookkeeping information
#'
#' @param models List of model objects
#' @keywords internal
build_param_info <- function(models) {
  n_eq <- length(models)
  coef_lengths <- vapply(models, function(m) m$n_coef, integer(1))
  coef_offsets <- c(0L, cumsum(coef_lengths))
  coef_names <- unlist(
    Map(
      function(m, idx) paste0(m$eq_name, ":", colnames(m$X)),
      models,
      seq_len(n_eq)
    ),
    use.names = FALSE
  )

  sigma_eq <- which(vapply(models, function(m) isTRUE(m$has_sigma), logical(1)))
  sigma_count <- length(sigma_eq)
  sigma_names <- if (sigma_count) {
    paste0("lnsigma_", vapply(models[sigma_eq], function(m) m$eq_name, character(1)))
  } else character()

  cut_counts <- vapply(models, function(m) m$n_cutpoints, integer(1))
  cut_offsets <- c(0L, cumsum(cut_counts))
  cut_names <- unlist(
    Map(
      function(m) if (m$n_cutpoints) paste0("cut_", m$eq_name, "_", seq_len(m$n_cutpoints)) else character(),
      models
    ),
    use.names = FALSE
  )

  if (n_eq > 1L) {
    pair_mat <- utils::combn(n_eq, 2L)
    corr_pairs <- t(pair_mat)
    corr_names <- apply(corr_pairs, 1L, function(idx) {
      paste0("atanhrho_", models[[idx[1]]]$eq_name, "_", models[[idx[2]]]$eq_name)
    })
  } else {
    corr_pairs <- matrix(integer(), ncol = 2L)
    corr_names <- character()
  }

  total_coef <- sum(coef_lengths)
  total_cut <- sum(cut_counts)
  total_sigma <- sigma_count
  total_corr <- length(corr_names)
  total_par <- total_coef + total_sigma + total_cut + total_corr

  list(
    n_eq = n_eq,
    coef = list(lengths = coef_lengths, offsets = coef_offsets),
    sigma = list(eq = sigma_eq, count = sigma_count),
    cut = list(counts = cut_counts, offsets = cut_offsets),
    corr = list(pairs = corr_pairs, count = total_corr),
    total_coef = total_coef,
    total_sigma = total_sigma,
    total_cut = total_cut,
    total_par = total_par,
    names = c(coef_names, sigma_names, cut_names, corr_names)
  )
}

#' Construct ordered probit cut points from unconstrained parameters
#'
#' @param par Vector of unconstrained cut parameters
#' @keywords internal
build_cutpoints <- function(par) {
  if (!length(par)) {
    return(numeric())
  }
  cuts <- numeric(length(par))
  cuts[1L] <- par[1L]
  if (length(par) > 1L) {
    diffs <- exp(par[-1L])
    cuts[-1L] <- cuts[-length(cuts)] + diffs
  }
  cuts
}

#' Safe log of determinant via Cholesky
#'
#' @param mat Symmetric matrix
#' @keywords internal
safe_chol_logdet <- function(mat) {
  R <- tryCatch(chol(mat), error = function(e) NULL)
  if (is.null(R)) {
    return(list(ok = FALSE, logdet = NA_real_))
  }
  logdet <- 2 * sum(log(diag(R)))
  list(ok = TRUE, logdet = logdet)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
