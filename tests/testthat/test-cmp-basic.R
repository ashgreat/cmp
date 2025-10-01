test_that("cmp function validates inputs correctly", {

  # Test invalid formula
  expect_error(cmp(indicators = "probit", data = mtcars),
               "argument \"formula\" is missing")

  # Test invalid indicators
  expect_error(cmp(formula = vs ~ mpg, data = mtcars, indicators = "invalid"),
               "Unsupported indicator type")

  # Test mismatched formula and indicators
  expect_error(cmp(formula = list(vs ~ mpg, am ~ hp), data = mtcars, indicators = "probit"),
               "Number of formulas \\(")
})

test_that("single equation probit model works", {
  skip_on_cran()

  data(mtcars)
  formula_list <- list(vs ~ mpg + hp)
  indicators <- "probit"

  model <- cmp(formula_list, data = mtcars, indicators = indicators,
              control = list(maxit = 10))

  expect_s3_class(model, "cmp")
  expect_true(length(model$coefficients) > 0)
  expect_true(is.finite(model$loglik))
})

test_that("bivariate probit model works", {
  skip_on_cran()

  data(mtcars)
  formula_list <- list(
    vs ~ mpg + hp,
    am ~ mpg + wt
  )
  indicators <- c("probit", "probit")

  model <- cmp(formula_list, data = mtcars, indicators = indicators,
              control = list(maxit = 10))

  expect_s3_class(model, "cmp")
  expect_equal(length(model$models), 2)
  expect_true(length(model$coefficients) > 0)
  expect_true(is.matrix(model$Sigma))
})

test_that("continuous model works", {
  skip_on_cran()

  data(mtcars)
  formula_list <- list(mpg ~ hp + wt)
  indicators <- "continuous"

  model <- cmp(formula_list, data = mtcars, indicators = indicators,
              control = list(maxit = 10))

  expect_s3_class(model, "cmp")
  expect_equal(length(model$models), 1)
  expect_true(all(is.finite(model$coefficients)))
})

test_that("print and summary methods work", {
  skip_on_cran()

  data(mtcars)
  formula_list <- list(vs ~ mpg)
  indicators <- "probit"

  model <- cmp(formula_list, data = mtcars, indicators = indicators,
              control = list(maxit = 5))

  expect_output(print(model), "Conditional Mixed Process Model")

  summ <- summary(model)
  expect_s3_class(summ, "summary.cmp")
  expect_output(print(summ), "Conditional Mixed Process Model")
})

test_that("predict method works", {
  skip_on_cran()

  data(mtcars)
  formula_list <- list(vs ~ mpg)
  indicators <- "probit"

  model <- cmp(formula_list, data = mtcars, indicators = indicators,
              control = list(maxit = 5))

  pred_link <- predict(model, type = "link")
  pred_response <- predict(model, type = "response")
  pred_prob <- predict(model, type = "prob")

  expect_true(is.matrix(pred_link))
  expect_true(is.matrix(pred_response))
  expect_true(is.list(pred_prob))
  expect_named(pred_prob, vapply(model$models, function(m) m$eq_name, character(1)))
})

test_that("utility helpers behave", {
  sigma <- c(1, 1.5)
  rho <- c(0.3)
  Sigma <- diag(sigma^2)
  Sigma[1, 2] <- Sigma[2, 1] <- rho * sigma[1] * sigma[2]
  chol_info <- safe_chol_logdet(Sigma)
  expect_true(chol_info$ok)
  expect_true(is.finite(chol_info$logdet))

  cuts <- build_cutpoints(c(0.2, log(0.5)))
  expect_length(cuts, 2)
  expect_true(cuts[2] > cuts[1])

  bounds <- resolve_bound(1, data.frame(x = 1:3))
  expect_equal(bounds, rep(1, 3))
})

test_that("model preparation flags missing data", {
  skip_on_cran()

  data(mtcars)
  mtcars_na <- mtcars
  mtcars_na$hp[1] <- NA

  spec <- list(mpg ~ hp)
  indicators <- "continuous"

  prep <- prepare_cmp_models(spec, mtcars_na, indicators)
  expect_false(prep$models[[1]]$mask[1])
})

test_that("validate_cmp_inputs catches dimension mismatch", {
  expect_error(
    validate_cmp_inputs(
      formula = list(vs ~ mpg, am ~ hp),
      indicators = "probit",
      data = mtcars
    ),
    "Number of formulas"
  )
})

test_that("left and right censoring work", {
  skip_on_cran()

  set.seed(123)
  n <- 80
  x <- rnorm(n)
  eps <- rnorm(n)
  y_star <- 1 + 0.5 * x + eps

  data_left <- data.frame(x = x, y = pmax(y_star, 0))
  indicators_left <- list(list(type = "left", lower = 0))
  model_left <- cmp(list(y ~ x), data = data_left, indicators = indicators_left,
                    control = list(maxit = 15))
  expect_s3_class(model_left, "cmp")
  expect_true(model_left$convergence %in% c(0, 1))

  y_star_right <- 0.5 - 0.4 * x + eps
  data_right <- data.frame(x = x, y = pmin(y_star_right, 1))
  indicators_right <- list(list(type = "right", upper = 1))
  model_right <- cmp(list(y ~ x), data = data_right, indicators = indicators_right,
                     control = list(maxit = 15))
  expect_s3_class(model_right, "cmp")
  expect_true(model_right$convergence %in% c(0, 1))
})

test_that("interval outcomes are supported", {
  skip_on_cran()

  set.seed(456)
  n <- 60
  x <- rnorm(n)
  eps <- rnorm(n)
  z <- 0.8 * x + eps
  lower <- z - 0.3
  upper <- z + 0.3
  exact_idx <- sample.int(n, size = 15)
  lower[exact_idx] <- upper[exact_idx] <- z[exact_idx]

  data_int <- data.frame(x = x, lower = lower, upper = upper)
  indicators <- list("interval")
  model <- cmp(list(cbind(lower, upper) ~ x), data = data_int, indicators = indicators,
               control = list(maxit = 20))

  expect_s3_class(model, "cmp")
  expect_true(is.finite(model$loglik))
})

test_that("ordered probit pathway runs", {
  skip_on_cran()

  set.seed(789)
  n <- 100
  x <- rnorm(n)
  latent <- 0.5 * x + rnorm(n)
  y <- cut(latent, breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("low", "mid", "high"), ordered_result = TRUE)

  data_op <- data.frame(y = y, x = x)
  indicators <- list("oprobit")
  model <- cmp(list(y ~ x), data = data_op, indicators = indicators,
               control = list(maxit = 20))

  expect_s3_class(model, "cmp")
  expect_true(length(model$cutpoints[[1]]) == 2)

  prob_pred <- predict(model, type = "prob")
  expect_true(is.list(prob_pred))
  expect_true(all(abs(rowSums(prob_pred[[1]]) - 1) < 1e-6, na.rm = TRUE))
})
