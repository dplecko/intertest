# small helpers -----------------------------------------------------------

make_data <- function(n = 200, seed = 1) {

  set.seed(seed)
  Z <- rnorm(n)
  X <- rbinom(n, 1, plogis(Z))
  W <- X + Z + rnorm(n)
  Y <- X * W + rnorm(n)
  data.frame(Z = Z, X = X, W = W, Y = Y)
}

make_binary_data <- function(n = 200, seed = 1) {

  set.seed(seed)
  Z <- rnorm(n)
  X <- rbinom(n, 1, plogis(Z))
  W <- X + Z + rnorm(n)
  Y <- rbinom(n, 1, plogis(X + W))
  data.frame(Z = Z, X = X, W = W, Y = Y)
}

# 1. return value structure -----------------------------------------------

test_that("inter_test() returns an intertest object with correct structure", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)

  expect_s3_class(fit, "intertest")
  expect_named(fit, c("results", "pso", "call", "graph", "n"))
  expect_true(is.data.frame(fit$results))
  expect_named(fit$results, c("measure", "value", "sd", "scale", "p_value"))
  expect_equal(fit$n, 200L)
  expect_null(fit$pso)
})

test_that("results contain all expected measures", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)

  expect_setequal(
    fit$results$measure,
    c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE",
      "tv", "ctfde", "ctfie", "ctfse", "ett")
  )
})

test_that("p_values are in [0, 1]", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)

  expect_true(all(fit$results$p_value >= 0 & fit$results$p_value <= 1,
                  na.rm = TRUE))
})

test_that("graph slot records the correct variable names", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)

  expect_equal(fit$graph$X, "X")
  expect_equal(fit$graph$Z, "Z")
  expect_equal(fit$graph$W, "W")
  expect_equal(fit$graph$Y, "Y")
})

# 2. S3 methods -----------------------------------------------------------

test_that("print.intertest() produces output", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)
  expect_output(print(fit), "TE x SE")
  expect_output(print(fit), "Difference scale")
})

test_that("summary.intertest() shows both ia and xspec sections", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)
  expect_output(summary(fit), "Interaction effects")
  expect_output(summary(fit), "Counterfactual effects")
  expect_output(summary(fit), "TV")
  expect_output(summary(fit), "Ctf-DE")
})

test_that("autoplot.intertest() returns a ggplot object", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y", nthread = 1)
  p <- ggplot2::autoplot(fit)

  expect_s3_class(p, "gg")
  expect_equal(length(p$layers), 3L)  # geom_col, geom_errorbar, geom_hline
})

# 3. multiple scales ------------------------------------------------------

test_that("two scales produce 2x the rows, one section per scale", {

  fit <- inter_test(make_binary_data(), "X", "Z", "W", "Y",
                    scale = c("difference", "log-risk"), nthread = 1)

  expect_equal(nrow(fit$results), 20L)  # 10 measures x 2 scales
  expect_setequal(unique(fit$results$scale), c("difference", "log-risk"))
})

test_that("autoplot facets when multiple scales requested", {

  fit <- inter_test(make_binary_data(), "X", "Z", "W", "Y",
                    scale = c("difference", "log-risk"), nthread = 1)
  p <- ggplot2::autoplot(fit)

  expect_s3_class(p$facet, "FacetWrap")
})

# 4. edge cases: graph structure ------------------------------------------

test_that("Z empty (no confounders)", {

  dat <- make_data()
  fit <- inter_test(dat, "X", W = "W", Y = "Y", nthread = 1)

  expect_s3_class(fit, "intertest")
  expect_length(fit$graph$Z, 0L)
})

test_that("W empty (no mediators)", {

  dat <- make_data()
  fit <- inter_test(dat, "X", Z = "Z", Y = "Y", nthread = 1)

  expect_s3_class(fit, "intertest")
  expect_length(fit$graph$W, 0L)
})

test_that("Z and W both empty", {

  dat <- make_data()
  fit <- inter_test(dat, "X", Y = "Y", nthread = 1)

  expect_s3_class(fit, "intertest")
})

# 5. preproc --------------------------------------------------------------

test_that("character X is auto-encoded with a message", {

  dat <- make_data()
  dat$X <- ifelse(dat$X == 1, "treated", "control")

  expect_message(
    fit <- inter_test(dat, "X", "Z", "W", "Y", nthread = 1),
    "Auto-encoded"
  )
  expect_s3_class(fit, "intertest")
})

test_that("factor Z with 3 levels gets dummy-expanded", {

  dat <- make_data()
  dat$Z <- cut(dat$Z, breaks = 3, labels = c("low", "mid", "high"))

  expect_message(
    fit <- inter_test(dat, "X", "Z", "W", "Y", nthread = 1),
    "Auto-encoded"
  )
  # dummy expansion: original Z replaced by Z1, Z2
  expect_false("Z" %in% fit$graph$Z)
  expect_true(any(grepl("^Z", fit$graph$Z)))
})

test_that("preproc = FALSE skips encoding on numeric data", {

  dat <- make_data()
  fit <- inter_test(dat, "X", "Z", "W", "Y", preproc = FALSE, nthread = 1)

  expect_s3_class(fit, "intertest")
})

# 6. save_pso -------------------------------------------------------------

test_that("save_pso = TRUE stores pseudo-outcomes keyed by scale", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y",
                    save_pso = TRUE, nthread = 1)

  expect_false(is.null(fit$pso))
  expect_true("difference" %in% names(fit$pso))
})

# 7. input guards ---------------------------------------------------------

test_that("log-odds requires Y in [0, 1]", {

  expect_error(
    inter_test(make_data(), "X", "Z", "W", "Y", scale = "log-odds", nthread = 1),
    regexp = "log-odds"
  )
})

test_that("log-risk requires Y >= 0", {

  dat <- make_data(); dat$Y <- dat$Y - 20
  expect_error(
    inter_test(dat, "X", "Z", "W", "Y", scale = "log-risk", nthread = 1),
    regexp = "log-risk"
  )
})

# 8. one_step_debias() directly -------------------------------------------

test_that("one_step_debias() returns a data frame with correct columns", {

  res <- one_step_debias(make_data(), "X", "Z", "W", "Y",
                         scale = "difference", nthread = 1)

  expect_true(is.data.frame(res))
  expect_named(res, c("measure", "value", "sd", "scale"))
  expect_equal(nrow(res), 10L)  # 5 ia + 5 xspec
})

test_that("one_step_debias() log-risk runs on binary Y", {

  res <- one_step_debias(make_binary_data(), "X", "Z", "W", "Y",
                         scale = "log-risk", nthread = 1)

  expect_true(is.data.frame(res))
  expect_true(all(res$scale == "log-risk"))
})

# 9. statistical tests (slow, skip on CRAN) -------------------------------

test_that("null SCM (class C) yields mostly non-significant IA measures", {

  skip_on_cran()
  skip_on_ci()

  set.seed(42)
  scm <- gen_from_scm("C", n = 2000)
  fit <- inter_test(scm$data, scm$mapping$X, scm$mapping$Z,
                    scm$mapping$W, scm$mapping$Y, nthread = 1)

  ia <- fit$results[fit$results$measure %in%
    c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE"), ]

  # under the null, expect < 40% rejection at alpha = 0.05
  expect_lt(mean(ia$p_value < 0.05, na.rm = TRUE), 0.4)
})

test_that("strong X*W interaction SCM yields significant 2nd-order IA measures", {

  skip_on_cran()
  skip_on_ci()

  # custom SCM with large X*W interaction in Y — designed for reliable detection
  set.seed(42)
  n  <- 1000
  Z  <- rnorm(n)
  X  <- rbinom(n, 1, plogis(Z))
  W  <- 2 * X + Z + rnorm(n, 0, 0.5)
  Y  <- 3 * X * W + rnorm(n)
  dat <- data.frame(Z = Z, X = X, W = W, Y = Y)

  fit <- inter_test(dat, "X", "Z", "W", "Y", nthread = 1)

  ia_2nd <- fit$results[fit$results$measure %in%
    c("TE x SE", "DE x IE", "DE x SE"), ]

  # all three 2nd-order measures with X*W signal should be highly significant
  expect_true(all(ia_2nd$p_value < 0.01))
})

# 10. log-odds scale ------------------------------------------------------

test_that("log-odds scale runs on binary Y and returns logit-scale estimates", {

  fit <- inter_test(make_binary_data(), "X", "Z", "W", "Y",
                    scale = "log-odds", nthread = 1)

  expect_s3_class(fit, "intertest")
  expect_true(all(fit$results$scale == "log-odds"))
  expect_true(all(is.finite(fit$results$value)))
})

# 11. structural edge cases -----------------------------------------------

test_that("ctfie is zero when W is empty (no mediator path)", {

  dat <- make_data()
  fit <- inter_test(dat, "X", Z = "Z", Y = "Y", nthread = 1)
  ctfie <- fit$results[fit$results$measure == "ctfie", "value"]

  expect_equal(ctfie, 0, tolerance = 1e-10)
})

test_that("ctfse is zero when Z is empty (no spurious path)", {

  dat <- make_data()
  fit <- inter_test(dat, "X", W = "W", Y = "Y", nthread = 1)
  ctfse <- fit$results[fit$results$measure == "ctfse", "value"]

  expect_equal(ctfse, 0, tolerance = 1e-10)
})

# 12. overlap message -----------------------------------------------------

test_that("overlap message triggered when extreme propensity weights present", {

  # mean shift of 3: Z|X=0 ~ N(0,1), Z|X=1 ~ N(3,1) — near-complete separation
  set.seed(1)
  n <- 500
  X <- rbinom(n, 1, 0.5)
  Z <- X * 3 + rnorm(n)
  W <- X + Z + rnorm(n)
  Y <- X + W + rnorm(n)
  dat <- data.frame(Z = Z, X = X, W = W, Y = Y)

  expect_message(
    inter_test(dat, "X", "Z", "W", "Y", eps_trim = 0.1, nthread = 1),
    "extreme propensity"
  )
})

# 13. nested_mean = "wregr" -----------------------------------------------

test_that("nested_mean = 'wregr' runs on difference scale with finite results", {

  fit <- inter_test(make_data(), "X", "Z", "W", "Y",
                    nested_mean = "wregr", nthread = 1)

  expect_s3_class(fit, "intertest")
  expect_true(all(is.finite(fit$results$value)))
})

test_that("nested_mean = 'wregr' errors on log-risk scale", {

  expect_error(
    inter_test(make_binary_data(), "X", "Z", "W", "Y",
               scale = "log-risk", nested_mean = "wregr", nthread = 1),
    "Weighted regression not available"
  )
})
