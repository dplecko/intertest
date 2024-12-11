
preproc_data <- function(data, X, Z, W, Y) {

  SFM <- list(X = X, Z = Z, W = W, Y = Y)
  for (cmp in c("X", "Z", "W", "Y")) {

    for (var in SFM[[cmp]]) {

      if (is.numeric(data[[var]]) || is.integer(data[[var]])) next
      if (is.logical(data[[var]])) {

        data[[var]] <- as.integer(data[[var]])
        next
      }

      if (is.character(data[[var]])) data[[var]] <- as.factor(data[[var]])

      # only factors at this point
      if (length(levels(data[[var]])) == 2) {

        data[[var]] <- as.integer(data[[var]] == levels(data[[var]])[1])
      } else {

        enc_mat <- model.matrix(~ ., data = data.frame(var = data[[var]]))[, -1]
        colnames(enc_mat) <- paste0(var, seq_len(ncol(enc_mat)))

        # remove the factor column
        SFM[[cmp]] <- setdiff(SFM[[cmp]], var)
        data[[var]] <- NULL

        # add the encoded data
        SFM[[cmp]] <- c(SFM[[cmp]], colnames(enc_mat))
        data <- cbind(data, enc_mat)
      }
    }
  }

  list(data = data, sfm = SFM)
}

#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train
cv_xgb <- function(df, y, weights = NULL, ...) {

  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y, weight = weights)

  binary <- all(y %in% c(0, 1))
  if (binary) {

    params <- list(objective = "binary:logistic", eval_metric = "logloss")
  } else {

    params <- list(objective = "reg:squarederror", eval_metric = "rmse")
  }

  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    prediction = TRUE,
    verbose = FALSE, ...
  )

  xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = cv$best_iteration,
    verbose = FALSE, ...
  )
  attr(xgb, "binary") <- binary

  xgb
}

pred_xgb <- function(xgb, df_test, intervention = NULL, X = "X") {

  if (!is.null(intervention)) {

    df_test[[X]] <- intervention
  }

  predict(xgb, as.matrix(df_test))
}

one_step_debias <- function(data, X, Z, W, Y, nested_mean = c("refit", "wregr"),
                            log_risk = FALSE, eps_trim = 0, ...) {

  n <- nrow(data)
  nested_mean <- match.arg(nested_mean, c("refit", "wregr"))

  # split into K folds
  K <- 10
  # folds <- sample(x = rep(1:K, each = n / K))
  folds <- sample(x = rep(1:K, each = ceiling(n / K)))[seq_len(n)]

  y <- data[[Y]]
  x <- data[[X]]
  pso <- list(list(list(list(), list()), list(list(), list())),
              list(list(list(), list()), list(list(), list())))

  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {

    pso[[xz+1]][[xw+1]][[xy+1]] <- rep(NA, n)
  }

  px_z <- px_zw <- rep(NA, n)

  if (length(Z) == 0 & length(W) == 0) {

    # get all potential outomes with no cross-fitting
    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {

      if (isFALSE(log_risk)) {

        pso[[xz+1]][[xw+1]][[xy+1]] <- (x == xy) / mean(x == xy) * y
      } else {

        pso[[xz+1]][[xw+1]][[xy+1]] <-
          log(mean(y[x==xy])) +
          (x == xy) / mean(x == xy) * (y / mean(y[x==xy]) - 1)
      }
    }

    # get the propensity weights
    px_z <- px_zw <- rep(mean(x), n)
  } else {

    # get the simplest pseudo-outcomes (no cross-fitting needed)
    for (xz in c(0, 1)) if (isFALSE(log_risk))
      pso[[xz+1]][[xz+1]][[xz+1]] <- (x == xz) / mean(x == xz) * y

    # cross-fitting procedure
    for (i in seq_len(K)) {

      # split into dev, val, tst
      tst <- folds == i
      dev <- folds %in% setdiff(seq_len(K), i)[1:6]
      val <- folds %in% setdiff(seq_len(K), i)[7:9]

      # develop models on dev
      if (length(Z) > 0) {

        mod_x_z <- cv_xgb(data[dev, Z], data[dev, X], ...)
        mod_y_xz <- cv_xgb(data[dev, c(X, Z)], data[dev, Y], ...)
      }

      if (length(W) > 0) {

        mod_x_zw <- cv_xgb(data[dev, c(Z, W)], data[dev, X], ...)
        mod_y_xzw <- cv_xgb(data[dev, c(X, Z, W)], data[dev, Y], ...)
      } else {

        # inherit from Z if W empty
        mod_x_zw <- mod_x_z
        mod_y_xzw <- mod_y_xz
      }

      # get the val set predictions
      px_zw_val <- pred_xgb(mod_x_zw, data[val, c(Z, W)])
      px_zw_val <- list(1 - px_zw_val, px_zw_val)

      if (length(Z) > 0) {

        px_z_val <- pred_xgb(mod_x_z, data[val, Z])
        px_z_val <- list(1 - px_z_val, px_z_val)
      }

      y_xzw_val <- list(
        pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], intervention = 0, X = X),
        pred_xgb(mod_y_xzw, data[val, c(X, Z, W)], intervention = 1, X = X)
      )

      # get the test set values
      px_zw_tst <- pred_xgb(mod_x_zw, data[tst, c(Z, W)])
      px_zw[tst] <- px_zw_tst
      px_zw_tst <- list(1 - px_zw_tst, px_zw_tst)

      if (length(Z) > 0) {

        px_z_tst <- pred_xgb(mod_x_z, data[tst, Z])
        px_z[tst] <- px_z_tst
        px_z_tst <- list(1 - px_z_tst, px_z_tst)
      }

      y_xzw_tst <- list(
        pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)], intervention = 0, X = X),
        pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)], intervention = 1, X = X)
      )

      if (length(Z) > 0) {

        y_xz_tst <- list(
          pred_xgb(mod_y_xz, data[tst, c(X, Z)], intervention = 0, X = X),
          pred_xgb(mod_y_xz, data[tst, c(X, Z)], intervention = 1, X = X)
        )
      }

      # get the total-effect pseudo-outcomes (ETT-like)
      if (length(Z) == 0) {

        for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {

          if (xz == xw & xz == xy) next

          if (isFALSE(log_risk)) {

            pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
              (x[tst] == xy) / mean(x == xw) * px_zw_tst[[xw+1]] / px_zw_tst[[xw+1]] *
              (y[tst] - y_xzw_tst[[xy+1]]) +
              (x[tst] == xw) / mean(x == xw) * y_xzw_tst[[xy+1]]
          } else {

            pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
              (x[tst] == xy) / mean(x == xw) * px_zw_tst[[xw+1]] / px_zw_tst[[xw+1]] *
              (y[tst] - y_xzw_tst[[xy+1]]) / y_xzw_tst[[xy+1]] +
              (x[tst] == xw) / mean(x == xw) * log(y_xzw_tst[[xy+1]])
          }
        }
      } else if (length(W) == 0) {

        for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {

          if (xz == xw & xz == xy) next

          if (isFALSE(log_risk)) {

            pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
              (x[tst] == xy) / mean(x == xz) * px_zw_tst[[xz+1]] / px_zw_tst[[xz+1]] *
              (y[tst] - y_xzw_tst[[xy+1]]) +
              (x[tst] == xz) / mean(x == xz) * y_xzw_tst[[xy+1]]
          } else {

            pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
              (x[tst] == xy) / mean(x == xz) * px_zw_tst[[xz+1]] / px_zw_tst[[xz+1]] *
              (y[tst] - y_xzw_tst[[xy+1]]) / y_xzw_tst[[xy+1]] +
              (x[tst] == xz) / mean(x == xz) * log(y_xzw_tst[[xy+1]])
          }
        }
      } else {
        for (xz in c(0, 1)) {

          xw <- xy <- 1 - xz
          if (isFALSE(log_risk)) {

            pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
              # Term I
              (x[tst] == xw) / mean(x == xz) * px_z_tst[[xz+1]] / px_z_tst[[xw+1]] *
              (y[tst] - y_xz_tst[[xw+1]]) +
              # Term II
              (x[tst] == xz) / mean(x == xz) * y_xz_tst[[xw+1]]
          }
        }
      }

      # nested means are not needed if either Z or W are empty
      if (length(Z) == 0 || length(W) == 0) next

      # get the nested mean
      ey_nest_tst <- list(NULL, NULL)
      for (xw in c(0, 1)) {

        xy <- 1 - xw
        if (nested_mean == "wregr") {

          assertthat::assert_that(
            isFALSE(log_risk),
            msg = "Weighted regression not available for log-risk scale."
          )

          weights <- ifelse(
            x[val] == xy,
            px_z_val[[xy+1]] / px_z_val[[xw+1]] *
              px_zw_val[[xw+1]] / px_zw_val[[xy+1]],
            px_z_val[[xy+1]] / px_z_val[[xw+1]]
          )

          mod_nested <- cv_xgb(data[val, Z], data[val, Y], weights = weights,
                               ...)
          ey_nest_tst[[xy+1]] <- pred_xgb(mod_nested, data[val, Z])
        } else if (nested_mean == "refit") {

          y_tilde <- pred_xgb(mod_y_xzw, data[val, c(X, Z, W)],
                              intervention = xy, X = X)
          if (log_risk) y_tilde <- log(y_tilde)
          mod_nested <- cv_xgb(data[val, c(X, Z)], y_tilde, ...)
          ey_nest_tst[[xy+1]] <- pred_xgb(mod_nested, data[tst, c(X, Z)],
                                          intervention = xw, X = X)
        }
      }

      for (xz in c(0, 1)) {

        for (xw in c(0, 1)) {

          for (xy in c(0, 1)) {

            if (xy == xw & isFALSE(log_risk)) next

            if (isFALSE(log_risk)) {

              pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
                # Term T1
                (x[tst] == xy) * (y[tst] - y_xzw_tst[[xy+1]]) *
                px_zw_tst[[xw+1]] / px_zw_tst[[xy+1]] *
                px_z_tst[[xz+1]] / px_z_tst[[xw+1]] * 1 / mean(x == xz) +
                # Term T2
                (x[tst] == xw) / mean(x == xz) * px_z_tst[[xz+1]] / px_z_tst[[xw+1]] *
                (y_xzw_tst[[xy+1]] - ey_nest_tst[[xy+1]]) +
                # Term T3
                (x[tst] == xz) / mean(x == xz) * ey_nest_tst[[xy+1]]
            } else {

              pso[[xz+1]][[xw+1]][[xy+1]][tst] <-
                # Term T1
                (x[tst] == xy) / mean(x == xz) *
                (y[tst] - y_xzw_tst[[xy+1]]) / y_xzw_tst[[xy+1]] *
                px_zw_tst[[xw+1]] / px_zw_tst[[xy+1]] *
                px_z_tst[[xz+1]] / px_z_tst[[xw+1]] +
                # Term T2
                (x[tst] == xw) / mean(x == xz) * px_z_tst[[xz+1]] / px_z_tst[[xw+1]] *
                (log(y_xzw_tst[[xy+1]]) - ey_nest_tst[[xy+1]]) +
                # Term T3
                (x[tst] == xz) / mean(x == xz) * ey_nest_tst[[xy+1]]
            }
          }
        }
      }
    }
  }

  if (all(is.na(px_z))) px_z <- rep(mean(x), n)
  extrm_pxz <- (px_z < eps_trim) | (1 - px_z < eps_trim)
  extrm_pxzw <-  (px_zw < eps_trim) | (1 - px_zw < eps_trim)
  extrm_idx <- extrm_pxz | extrm_pxzw

  if (mean(extrm_idx) > 0.02) {
    message(round(100 * mean(extrm_idx), 2),
            "% of extreme P(x | z) or P(x | z, w) probabilities at threshold",
            " = ", eps_trim, ".\n",
            "Reported results are for the overlap population. ",
            "Consider investigating overlap issues.")
  }

  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1))
    pso[[xz+1]][[xw+1]][[xy+1]][extrm_idx] <- NA

  ias <- measure_spec()

  res <- c()
  scale <- if (isFALSE(log_risk)) "difference" else "log-risk"
  for (i in seq_along(ias)) {

    pseudo_out <- 0
    for (j in seq_along(ias[[i]]$sgn)) {

      xz <- ias[[i]]$spc[[j]][1]
      xw <- ias[[i]]$spc[[j]][2]
      xy <- ias[[i]]$spc[[j]][3]
      pseudo_out <- pseudo_out + ias[[i]]$sgn[j] * pso[[xz+1]][[xw+1]][[xy+1]]
    }
    psi_osd <- mean(pseudo_out, na.rm = TRUE)
    dev <- sqrt(var(pseudo_out, na.rm = TRUE) / sum(!is.na(pseudo_out)))

    res <- rbind(
      res,
      data.frame(measure = ias[[i]]$ia, value = psi_osd, sd = dev, scale = scale)
    )
  }

  pw <- list(px_z = px_z, px_zw = px_zw)
  attr(res, "pw") <- pw
  res
}

measure_spec <- function(spec = c("xspec", "ia")) {

  meas <- list()
  if (is.element("ia", spec)) {

    ia <- list(
      tese = list(
        sgn = c(1, -1, -1, 1),
        spc = list(
          c(0, 1, 1), c(0, 0, 0), c(1, 1, 1), c(1, 0, 0)
        ),
        ia = "TE x SE"
      ),
      deie0 = list(
        sgn = c(1, -1, -1, 1),
        spc = list(
          c(0, 0, 1), c(0, 0, 0), c(0, 1, 1), c(0, 1, 0)
        ),
        ia = "DE x IE"
      ),
      # deie1 = list(
      #   sgn = c(1, -1, -1, 1),
      #   spc = list(
      #     c(1, 0, 1), c(1, 0, 0), c(1, 1, 1), c(1, 1, 0)
      #   ),
      #   ia = "DE x IE"
      # ),
      dese = list(
        sgn = c(1, -1, -1, 1),
        spc = list(
          c(0, 0, 1), c(0, 0, 0), c(1, 0, 1), c(1, 0, 0)
        ),
        ia = "DE x SE"
      ),
      iese = list(
        sgn = c(1, -1, -1, 1),
        spc = list(
          c(0, 1, 0), c(0, 0, 0), c(1, 1, 0), c(1, 0, 0)
        ),
        ia = "IE x SE"
      ),
      deiese = list(
        sgn = c(c(1, -1, -1, 1), -c(1, -1, -1, 1)),
        spc = list(
          c(0, 0, 1), c(0, 0, 0), c(0, 1, 1), c(0, 1, 0),
          c(1, 0, 1), c(1, 0, 0), c(1, 1, 1), c(1, 1, 0)
        ),
        ia = "DE x IE x SE"
      )
    )
    meas <- c(meas, ia)
  }

  if (is.element("xspec", spec)) {

    xspec <- list(
      tv = list(
        sgn = c(1, -1),
        spc = list(c(1, 1, 1), c(0, 0, 0)),
        ia = "tv"
      ),
      ctfde = list(
        sgn = c(1, -1),
        spc = list(c(0, 0, 1), c(0, 0, 0)),
        ia = "ctfde"
      ),
      ctfie = list(
        sgn = c(1, -1),
        spc = list(c(0, 0, 1), c(0, 1, 1)),
        ia = "ctfie"
      ),
      ctfse = list(
        sgn = c(1, -1),
        spc = list(c(0, 1, 1), c(1, 1, 1)),
        ia = "ctfse"
      ),
      ett = list(
        sgn = c(1, -1),
        spc = list(c(0, 1, 1), c(0, 0, 0)),
        ia = "ett"
      )
    )
    meas <- c(meas, xspec)
  }

  meas
}
