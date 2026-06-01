#' @importFrom stats predict
#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train
cv_xgb <- function(df, y, weights = NULL, ...) {

  dtrain <- xgb.DMatrix(data = as.matrix(df), label = y, weight = weights)

  binary <- all(y %in% c(0, 1))

  # merge ... into params so xgboost-level args (e.g. nthread) are passed correctly
  base_params <- if (binary) {
    list(objective = "binary:logistic", eval_metric = "logloss")
  } else {
    list(objective = "reg:squarederror", eval_metric = "rmse")
  }
  params <- c(base_params, list(...))

  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 10,
    prediction = TRUE,
    verbose = FALSE
  )

  xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = cv$early_stop$best_iteration,
    verbose = FALSE
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

cross_fit <- function(data, X, Z, W, Y, nested_mean, scale = "diff", ...) {

  if (length(Z) == 0 & length(W) == 0) {

    px <- rep(mean(data[[X]]), nrow(data))
    return(list(px_z = list(1-px, px), px_zw = list(1-px, px)))
  }

  # split into K folds
  n <- nrow(data)
  K <- 10
  folds <- sample(x = rep(1:K, each = ceiling(n / K)))[seq_len(n)]

  # elements to be filled
  y_xzw <- y_xz <- px_z <- px_zw <- ey_nest <- list(rep(NA, n), rep(NA, n))

  # x, y data
  y <- data[[Y]]
  x <- data[[X]]

  # cross-fit
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

    # get the val set predictions (needed for nested means)
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
    px_zw[[1 + 0]][tst] <- 1 - px_zw_tst
    px_zw[[1 + 1]][tst] <- px_zw_tst

    if (length(Z) > 0) {

      px_z_tst <- pred_xgb(mod_x_z, data[tst, Z])
      px_z[[1 + 0]][tst] <- 1 - px_z_tst
      px_z[[1 + 1]][tst] <- px_z_tst
    } else {

      px_z[[1 + 0]][tst] <- 1 - mean(x)
      px_z[[1 + 1]][tst] <- mean(x)
    }

    y_xzw[[1 + 0]][tst] <- pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)],
                                    intervention = 0, X = X)
    y_xzw[[1 + 1]][tst] <- pred_xgb(mod_y_xzw, data[tst, c(X, Z, W)],
                                    intervention = 1, X = X)

    if (length(Z) > 0) {

      y_xz[[1 + 0]][tst] <- pred_xgb(mod_y_xz, data[tst, c(X, Z)],
                                     intervention = 0, X = X)
      y_xz[[1 + 1]][tst] <- pred_xgb(mod_y_xz, data[tst, c(X, Z)],
                                     intervention = 1, X = X)
    }

    # nested means are not needed if either Z or W are empty
    if (length(Z) == 0 || length(W) == 0) next

    for (xw in c(0, 1)) {

      xy <- 1 - xw
      if (nested_mean == "wregr") {

        if (scale != "diff") stop("Weighted regression not available for log-risk or log-odds scale.")

        weights <- ifelse(
          x[val] == xy,
          px_z_val[[xy+1]] / px_z_val[[xw+1]] *
            px_zw_val[[xw+1]] / px_zw_val[[xy+1]],
          px_z_val[[xy+1]] / px_z_val[[xw+1]]
        )

        mod_nested <- cv_xgb(data[val, Z], data[val, Y], weights = weights,
                             ...)
        ey_nest[[xy+1]][tst] <- pred_xgb(mod_nested, data[tst, Z])
      } else if (nested_mean == "refit") {

        y_tilde <- pred_xgb(mod_y_xzw, data[val, c(X, Z, W)],
                            intervention = xy, X = X)
        if (scale == "logr") y_tilde <- log(y_tilde)
        else if (scale == "logo") y_tilde <- log(y_tilde / (1 - y_tilde))
        mod_nested <- cv_xgb(data[val, c(X, Z)], y_tilde, ...)
        ey_nest[[xy+1]][tst] <- pred_xgb(mod_nested, data[tst, c(X, Z)],
                                         intervention = xw, X = X)
      }
    }
  }

  list(
    y_xzw = y_xzw,
    y_xz = y_xz,
    px_z = px_z,
    px_zw = px_zw,
    ey_nest = ey_nest
  )
}
