# gen_from_scm: generate synthetic data from a named SCM class.
# Used in statistical tests. Classes A-G, I-K are included; H is omitted
# (has a coding issue with null Z references).

gen_from_scm <- function(sclass, n, X_to_W = "random", X_to_Y = "random") {

  gen_X <- function(Z, option) {
    if (option == "random") {
      return(rbinom(n, 1, prob = plogis(Z %*% c(0.3, -0.2, 0.5) + 0.2 * (Z[,1]^2))))
    } else if (option == 0) {
      return(rep(0, n))
    } else if (option == 1) {
      return(rep(1, n))
    }
  }

  prob_y <- NULL

  if (sclass == "A") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- Z %*% matrix(c(0.4, 0.1, -0.3, 0.2, -0.1, 0.3, 0.3, -0.2, 0.1), 3, 3) +
      X_W * matrix(c(0.5, 0.4, 0.3), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    Y <- W %*% c(0.5, 0.4, 0.3) + Z %*% c(0.2, 0.1, 0.4) + X_Y * 0.7 +
      X_Y * W[,1] * 0.2 + rnorm(n)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "B") {

    Z <- cbind(rexp(n, rate = 1), rnorm(n, 5, 1), runif(n, -2, 2))
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- Z %*% matrix(c(0.3, -0.5, 0.2, -0.1, 0.3, 0.1, 0.2, 0.2, -0.3), 3, 3) +
      X_W * matrix(c(0.2, 0.1, 0.4), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    Y <- W %*% c(0.4, 0.3, 0.2) + Z %*% c(0.1, 0.3, 0.2) + X_Y * 0.4 + rnorm(n) +
      (W %*% c(0.1, -0.3, -0.3)) * (Z %*% c(0.1, -0.2, 0.2))
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "C") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- matrix(c((Z[,1]^2) * 0.3, Z[,2] * 0.5, X_W * 0.4), n, 3) +
      matrix(rnorm(n * 3), n, 3)
    Y <- W %*% c(0.3, 0.2, 0.1) + Z %*% c(0.2, 0.1, 0.3) + rnorm(n)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = FALSE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "D") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- Z %*% matrix(c(0.4, 0.1, -0.3, 0.2, -0.1, 0.3, 0.3, -0.2, 0.1), 3, 3) +
      X_W * matrix(c(0.5, 0.4, 0.3), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    Y <- W %*% c(0.5, 0.4, 0.3) + Z %*% c(0.2, 0.1, 0.4) +
      X_Y * Z[,1] * 0.3 + rnorm(n)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "E") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- matrix(c((Z[,1]^2) * 0.3, Z[,2] * 0.5, X_W * 0.4), n, 3) +
      matrix(rnorm(n * 3), n, 3)
    Y <- W %*% c(0.4, 0.3, 0.2) + Z %*% c(0.2, 0.1, 0.3) +
      X_Y * Z[,1] * W[,3] * 0.5 + Z[,2] * W[,3] * (-0.4) +
      X_Y * Z[, 3] + rnorm(n)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = TRUE)
    )

  } else if (sclass == "F") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- Z %*% matrix(c(0.3, -0.2, 0.4, 0.5, -0.3, 0.1, 0.2, 0.1, -0.2), 3, 3) +
      X_W * matrix(c(0.4, 0.3, 0.5), n, 3, byrow = TRUE) + matrix(rnorm(n * 3), n, 3)
    logit_Y <- W %*% c(0.3, 0.4, 0.2) + Z %*% c(0.1, 0.3, 0.2) + X_Y * 0.6
    prob_y <- plogis(logit_Y)
    Y <- rbinom(n, 1, prob = prob_y)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = TRUE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "G") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- matrix(c((Z[,1]^2) * 0.2, Z[,2] * 0.4, X_W * 0.5), n, 3) +
      matrix(rnorm(n * 3), n, 3)
    risk_ratio <- exp(W %*% c(0.2, 0.3, 0.1) + Z %*% c(0.2, -0.1, 0.3) + X_Y * 0.5)
    prob_y <- risk_ratio / (1 + risk_ratio)
    Y <- rbinom(n, 1, prob = prob_y)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = FALSE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = TRUE)
    )

  } else if (sclass == "I") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- NULL
    Y <- Z %*% c(0.2, 0.1, 0.3) +
      X_Y * Z[,1] * 0.5 + Z[,2] * (-0.4) +
      X_Y * Z[, 3]^2 + rnorm(n)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "J") {

    Z <- NULL
    X <- rbinom(n, 1, 0.5)
    X_W <- if (X_to_W == "random") X else gen_X(Z, X_to_W)
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- X_W * matrix(c(0.4, 0.3, 0.5), n, 3, byrow = TRUE) +
      matrix(rnorm(n * 3), n, 3)
    logit_Y <- W %*% c(0.3, 0.4, 0.2) + W[, 1]^2 * X_Y + X_Y * 0.6
    prob_y <- plogis(logit_Y)
    Y <- rbinom(n, 1, prob = prob_y)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = FALSE),
      data.frame(measure = "DE x IE", gt = TRUE),
      data.frame(measure = "DE x SE", gt = FALSE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )

  } else if (sclass == "K") {

    Z <- matrix(rnorm(n * 3), n, 3)
    X <- gen_X(Z, "random")
    X_Y <- if (X_to_Y == "random") X else gen_X(Z, X_to_Y)
    W <- NULL
    logit_Y <- Z %*% c(0.1, 0.3, 0.2) + Z[, 2]^2 * X_Y + X_Y * 0.6
    prob_y <- plogis(logit_Y)
    Y <- rbinom(n, 1, prob = prob_y)
    gt <- rbind(
      data.frame(measure = "TE x SE", gt = TRUE),
      data.frame(measure = "DE x IE", gt = FALSE),
      data.frame(measure = "DE x SE", gt = TRUE),
      data.frame(measure = "IE x SE", gt = FALSE),
      data.frame(measure = "DE x IE x SE", gt = FALSE)
    )
  }

  col_nms <- function(lab, A)
    if (is.null(A)) return(NULL) else paste0(lab, seq_len(ncol(A)))

  nms <- c("X", col_nms("Z", Z), col_nms("W", W), "Y")

  data <- as.data.frame(cbind(X, Z, W, Y))
  names(data) <- nms

  ret <- list(
    data    = data,
    mapping = list(X = "X", Z = col_nms("Z", Z), W = col_nms("W", W), Y = "Y"),
    gt      = gt
  )

  if (!is.null(prob_y)) ret$prob_y <- prob_y
  ret
}
