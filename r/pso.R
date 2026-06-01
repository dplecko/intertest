pso <- function(cfit, data, X, Z, W, Y, scale = c("diff", "logr", "logo"), ...) {

  scale  <- match.arg(scale)
  link   <- switch(scale,
                   diff = function(m) m,
                   logr = function(m) log(m),
                   logo = function(m) log(m / (1 - m)))
  dlink  <- switch(scale,
                   diff = function(m) rep(1, length(m)),
                   logr = function(m) 1 / m,
                   logo = function(m) 1 / (m * (1 - m)))

  n      <- nrow(data)
  y_xzw  <- cfit$y_xzw
  y_xz   <- cfit$y_xz
  px_z   <- cfit$px_z
  px_zw  <- cfit$px_zw
  ey_nest <- cfit$ey_nest
  y      <- data[[Y]]
  x      <- data[[X]]

  phi <- list(list(list(list(), list()), list(list(), list())),
              list(list(list(), list()), list(list(), list())))
  for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1))
    phi[[xz+1]][[xw+1]][[xy+1]] <- rep(NA, n)

  if (length(Z) == 0 & length(W) == 0) {

    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      mu <- mean(y[x == xy])
      phi[[xz+1]][[xw+1]][[xy+1]] <- link(mu) +
        (x == xy) / mean(x == xy) * dlink(mu) * (y - mu)
    }

  } else if (length(Z) == 0) {

    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      mu_hat <- y_xzw[[xy+1]]
      phi[[xz+1]][[xw+1]][[xy+1]] <-
        (x == xy) / mean(x == xw) *
        dlink(mu_hat) * (y - mu_hat) +
        (x == xw) / mean(x == xw) * link(mu_hat)
    }

  } else if (length(W) == 0) {

    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      mu_hat <- y_xzw[[xy+1]]
      phi[[xz+1]][[xw+1]][[xy+1]] <-
        (x == xy) / mean(x == xz) *
        dlink(mu_hat) * (y - mu_hat) +
        (x == xz) / mean(x == xz) * link(mu_hat)
    }

  } else {

    for (xz in c(0, 1)) for (xw in c(0, 1)) for (xy in c(0, 1)) {
      if (xw == xy) {
        mu_hat <- y_xz[[xy + 1]]
        phi[[xz+1]][[xw+1]][[xy+1]] <-
          (x == xy) / mean(x == xz) * px_z[[xz+1]] / px_z[[xw+1]] *
          dlink(mu_hat) * (y - mu_hat) +
          (x == xz) / mean(x == xz) * link(mu_hat)
      } else {
        mu_hat <- y_xzw[[xy + 1]]
        phi[[xz+1]][[xw+1]][[xy+1]] <-
          (x == xy) / mean(x == xz) *
          dlink(mu_hat) * (y - mu_hat) *
          px_zw[[xw+1]] / px_zw[[xy+1]] *
          px_z[[xz+1]] / px_z[[xw+1]] +
          (x == xw) / mean(x == xz) * px_z[[xz+1]] / px_z[[xw+1]] *
          (link(mu_hat) - ey_nest[[xy+1]]) +
          (x == xz) / mean(x == xz) * ey_nest[[xy+1]]
      }
    }
  }

  phi
}
