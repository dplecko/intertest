measure_spec <- function(spec = c("ia", "xspec")) {

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
