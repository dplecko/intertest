#' @importFrom stats model.matrix
preproc_data <- function(data, X, Z, W, Y) {

  SFM <- list(X = X, Z = Z, W = W, Y = Y)
  encoded <- character(0)

  for (cmp in c("X", "Z", "W", "Y")) {

    for (var in SFM[[cmp]]) {

      if (is.numeric(data[[var]]) || is.integer(data[[var]])) next
      if (is.logical(data[[var]])) {

        data[[var]] <- as.integer(data[[var]])
        encoded <- c(encoded, var)
        next
      }

      if (is.character(data[[var]])) data[[var]] <- as.factor(data[[var]])

      # only factors at this point
      if (length(levels(data[[var]])) == 2) {

        data[[var]] <- as.integer(data[[var]] == levels(data[[var]])[1])
        encoded <- c(encoded, var)
      } else {

        enc_mat <- model.matrix(~ ., data = data.frame(var = data[[var]]))[, -1]
        colnames(enc_mat) <- paste0(var, seq_len(ncol(enc_mat)))

        # remove the factor column
        SFM[[cmp]] <- setdiff(SFM[[cmp]], var)
        data[[var]] <- NULL

        # add the encoded data
        SFM[[cmp]] <- c(SFM[[cmp]], colnames(enc_mat))
        data <- cbind(data, enc_mat)
        encoded <- c(encoded, var)
      }
    }
  }

  list(data = data, sfm = SFM, encoded = encoded)
}
