
library(intertest)
library(data.table)

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
source(file.path(root, "scripts", "helpers-load-data.R"))
set.seed(2024)

ia_inspect <- function(dataset) {

  if (length(dataset) > 1) {

    return(
      do.call(rbind, lapply(dataset, ia_inspect))
    )
  }

  dat <- load_data(dataset)
  dat$data <- as.data.frame(dat$data)
  dat$data <- dat$data[, unlist(dat$sfm[c("X", "Z", "W", "Y")])]
  dat <- intertest:::preproc_data(dat$data, dat$sfm$X, dat$sfm$Z, dat$sfm$W, dat$sfm$Y)
  data <- dat$data
  sfm <- dat$sfm

  # narrow down the quantities
  target_ia <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
  # check if the outcome is binary
  binary <- length(unique(data[[sfm$Y]])) == 2

  itw <- NULL
  for (scale in unique(c("difference", if (binary) "log-risk"))) {

    iat <- inter_test(data, sfm$X, sfm$Z, sfm$W, sfm$Y, scale = scale)
    iat <- as.data.table(iat$results)
    iat <- iat[measure %in% target_ia]
    itw <- rbind(itw, iat)
  }

  itw[, dataset := dataset]
  return(itw)
}

ia_itw_latex <- function(itw) {
  
  measures <- c("TE x SE", "DE x IE", "DE x SE", "IE x SE", "DE x IE x SE")
  
  ds_names <- list(
    compas = "COMPAS", census = "Census 2018", credit = "UCI Credit",
    mimic = "MIMIC-IV", heloc = "HELOC", adult = "UCI Adult",
    wine = "UCI Wine Quality", bank = "UCI Bank Marketing",
    obesity = "UCI Obesity", diabetes = "BRFSS Diabetes"
  )
  glyph <- function(p) if (!length(p)) "--" else if (p < 0.05) "$\\bullet$" else "$\\circ$"
  
  for (dts in unique(itw$dataset)) {
    d <- itw[dataset == dts]
    n_samp <- nrow(load_data(dts)$data)
    cells <- vapply(measures, function(m) {
      ps <- d[measure == m]$p_value      # 1 (cts) or 2 (binary) values
      paste(vapply(ps, glyph, character(1)), collapse = "/")
    }, character(1))
    cat(paste(c(ds_names[[dts]], n_samp, cells), collapse = " & "),
        "\\\\ \\hline \n")
  }
}

datasets <- c("compas", "census", "credit", "mimic", "heloc", "adult", "wine",
              "bank", "obesity", "diabetes")

# run the interaction testing on the datasets
itw <- ia_inspect(datasets)

# get the latex table rows
ia_itw_latex(itw)

# number of rejected hypotheses with Benjamini-Hochberg
sum(p.adjust(itw$p_value, method = "BH") < 0.05) # 18
