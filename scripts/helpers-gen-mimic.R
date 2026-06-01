# Generate data/dat-miiv-death.RData from MIMIC-IV via ricu.
#
# Prerequisites:
#   - ricu installed and configured with MIMIC-IV access
#   - config/concept-dict.json in the project root (already present)
#
# Run once from the project root:
#   Rscript scripts/helpers-gen-mimic.R

root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))
Sys.setenv("RICU_CONFIG_PATH" = file.path(root, "config"))
library(ricu)
library(data.table)

# ── ricu callbacks (referenced in config/concept-dict.json) ──────────────────

miiv_adm_epi_cb <- function(x, val_var, extra_var, env) {

  x <- merge(x, env$icustays[, c("stay_id", "intime")], by = "stay_id")
  x <- as_id_tbl(x, id_vars = "subject_id")
  x <- data.table::setorderv(x, cols = c("subject_id", "intime"))
  x[, adm_episode := seq_along(stay_id), by = "subject_id"]
  x <- as_id_tbl(x, id_vars = "stay_id")
  x[, subject_id := adm_episode]
}

acute_dayone <- function(sofa, ...) {

  ind_var <- index_var(sofa)
  sofa    <- sofa[get(ind_var) == hours(24L)]
  sofa[, acu_24 := sofa]
  sofa[, c(ind_var, "sofa") := NULL]
  sofa
}

lact_dayone <- function(...) {

  tbl     <- list(...)[[1]]
  ind_var <- index_var(tbl)
  tbl     <- tbl[get(ind_var) >= hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(lact_24 = max(lact)), by = c(id_vars(tbl))]
}

ast_dayone <- function(...) {

  tbl     <- list(...)[[1]]
  ind_var <- index_var(tbl)
  tbl     <- tbl[get(ind_var) >= hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(ast_24 = max(ast)), by = c(id_vars(tbl))]
}

pafi_dayone <- function(...) {

  tbl     <- list(...)[[1]]
  ind_var <- index_var(tbl)
  tbl     <- tbl[get(ind_var) >= hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(pafi_24 = min(pafi)), by = c(id_vars(tbl))]
}

# ── load and process MIMIC-IV ─────────────────────────────────────────────────

sel_coh     <- load_concepts(c("adm_episode", "age", "death"), "miiv", verbose = FALSE)
patient_ids <- id_col(sel_coh[adm_episode == 1 & age >= 18])

dat <- load_concepts(
  c("death", "age", "acu_24", "charlson", "adm_diag", "elective", "sex", "race"),
  "miiv", patient_ids = patient_ids, verbose = FALSE
)

dat <- dat[race %in% c("African American", "Caucasian")]
dat[, sex      := as.integer(sex == "Male")]
dat[, majority := as.integer(race == "Caucasian")]
dat[, race     := NULL]
dat[, c(index_var(dat)) := NULL]

dat <- dat[majority %in% c(0L, 1L)]  # Caucasian vs African American only

# map adm_diag (string) to integer diag_index
diag_dt <- structure(
  list(diag_index = 0:17,
       adm_diag   = c("MED", "CMED", "NMED", "OMED", "PSYCH", "GU", "TRAUM",
                      "ENT", "CSURG", "NSURG", "ORTHO", "PSURG", "SURG",
                      "TSURG", "VSURG", "GYN", "OBS", "DENT")),
  row.names = c(NA, -18L), class = c("data.table", "data.frame")
)

dat <- merge(dat, diag_dt, by = "adm_diag", all.x = TRUE)
dat[, adm_diag := NULL]

imp_lst <- list(age = 65, acu_24 = median(dat$acu_24, na.rm = TRUE),
                death = FALSE, elective = NA, charlson = 0)
for (i in seq_len(ncol(dat))) {
  var <- names(dat)[i]
  if (any(is.na(dat[[var]])) && !is.null(imp_lst[[var]]))
    dat[is.na(get(var)), c(var) := imp_lst[[var]]]
}

dat <- dat[complete.cases(dat)]

out_path <- file.path(root, "data", "dat-miiv-death.RData")
save(dat, file = out_path)
cat("Saved", nrow(dat), "rows to", out_path, "\n")
