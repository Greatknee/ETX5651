# ETX5651 Layer 1
# 1. Settings fixed before examining significance ------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[1L] else "data.xlsx"
output_dir <- if (length(args) >= 2L) args[2L] else "results_layer1"
L <- 8L                         # two years of quarterly dynamics
selected_lags <- c(1L, 8L)
alpha <- 0.05
lambda <- 0                     # unthresholded GLP Q; see numerical note below
install_missing <- TRUE

packages <- c("testcorr", "urca", "ggplot2", "readxl")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) && install_missing)
  install.packages(missing, repos = "https://cloud.r-project.org")
if (any(!vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
  stop("Install missing packages: ", paste(packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)], collapse = ", "))
if (utils::packageVersion("testcorr") < "0.3.0")
  stop("Please update testcorr to version 0.3.0 or later.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name, directory = output_dir, row.names = FALSE) {
  write.csv(x, file.path(directory, name), row.names = row.names, na = "")
}

# 2. Read, validate and stack the four variables -------------------------------
if (grepl("\\.xlsx$", input_file, ignore.case = TRUE)) {
  d <- as.data.frame(readxl::read_excel(
    input_file, sheet = "Data", na = c("", "NA", "NaN", "."),
    .name_repair = "minimal"
  ))
} else {
  d <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA", "NaN", "."))
}
source_columns <- c("us_gdp_growth_qoq_pct", "us_nfci",
                    "aus_gdp_growth_qoq_pct", "aus_fci")

d$quarter <- toupper(trimws(d$quarter))
if (anyNA(d$quarter) || any(!grepl("^[0-9]{4}Q[1-4]$", d$quarter)))
  stop("quarter must have format YYYYQ1, ..., YYYYQ4.")
year <- as.integer(substr(d$quarter, 1L, 4L))
quarter_number <- as.integer(substr(d$quarter, 6L, 6L))
quarter_id <- 4L * year + quarter_number - 1L
o <- order(quarter_id)
d <- d[o, , drop = FALSE]
quarter_id <- quarter_id[o]
if (anyDuplicated(quarter_id)) stop("Duplicate quarters: resolve them first.")
if (any(diff(quarter_id) != 1L))
  stop("Missing calendar quarters: do not compress gaps before computing lags.")
if (!all(vapply(d[source_columns], is.numeric, logical(1))))
  stop("The four data columns must be numeric.")
missing_table <- data.frame(variable = source_columns,
                            missing = colSums(is.na(d[source_columns])))
write_csv(missing_table, "missing_observations.csv")
if (any(!is.finite(as.matrix(d[source_columns]))))
  stop("Missing/nonfinite values found. Document a treatment; no silent deletion.")

# y_t = (g_US,t, NFCI_US,t, g_AUS,t, FCI_AUS,t)'.
# Y is T x 4: row t is y_t'. Country blocks occupy consecutive columns.
# GDP columns are ALREADY quarter-on-quarter percentages: do not difference,
# log, multiply by 100, or annualize them again. FCI columns remain in levels.
Y <- as.matrix(d[source_columns])
colnames(Y) <- c("US_GDP", "US_NFCI", "AUS_GDP", "AUS_FCI")
if (any(apply(Y, 2L, sd) == 0)) stop("A variable is constant.")
Y_ts <- ts(Y, start = c(quarter_id[1L] %/% 4L, quarter_id[1L] %% 4L + 1L),
           frequency = 4)
dates <- as.Date(sprintf("%04d-%02d-01", quarter_id %/% 4L,
                         3L * (quarter_id %% 4L) + 1L))
write_csv(data.frame(quarter = d$quarter, Y), "stacked_data.csv")

# Source metadata are NOT embedded in the CSV. Complete these before submission.
metadata <- data.frame(
  variable = colnames(Y), original_column = source_columns,
  unit = c("percent q/q, not annualized", "index", "percent q/q, not annualized", "index"),
  transformation = c("as supplied", "level as supplied", "as supplied", "level as supplied"),
  source = c("TODO: verify GDP provider and series ID",
             "US NFCI: verify provider, vintage and quarterly aggregation",
             "TODO: verify GDP provider and series ID",
             "TODO: identify Australian FCI construction and source"),
  definition_to_verify = c("real/nominal; seasonal adjustment; vintage",
                           "quarterly mean/end-quarter; sign convention",
                           "real/nominal; seasonal adjustment; vintage",
                           "quarterly aggregation; whether higher means tighter"))
write_csv(metadata, "data_dictionary_COMPLETE_SOURCES.csv")


# 3. Descriptive statistics, outliers and stationarity --------------------------
summary_table <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
  x <- Y[, j]; z <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  data.frame(variable = colnames(Y)[j], n = length(x), mean = mean(x), sd = sd(x),
             min = min(x), q25 = unname(quantile(x, .25)), median = median(x),
             q75 = unname(quantile(x, .75)), max = max(x),
             skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3)
}))
write_csv(summary_table, "summary_statistics.csv")
outliers <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
  x <- Y[, j]; s <- mad(x)
  score <- if (s > 0) abs(x - median(x)) / s else rep(NA_real_, length(x))
  ix <- which(score > 3.5)
  data.frame(variable = rep(colnames(Y)[j], length(ix)), quarter = d$quarter[ix],
             value = x[ix], robust_z = score[ix])
}))
write_csv(outliers, "outliers_flagged_NOT_REMOVED.csv")

# ADF: intercept, no deterministic trend, BIC selection from up to 8 lags.
# KPSS: null of level stationarity; the truncation lag is fixed at 8 quarters.
# Neither test is a mechanical transformation rule or proof of GLP assumptions.
stationarity <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
  a <- urca::ur.df(Y[, j], type = "drift", lags = L, selectlags = "BIC")
  k <- urca::ur.kpss(Y[, j], type = "mu", use.lag = L)
  astat <- as.numeric(a@teststat[1L, "tau2"])
  acv <- as.numeric(a@cval["tau2", "5pct"])
  kstat <- as.numeric(k@teststat)
  kcv <- as.numeric(k@cval[1L, "5pct"])
  writeLines(capture.output(summary(a)), file.path(output_dir,
             paste0("ADF_", colnames(Y)[j], ".txt")))
  data.frame(variable = colnames(Y)[j], ADF_statistic = astat,
             ADF_cv_5pct = acv, ADF_reject_unit_root = astat < acv,
             KPSS_statistic = kstat, KPSS_cv_5pct = kcv,
             KPSS_reject_level_stationarity = kstat > kcv)
}))
write_csv(stationarity, "stationarity_diagnostics.csv")

pdf(file.path(output_dir, "series.pdf"), width = 11, height = 5,
    useDingbats = FALSE)
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1), las = 1, bty = "l")
for (j in seq_len(ncol(Y))) {
  plot(dates, Y[, j], type = "l", col = "#176caa", lwd = 1.5,
       main = colnames(Y)[j], xlab = "", ylab = "")
  abline(h = 0, col = "grey75", lty = 2)
}
dev.off()

# 4. Helpers for supplementary HAC inference and result tables ----------------
# HAC tests the mean of the lagged cross-product, using Bartlett weights.
# Under zero covariance and constant means, estimating means does not alter
# its first-order null limit. This is a sensitivity analysis, NOT the GLP test.
# It assumes appropriate finite moments and stationary short-memory products;
# it is not a remedy for unit roots, mean breaks or infinite moments.
lag_products <- function(x, y, k) {
  n <- length(x); x <- x - mean(x); y <- y - mean(y)
  if (k >= 0L) x[seq.int(k + 1L, n)] * y[seq_len(n - k)] else
    x[seq_len(n + k)] * y[seq.int(1L - k, n)]
}
decide <- function(p) ifelse(is.na(p), "Not available",
                            ifelse(p < alpha, "Reject H0", "Do not reject H0"))
add_decisions <- function(tab) {
  tab$conclusion_SN <- decide(tab$p_SN)
  tab
}
plot_corr <- function(lags, values, bands, title, subtitle = NULL) {
  yr <- range(c(values, bands, 0), finite = TRUE)
  plot(lags, values, type = "h", lwd = 2, col = "#176caa", ylim = yr,
       main = title, xlab = "Lag (quarters)", ylab = "Correlation", xaxt = "n")
  axis(1, at = lags)
  abline(h = 0, col = "grey60")
  lines(lags, bands[, 1L], col = "#ba5b2d", lty = 2)
  lines(lags, bands[, 2L], col = "#ba5b2d", lty = 2)
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = .15, cex = .7)
  }
}

# Check positive definiteness, including cases where an indefinite matrix
# happens to produce a positive quadratic form. Match GLP overlap normalization.
sn_min_eigen <- function(x, y, endpoint, threshold) {
  if (endpoint < 0L) { tmp <- x; x <- y; y <- tmp }
  m <- abs(endpoint); n <- length(x)
  x <- x - mean(x); y <- y - mean(y)
  E <- matrix(0, n, m)
  for (k in seq_len(m)) E[(k + 1L):n, k] <- x[(k + 1L):n] * y[seq_len(n - k)]
  R <- diag(m)
  if (m > 1L) for (j in seq_len(m - 1L)) for (k in (j + 1L):m) {
    ix <- (k + 1L):n; a <- E[ix, j]; b <- E[ix, k]; ab <- a * b
    den <- sqrt(sum(ab^2))
    if (den > 0 && abs(sum(ab) / den) > threshold)
      R[j, k] <- R[k, j] <- sum(ab) / sqrt(sum(a^2) * sum(b^2))
  }
  min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
}

# 5. All four ACFs, all six pairs, both CCF directions, GLP tests ---------------
# Convention: rho_xy(k) = Corr(x_t, y_{t-k}), as used by R's ccf/testcorr.
# k > 0: y precedes x; k < 0: x precedes y. Lead-lag association alone is
# neither a conditional Granger-causality test nor structural causality.
#
# GLP individual statistic: t_tilde(k) = sum(e_tk)/sqrt(sum(e_tk^2)).
# e_tk = (x_t - xbar)(y_{t-k} - ybar), with overlapping observations only.
# GLP cumulative statistic: Q = t_tilde' R^(-1) t_tilde, GLP eq. (37).
# lambda > 0 gives the thresholded version in GLP eqs. (38)-(39).
# It is NOT generally sum(t_tilde^2), nor the ordinary Ljung-Box statistic.
# We use the authors' testcorr implementation, normal/chi-square null laws.
# m0=1 EXCLUDES lag zero from directional joint tests; lag zero is tested alone.
analyse <- function(X, quarter, directory, draw = TRUE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  n <- nrow(X); nm <- colnames(X)
  if (n <= 3L * L) stop("Insufficient sample for selected lags.")
  C <- cor(X)
  write_csv(C, "correlation_matrix.csv", directory, row.names = TRUE)
  ac_objects <- setNames(lapply(seq_len(ncol(X)), function(j)
    testcorr::ac.test(as.numeric(X[, j]), max.lag = L, m0 = 1,
                      alpha = alpha, lambda = lambda, plot = FALSE)), nm)
  pairs <- combn(nm, 2, simplify = FALSE)
  cc_objects <- setNames(lapply(pairs, function(v)
    testcorr::cc.test(as.numeric(X[, v[1L]]), as.numeric(X[, v[2L]]),
                      max.lag = L, m0 = 1, alpha = alpha,
                      lambda = lambda, plot = FALSE)),
    vapply(pairs, paste, character(1), collapse = "__"))

  individual <- list(); cumulative <- list()
  collect <- function(a, xname, yname, kind) {
    lags <- a$lag
    rho <- if (kind == "ACF") a$ac else a$cc
    ind <- data.frame(type = kind, x = xname, y = yname, lag = lags,
                      n_pairs = n - abs(lags), correlation = rho,
                      statistic_SN = a$ttilde, p_SN = 2 * pnorm(-abs(a$ttilde)),
                      abs_critical_SN = qnorm(1 - alpha / 2),
                      null_band_low = a$rcb[, 1L], null_band_high = a$rcb[, 2L]
                      )
    ind$null <- sprintf("Corr(%s_t,%s_{t-(%d)})=0", xname, yname, lags)
    ind$timing <- if (kind == "ACF") "Own-series persistence" else
      ifelse(lags > 0, paste(yname, "precedes", xname),
             ifelse(lags < 0, paste(xname, "precedes", yname), "Contemporaneous"))
    ind$primary_lag <- abs(lags) %in% c(0L, selected_lags)
    qc <- a$qtilde
    # Thresholding can yield an indefinite matrix in a finite sample. Such a
    # negative quadratic form is not a valid chi-square statistic; do not clip it.
    eigen_min <- vapply(a$lagc, function(k)
      sn_min_eigen(X[, xname], X[, yname], k, lambda), numeric(1))
    valid <- is.finite(qc) & qc >= 0 & is.finite(eigen_min) & eigen_min > 1e-8
    if (any(!valid)) warning("Invalid cumulative statistic for ", xname, "/", yname)
    df <- abs(a$lagc)   # m0=1 -> df=m, separately for positive and negative lags
    joint <- data.frame(type = kind, x = xname, y = yname,
                        endpoint_lag = a$lagc, df = df, statistic_SN = qc,
                        matrix_min_eigenvalue = eigen_min, valid_matrix = valid,
                        critical_SN = qchisq(1 - alpha, df),
                        p_SN = ifelse(valid, pchisq(qc, df, lower.tail = FALSE), NA_real_))
    joint$null <- vapply(a$lagc, function(k) {
      tested <- if (k > 0) paste0("1:", k) else paste0(k, ":-1")
      paste0("Corr(", xname, "_t,", yname, "_{t-k})=0 for all k in ", tested)
    }, character(1))
    joint$timing <- if (kind == "ACF") "Own-series persistence" else
      ifelse(a$lagc > 0, paste(yname, "precedes", xname), paste(xname, "precedes", yname))
    joint$primary_range <- abs(a$lagc) == L
    list(individual = ind, cumulative = joint)
  }
  for (j in nm) {
    a <- collect(ac_objects[[j]], j, j, "ACF")
    individual[[length(individual) + 1L]] <- a$individual
    cumulative[[length(cumulative) + 1L]] <- a$cumulative
  }
  for (i in seq_along(pairs)) {
    v <- pairs[[i]]; a <- collect(cc_objects[[i]], v[1L], v[2L], "CCF")
    individual[[length(individual) + 1L]] <- a$individual
    cumulative[[length(cumulative) + 1L]] <- a$cumulative
  }
  I <- add_decisions(do.call(rbind, individual))
  Q <- add_decisions(do.call(rbind, cumulative))


  write_csv(I, "individual_tests_all_lags.csv", directory)
  write_csv(I[I$primary_lag, ], "individual_tests_selected_lags.csv", directory)
  write_csv(I[I$lag == 0L, ], "contemporaneous_tests.csv", directory)
  write_csv(Q, "cumulative_SN_tests_all_ranges.csv", directory)
  write_csv(Q[Q$primary_range, ], "cumulative_SN_tests_lags1to8.csv", directory)
  saveRDS(list(Y = X, quarter = quarter, ACF = ac_objects, CCF = cc_objects,
               individual = I, cumulative = Q), file.path(directory, "analysis.rds"))

  if (draw) {
    pdf(file.path(directory, "ACF_SN.pdf"), width = 11, height = 6, useDingbats = FALSE)
    par(mfrow = c(2, 2), mar = c(4, 4, 4, 1), bty = "l", las = 1)
    for (j in nm) {
      a <- ac_objects[[j]]
      plot_corr(a$lag, a$ac, a$rcb, j)
    }
    dev.off()
    pdf(file.path(directory, "CCF_both_directions_SN.pdf"), width = 12, height = 9,
        useDingbats = FALSE)
    par(mfrow = c(3, 2), mar = c(4, 4, 4, 1), bty = "l", las = 1)
    for (i in seq_along(pairs)) {
      v <- pairs[[i]]; a <- cc_objects[[i]]
      plot_corr(a$lag, a$cc, a$rcb, paste(v, collapse = " / "),
                paste0("k > 0: ", v[2L], " precedes ", v[1L], "; k < 0: reverse"))
    }
    dev.off()
  }
  list(correlation = C, individual = I, cumulative = Q)
}

full <- analyse(Y, d$quarter, file.path(output_dir, "full_sample"))

