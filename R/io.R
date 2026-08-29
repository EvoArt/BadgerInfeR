# Turning R data into the binary's CSV schema, and its CSVs back into R.

# The binary's reader is deliberately minimal: a header row, then numeric
# fields, with empty or unparseable cells meaning "not tested". It does not
# handle quoting, so everything written here must be plain numeric.
write_test_matrix <- function(test_mat, path) {
  if (is.data.frame(test_mat)) test_mat <- as.matrix(test_mat)
  if (!is.matrix(test_mat)) {
    stop("`test_mat` must be a matrix or data frame.")
  }
  if (ncol(test_mat) < 9L) {
    stop("`test_mat` needs at least 9 columns ",
         "(time, id, captured, test1..test6); got ", ncol(test_mat), ".")
  }
  m <- test_mat[, 1:9, drop = FALSE]
  storage.mode(m) <- "double"

  if (any(is.na(m[, 1:3]))) {
    stop("`time`, `id` and `captured` must not contain NA.")
  }
  if (any(m[, 1] < 1)) stop("`time` must be a positive integer timestep.")
  if (!all(m[, 3] %in% c(0, 1))) stop("`captured` must be 0 or 1.")

  df <- as.data.frame(m)
  names(df) <- c("time", "id", "captured", paste0("test", 1:6))
  # NA -> empty field: that is what the binary reads as "not tested".
  utils::write.csv(df, path, row.names = FALSE, na = "", quote = FALSE)
  invisible(path)
}

# Prevalence is reported per timestep; the bundles want calendar years too.
add_calendar_year <- function(df, start_year, seasons) {
  if (is.null(df) || start_year <= 0) return(df)
  df$calendar_year <- start_year + (df$time - 1L) %/% seasons
  df
}

add_year_calendar <- function(df, start_year) {
  if (is.null(df) || start_year <= 0) return(df)
  df$calendar_year <- start_year + (df$year_index - 1L)
  df
}

# Se/Sp as a tidy frame, matching the sensitivity_specificity.csv the bundles
# write.
#
# `estimated` comes from the MODEL SPEC, not from whether Se/Sp rows are
# present: the binary writes Se[k]/Sp[k] rows in both cases -- under a fixed
# model they are simply the Table 1 constants echoed back. Inferring
# "estimated" from their presence labels every fixed fit as estimated, which
# is exactly the bug this comment exists to prevent recurring.
extract_sesp <- function(params, infer_sesp) {
  if (is.null(params)) return(NULL)
  valcol <- if ("estimate" %in% names(params)) "estimate" else "mean"
  pick <- function(pat) {
    i <- grep(pat, params$parameter)
    if (length(i) == 0L) return(numeric(0))
    params[[valcol]][i]
  }
  se <- pick("^Se\\[")
  sp <- pick("^Sp\\[")
  if (length(se) == 0L || length(sp) == 0L) {
    se <- BADGER_SE
    sp <- BADGER_SP
  }
  estimated <- isTRUE(infer_sesp)
  data.frame(
    test      = paste0("test_", seq_along(se)),
    test_name = BADGER_TESTS[seq_along(se)],
    Se        = se,
    Sp        = sp,
    Youden    = se + sp - 1,
    estimated = estimated,
    stringsAsFactors = FALSE
  )
}
