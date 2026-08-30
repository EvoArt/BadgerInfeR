# Running several chains concurrently.
#
# The binary is single-threaded and holds no shared state, so process-level
# fan-out is the whole mechanism -- there is no threading to do inside it. On
# this project's 8-core box four concurrent processes measured 3.96x over four
# serial runs (`sql e2e badgerinfer/results/timings_4chain_parallel_cli.csv`),
# which is close to ideal.
#
# PSOCK, not `parallel::mclapply`: mclapply forks, and there is no fork on
# Windows, which is where this package is developed and run.

# One chain, in a worker. Returns the fit plus its own elapsed time.
#
# `seed` is passed explicitly rather than derived inside the worker: four chains
# at the same seed are four IDENTICAL chains, which produces an Rhat of exactly
# 1 and looks like perfect convergence while meaning nothing.
.run_one_chain <- function(seed, args) {
  t0 <- Sys.time()
  args$seed <- seed
  fit <- do.call(badger_hmm_one, args)
  fit$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit
}

# Stack a per-chain data frame, tagging each row with its chain index.
#
# Returns NULL when no chain produced the element (e.g. `draws` under
# method = "map"), so the combined fit keeps the same NULL-means-absent contract
# a single-chain fit has.
.stack_chain_frames <- function(fits, element) {
  parts <- lapply(seq_along(fits), function(i) {
    d <- fits[[i]][[element]]
    if (is.null(d) || nrow(d) == 0L) return(NULL)
    d$chain <- i
    d
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(NULL)
  do.call(rbind, parts)
}

# Split-Rhat and bulk ESS, computed from the pooled draws.
#
# Worth having because the reason to run four chains at all is to find out
# whether they agree; returning the draws without a convergence measure leaves
# the caller to rediscover this. Uses the rank-normalised split-Rhat of
# Vehtari et al. (2021), which is what `posterior::rhat` reports.
.chain_diagnostics <- function(draws) {
  if (is.null(draws) || !"chain" %in% names(draws)) return(NULL)
  pars <- setdiff(names(draws), c("chain", "draw", "iteration"))
  pars <- pars[vapply(draws[pars], is.numeric, logical(1))]
  if (length(pars) == 0L) return(NULL)

  chains <- split(seq_len(nrow(draws)), draws$chain)
  m <- length(chains)
  if (m < 2L) return(NULL)
  n <- min(vapply(chains, length, integer(1)))
  if (n < 4L) return(NULL)

  stats <- lapply(pars, function(p) {
    # Columns with no variance (a fixed Se/Sp echoed back) give 0/0. Report NA
    # rather than NaN: "not applicable" is the honest answer for a constant,
    # and it keeps a fixed-Se/Sp fit from looking like a convergence failure.
    x <- vapply(chains, function(idx) draws[[p]][idx[seq_len(n)]], numeric(n))
    if (stats::sd(as.vector(x)) < .Machine$double.eps^0.5) {
      return(c(rhat = NA_real_, ess_bulk = NA_real_))
    }
    # Split each chain in half: split-Rhat catches a chain that is drifting
    # within itself, which plain Rhat does not.
    half <- n %/% 2L
    sp <- cbind(x[seq_len(half), , drop = FALSE],
                x[(n - half + 1L):n, , drop = FALSE])
    # Rank-normalise, then Rhat on the normal scores.
    r <- matrix(rank(as.vector(sp), ties.method = "average"), nrow = half)
    z <- stats::qnorm((r - 0.375) / (length(r) + 0.25))
    c(rhat = .rhat_from_matrix(z), ess_bulk = .ess_from_matrix(z))
  })
  out <- as.data.frame(do.call(rbind, stats))
  data.frame(parameter = pars, rhat = out$rhat, ess_bulk = out$ess_bulk,
             stringsAsFactors = FALSE)
}

# Rhat for an n x m matrix of draws (columns are chains).
.rhat_from_matrix <- function(x) {
  n <- nrow(x); m <- ncol(x)
  cm <- colMeans(x)
  B <- n * stats::var(cm)
  W <- mean(apply(x, 2, stats::var))
  if (W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

# Bulk ESS for an n x m matrix, from the summed autocorrelations. Geyer's
# initial positive sequence: truncate at the first negative pair, which is what
# keeps the estimate from running away on a noisy tail.
.ess_from_matrix <- function(x) {
  n <- nrow(x); m <- ncol(x)
  W <- mean(apply(x, 2, stats::var))
  if (W <= 0) return(NA_real_)
  acov <- vapply(seq_len(m), function(j) {
    v <- x[, j] - mean(x[, j])
    a <- stats::acf(v, lag.max = n - 1L, type = "covariance",
                    plot = FALSE, demean = FALSE)$acf[, , 1]
    a
  }, numeric(n))
  rho <- 1 - (W - rowMeans(acov)) / W
  # rho is indexed from lag 0, so rho[k + 1] is lag k. Pairs are (1,2), (3,4),
  # ... and the bound must leave BOTH elements of the pair in range: reading
  # rho[n + 1] would silently return NA and end the sum one pair early.
  t <- 1L
  s <- 0
  while (t + 2L <= n) {
    p <- rho[t + 1L] + rho[t + 2L]
    if (is.na(p) || p < 0) break
    s <- s + p
    t <- t + 2L
  }
  ess <- n * m / (1 + 2 * s)
  min(ess, n * m)
}

# Merge per-chain fits into one badger_fit.
#
# The combined object keeps the class and element names of a single-chain fit,
# so `print()`, the sql-e2e scripts and anything else downstream keep working;
# what changes is that the row-wise elements gain a `chain` column and
# `$parameters` is recomputed from the pooled draws.
.combine_chain_fits <- function(fits, settings) {
  stacked <- c("draws", "prevalence_draws", "trajectories")
  out <- list()

  # Non-draw outputs are per-chain summaries of the same posterior. Take chain
  # 1's, which is a valid summary, rather than averaging point estimates across
  # chains -- averaging quantiles is not a quantile of anything.
  for (el in c("prevalence", "p_infected", "year_effects", "sesp")) {
    out[[el]] <- fits[[1]][[el]]
  }
  for (el in stacked) out[[el]] <- .stack_chain_frames(fits, el)

  out$draws <- .stack_chain_frames(fits, "draws")
  diag <- .chain_diagnostics(out$draws)

  # Recompute the parameter table from the pooled draws so the means and
  # intervals reflect all chains, not just the first.
  params <- fits[[1]]$parameters
  if (!is.null(out$draws) && !is.null(params) && "parameter" %in% names(params)) {
    pooled <- .pool_parameters(out$draws, params)
    if (!is.null(pooled)) params <- pooled
  }
  if (!is.null(params) && !is.null(diag)) {
    params <- merge(params, diag, by = "parameter", all.x = TRUE,
                    sort = FALSE)
  }
  out$parameters <- params
  out$diagnostics <- diag

  out$chain_times <- vapply(fits, function(f) f$elapsed %||% NA_real_,
                            numeric(1))
  out$log <- lapply(fits, `[[`, "log")
  out$settings <- settings
  class(out) <- "badger_fit"
  out
}

# Recompute mean/sd/quantile columns from the pooled draws, preserving whatever
# column names the binary's own parameters.csv used.
.pool_parameters <- function(draws, template) {
  pars <- intersect(template$parameter, names(draws))
  if (length(pars) == 0L) return(NULL)
  cols <- names(template)
  res <- template[match(pars, template$parameter), , drop = FALSE]
  for (p in seq_along(pars)) {
    x <- draws[[pars[p]]]
    if (!is.numeric(x)) next
    if ("mean" %in% cols)     res$mean[p]     <- mean(x)
    if ("estimate" %in% cols) res$estimate[p] <- mean(x)
    if ("sd" %in% cols)       res$sd[p]       <- stats::sd(x)
    q <- stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE, na.rm = TRUE)
    if ("q2.5" %in% cols)   res$`q2.5`[p]  <- q[1]
    if ("median" %in% cols) res$median[p]  <- q[2]
    if ("q97.5" %in% cols)  res$`q97.5`[p] <- q[3]
  }
  res
}

`%||%` <- function(a, b) if (is.null(a)) b else a
