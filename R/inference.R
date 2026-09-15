# The one function scripts call. Everything else in this package supports it.

#' Fit the badger diagnostic HMM, on one or several chains
#'
#' The user-facing entry point. With `chains = 1` (the default) this is exactly
#' the single-chain fit the package has always done; with `chains > 1` it fans
#' the chains out across concurrent worker processes and combines them into one
#' `badger_fit`.
#'
#' @inheritParams badger_hmm_one
#' @param seed RNG seed. With `chains > 1` this is the seed of the FIRST chain;
#'   chain `k` uses `seed + k - 1`. The chains must differ: four chains at one
#'   seed are four identical chains, which gives an Rhat of exactly 1 and tells
#'   you nothing.
#' @param chains Number of chains to run. Defaults to 1.
#'
#'   The binary is single-threaded and holds no shared state, so chains are run
#'   as concurrent PSOCK worker processes -- that is the whole mechanism, and
#'   there is no threading to do inside the binary. On 8 cores, 4 chains
#'   measured 3.96x over running them one after another. Each worker loads the
#'   cohort, so memory scales with `chains`; on a machine short of RAM keep it
#'   at 1 or 2, or use `cores` to stagger them.
#' @param cores Workers to run concurrently. Defaults to `chains`. Lower it to
#'   run, say, 4 chains two at a time on a memory-constrained box.
#'
#' @return A `badger_fit`, as [badger_hmm_one()] returns.
#'
#'   With `chains > 1` the shape is unchanged -- it is still one `badger_fit`,
#'   so `print()` and existing scripts keep working -- but the per-draw elements
#'   (`draws`, `prevalence_draws`, `trajectories`) gain a `chain` column and
#'   hold every chain's draws; `parameters` is recomputed from the pooled draws
#'   and carries `rhat` and `ess_bulk`; `diagnostics` holds those per parameter;
#'   `chain_times` holds each chain's wall time; and `log` becomes a list of the
#'   per-chain logs. `chains = 1` returns exactly what it always did.
#' @export
badger_hmm <- function(test_mat,
                       model = "fixed",
                       method = c("nutsthenhmc", "nuts", "map", "hmc"),
                       draws = 1000L,
                       warmup = 1000L,
                       traj_draws = 0L,
                       metric = NULL,
                       write_metric = NULL,
                       target_acc = 0.8,
                       seed = 1L,
                       chains = 1L,
                       cores = NULL,
                       seasons = 4L,
                       repeat_captures = c("stack", "pool", "last"),
                       year_process = c("rw1c", "rw1", "iid", "rw2", "none",
                                        "shrunk"),
                       year_weight = TRUE,
                       accuracy_prior = c("logitnormal", "beta"),
                       penalty = TRUE,
                       start_year = 0L,
                       workdir = tempfile("badgerinfer"),
                       keep_workdir = FALSE) {
  method <- match.arg(method)
  repeat_captures <- match.arg(repeat_captures)
  year_process <- match.arg(year_process)
  accuracy_prior <- match.arg(accuracy_prior)
  chains <- as.integer(chains)
  if (is.na(chains) || chains < 1L) stop("`chains` must be a positive integer.")

  # MAP is an optimisation, not a sampler: every chain would start from the
  # same place and return the same mode. Refuse rather than silently spending
  # four processes on identical work.
  if (method == "map" && chains > 1L) {
    stop("`chains > 1` is meaningless for method = \"map\": it is an ",
         "optimisation, not a sampler, and every chain would return the same ",
         "mode. Use chains = 1.")
  }

  args <- list(test_mat = test_mat, model = model, method = method,
               draws = draws, warmup = warmup, traj_draws = traj_draws,
               metric = metric, write_metric = write_metric,
               target_acc = target_acc, seasons = seasons,
               repeat_captures = repeat_captures, year_process = year_process,
               year_weight = year_weight, accuracy_prior = accuracy_prior,
               penalty = penalty, start_year = start_year,
               keep_workdir = keep_workdir)

  if (chains == 1L) {
    # The old path, unchanged, including the caller's own `workdir`. Only the
    # multi-chain branch needs a workdir per chain.
    fit <- do.call(badger_hmm_one,
                   c(args, list(seed = seed, workdir = workdir)))
    fit$settings$chains <- 1L
    return(fit)
  }

  if (is.null(cores)) cores <- chains
  cores <- max(1L, min(as.integer(cores), chains))

  # Each chain needs its OWN workdir: the binary writes data.csv and out/ under
  # fixed names, so chains sharing a directory would overwrite each other's
  # input and read each other's output.
  base <- workdir
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  if (!keep_workdir) on.exit(unlink(base, recursive = TRUE), add = TRUE)

  seeds <- as.integer(seed) + seq_len(chains) - 1L
  chain_args <- lapply(seq_len(chains), function(k) {
    c(args, list(workdir = file.path(base, paste0("chain", k))))
  })

  # PSOCK, not fork: `parallel::mclapply` does not work on Windows, which is
  # where this package is developed and run. Workers are given the binary path
  # explicitly, because badger_binary() resolves it from an environment
  # variable a fresh worker does not inherit.
  bin <- badger_binary()
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, "bin", envir = environment())
  invisible(parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(BadgerInfeR))
    Sys.setenv(TEST2INF_BIN = bin)
    NULL
  }))

  fits <- parallel::clusterMap(cl, .run_one_chain,
                               seed = seeds, args = chain_args,
                               SIMPLIFY = FALSE, USE.NAMES = FALSE)

  settings <- fits[[1]]$settings
  settings$seed <- seed
  settings$chains <- chains
  settings$cores <- cores
  settings$seeds <- seeds
  .combine_chain_fits(fits, settings)
}

#' Fit the badger diagnostic HMM (one chain)
#'
#' Runs one of the five supported model variants through the standalone
#' `test2inf` binary and returns its outputs as data frames.
#'
#' @param test_mat A matrix or data frame with columns
#'   `time, id, captured, test1..test6`. Extra columns are ignored; the six test
#'   columns must be in [BADGER_TESTS] order. `NA` means "not tested".
#'   Rows with `captured == 0` are optional: each badger is expanded onto the
#'   full timestep grid between its first and last capture either way.
#' @param model One of [badger_models()]. Defaults to `"fixed"` -- accuracy
#'   held at the Table 1 values, which is the main comparator.
#' @param method `"nutsthenhmc"` (default), `"nuts"`, `"map"` or `"hmc"`.
#'
#'   `"nutsthenhmc"` is the default for full posterior sampling. NUTS adapts
#'   the step size and metric during warmup; fixed-length HMC then draws with
#'   them, skipping the tree building NUTS pays on every draw. The trajectory
#'   length comes from NUTS's own mean tree depth. Nothing is pre-tuned, so it
#'   is safe on any data. On the badger cohort, 1000 draws after 1000 warmup
#'   iterations takes about 2m20s per chain.
#'
#'   `"nuts"` is pure NUTS throughout: slower, and the one to compare against
#'   when a sampler looks wrong.
#'
#'   `"map"` is a point estimate by L-BFGS. It maximises the density on the
#'   **constrained** scale, matching Turing's `maximum_a_posteriori`, so it
#'   reproduces the reference fits. Note that a MAP estimate is not a posterior
#'   summary: it reports no uncertainty, and for a positive-constrained
#'   variance such as `sigma_g` the mode and the posterior mean differ
#'   materially (0.244 against 0.172 on this cohort).
#'
#'   `"hmc"` is fixed-length HMC and now **requires** an explicit `metric`.
#'   The metrics compiled into the binary were tuned for the superseded model
#'   variants -- plain `rw1`, Beta priors, no increment prior -- so they do not
#'   describe the current posterior. To use it, tune once with `write_metric =`
#'   under `"nuts"` or `"nutsthenhmc"`, then pass that file back as `metric =`.
#' @param metric Optional path to an inverse mass matrix CSV, overriding the
#'   compiled-in metric for `model`. Headerless: `n` rows of `n` values, or a
#'   single row/column of `n` values for a diagonal metric.
#' @param write_metric Optional path to write the metric actually used, so a
#'   run can be replayed or a NUTS-adapted metric reused by HMC.
#' @param draws,warmup Posterior draws and warmup iterations (sampling methods).
#' @param traj_draws Number of evenly spaced draws to reconstruct per-draw
#'   prevalence and infection times from. 0 skips that work.
#' @param target_acc NUTS target acceptance probability.
#' @param seed RNG seed.
#' @param seasons Timesteps per year.
#' @param repeat_captures How to treat multiple captures in one timestep:
#'   `"stack"`, `"pool"` or `"last"`.
#' @param year_process `"rw1c"` (the default), `"rw1"`, `"iid"`, `"rw2"`,
#'   `"none"` or `"shrunk"`. These are different models, not tuning options:
#'   each carries its own `sigma_g` prior (rw1 and rw1c 0.05, rw2 0.01, shrunk
#'   0.10, iid and none 0.30).
#'
#'   `"rw1c"` is `rw1` constrained to sum to zero. Plain `rw1` leaves the
#'   **level** of the year series unidentified -- `gamma[1]` is free and every
#'   later year is an increment from it -- so years with no test data drift,
#'   and the walk extrapolates backwards from the first year that does have
#'   data. Centring makes `gamma` a pure deviation series and hands the level
#'   to `alpha`, which its own prior and every capture in the study anchor.
#'
#'   `ar1` is not available here: it adds a `rho` parameter that changes the
#'   parameter-vector layout.
#' @param year_weight Apply the data-adaptive year-increment prior (default
#'   `TRUE`). Each year-to-year step is shrunk in proportion to how little test
#'   data informs it, so a year with no tests is pinned to its neighbour while
#'   well-tested years are left to the likelihood. Together with `"rw1c"` this
#'   is what stops sparsely tested early years inheriting the level of the
#'   better-tested ones. It does not pull hazards toward zero, and it should
#'   not: 1976-77 were tested by culture alone, whose sensitivity is 0.10, so
#'   some residual early hazard is honest rather than an artefact.
#' @param accuracy_prior `"logitnormal"` (default) or `"beta"`. Only used when
#'   the model infers Se/Sp.
#'
#'   Under `"logitnormal"` each Table 1 value is the prior **median** and the
#'   density vanishes at both probability boundaries. Under `"beta"` the shape
#'   parameters are floored at 0.5, which both distorts the requested moments
#'   -- culture specificity asks for a mean of 0.999 and gets `Beta(3.26, 0.5)`,
#'   whose mean is 0.867 -- and makes the density diverge at the boundary, so
#'   there is no interior maximum. `method = "map"` refuses a free Beta prior
#'   for that reason.
#' @param penalty Apply the Se+Sp>1 identifiability penalty.
#' @param start_year Calendar year of `time == 1`, used to label year outputs.
#'   0 means "unknown", and calendar columns are omitted.
#' @param workdir Directory for the binary's CSV input and output. Deleted on
#'   exit unless `keep_workdir = TRUE`.
#' @param keep_workdir Keep `workdir` for inspection.
#'
#' @return A list with elements `parameters`, `prevalence`, `p_infected`,
#'   `year_effects`, `draws`, `prevalence_draws`, `trajectories`, `sesp`,
#'   `settings` and `log`. Sampling-only outputs are `NULL` under `method = "map"`.
#' @keywords internal
badger_hmm_one <- function(test_mat,
                       model = "fixed",
                       method = c("nutsthenhmc", "nuts", "map", "hmc"),
                       draws = 1000L,
                       warmup = 1000L,
                       traj_draws = 0L,
                       metric = NULL,
                       write_metric = NULL,
                       target_acc = 0.8,
                       seed = 1L,
                       seasons = 4L,
                       repeat_captures = c("stack", "pool", "last"),
                       year_process = c("rw1c", "rw1", "iid", "rw2", "none",
                                        "shrunk"),
                       year_weight = TRUE,
                       accuracy_prior = c("logitnormal", "beta"),
                       penalty = TRUE,
                       start_year = 0L,
                       workdir = tempfile("badgerinfer"),
                       keep_workdir = FALSE) {
  method <- match.arg(method)
  repeat_captures <- match.arg(repeat_captures)
  year_process <- match.arg(year_process)
  accuracy_prior <- match.arg(accuracy_prior)
  spec <- resolve_model(model)

  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (!keep_workdir) on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  csv <- file.path(workdir, "data.csv")
  write_test_matrix(test_mat, csv)
  outdir <- file.path(workdir, "out")

  args <- c("--data", csv, "--out", outdir,
            "--method", method,
            "--seasons", seasons,
            "--repeat", repeat_captures,
            "--year-process", year_process,
            "--seed", seed)
  if (!isTRUE(year_weight))    args <- c(args, "--no-year-weight")
  if (isTRUE(spec$infer_sesp)) {
    args <- c(args, "--infer-sesp",
              "--accuracy-prior", accuracy_prior)
  }
  # A pinned specificity works in both accuracy modes: with Se/Sp inferred the
  # assay stops being a parameter, and with them fixed the pin overrides that
  # assay's Table 1 value.
  if (!is.null(spec$sp_pin)) {
    args <- c(args, "--pin-sp",
              sprintf("%d=%s", as.integer(spec$sp_pin[1]),
                      format(spec$sp_pin[2], digits = 15)))
  }
  if (!isTRUE(penalty))        args <- c(args, "--no-penalty")
  if (method %in% c("nuts", "hmc", "nutsthenhmc")) {
    args <- c(args, "--draws", draws, "--warmup", warmup)
    if (traj_draws > 0L) args <- c(args, "--traj-draws", traj_draws)
  }
  if (method %in% c("nuts", "nutsthenhmc")) {
    args <- c(args, "--accept", target_acc)
  }
  # `--model` selects a compiled-in pre-tuned metric, and the stored sets were
  # tuned for the superseded five-variant models. Passing a current model name
  # would either fail the lookup or, worse, match a set tuned under different
  # priors. The binary refuses the mismatched combinations; this keeps the
  # package from asking for them at all.
  if (method == "hmc" && is.null(metric)) {
    stop("method = \"hmc\" needs a `metric`: the binary's compiled-in ",
         "pre-tuned metrics were tuned for the superseded model variants ",
         "(plain rw1, Beta priors, no increment prior) and are not valid for ",
         "the current ones.\nUse method = \"nutsthenhmc\", which adapts from ",
         "the data, or tune once with write_metric = under \"nuts\" and pass ",
         "that file back as metric =.")
  }
  if (!is.null(metric)) {
    if (!file.exists(metric)) stop("metric file not found: ", metric)
    args <- c(args, "--metric", normalizePath(metric))
  }
  if (!is.null(write_metric)) {
    args <- c(args, "--write-metric", write_metric)
  }

  # No fallback path here any more. The package never passes `--model`, because
  # the compiled-in pre-tuned metrics belong to the superseded variants, and
  # `method = "hmc"` without an explicit `metric` is refused above. The default
  # method adapts to whatever data it is given, so there is nothing left to
  # fall back FROM.
  log <- run_test2inf(args)

  rd <- function(f) {
    p <- file.path(outdir, f)
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
  }
  params <- rd("parameters.csv")

  out <- list(
    parameters       = params,
    prevalence       = add_calendar_year(rd("prevalence.csv"), start_year, seasons),
    p_infected       = rd("p_infected.csv"),
    year_effects     = add_year_calendar(rd("year_effects.csv"), start_year),
    draws            = rd("draws.csv"),
    prevalence_draws = rd("prevalence_draws.csv"),
    trajectories     = rd("trajectories.csv"),
    sesp             = extract_sesp(params, spec$infer_sesp, spec$sp_pin),
    settings = list(model = model, label = badger_label(model), method = method,
                    infer_sesp = spec$infer_sesp, sp_pin = spec$sp_pin,
                    accuracy_prior = accuracy_prior, year_weight = year_weight,
                    draws = draws, warmup = warmup, traj_draws = traj_draws,
                    metric = metric, write_metric = write_metric,
                    seed = seed, seasons = seasons,
                    repeat_captures = repeat_captures,
                    year_process = year_process, penalty = penalty,
                    start_year = start_year),
    log = log
  )
  class(out) <- "badger_fit"
  out
}

#' @export
print.badger_fit <- function(x, ...) {
  s <- x$settings
  cat("<badger_fit>", s$label, "/", s$method, "\n")
  cat("  year process: ", s$year_process, "   seed: ", s$seed, "\n", sep = "")
  if (!is.null(x$parameters)) {
    cat("  parameters:   ", nrow(x$parameters), "\n", sep = "")
  }
  if (!is.null(x$draws)) {
    cat("  draws:        ", nrow(x$draws), "\n", sep = "")
  }
  if (!is.null(x$trajectories)) {
    cat("  trajectories: ", nrow(x$trajectories), " rows\n", sep = "")
  }
  # The binary reports divergences and step size on stdout; surface them,
  # because a run with post-warmup divergences should not pass unnoticed.
  diag <- grep("^(step size|warmup|sampling|WARNING)", x$log, value = TRUE)
  for (l in diag) cat("  ", l, "\n", sep = "")
  invisible(x)
}
