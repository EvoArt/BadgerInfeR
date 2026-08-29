# The one function scripts call. Everything else in this package supports it.

#' Fit the badger diagnostic HMM
#'
#' Runs one of the five supported model variants through the standalone
#' `test2inf` binary and returns its outputs as data frames.
#'
#' @param test_mat A matrix or data frame with columns
#'   `time, id, captured, test1..test6`. Extra columns are ignored; the six test
#'   columns must be in [BADGER_TESTS] order. `NA` means "not tested".
#'   Rows with `captured == 0` are optional: each badger is expanded onto the
#'   full timestep grid between its first and last capture either way.
#' @param model One of [badger_models()].
#' @param method `"hmc"` (default), `"nuts"` or `"map"`.
#'
#'   `"hmc"` runs fixed-length HMC with the step size, trajectory length and
#'   **dense** mass matrix tuned offline for this exact model variant and
#'   compiled into the binary. On the badger cohort that samples ~8x more
#'   effective draws per second than adapting a diagonal metric from scratch,
#'   because the dense metric removes the funnel that rw1's `cumsum` creates:
#'   it steps at eps ~0.3 rather than ~0.047. Adaptation is not repeated.
#'
#'   The tuned values are specific to the cohort they were tuned on. Only the
#'   parameter count is checked, so pointing them at a very different dataset
#'   samples badly rather than failing. For other data either use `"nuts"`,
#'   which adapts from scratch and is always safe, or tune once and reuse:
#'   run with `write_metric =` under `"nuts"`, then pass that file as `metric =`
#'   to subsequent `"hmc"` runs. No rebuild and no Julia required.
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
#' @param year_process `"rw1"` (the default, and what every sql-e2e fit used) or
#'   `"iid"`. These are different models, not tuning options.
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
#' @export
badger_hmm <- function(test_mat,
                       model = "all_fixed",
                       method = c("hmc", "nuts", "map"),
                       draws = 1000L,
                       warmup = 1000L,
                       traj_draws = 0L,
                       metric = NULL,
                       write_metric = NULL,
                       target_acc = 0.8,
                       seed = 1L,
                       seasons = 4L,
                       repeat_captures = c("stack", "pool", "last"),
                       year_process = c("rw1", "iid"),
                       penalty = TRUE,
                       start_year = 0L,
                       workdir = tempfile("badgerinfer"),
                       keep_workdir = FALSE) {
  method <- match.arg(method)
  repeat_captures <- match.arg(repeat_captures)
  year_process <- match.arg(year_process)
  spec <- resolve_model(model)

  if (identical(model, "culture_only_fixed") && isTRUE(spec$infer_sesp)) {
    stop("culture-only with inferred Se/Sp is not identifiable.")  # defensive
  }

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
  if (!is.null(spec$tests)) {
    args <- c(args, "--tests", paste(spec$tests, collapse = ","))
  }
  if (isTRUE(spec$infer_sesp)) args <- c(args, "--infer-sesp")
  if (!isTRUE(penalty))        args <- c(args, "--no-penalty")
  if (method %in% c("nuts", "hmc")) {
    args <- c(args, "--draws", draws, "--warmup", warmup)
    if (traj_draws > 0L) args <- c(args, "--traj-draws", traj_draws)
  }
  if (method == "nuts") args <- c(args, "--accept", target_acc)
  if (method == "hmc")  args <- c(args, "--model", model)
  if (!is.null(metric)) {
    if (!file.exists(metric)) stop("metric file not found: ", metric)
    args <- c(args, "--metric", normalizePath(metric))
  }
  if (!is.null(write_metric)) {
    args <- c(args, "--write-metric", write_metric)
  }

  # The compiled-in metrics are sized for the badger cohort (58 or 70
  # parameters). On data with a different number of years the binary refuses
  # them, which would make the DEFAULT method fail on perfectly reasonable
  # input. Fall back to NUTS with a warning rather than erroring: adapting from
  # scratch is always correct, just slower.
  fell_back <- FALSE
  log <- tryCatch(run_test2inf(args), error = function(e) {
    if (method == "hmc" && is.null(metric) &&
        grepl("was tuned for", conditionMessage(e), fixed = TRUE)) {
      warning("The pre-tuned metric for '", model, "' does not fit this data ",
              "(it is tuned for the badger cohort). Falling back to method = ",
              "\"nuts\", which adapts from scratch. To use HMC here, tune once ",
              "with write_metric = and pass it back as metric =.",
              call. = FALSE)
      fell_back <<- TRUE
      nuts_args <- args[!args %in% c("--model", model)]
      nuts_args[nuts_args == "hmc"] <- "nuts"
      return(run_test2inf(c(nuts_args, "--accept", target_acc)))
    }
    stop(e)
  })
  # `fell_back` is assigned with <<- from inside the handler, so it must exist
  # in this frame first -- otherwise <<- would walk past it to the global env.
  if (fell_back) method <- "nuts"

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
    sesp             = extract_sesp(params, spec$infer_sesp),
    settings = list(model = model, label = badger_label(model), method = method,
                    tests = spec$tests, infer_sesp = spec$infer_sesp,
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
