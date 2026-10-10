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
                       method = c("nuts", "nutsthenhmc", "map", "hmc"),
                       draws = 1000L,
                       warmup = 1000L,
                       traj_draws = 0L,
                       traj_mode = c("threshold", "sample"),
                       states = 2L,
                       e_mult = NULL,
                       birth = NULL,
                       init_starts = 1L,
                       init_jitter = 0,
                       metric = NULL,
                       write_metric = NULL,
                       step_size = NULL,
                       n_leapfrog = NULL,
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
  traj_mode <- match.arg(traj_mode)
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
               traj_mode = traj_mode, states = states, e_mult = e_mult,
               birth = birth, init_starts = init_starts,
               init_jitter = init_jitter,
               metric = metric, write_metric = write_metric,
               step_size = step_size, n_leapfrog = n_leapfrog,
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
#' Runs one of the supported model variants through the standalone `test2inf`
#' binary and returns its outputs as data frames.
#'
#' @param test_mat A matrix or data frame with columns
#'   `time, id, captured, test1..test6`. Extra columns are ignored; the six test
#'   columns must be in [BADGER_TESTS] order. `NA` means "not tested".
#'   Rows with `captured == 0` are optional: each badger is expanded onto the
#'   full timestep grid between its first and last capture either way.
#' @param model One of [badger_models()].
#' @param method `"nuts"` (default), `"nutsthenhmc"`, `"map"` or `"hmc"`.
#'
#'   `"nutsthenhmc"` is quicker per iteration but mixes far worse here (min ESS
#'   46 against 843); check ESS if you use it. `"map"` maximises the
#'   constrained-scale density, as Turing does, and reports no uncertainty.
#'   `"hmc"` adapts nothing, so it needs `metric`, `step_size` and
#'   `n_leapfrog`: the compiled-in values were tuned for the superseded
#'   variants.
#' @param metric Optional path to an inverse mass matrix CSV. Headerless: `n`
#'   rows of `n` values, or a single row/column of `n` for a diagonal metric.
#' @param write_metric Optional path to write the metric actually used.
#' @param step_size,n_leapfrog Step size and trajectory length for
#'   `method = "hmc"`, which adapts neither. Under `"nutsthenhmc"`,
#'   `n_leapfrog` overrides the length inferred from NUTS's tree depth.
#' @param draws,warmup Posterior draws and warmup iterations (sampling methods).
#' @param traj_draws Number of evenly spaced draws to reconstruct per-draw
#'   prevalence and infection times from. 0 skips that work.
#' @param traj_mode How the infection time in `trajectories` is found for each
#'   kept draw. `"threshold"` (default): first timestep where smoothed
#'   P(infected) reaches 0.5. `"sample"`: sample a state path given that draw's
#'   parameters, so a badger with a 20% chance of infection is infected in
#'   about 20% of draws.
#' @param states `2` (default) for uninfected -> infected, or `3` for
#'   uninfected -> exposed -> infectious. The three-state model adds `tau`
#'   (per-timestep exposed -> infectious hazard), `rho` (share of badgers
#'   infected at first capture that are already infectious) and a multiplier
#'   per assay with P(+ | exposed) = `e_mult[k]` * `Se[k]`. `p_infected` is then
#'   P(exposed or infectious), with `p_exposed` and `p_infectious` alongside,
#'   and `trajectories` has `infectious_time`.
#' @param e_mult Multipliers to fix, for `states = 3`: length six in
#'   [BADGER_TESTS] order, values in `[0, 1]`, `NA` to estimate. `NULL`
#'   (default) estimates all six with a Uniform(0, 1) prior.
#' @param birth Birth timesteps as a numeric vector named by badger `id`, on
#'   the same grid as `test_mat`. These badgers start susceptible at birth;
#'   badgers that are missing or `NA` keep the estimated probability of being
#'   infected at first capture. `NULL` (default) uses that for everyone.
#'   Births before timestep 1 or after the first capture are ignored. Rows in
#'   `p_infected` then start at birth.
#' @param init_starts Number of starting points for the optimiser that
#'   initialises each chain; chains start from the best optimum found. The
#'   starts are the same for every chain. All optima are printed in `log`.
#'   Useful when chains can get stuck in a minor mode (the three-state model on
#'   the Woodchester data has one).
#' @param init_jitter Sd of noise added to each chain's starting point
#'   (unconstrained scale, seeded by the chain's seed). Chains then start
#'   apart, so separated modes show up as a large `rhat`.
#' @param target_acc NUTS target acceptance probability.
#' @param seed RNG seed.
#' @param seasons Timesteps per year.
#' @param repeat_captures How to treat multiple captures in one timestep:
#'   `"stack"`, `"pool"` or `"last"`.
#' @param year_process `"rw1c"` (default), `"rw1"`, `"iid"`, `"rw2"`, `"none"`
#'   or `"shrunk"`. Different models, not tuning options: each carries its own
#'   `sigma_g` prior. `"rw1c"` is `rw1` constrained to sum to zero, which
#'   identifies the level of the year series that plain `rw1` leaves free.
#'   `ar1` is unavailable (it adds a parameter).
#' @param year_weight Data-adaptive year-increment prior (default `TRUE`).
#'   Shrinks each year-to-year step by how little test data informs it, so an
#'   untested year is pinned to its neighbour.
#' @param accuracy_prior `"logitnormal"` (default) or `"beta"`, used only when
#'   the model infers Se/Sp. Beta's shape floor distorts the requested moments
#'   and diverges at the boundary, so `method = "map"` refuses a free one.
#' @param penalty Apply the Se+Sp>1 identifiability penalty.
#' @param start_year Calendar year of `time == 1`, used to label year outputs.
#'   0 means "unknown", and calendar columns are omitted.
#' @param workdir Directory for the binary's CSV input and output. Deleted on
#'   exit unless `keep_workdir = TRUE`.
#' @param keep_workdir Keep `workdir` for inspection.
#'
#' @return A list with elements `parameters`, `prevalence`, `p_infected`,
#'   `year_effects`, `draws`, `prevalence_draws`, `trajectories`, `sesp`,
#'   `settings` and `log` (and `prevalence_infectious_draws` if `states = 3`). Sampling-only outputs are `NULL` under `method = "map"`.
#' @keywords internal
badger_hmm_one <- function(test_mat,
                       model = "fixed",
                       method = c("nuts", "nutsthenhmc", "map", "hmc"),
                       draws = 1000L,
                       warmup = 1000L,
                       traj_draws = 0L,
                       traj_mode = c("threshold", "sample"),
                       states = 2L,
                       e_mult = NULL,
                       birth = NULL,
                       init_starts = 1L,
                       init_jitter = 0,
                       metric = NULL,
                       write_metric = NULL,
                       step_size = NULL,
                       n_leapfrog = NULL,
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
  traj_mode <- match.arg(traj_mode)
  states <- as.integer(states)
  if (length(states) != 1L || is.na(states) || !states %in% 2:3) {
    stop("`states` must be 2 or 3.")
  }
  if (!is.null(e_mult)) {
    if (!is.numeric(e_mult) || length(e_mult) != length(BADGER_TESTS) ||
        any(!is.na(e_mult) & (e_mult < 0 | e_mult > 1))) {
      stop("`e_mult` must be NULL or ", length(BADGER_TESTS),
           " values in [0, 1], with NA for a multiplier to estimate.")
    }
  }
  if (!is.null(birth)) {
    bid <- suppressWarnings(as.integer(names(birth)))
    if (!is.numeric(birth) || is.null(names(birth)) || anyNA(bid) ||
        anyDuplicated(bid)) {
      stop("`birth` must be a numeric vector of birth timesteps named by ",
           "badger id, one entry per badger.")
    }
    known <- !is.na(birth)
    if (any(birth[known] != round(birth[known]))) {
      stop("`birth` must hold whole timesteps on the same grid as `test_mat`.")
    }
  }
  if (!is.numeric(init_jitter) || length(init_jitter) != 1L ||
      is.na(init_jitter) || init_jitter < 0) {
    stop("`init_jitter` must be a single non-negative number.")
  }
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
  if (states == 3L) {
    args <- c(args, "--states", 3L)
    for (k in which(!is.na(e_mult))) {
      args <- c(args, "--e-mult",
                sprintf("%d=%s", k, format(e_mult[k], digits = 15)))
    }
  }
  if (init_jitter > 0) args <- c(args, "--init-jitter", init_jitter)
  init_starts <- as.integer(init_starts)
  if (length(init_starts) != 1L || is.na(init_starts) || init_starts < 1L) {
    stop("`init_starts` must be a positive integer.")
  }
  if (init_starts > 1L) args <- c(args, "--init-starts", init_starts)
  if (!is.null(birth) && any(!is.na(birth))) {
    bcsv <- file.path(workdir, "birth.csv")
    keep <- !is.na(birth)
    utils::write.csv(data.frame(id = as.integer(names(birth))[keep],
                                birth = as.integer(birth[keep])),
                     bcsv, row.names = FALSE, quote = FALSE)
    args <- c(args, "--birth", bcsv)
  }
  if (method %in% c("nuts", "hmc", "nutsthenhmc")) {
    args <- c(args, "--draws", draws, "--warmup", warmup)
    if (traj_draws > 0L) {
      args <- c(args, "--traj-draws", traj_draws, "--traj-mode", traj_mode)
    }
  }
  if (method %in% c("nuts", "nutsthenhmc")) {
    args <- c(args, "--accept", target_acc)
  }
  # `--model` selects a compiled-in pre-tuned metric, and the stored sets were
  # tuned for the superseded five-variant models. Passing a current model name
  # would either fail the lookup or, worse, match a set tuned under different
  # priors. The binary refuses the mismatched combinations; this keeps the
  # package from asking for them at all.
  if (method == "hmc" &&
      (is.null(metric) || is.null(step_size) || is.null(n_leapfrog))) {
    stop("method = \"hmc\" adapts nothing, so it needs `metric`, `step_size` ",
         "and `n_leapfrog` together. The binary's compiled-in pre-tuned ",
         "values belong to the superseded variants (plain rw1, Beta priors, ",
         "no increment prior) and are not valid here.\n",
         "Use method = \"nuts\", or adapt once with write_metric = and reuse ",
         "the exported metric with the step size that run reports.")
  }
  if (!is.null(metric)) {
    if (!file.exists(metric)) stop("metric file not found: ", metric)
    args <- c(args, "--metric", normalizePath(metric))
  }
  if (!is.null(write_metric)) {
    args <- c(args, "--write-metric", write_metric)
  }
  if (!is.null(step_size))  args <- c(args, "--hmc-eps", step_size)
  if (!is.null(n_leapfrog)) args <- c(args, "--hmc-L", n_leapfrog)

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
    prevalence_infectious_draws = rd("prevalence_infectious_draws.csv"),
    trajectories     = rd("trajectories.csv"),
    sesp             = extract_sesp(params, spec$infer_sesp, spec$sp_pin),
    settings = list(model = model, label = badger_label(model), method = method,
                    infer_sesp = spec$infer_sesp, sp_pin = spec$sp_pin,
                    accuracy_prior = accuracy_prior, year_weight = year_weight,
                    draws = draws, warmup = warmup, traj_draws = traj_draws,
                    traj_mode = traj_mode, states = states,
                    e_mult = if (states == 3L) e_mult else NULL,
                    init_jitter = init_jitter, init_starts = init_starts,
                    n_birth = if (is.null(birth)) 0L else sum(!is.na(birth)),
                    metric = metric, write_metric = write_metric,
                    step_size = step_size, n_leapfrog = n_leapfrog,
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
  if (identical(s$states, 3L)) {
    cat("  states:       3 (uninfected, exposed, infectious)\n")
  }
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
  diag <- grep("^(births|init:|step size|warmup|sampling|WARNING)", x$log, value = TRUE)
  for (l in diag) cat("  ", l, "\n", sep = "")
  invisible(x)
}
