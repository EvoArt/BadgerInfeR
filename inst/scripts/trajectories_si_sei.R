# Posterior infection trajectories for the SI and SEI models, written in the
# layout of the `sql e2e` results folder.
#
#   Rscript trajectories_si_sei.R <input.csv> <key.csv> <outdir> [chains] [cores] [birth.csv]
#
#   input.csv  time, id, captured, test1..test6
#   key.csv    id, tattoo
#   birth.csv  id, birth   (optional; birth timestep on the model grid, NA = unknown)
#
# Each fit is `chains` NUTS chains of 1000 draws after 1000 warmup. Trajectories
# are sampled state paths, 250 per chain. If birth.csv is given, badgers with a
# birth quarter start susceptible there in all fits.

suppressPackageStartupMessages(library(BadgerInfeR))

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 3L) {
  stop("usage: trajectories_si_sei.R <input.csv> <key.csv> <outdir> ",
       "[chains] [cores]")
}
input_csv <- argv[1]
key_csv   <- argv[2]
outdir    <- argv[3]
CHAINS    <- if (length(argv) >= 4L) as.integer(argv[4]) else 4L
CORES     <- if (length(argv) >= 5L) as.integer(argv[5]) else CHAINS
birth_csv <- if (length(argv) >= 6L) argv[6] else ""

START_YEAR <- 1976L      # time == 1 is the first quarter of 1976
SEASONS    <- 4L
DRAWS      <- 1000L
WARMUP     <- 1000L
TRAJ_DRAWS <- 250L       # per chain
SEED       <- 1L
# Chain starts. The SI fits use jittered starts. The SEI posterior has a minor
# second mode that chains get stuck in (about 1e-8 of the mass by a Laplace
# estimate), so its chains start from the best of several optimiser runs. The
# optima found are written to init_optima.csv.
INIT_JITTER <- 1
INIT_STARTS <- 8L

# label -> model settings
FITS <- list(
  "SI (fixed accuracy)"     = list(model = "fixed",    states = 2L,
                                   jitter = INIT_JITTER, starts = 1L),
  "SI (inferred accuracy)"  = list(model = "inferred", states = 2L,
                                   jitter = INIT_JITTER, starts = 1L),
  "SEI (inferred accuracy)" = list(model = "inferred", states = 3L,
                                   jitter = 0, starts = INIT_STARTS)
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
d <- utils::read.csv(input_csv)
test_matrix <- as.matrix(d[, c("time", "id", "captured", paste0("test", 1:6))])
key <- unique(utils::read.csv(key_csv, stringsAsFactors = FALSE)[, c("id", "tattoo")])
stopifnot(!anyDuplicated(key$id), all(d$id %in% key$id))

birth <- NULL
if (nzchar(birth_csv)) {
  b <- utils::read.csv(birth_csv)
  birth <- stats::setNames(b$birth, b$id)
  message(sum(!is.na(birth)), " badgers with a birth quarter, ",
          sum(is.na(birth)), " without")
}

year_of_time <- function(time) START_YEAR + (time - 1L) %/% SEASONS

# Fits are cached in <outdir>/fits and reused. Delete that folder to refit.
fitdir <- file.path(outdir, "fits")
dir.create(fitdir, showWarnings = FALSE)
fits <- list()
for (lab in names(FITS)) {
  spec <- FITS[[lab]]
  rds <- file.path(fitdir, paste0(gsub("[^a-z]+", "_", tolower(lab)), ".rds"))
  if (file.exists(rds)) {
    message("reusing: ", lab)
    fits[[lab]] <- readRDS(rds)
    next
  }
  message("fitting: ", lab)
  t0 <- Sys.time()
  fits[[lab]] <- badger_hmm(test_matrix, model = spec$model, method = "nuts",
                            states = spec$states, draws = DRAWS, warmup = WARMUP,
                            traj_draws = TRAJ_DRAWS, traj_mode = "sample",
                            chains = CHAINS, cores = CORES, seed = SEED,
                            init_jitter = spec$jitter, init_starts = spec$starts,
                            birth = birth,
                            seasons = SEASONS, start_year = START_YEAR)
  fits[[lab]]$wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  %.0f s", fits[[lab]]$wall))
  saveRDS(fits[[lab]], rds)
}

# ---- drop chains that ended up in another mode --------------------------------
# If the chains disagree, keep the largest subset that agrees and record which
# were dropped in diagnostics.csv. Stops if it is a tie. Check separately that
# the kept mode is the one with the posterior mass.
RHAT_OK <- 1.05
max_rhat <- function(f) max(f$diagnostics$rhat, na.rm = TRUE)
chain_note <- list()
for (lab in names(fits)) {
  f <- fits[[lab]]
  chain_note[[lab]] <- list(kept = seq_len(CHAINS), dropped = integer(0),
                            rhat_all = max_rhat(f))
  if (max_rhat(f) <= RHAT_OK) next
  best <- NULL
  for (size in (CHAINS - 1L):2L) {
    sets <- utils::combn(CHAINS, size, simplify = FALSE)
    ok <- Filter(function(k) max_rhat(badger_keep_chains(f, k)) <= RHAT_OK, sets)
    if (length(ok) == 1L) { best <- ok[[1]]; break }
    if (length(ok) > 1L) {
      stop(lab, ": chains disagree (max Rhat ", round(max_rhat(f), 2),
           ") and more than one set of ", size, " agrees among themselves. ",
           "Decide which mode to report before writing results.")
    }
  }
  if (is.null(best)) {
    stop(lab, ": chains disagree (max Rhat ", round(max_rhat(f), 2),
         ") and no subset of them agrees.")
  }
  message(lab, ": max Rhat ", round(max_rhat(f), 2), " over all chains; keeping ",
          paste(best, collapse = ","), ", which agree (max Rhat ",
          round(max_rhat(badger_keep_chains(f, best)), 4), ")")
  chain_note[[lab]]$kept <- best
  chain_note[[lab]]$dropped <- setdiff(seq_len(CHAINS), best)
  fits[[lab]] <- badger_keep_chains(f, best)
}

# ---- optima from the multi-start initialisation -------------------------------
optima <- do.call(rbind, lapply(names(fits), function(lab) {
  logs <- fits[[lab]]$log
  l <- if (is.list(logs)) logs[[1]] else logs      # the same for every chain
  ln <- grep("^init start ", l, value = TRUE)
  if (length(ln) == 0L) return(NULL)
  lp <- as.numeric(sub("^init start [0-9]+: log posterior +", "", ln))
  data.frame(model = lab, start = seq_along(lp), log_posterior = lp,
             below_best = max(lp) - lp)
}))
if (!is.null(optima)) {
  utils::write.csv(optima, file.path(outdir, "init_optima.csv"), row.names = FALSE)
  print(optima)
}

# ---- diagnostics ------------------------------------------------------------
# Divergences come from each chain's log; Rhat and ESS from the pooled draws.
divergences <- function(fit) {
  logs <- if (is.list(fit$log)) fit$log else list(fit$log)
  sum(vapply(logs, function(l) {
    ln <- grep("^sampling: +divergences:", l, value = TRUE)
    if (length(ln) == 0L) return(NA_real_)
    as.numeric(sub("^sampling: +divergences: +([0-9]+).*$", "\\1", ln[1]))
  }, numeric(1)))
}
diag <- do.call(rbind, lapply(names(fits), function(lab) {
  f <- fits[[lab]]
  dg <- f$diagnostics
  data.frame(model = lab, chains = CHAINS, draws_per_chain = DRAWS,
             chains_kept = paste(chain_note[[lab]]$kept, collapse = " "),
             chains_dropped = paste(chain_note[[lab]]$dropped, collapse = " "),
             max_rhat_all_chains = chain_note[[lab]]$rhat_all,
             parameters = nrow(dg),
             max_rhat = max(dg$rhat, na.rm = TRUE),
             min_ess_bulk = min(dg$ess_bulk, na.rm = TRUE),
             median_ess_bulk = stats::median(dg$ess_bulk, na.rm = TRUE),
             worst_parameter = dg$parameter[which.min(dg$ess_bulk)],
             divergences = divergences(f),
             wall_seconds = round(f$wall),
             stringsAsFactors = FALSE)
}))
utils::write.csv(diag, file.path(outdir, "diagnostics.csv"), row.names = FALSE)
print(diag)

utils::write.csv(
  do.call(rbind, lapply(names(fits), function(lab) {
    cbind(model = lab, fits[[lab]]$diagnostics)
  })),
  file.path(outdir, "diagnostics_by_parameter.csv"), row.names = FALSE)

# ---- per-capture probabilities, wide by model -------------------------------
per_capture <- merge(d[, c("id", "time")], key, by = "id")
per_capture <- unique(per_capture)
per_capture$year <- year_of_time(per_capture$time)
per_capture$quarter <- (per_capture$time - 1L) %% SEASONS + 1L
for (lab in names(fits)) {
  p <- fits[[lab]]$p_infected
  i <- match(paste(per_capture$id, per_capture$time), paste(p$id, p$time))
  per_capture[[lab]] <- p$p_infected[i]
  if (FITS[[lab]]$states == 3L) {
    per_capture[[paste(lab, "exposed")]]    <- p$p_exposed[i]
    per_capture[[paste(lab, "infectious")]] <- p$p_infectious[i]
  }
}
per_capture <- per_capture[order(per_capture$tattoo, per_capture$time), ]
per_capture$id <- NULL
utils::write.csv(per_capture, file.path(outdir, "p_infected_per_capture.csv"),
                 row.names = FALSE)

# ---- sampled trajectories ---------------------------------------------------
# One row per (model, badger, draw): infection_time and, for SEI,
# infectious_time. 0 = not within the record. `draw` runs across chains, so
# rows with the same draw (within a model) share one set of parameters.
traj <- do.call(rbind, lapply(names(fits), function(lab) {
  tj <- fits[[lab]]$trajectories
  data.frame(model = lab,
             tattoo = key$tattoo[match(tj$id, key$id)],
             draw = (tj$chain - 1L) * TRAJ_DRAWS + tj$draw,
             infection_time = as.integer(tj$infection_time),
             infectious_time = if ("infectious_time" %in% names(tj))
               as.integer(tj$infectious_time) else NA_integer_,
             stringsAsFactors = FALSE)
}))
utils::write.csv(traj, file.path(outdir, "trajectories_by_draw.csv"),
                 row.names = FALSE)

# Compact form: how many draws put the event at each timestep.
compact <- function(col, yearcol) {
  x <- traj[!is.na(traj[[col]]), ]
  out <- stats::aggregate(list(n_draws = rep(1L, nrow(x))),
                          by = list(model = x$model, tattoo = x$tattoo,
                                    time = x[[col]]), FUN = sum)
  out[[yearcol]] <- ifelse(out$time == 0L, NA_integer_, year_of_time(out$time))
  names(out)[names(out) == "time"] <- col
  out[order(out$model, out$tattoo, out[[col]]), ]
}
utils::write.csv(compact("infection_time", "infection_year"),
                 file.path(outdir, "trajectories.csv"), row.names = FALSE)
utils::write.csv(compact("infectious_time", "infectious_year"),
                 file.path(outdir, "trajectories_infectious.csv"),
                 row.names = FALSE)

# ---- annual prevalence, with intervals from the draws -----------------------
annual_draws <- function(pd) {
  # Average within a year for every draw; one row per (year, draw).
  dcols <- grep("^draw", names(pd), value = TRUE)
  yr <- year_of_time(pd$time)
  ok <- !is.na(pd[[dcols[1]]])
  parts <- lapply(split(pd[ok, ], pd$chain[ok]), function(ch) {
    y <- year_of_time(ch$time)
    m <- rowsum(as.matrix(ch[, dcols]), y) / as.vector(table(y))
    data.frame(calendar_year = rep(as.integer(rownames(m)), ncol(m)),
               draw = (ch$chain[1] - 1L) * TRAJ_DRAWS +
                 rep(seq_len(ncol(m)), each = nrow(m)),
               prevalence = as.vector(m))
  })
  do.call(rbind, parts)
}
prev_traj <- do.call(rbind, lapply(names(fits), function(lab) {
  f <- fits[[lab]]
  out <- cbind(model = lab, state = "infected", annual_draws(f$prevalence_draws))
  if (!is.null(f$prevalence_infectious_draws)) {
    out <- rbind(out, cbind(model = lab, state = "infectious",
                            annual_draws(f$prevalence_infectious_draws)))
  }
  out
}))
utils::write.csv(prev_traj, file.path(outdir, "prevalence_trajectories.csv"),
                 row.names = FALSE)

q <- function(p) function(x) stats::quantile(x, p, names = FALSE)
grp <- list(model = prev_traj$model, state = prev_traj$state,
            calendar_year = prev_traj$calendar_year)
by_year <- stats::aggregate(list(prevalence = prev_traj$prevalence), grp, mean)
by_year$lower <- stats::aggregate(prev_traj$prevalence, grp, q(0.025))$x
by_year$upper <- stats::aggregate(prev_traj$prevalence, grp, q(0.975))$x
by_year <- by_year[order(by_year$model, by_year$state, by_year$calendar_year), ]
utils::write.csv(by_year, file.path(outdir, "prevalence_by_year.csv"),
                 row.names = FALSE)

# ---- parameters -------------------------------------------------------------
utils::write.csv(
  do.call(rbind, lapply(names(fits), function(lab) {
    cbind(model = lab, fits[[lab]]$sesp)
  })),
  file.path(outdir, "sensitivity_specificity.csv"), row.names = FALSE)

message("wrote ", normalizePath(outdir))
