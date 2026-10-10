# These tests need the test2inf binary. Set TEST2INF_BIN, or install the
# bundle into inst/bin/. Without it every test skips rather than fails, so the
# suite is still runnable on a machine that has not built the CLI.
have_binary <- function() {
  ok <- tryCatch({ check_badger_binary(); TRUE }, error = function(e) FALSE)
  if (!ok) skip("test2inf binary not available (set TEST2INF_BIN)")
  invisible(TRUE)
}

test_data <- function() {
  p <- Sys.getenv("BADGER_TEST_CSV", "")
  if (!nzchar(p) || !file.exists(p)) skip("set BADGER_TEST_CSV to a sim.csv")
  as.matrix(utils::read.csv(p))
}

test_that("the model list is exactly the four current variants", {
  expect_setequal(badger_models(),
                  c("fixed", "fixed_cult1", "inferred", "inferred_cult1"))
  expect_equal(badger_label("fixed_cult1"), "Fixed accuracy, culture Sp = 1")
  expect_error(badger_hmm(matrix(0, 1, 9), model = "nonsense"), "Unknown model")
  # The superseded names must fail with an explanation, not a bare "unknown":
  # a script still asking for one wants archived output, not a fresh fit under
  # a different model that happens to share the name.
  expect_error(badger_hmm(matrix(0, 1, 9), model = "all_fixed"), "superseded")
})

test_that("the defaults are the current model and sampler", {
  fml <- formals(badger_hmm)
  expect_equal(fml$model, "fixed")
  expect_equal(eval(fml$method)[1], "nuts")
  expect_equal(eval(fml$year_process)[1], "rw1c")
  expect_true(fml$year_weight)
  expect_equal(eval(fml$accuracy_prior)[1], "logitnormal")
})

test_that("hmc refuses to run on stale compiled-in constants", {
  # The pre-tuned values belong to the superseded variants. Sampling with one
  # would not fail, it would sample badly and quietly, so the package declines
  # to ask for it -- and hmc adapts nothing, so a metric alone is not enough.
  # A VALID matrix: the refusal has to come from the method check, not from
  # input validation, or the test would pass for the wrong reason.
  m <- matrix(1, 2, 9); m[, 1] <- c(1, 2)
  expect_error(badger_hmm(m, "fixed", "hmc"), "adapts nothing")
  expect_error(badger_hmm(m, "fixed", "hmc", metric = "m.csv"),
               "adapts nothing")
})

test_that("the test matrix is validated before the binary is invoked", {
  wtm <- BadgerInfeR:::write_test_matrix
  expect_error(wtm(matrix(0, 2, 3), tempfile()), "at least 9")
  bad <- matrix(1, 2, 9); bad[1, 3] <- 7
  expect_error(wtm(bad, tempfile()), "captured")
  bad2 <- matrix(1, 2, 9); bad2[1, 1] <- 0
  expect_error(wtm(bad2, tempfile()), "positive integer")
})

test_that("NA test results survive the round trip as empty fields", {
  m <- matrix(1, 2, 9)
  m[, 1] <- c(1, 2); m[, 2] <- 1; m[, 3] <- 1
  m[1, 4] <- NA
  p <- tempfile(fileext = ".csv")
  BadgerInfeR:::write_test_matrix(m, p)
  txt <- readLines(p)
  expect_true(grepl(",,", txt[2]) || grepl(",$", txt[2]))
})

test_that("MAP runs and returns the documented fields", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "fixed", "map", start_year = 1976)
  expect_s3_class(fit, "badger_fit")
  expect_true(all(c("parameters", "prevalence", "p_infected",
                    "year_effects", "sesp", "settings") %in% names(fit)))
  expect_true(nrow(fit$parameters) > 0)
  # MAP has no draws, and must say so rather than inventing them.
  expect_null(fit$draws)
  expect_null(fit$prevalence_draws)
  expect_true("calendar_year" %in% names(fit$prevalence))
})

test_that("Se/Sp are flagged estimated only when the model infers them", {
  have_binary(); m <- test_data()
  fixed <- badger_hmm(m, "fixed", "map")
  inferred <- badger_hmm(m, "inferred", "map")
  expect_false(any(fixed$sesp$estimated))
  expect_true(all(inferred$sesp$estimated))
  # A fixed fit must report exactly the Table 1 constants.
  expect_equal(fixed$sesp$Se, BADGER_SE, tolerance = 1e-8)
  expect_equal(fixed$sesp$Sp, BADGER_SP, tolerance = 1e-8)
  # An inferred fit should have moved off them.
  expect_false(isTRUE(all.equal(inferred$sesp$Se, BADGER_SE, tolerance = 1e-6)))
})

test_that("a pinned culture specificity is reported as a fixed 1", {
  have_binary(); m <- test_data()
  k <- BADGER_CULTURE

  # Inferred: culture stops being a parameter, so it has no row in the output
  # and must be filled from the pin rather than read positionally -- reading by
  # position would slide every later assay onto its neighbour's value.
  pin <- badger_hmm(m, "inferred_cult1", "map")
  expect_equal(nrow(pin$sesp), length(BADGER_TESTS))
  expect_equal(pin$sesp$Sp[k], 1)
  expect_false(pin$sesp$estimated[k])
  expect_true(all(pin$sesp$estimated[-k]))
  expect_false(any(is.na(pin$sesp$Sp)))

  # Fixed: the pin OVERRIDES the Table 1 value rather than removing a
  # parameter, so it is still a different model from plain "fixed".
  fx <- badger_hmm(m, "fixed_cult1", "map")
  expect_equal(fx$sesp$Sp[k], 1)
  expect_equal(fx$sesp$Sp[-k], BADGER_SP[-k], tolerance = 1e-8)
  expect_false(any(fx$sesp$estimated))
})

test_that("NUTS produces draws, and traj_draws produces trajectory outputs", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "fixed", "nuts", draws = 100, warmup = 100,
                    traj_draws = 5, seed = 2)
  expect_equal(nrow(fit$draws), 100)
  expect_equal(ncol(fit$prevalence_draws), 6)          # time + 5 draws
  expect_equal(length(unique(fit$trajectories$draw)), 5)
  # infection_time is 0 (never) or a real timestep.
  expect_true(all(fit$trajectories$infection_time >= 0))
  # Without traj_draws those outputs are absent, not empty.
  bare <- badger_hmm(m, "fixed", "nuts", draws = 50, warmup = 50, seed = 2)
  expect_null(bare$prevalence_draws)
  expect_null(bare$trajectories)
})

test_that("the year process is a model choice and changes the answer", {
  have_binary(); m <- test_data()
  a <- badger_hmm(m, "fixed", "map", year_process = "rw1")
  b <- badger_hmm(m, "fixed", "map", year_process = "iid")
  expect_equal(a$settings$year_process, "rw1")
  # Same parameter count, different fit: rw1 and iid are different models.
  expect_equal(nrow(a$parameters), nrow(b$parameters))
  va <- a$parameters$estimate; vb <- b$parameters$estimate
  expect_false(isTRUE(all.equal(va, vb, tolerance = 1e-6)))
})

test_that("a metric can be exported and reused, without a rebuild", {
  have_binary(); m <- test_data()
  # `hmc` now always needs a metric, because the compiled-in ones belong to the
  # superseded variants. The workflow is: adapt once with nuts, export, then
  # reuse under hmc. This is the whole reason the metric is not locked inside
  # the binary.
  mf <- tempfile(fileext = ".csv")
  a <- badger_hmm(m, "fixed", "nuts", draws = 100, warmup = 100,
                  write_metric = mf, seed = 4)
  expect_true(file.exists(mf))
  metric <- as.matrix(utils::read.csv(mf, header = FALSE))
  # Square, and sized to the sampled parameter vector (ncol of draws, not
  # nrow of the parameters table -- that lists Se/Sp rows too).
  expect_equal(nrow(metric), ncol(metric))
  expect_equal(nrow(metric), ncol(a$draws))
  expect_true(all(diag(metric) > 0))

  # Replayed under hmc, with the step size and length supplied alongside.
  b <- badger_hmm(m, "fixed", "hmc", draws = 100, warmup = 100,
                  metric = mf, step_size = 0.1, n_leapfrog = 10, seed = 4)
  expect_equal(nrow(b$draws), 100)
  expect_equal(b$settings$metric, mf)
  expect_equal(b$settings$method, "hmc")

  expect_error(badger_hmm(m, "fixed", "hmc", metric = "no_such_file.csv",
                          step_size = 0.1, n_leapfrog = 10),
               "metric file not found")
})

test_that("the default method samples without any pre-tuned constants", {
  # nuts adapts from scratch, so it works on any data with nothing compiled in.
  # There is no fallback path any more, and no warning: the default just runs.
  have_binary(); m <- test_data()
  f <- expect_no_warning(badger_hmm(m, "fixed", draws = 50, warmup = 50))
  expect_equal(f$settings$method, "nuts")
  expect_equal(nrow(f$draws), 50)
})

test_that("nutsthenhmc and nuts agree on the posterior", {
  # Two samplers on one target must give one answer. A bad hand-off from the
  # adapting phase to the fixed-length phase would still produce draws -- of
  # the wrong distribution.
  have_binary(); m <- test_data()
  a <- badger_hmm(m, "fixed", "nuts", draws = 600, warmup = 400, seed = 11)
  b <- badger_hmm(m, "fixed", "nutsthenhmc", draws = 600, warmup = 400, seed = 11)
  for (j in seq_len(ncol(a$draws))) {
    ma <- mean(a$draws[[j]]); mb <- mean(b$draws[[j]])
    sa <- stats::sd(a$draws[[j]])
    expect_lt(abs(ma - mb), 0.4 * max(sa, 1e-6))
  }
})

test_that("a fit is reproducible from its seed", {
  have_binary(); m <- test_data()
  f1 <- badger_hmm(m, "fixed", "nuts", draws = 50, warmup = 50, seed = 7)
  f2 <- badger_hmm(m, "fixed", "nuts", draws = 50, warmup = 50, seed = 7)
  expect_equal(f1$draws, f2$draws)
})

test_that("the three-state model is opt-in and validated up front", {
  fml <- formals(badger_hmm)
  expect_equal(fml$states, 2L)
  expect_equal(eval(fml$traj_mode)[1], "threshold")
  expect_null(fml$e_mult)
  expect_equal(fml$init_jitter, 0)
  m <- matrix(1, 2, 9); m[, 1] <- c(1, 2)
  expect_error(badger_hmm(m, "inferred", "map", states = 4L), "2 or 3")
  expect_error(badger_hmm(m, "inferred", "map", states = 3L, e_mult = 0.5),
               "values in")
  expect_error(badger_hmm(m, "inferred", "map", states = 3L,
                          e_mult = c(NA, NA, 1.5, NA, NA, NA)), "values in")
  expect_error(badger_hmm(m, "inferred", "map", init_jitter = -1),
               "non-negative")
})

test_that("a three-state fit splits infected into exposed and infectious", {
  have_binary()
  fit <- badger_hmm(test_data(), "inferred", "nuts", states = 3L,
                    draws = 40L, warmup = 60L, traj_draws = 8L,
                    traj_mode = "sample", seed = 3L)
  expect_true(all(c("tau", "rho", paste0("e_mult[", 1:6, "]")) %in%
                    fit$parameters$parameter))
  expect_true(all(c("p_exposed", "p_infectious") %in% names(fit$p_infected)))
  expect_equal(fit$p_infected$p_infected,
               fit$p_infected$p_exposed + fit$p_infected$p_infectious,
               tolerance = 1e-6)
  tj <- fit$trajectories
  expect_true("infectious_time" %in% names(tj))
  # infectious_time needs an infection_time and can't precede it
  inf <- tj$infectious_time > 0
  expect_true(all(tj$infection_time[inf] > 0))
  expect_true(all(tj$infectious_time[inf] >= tj$infection_time[inf]))
  expect_false(is.null(fit$prevalence_infectious_draws))
  expect_equal(fit$settings$states, 3L)
})

test_that("a fixed multiplier is a constant, not a parameter", {
  have_binary()
  fit <- badger_hmm(test_data(), "inferred", "map", states = 3L,
                    e_mult = c(NA, NA, 0, NA, NA, 1))
  p <- fit$parameters
  expect_equal(p$estimate[p$parameter == "e_mult[3]"], 0)
  expect_equal(p$estimate[p$parameter == "e_mult[6]"], 1)
  free <- p$estimate[p$parameter %in% paste0("e_mult[", c(1, 2, 4, 5), "]")]
  expect_true(all(free >= 0 & free <= 1))
})

test_that("jittered starts differ between seeds and repeat within one", {
  have_binary()
  f <- function(seed) badger_hmm(test_data(), "inferred", "nuts", states = 3L,
                                 draws = 20L, warmup = 20L, seed = seed,
                                 init_jitter = 1)$draws
  expect_equal(f(1L), f(1L))
  expect_false(isTRUE(all.equal(f(1L), f(2L))))
})

test_that("birth is validated before the binary is invoked", {
  m <- matrix(1, 2, 9); m[, 1] <- c(1, 2)
  expect_null(formals(badger_hmm)$birth)
  expect_error(badger_hmm(m, "inferred", "map", birth = c(3, 4)), "named by")
  expect_error(badger_hmm(m, "inferred", "map", birth = c(a = 3)), "named by")
  expect_error(badger_hmm(m, "inferred", "map", birth = c(`1` = 2.5)),
               "whole timesteps")
})

test_that("a badger with a known birth starts susceptible there", {
  have_binary()
  m <- test_data()
  first <- tapply(m[, "time"], m[, "id"], min)
  ok <- first > 4
  birth <- first[ok] - 4                 # named by id already
  birth[1] <- NA                         # one explicitly unknown
  fit0 <- badger_hmm(m, "inferred", "map", states = 3L)
  fit  <- badger_hmm(m, "inferred", "map", states = 3L, birth = birth)
  expect_equal(fit$settings$n_birth, sum(!is.na(birth)))
  expect_true(any(grepl("^births: ", fit$log)))
  # four extra rows per born badger, p = 0 at birth
  p <- fit$p_infected
  expect_equal(nrow(p), nrow(fit0$p_infected) + 4L * sum(!is.na(birth)))
  born <- names(birth)[!is.na(birth)]
  at_birth <- p[paste(p$id, p$time) %in% paste(born, birth[born]), ]
  expect_equal(nrow(at_birth), length(born))
  expect_true(all(at_birth$p_infected == 0))
  # other badgers unchanged
  others <- setdiff(unique(p$id), as.integer(born))
  expect_equal(sum(p$id %in% others), sum(fit0$p_infected$id %in% others))
})

test_that("init_starts picks the best of several optimisations", {
  expect_equal(formals(badger_hmm)$init_starts, 1L)
  m <- matrix(1, 2, 9); m[, 1] <- c(1, 2)
  expect_error(badger_hmm(m, "inferred", "map", init_starts = 0L),
               "positive integer")
  have_binary()
  lp <- function(fit) as.numeric(sub("^log posterior: ([-0-9.eE+]+).*$", "\\1",
                                     grep("^log posterior:", fit$log, value = TRUE)))
  one  <- badger_hmm(test_data(), "inferred", "map", states = 3L)
  many <- badger_hmm(test_data(), "inferred", "map", states = 3L,
                     init_starts = 6L)
  expect_length(grep("^init start ", many$log), 6L)
  expect_true(any(grepl("^init: best of 6 starts", many$log)))
  expect_gte(lp(many), lp(one) - 1e-6)
  # starts don't depend on the seed
  again <- badger_hmm(test_data(), "inferred", "map", states = 3L,
                      init_starts = 6L, seed = 99L)
  expect_equal(again$parameters, many$parameters)
})
