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

test_that("the model list is exactly the five sql-e2e variants", {
  expect_setequal(badger_models(),
                  c("all_fixed", "all_inferred", "culture_only_fixed",
                    "no_dpp_fixed", "no_dpp_inferred"))
  # culture-only + inferred is unidentifiable and must not be offered at all.
  expect_false("culture_only_inferred" %in% badger_models())
  expect_equal(badger_label("no_dpp_inferred"), "No DPP (inferred)")
  expect_error(badger_hmm(matrix(0, 1, 9), model = "nonsense"), "Unknown model")
})

test_that("the test matrix is validated before the binary is invoked", {
  expect_error(write_test_matrix(matrix(0, 2, 3), tempfile()), "at least 9")
  bad <- matrix(1, 2, 9); bad[1, 3] <- 7
  expect_error(write_test_matrix(bad, tempfile()), "captured")
  bad2 <- matrix(1, 2, 9); bad2[1, 1] <- 0
  expect_error(write_test_matrix(bad2, tempfile()), "positive integer")
})

test_that("NA test results survive the round trip as empty fields", {
  m <- matrix(1, 2, 9)
  m[, 1] <- c(1, 2); m[, 2] <- 1; m[, 3] <- 1
  m[1, 4] <- NA
  p <- tempfile(fileext = ".csv")
  write_test_matrix(m, p)
  txt <- readLines(p)
  expect_true(grepl(",,", txt[2]) || grepl(",$", txt[2]))
})

test_that("MAP runs and returns the documented fields", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "all_fixed", "map", start_year = 1976)
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
  fixed <- badger_hmm(m, "all_fixed", "map")
  inferred <- badger_hmm(m, "all_inferred", "map")
  expect_false(any(fixed$sesp$estimated))
  expect_true(all(inferred$sesp$estimated))
  # A fixed fit must report exactly the Table 1 constants.
  expect_equal(fixed$sesp$Se, BADGER_SE, tolerance = 1e-8)
  expect_equal(fixed$sesp$Sp, BADGER_SP, tolerance = 1e-8)
  # An inferred fit should have moved off them.
  expect_false(isTRUE(all.equal(inferred$sesp$Se, BADGER_SE, tolerance = 1e-6)))
})

test_that("NUTS produces draws, and traj_draws produces trajectory outputs", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "all_fixed", "nuts", draws = 100, warmup = 100,
                    traj_draws = 5, seed = 2)
  expect_equal(nrow(fit$draws), 100)
  expect_equal(ncol(fit$prevalence_draws), 6)          # time + 5 draws
  expect_equal(length(unique(fit$trajectories$draw)), 5)
  # infection_time is 0 (never) or a real timestep.
  expect_true(all(fit$trajectories$infection_time >= 0))
  # Without traj_draws those outputs are absent, not empty.
  bare <- badger_hmm(m, "all_fixed", "nuts", draws = 50, warmup = 50, seed = 2)
  expect_null(bare$prevalence_draws)
  expect_null(bare$trajectories)
})

test_that("the year process is a model choice and changes the answer", {
  have_binary(); m <- test_data()
  a <- badger_hmm(m, "all_fixed", "map", year_process = "rw1")
  b <- badger_hmm(m, "all_fixed", "map", year_process = "iid")
  expect_equal(a$settings$year_process, "rw1")
  # Same parameter count, different fit: rw1 and iid are different models.
  expect_equal(nrow(a$parameters), nrow(b$parameters))
  va <- a$parameters$estimate; vb <- b$parameters$estimate
  expect_false(isTRUE(all.equal(va, vb, tolerance = 1e-6)))
})

test_that("a metric can be exported and reused, without a rebuild", {
  have_binary(); m <- test_data()
  # The compiled-in metrics are tuned for the badger cohort. For any other
  # data the workflow is: adapt once with nuts, export, then reuse under hmc.
  # This is the whole reason the metric is not locked inside the binary.
  mf <- tempfile(fileext = ".csv")
  a <- badger_hmm(m, "all_fixed", "nuts", draws = 100, warmup = 100,
                  write_metric = mf, seed = 4)
  expect_true(file.exists(mf))
  metric <- as.matrix(utils::read.csv(mf, header = FALSE))
  # Square, and sized to the sampled parameter vector (ncol of draws, not
  # nrow of the parameters table -- that lists Se/Sp rows too).
  expect_equal(nrow(metric), ncol(metric))
  expect_equal(nrow(metric), ncol(a$draws))
  expect_true(all(diag(metric) > 0))

  # With a metric that DOES fit, hmc runs and no fallback happens.
  b <- badger_hmm(m, "all_fixed", "hmc", draws = 100, warmup = 100,
                  metric = mf, seed = 4)
  expect_equal(nrow(b$draws), 100)
  expect_equal(b$settings$metric, mf)
  expect_equal(b$settings$method, "hmc")

  expect_error(badger_hmm(m, "all_fixed", "hmc", metric = "no_such_file.csv"),
               "metric file not found")
})

test_that("hmc is the default, and falls back when its metric does not fit", {
  # Pre-tuned HMC is ~14x faster in wall time than adapting from scratch on the
  # badger cohort, which is why it is the default. But the compiled metrics are
  # sized for that cohort (58/70 parameters), so on other data the default must
  # DEGRADE to nuts rather than error -- a default that fails on reasonable
  # input is not a default.
  have_binary(); m <- test_data()
  expect_warning(f <- badger_hmm(m, "all_fixed", draws = 50, warmup = 50),
                 "does not fit this data")
  expect_equal(f$settings$method, "nuts")   # reports what it actually ran
  expect_equal(nrow(f$draws), 50)
})

test_that("a fit is reproducible from its seed", {
  have_binary(); m <- test_data()
  f1 <- badger_hmm(m, "all_fixed", "nuts", draws = 50, warmup = 50, seed = 7)
  f2 <- badger_hmm(m, "all_fixed", "nuts", draws = 50, warmup = 50, seed = 7)
  expect_equal(f1$draws, f2$draws)
})
