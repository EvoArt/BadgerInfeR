# Multi-chain fits. The diagnostics tests are pure R and run anywhere; the
# fitting tests need the binary and skip without it.

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

test_that("Rhat is near 1 for agreeing chains and large for disagreeing ones", {
  set.seed(1)
  n <- 400
  # Four chains from the same distribution: Rhat should be ~1.
  same <- data.frame(
    chain = rep(1:4, each = n),
    x = c(rnorm(n), rnorm(n), rnorm(n), rnorm(n))
  )
  d <- BadgerInfeR:::.chain_diagnostics(same)
  expect_equal(d$rhat[d$parameter == "x"], 1, tolerance = 0.05)
  expect_true(d$ess_bulk[d$parameter == "x"] > 100)

  # Four chains centred in different places: Rhat must be clearly above 1.
  # This is the failure the diagnostic exists to catch, so it has to be loud.
  diff <- data.frame(
    chain = rep(1:4, each = n),
    x = c(rnorm(n, 0), rnorm(n, 5), rnorm(n, 10), rnorm(n, 15))
  )
  d2 <- BadgerInfeR:::.chain_diagnostics(diff)
  expect_true(d2$rhat[d2$parameter == "x"] > 1.5)
})

test_that("a constant parameter reports NA rather than NaN", {
  # A fixed-Se/Sp fit echoes constants back. 0/0 is not a convergence failure,
  # and reporting NaN there would make every fixed fit look broken.
  df <- data.frame(chain = rep(1:4, each = 50), k = 0.407)
  d <- BadgerInfeR:::.chain_diagnostics(df)
  expect_true(is.na(d$rhat[d$parameter == "k"]))
  expect_false(is.nan(d$rhat[d$parameter == "k"]))
})

test_that("diagnostics need at least two chains", {
  df <- data.frame(chain = rep(1L, 100), x = rnorm(100))
  expect_null(BadgerInfeR:::.chain_diagnostics(df))
})

test_that("stacking tags each chain and drops absent elements", {
  fits <- list(
    list(draws = data.frame(a = 1:3), trajectories = NULL),
    list(draws = data.frame(a = 4:6), trajectories = NULL)
  )
  st <- BadgerInfeR:::.stack_chain_frames(fits, "draws")
  expect_equal(nrow(st), 6)
  expect_equal(st$chain, rep(1:2, each = 3))
  # Absent everywhere means NULL, not a zero-row frame: that is the contract a
  # single-chain fit has for `method = "map"`.
  expect_null(BadgerInfeR:::.stack_chain_frames(fits, "trajectories"))
})

test_that("chains must be a positive integer and map refuses to fan out", {
  m <- matrix(0, 1, 9)
  expect_error(badger_hmm(m, chains = 0), "positive integer")
  expect_error(badger_hmm(m, chains = -2), "positive integer")
  # MAP is an optimisation: four chains would return the same mode four times.
  expect_error(badger_hmm(m, method = "map", chains = 4), "meaningless")
})

test_that("chains = 1 returns exactly the single-chain object", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "all_fixed", "nuts", draws = 50, warmup = 50, seed = 3)
  expect_s3_class(fit, "badger_fit")
  expect_equal(fit$settings$chains, 1L)
  # No chain column and no pooling machinery on the single-chain path.
  expect_false("chain" %in% names(fit$draws))
  expect_null(fit$diagnostics)
})

test_that("four chains pool their draws and carry convergence diagnostics", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "all_fixed", "nuts", draws = 100, warmup = 100,
                    seed = 1, chains = 4)
  expect_s3_class(fit, "badger_fit")
  expect_equal(nrow(fit$draws), 400)
  expect_equal(sort(unique(fit$draws$chain)), 1:4)
  expect_equal(fit$settings$chains, 4L)
  # Seeds must differ, or the four chains are one chain copied four times.
  expect_equal(fit$settings$seeds, 1:4)
  expect_equal(length(unique(fit$settings$seeds)), 4L)
  expect_true(!is.null(fit$diagnostics))
  expect_true(all(c("rhat", "ess_bulk") %in% names(fit$parameters)))
  expect_equal(length(fit$chain_times), 4L)
  expect_equal(length(fit$log), 4L)
})

test_that("chains differ from one another but agree on the posterior", {
  have_binary(); m <- test_data()
  fit <- badger_hmm(m, "all_fixed", "nuts", draws = 200, warmup = 200,
                    seed = 1, chains = 4)
  pars <- setdiff(names(fit$draws), "chain")
  p <- pars[which.max(vapply(fit$draws[pars], stats::sd, numeric(1)))]
  by_chain <- split(fit$draws[[p]], fit$draws$chain)
  # Not identical: that would mean the seeds did not take effect.
  expect_false(isTRUE(all.equal(by_chain[[1]], by_chain[[2]])))
  # But converged: Rhat near 1 on a parameter that actually varies.
  r <- fit$diagnostics$rhat[fit$diagnostics$parameter == p]
  expect_true(is.na(r) || r < 1.1)
})

test_that("running four chains beats running them one after another", {
  have_binary(); m <- test_data()
  # The binary is single-threaded, so concurrency is the only mechanism; if
  # this ever stops holding, the feature has silently become pointless.
  t_par <- system.time(
    badger_hmm(m, "all_fixed", "nuts", draws = 200, warmup = 200,
               seed = 1, chains = 4))[["elapsed"]]
  t_ser <- system.time(
    badger_hmm(m, "all_fixed", "nuts", draws = 200, warmup = 200,
               seed = 1, chains = 4, cores = 1))[["elapsed"]]
  expect_lt(t_par, t_ser)
})
