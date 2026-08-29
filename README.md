# BadgerInfeR

Badger diagnostic HMM inference, wrapping the standalone `test2inf` binary
instead of a Julia session. No Julia, no JuliaCall, no compilation at run time.

Restricted on purpose to the five model variants used in the badger project's
sql-e2e bundles: all tests, no-DPP and culture-only panels, each with Se/Sp
fixed at Table 1 values or inferred.

```r
library(BadgerInfeR)
badger_platform()               # is a usable binary present?

fit <- badger_hmm(test_matrix, model = "all_fixed", start_year = 1976,
                  draws = 1000, warmup = 200, traj_draws = 200)
fit$prevalence
```

## Portability: read this before assuming it will work

**The package is only as portable as the binary it wraps, and it does not ship
one.** `install_github()` alone gives you a package that errors on first use
with instructions, on every platform including Windows.

Why: the bundle is ~23 MB per platform, and JuliaC has no cross-compilation, so
a bundle must be **built on the platform it will run on**. As of now only a
Windows x86-64 bundle has been built and tested.

| platform | status |
|---|---|
| Windows x86-64 | built, tested, all results in this project produced with it |
| Linux x86-64 | not built — `build.sh` should work there, untested |
| macOS (Intel or ARM) | not built; also needs codesigning/notarisation, or Gatekeeper quarantines a downloaded binary |

To get a working install on a platform with no bundle:

1. clone the CLI source, install Julia 1.12+, run `test2inf-cli/build.sh`
2. set `TEST2INF_BIN` to `dist/bin/test2inf` (or `test2inf.exe`)

Ship the **whole `dist/` directory**, not just the executable: the runtime
libraries beside it are required.

### Making it genuinely installable

Not done yet. The realistic route, per the CLI's handover:

- build bundles on each target platform (CI with Windows/Linux/macOS runners)
- attach them to a GitHub release
- have the package download the right one on first use into a user cache
  (`tools::R_user_dir`), with a checksum check
- codesign and notarise the macOS bundle

Until that exists, treat this as an in-house package for machines that already
have, or can build, a bundle.

## Why it is fast

Two things, both measured on the real cohort (3,224 badgers, 29,737 timepoints,
58 parameters):

- **an analytic gradient** compiled into the binary: 1.63 ms against Mooncake's
  6.02 ms through Turing/DynamicPPL, and 29.5 ms through ForwardDiff
- **pre-tuned dense mass matrices** per model variant, so no adaptation happens
  at run time: 30 s for 1000 draws against 431 s for NUTS adapting from scratch

End to end, the five-model sql-e2e bundle went from about an hour per model to
3m11s for all five.

The pre-tuned metrics are specific to the badger cohort. On other data either
use `method = "nuts"` (adapts from scratch, always correct) or tune once and
reuse:

```r
a <- badger_hmm(m, "all_fixed", "nuts", write_metric = "metric.csv")
b <- badger_hmm(m, "all_fixed", "hmc",  metric = "metric.csv")
```

If a pre-tuned metric does not fit your data, `method = "hmc"` warns and falls
back to NUTS rather than failing.

## Model variants

`all_fixed`, `all_inferred`, `culture_only_fixed`, `no_dpp_fixed`,
`no_dpp_inferred`. Culture-only with inferred Se/Sp is deliberately absent: it
is not identifiable, and which of the two modes the optimiser finds depends only
on the seed.

The year-effect process defaults to `"rw1"`, matching `test2infeR`'s default and
every published sql-e2e fit. `"iid"` is a different model, not a tuning option.
