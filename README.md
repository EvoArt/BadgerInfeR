# BadgerInfeR

R interface to the badger diagnostic HMM. Wraps the standalone `test2inf`
binary; no Julia or compilation at run time.

## Install

```r
remotes::install_github("EvoArt/BadgerInfeR")
BadgerInfeR::install_test2inf()   # downloads the binary for this platform
```

`install_test2inf()` is called automatically on first use unless
`BADGERINFER_NO_DOWNLOAD` is set. Prebuilt bundles exist for Windows x86-64,
Linux x86-64 and macOS arm64. On any other platform — including Intel macOS,
whose CI runner is being retired — build from
[test2inf-cli](https://github.com/EvoArt/test2inf-cli) and set `TEST2INF_BIN`.

## badger_hmm()

```r
badger_hmm(test_mat,
           model = "fixed",
           method = c("nuts", "nutsthenhmc", "map", "hmc"),
           draws = 1000L, warmup = 1000L, traj_draws = 0L,
           metric = NULL, write_metric = NULL,
           target_acc = 0.8, seed = 1L, chains = 1L, cores = NULL,
           seasons = 4L,
           repeat_captures = c("stack", "pool", "last"),
           year_process = c("rw1c", "rw1", "iid", "rw2", "none", "shrunk"),
           year_weight = TRUE,
           accuracy_prior = c("logitnormal", "beta"),
           penalty = TRUE, start_year = 0L,
           workdir = tempfile("badgerinfer"), keep_workdir = FALSE)
```

| arg | |
|---|---|
| `test_mat` | matrix/data frame: `time, id, captured, test1..test6`. `NA` = not tested. Rows with `captured == 0` optional. |
| `model` | one of `badger_models()` |
| `method` | see Samplers below |
| `draws`, `warmup` | sampling iterations |
| `traj_draws` | draws used for per-draw prevalence and infection times; `0` skips |
| `metric` | path to an inverse mass matrix CSV |
| `write_metric` | path to write the metric used |
| `target_acc` | target acceptance during adaptation |
| `seed` | RNG seed; chain `k` uses `seed + k - 1` |
| `chains`, `cores` | chains to run, and how many at once |
| `seasons` | timesteps per year |
| `repeat_captures` | how to treat multiple captures in one timestep |
| `year_process` | `rw1c`, `rw1`, `iid`, `rw2`, `none`, `shrunk` — different models, not tuning options. `ar1` is unavailable (it adds a parameter). |
| `year_weight` | data-adaptive year-increment prior |
| `accuracy_prior` | `logitnormal` or `beta`; only used when Se/Sp are inferred |
| `penalty` | Se+Sp>1 identifiability penalty |
| `start_year` | calendar year of `time == 1`; `0` omits calendar columns |
| `workdir`, `keep_workdir` | scratch directory for the binary's I/O |

Returns an object of class `badger_fit`, a list of:

| element | |
|---|---|
| `parameters` | data frame; `estimate` (map) or `mean`/`sd` (sampling) |
| `prevalence` | prevalence per timestep, `n_captures` |
| `p_infected` | `id, time, p_infected` |
| `year_effects` | `year_index, gamma, annual_hazard` |
| `draws` | posterior draws, one column per parameter; `NULL` for map |
| `prevalence_draws` | `time` + one column per draw; `NULL` unless `traj_draws > 0` |
| `trajectories` | `id, draw, infection_time` (`0` = never); `NULL` unless `traj_draws > 0` |
| `sesp` | `test, test_name, Se, Sp, Youden, estimated` |
| `settings` | the arguments used |
| `log` | binary stdout, including divergences and step size |

`method = "hmc"` requires an explicit `metric`: the metrics compiled into the
binary were tuned for the superseded variants (plain `rw1`, Beta priors, no
increment prior), so they do not describe the current posterior. Tune once and
reuse:

```r
a <- badger_hmm(m, "fixed", "nuts", write_metric = "metric.csv")
b <- badger_hmm(m, "fixed", "hmc",  metric = "metric.csv")
```

The default `nuts` needs none of this — it adapts to whatever data it is
given.

## Other functions

| | |
|---|---|
| `badger_models(labels = FALSE)` | the four model names; `labels = TRUE` gives display names |
| `badger_label(model)` | display name for a model key |
| `badger_binary()` | path to the executable |
| `check_badger_binary()` | error if it does not run |
| `badger_platform()` | print OS, arch, and whether a binary was found |
| `badger_platform_tag()` | platform string used for release assets |
| `badger_cache_dir()` | where downloaded bundles are kept |
| `install_test2inf(tag = "latest", force = FALSE, quiet = FALSE)` | download and install the binary |

`BADGER_TESTS`, `BADGER_SE`, `BADGER_SP`: test names and Table 1 values, in
column order. `BADGER_CULTURE` is culture's index within them.

## Models

Four fits, differing only in how test accuracy is treated. Everything else —
the quarterly grid, stacked within-quarter captures, `rw1c` + the increment
prior — is held constant, so what separates them is the accuracy assumption
and nothing else.

| key | Se/Sp | culture Sp |
|---|---|---|
| `fixed` | Table 1 values | 0.999 (Table 1) |
| `fixed_cult1` | Table 1 values | pinned to 1 |
| `inferred` | estimated | estimated |
| `inferred_cult1` | estimated | pinned to 1 |

Pinning culture specificity to one states that a culture positive is never a
false positive — the carry-forward assumption, written into the model rather
than imposed on the data. With `inferred` it removes a parameter; with `fixed`
it overrides the tabulated 0.999.

The earlier five variants (`all_fixed`, `all_inferred`, `culture_only_fixed`,
`no_dpp_fixed`, `no_dpp_inferred`) are **superseded and removed**. They used
plain `rw1`, which leaves the level of the year series unidentified, and Beta
accuracy priors, whose shape floor at 0.5 distorts the requested moments
(culture Sp: asked for a mean of 0.999, got `Beta(3.26, 0.5)`, mean 0.867).
Asking for one by name is an error, not a silent substitution.

## The year model

`rw1c` plus the data-adaptive increment prior, both on by default.

Plain `rw1` leaves the **level** of the year series unidentified, so untested
years drift and the walk extrapolates backwards from the first year with data.
`rw1c` centres the series, handing the level to the seasonal intercepts; the
increment prior then pins an untested year to its neighbour while leaving
well-tested years to the likelihood.

This fixes identifiability; it does not push early hazards toward zero, and
should not — 1976–77 were tested by culture alone (Se = 0.10), so some
residual early hazard is honest.

## Samplers

| method | | |
|---|---|---|
| `nuts` | default | 142 s, min ESS 843 / 1000 draws |
| `nutsthenhmc` | NUTS adapts, jittered-L HMC draws | 99 s, min ESS 46 |
| `map` | constrained-scale mode, as Turing | 9 s, no uncertainty |
| `hmc` | fixed L, needs `metric` | |

Measured on the badger cohort, 1000 draws after 1000 warmup. `nutsthenhmc` is
quicker per iteration and about a thirteenth the effective draws per second;
acceptance and divergences look healthy either way, so only ESS shows it.

## Environment variables

| | |
|---|---|
| `TEST2INF_BIN` | path to the executable, overrides all lookup |
| `BADGERINFER_CACHE` | where to keep downloaded bundles |
| `BADGERINFER_NO_DOWNLOAD` | never download automatically |
