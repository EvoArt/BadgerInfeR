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
           model = "all_fixed",
           method = c("hmc", "nuts", "map"),
           draws = 1000L, warmup = 1000L, traj_draws = 0L,
           metric = NULL, write_metric = NULL,
           target_acc = 0.8, seed = 1L, seasons = 4L,
           repeat_captures = c("stack", "pool", "last"),
           year_process = c("rw1", "iid"),
           penalty = TRUE, start_year = 0L,
           workdir = tempfile("badgerinfer"), keep_workdir = FALSE)
```

| arg | |
|---|---|
| `test_mat` | matrix/data frame: `time, id, captured, test1..test6`. `NA` = not tested. Rows with `captured == 0` optional. |
| `model` | one of `badger_models()` |
| `method` | `hmc` uses a pre-tuned dense metric; `nuts` adapts from scratch; `map` is the posterior mode |
| `draws`, `warmup` | sampling iterations |
| `traj_draws` | draws used for per-draw prevalence and infection times; `0` skips |
| `metric` | path to an inverse mass matrix CSV, overrides the compiled-in one |
| `write_metric` | path to write the metric used |
| `target_acc` | NUTS target acceptance |
| `seed` | RNG seed |
| `seasons` | timesteps per year |
| `repeat_captures` | how to treat multiple captures in one timestep |
| `year_process` | `rw1` and `iid` are different models, not tuning options |
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

`method = "hmc"` falls back to `nuts` with a warning when the compiled metric
does not match the supplied data. To use HMC on other data:

```r
a <- badger_hmm(m, "all_fixed", "nuts", write_metric = "metric.csv")
b <- badger_hmm(m, "all_fixed", "hmc",  metric = "metric.csv")
```

## Other functions

| | |
|---|---|
| `badger_models(labels = FALSE)` | the five model names; `labels = TRUE` gives display names |
| `badger_label(model)` | display name for a model key |
| `badger_binary()` | path to the executable |
| `check_badger_binary()` | error if it does not run |
| `badger_platform()` | print OS, arch, and whether a binary was found |
| `badger_platform_tag()` | platform string used for release assets |
| `badger_cache_dir()` | where downloaded bundles are kept |
| `install_test2inf(tag = "latest", force = FALSE, quiet = FALSE)` | download and install the binary |

`BADGER_TESTS`, `BADGER_SE`, `BADGER_SP`: test names and Table 1 values, in
column order.

## Models

| key | tests | Se/Sp |
|---|---|---|
| `all_fixed` | all six | fixed |
| `all_inferred` | all six | inferred |
| `culture_only_fixed` | Culture | fixed |
| `no_dpp_fixed` | not DPP | fixed |
| `no_dpp_inferred` | not DPP | inferred |

Culture-only with inferred Se/Sp is not offered: it is not identifiable.

## Environment variables

| | |
|---|---|
| `TEST2INF_BIN` | path to the executable, overrides all lookup |
| `BADGERINFER_CACHE` | where to keep downloaded bundles |
| `BADGERINFER_NO_DOWNLOAD` | never download automatically |
