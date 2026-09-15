# The four model variants, named once so scripts cannot drift apart.
#
# These are the fits in `identify differences/the case`: fixed or inferred
# test accuracy, each with culture specificity either free or pinned to one.
# All four share the same quarterly grid, stacked within-quarter captures and
# `rw1c` + data-adaptive year model; those were held constant deliberately, so
# that what differs between the four is the accuracy assumption and nothing
# else.
#
# The earlier five-variant set (all_fixed, all_inferred, culture_only_fixed,
# no_dpp_fixed, no_dpp_inferred) is superseded and has been removed. Those
# fits used plain `rw1`, which leaves the level of the year series
# unidentified, and Beta accuracy priors, whose shape floor distorts the
# requested moments. Anything still reading those names wants the archived
# outputs, not a fresh fit.

#' Test names, in the column order the model expects
#' @export
BADGER_TESTS <- c("ELISA_BELT", "ELISA_INDIRECT", "Culture", "DPP",
                  "IGRA", "StatPak")

#' Table 1 sensitivity values, in `BADGER_TESTS` order
#'
#' Used directly when accuracy is fixed, and as the prior MEDIAN when it is
#' inferred.
#' @export
BADGER_SE <- c(0.407, 0.407, 0.100, 0.650, 0.809, 0.492)

#' Table 1 specificity values, in `BADGER_TESTS` order
#' @export
BADGER_SP <- c(0.943, 0.943, 0.999, 0.943, 0.936, 0.931)

#' Index of the culture assay within [BADGER_TESTS]
#'
#' Named rather than written as a bare `3`, because pinning the wrong assay's
#' specificity would be a silent error: the fit would succeed and mean
#' something different.
#' @export
BADGER_CULTURE <- 3L

# name -> (infer_sesp, sp_pin). sp_pin is NA for "leave it to the model".
.BADGER_MODELS <- list(
  fixed          = list(infer_sesp = FALSE, sp_pin = NULL),
  fixed_cult1    = list(infer_sesp = FALSE, sp_pin = c(BADGER_CULTURE, 1)),
  inferred       = list(infer_sesp = TRUE,  sp_pin = NULL),
  inferred_cult1 = list(infer_sesp = TRUE,  sp_pin = c(BADGER_CULTURE, 1))
)

.BADGER_LABELS <- c(
  fixed          = "Fixed accuracy",
  fixed_cult1    = "Fixed accuracy, culture Sp = 1",
  inferred       = "Inferred accuracy",
  inferred_cult1 = "Inferred accuracy, culture Sp = 1"
)

#' The supported model variants
#'
#' Four fits differing only in how test accuracy is treated:
#'
#' \describe{
#'   \item{`fixed`}{Sensitivity and specificity held at the Table 1 values.
#'     The main comparator, and the default.}
#'   \item{`fixed_cult1`}{As `fixed`, but culture specificity is one. This is
#'     the carry-forward assumption -- that a culture positive is never a false
#'     positive -- stated inside the model rather than imposed on the data.}
#'   \item{`inferred`}{Accuracy estimated, with logit-normal priors centred on
#'     the Table 1 values as medians.}
#'   \item{`inferred_cult1`}{As `inferred`, but culture specificity pinned to
#'     one, so it is no longer a parameter at all.}
#' }
#'
#' @param labels If `TRUE`, return the display names rather than the keys.
#' @return A character vector of model names.
#' @export
badger_models <- function(labels = FALSE) {
  if (labels) unname(.BADGER_LABELS[names(.BADGER_MODELS)])
  else names(.BADGER_MODELS)
}

#' Display label for a model key
#' @param model A model key from [badger_models()].
#' @return The display name for that variant.
#' @export
badger_label <- function(model) {
  model <- match.arg(model, names(.BADGER_MODELS))
  unname(.BADGER_LABELS[[model]])
}

resolve_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    stop("`model` must be a single string; one of: ",
         paste(names(.BADGER_MODELS), collapse = ", "))
  }
  if (!model %in% names(.BADGER_MODELS)) {
    old <- c("all_fixed", "all_inferred", "culture_only_fixed",
             "no_dpp_fixed", "no_dpp_inferred")
    hint <- if (model %in% old) {
      paste0("\n'", model, "' was one of the earlier five variants. Those are ",
             "superseded: they used plain 'rw1', which leaves the year ",
             "series' level unidentified, and Beta accuracy priors whose ",
             "shape floor distorts the requested moments. Use 'fixed' for the ",
             "closest current equivalent.")
    } else ""
    stop("Unknown model '", model, "'. Supported models are: ",
         paste(names(.BADGER_MODELS), collapse = ", "), ".", hint)
  }
  .BADGER_MODELS[[model]]
}
