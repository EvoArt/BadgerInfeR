# The five model variants, named once so scripts cannot drift apart.
#
# These are exactly the fits in `sql e2e/scripts/02_fit_and_plot.R` and
# `sql e2e nuts/scripts/02_fit.R`. Restricting the package to them is deliberate:
# the pre-tuned HMC parameters compiled into the binary are specific to these
# variants, and culture-only + inferred Se/Sp is not identifiable at all.

#' Test names, in the column order the model expects
#' @export
BADGER_TESTS <- c("ELISA_BELT", "ELISA_INDIRECT", "Culture", "DPP",
                  "IGRA", "StatPak")

#' Table 1 sensitivity values, in `BADGER_TESTS` order
#' @export
BADGER_SE <- c(0.407, 0.407, 0.100, 0.650, 0.809, 0.492)

#' Table 1 specificity values, in `BADGER_TESTS` order
#' @export
BADGER_SP <- c(0.943, 0.943, 0.999, 0.943, 0.936, 0.931)

# name -> (test indices, infer_sesp). NULL tests means all six.
.BADGER_MODELS <- list(
  all_fixed          = list(tests = NULL,             infer_sesp = FALSE),
  all_inferred       = list(tests = NULL,             infer_sesp = TRUE),
  culture_only_fixed = list(tests = 3L,               infer_sesp = FALSE),
  no_dpp_fixed       = list(tests = c(1L,2L,3L,5L,6L), infer_sesp = FALSE),
  no_dpp_inferred    = list(tests = c(1L,2L,3L,5L,6L), infer_sesp = TRUE)
)

# The display names used in the sql-e2e CSV outputs, so a bundle written with
# this package is column-compatible with one written by test2infeR.
.BADGER_LABELS <- c(
  all_fixed          = "All tests (fixed)",
  all_inferred       = "All tests (inferred)",
  culture_only_fixed = "Culture only (fixed)",
  no_dpp_fixed       = "No DPP (fixed)",
  no_dpp_inferred    = "No DPP (inferred)"
)

#' The supported model variants
#'
#' @param labels If `TRUE`, return the display names used in bundle outputs
#'   ("All tests (fixed)" and so on) rather than the internal keys.
#' @return A character vector of model names.
#' @export
badger_models <- function(labels = FALSE) {
  if (labels) unname(.BADGER_LABELS[names(.BADGER_MODELS)])
  else names(.BADGER_MODELS)
}

#' Display label for a model key
#' @param model A model key from [badger_models()].
#' @return The display name used in bundle CSV columns.
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
    stop("Unknown model '", model, "'. Supported models are: ",
         paste(names(.BADGER_MODELS), collapse = ", "),
         ".\nThis package is deliberately restricted to the five variants used ",
         "in the sql-e2e bundles.")
  }
  .BADGER_MODELS[[model]]
}
