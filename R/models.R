# The four model variants, named once so scripts cannot drift apart. These are
# the fits in `identify differences/the case`. All four share the same
# quarterly grid, stacked captures and rw1c+yw year model, so what differs
# between them is the accuracy assumption and nothing else.
#
# The earlier five (all_fixed, all_inferred, culture_only_fixed, no_dpp_fixed,
# no_dpp_inferred) are superseded: plain rw1 leaves the year series' level
# unidentified, and Beta accuracy priors distort the requested moments.

#' Test names, in the column order the model expects
#' @export
BADGER_TESTS <- c("ELISA_BELT", "ELISA_INDIRECT", "Culture", "DPP",
                  "IGRA", "StatPak")

#' Table 1 sensitivity values, in `BADGER_TESTS` order
#'
#' Used directly when accuracy is fixed, and as the prior median when inferred.
#' @export
BADGER_SE <- c(0.407, 0.407, 0.100, 0.650, 0.809, 0.492)

#' Table 1 specificity values, in `BADGER_TESTS` order
#' @export
BADGER_SP <- c(0.943, 0.943, 0.999, 0.943, 0.936, 0.931)

#' Index of the culture assay within [BADGER_TESTS]
#'
#' Named rather than a bare `3`: pinning the wrong assay would fit silently.
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
#' Four fits differing only in how test accuracy is treated: `fixed`,
#' `fixed_cult1`, `inferred`, `inferred_cult1`. The `cult1` pair pin culture
#' specificity to one.
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
