# Locating and running the trimmed test2inf executable.
#
# The whole point of this package is that the target machine needs neither Julia
# nor R-side compilation: `test2inf` is a self-contained ~23 MB bundle. That
# bundle is too large to ship inside the package, so it is located at runtime.

#' Locate the test2inf binary
#'
#' Resolution order:
#' 1. the `TEST2INF_BIN` environment variable, if it points at a real file;
#' 2. `inst/bin/` inside the installed package;
#' 3. a `dist/bin/` directory beside the working directory (development layout).
#'
#' @return Absolute path to the executable.
#' @export
badger_binary <- function() {
  exe_name <- if (.Platform$OS.type == "windows") "test2inf.exe" else "test2inf"

  env <- Sys.getenv("TEST2INF_BIN", "")
  if (nzchar(env)) {
    if (!file.exists(env)) {
      stop("TEST2INF_BIN is set to '", env, "' but no such file exists.")
    }
    return(normalizePath(env))
  }

  pkg <- system.file("bin", exe_name, package = "BadgerInfeR")
  if (nzchar(pkg) && file.exists(pkg)) return(normalizePath(pkg))

  dev <- file.path(getwd(), "dist", "bin", exe_name)
  if (file.exists(dev)) return(normalizePath(dev))

  # A bundle downloaded by install_test2inf() on a previous run. Cache
  # directories are named "<platform>-<tag>"; prefer the most recently written,
  # so an upgrade takes effect without the user clearing anything.
  plat <- badger_platform_tag()
  if (!is.na(plat)) {
    root <- badger_cache_dir()
    hits <- list.files(root, pattern = paste0("^", plat, "-"), full.names = TRUE)
    hits <- file.path(hits, "bin", exe_name)
    hits <- hits[file.exists(hits)]
    if (length(hits) > 0) {
      newest <- hits[which.max(file.mtime(hits))]
      return(normalizePath(newest))
    }
  }

  # Nothing local: fetch the bundle built for this platform, unless the user
  # has opted out of network access.
  if (!is.na(plat) && !nzchar(Sys.getenv("BADGERINFER_NO_DOWNLOAD", ""))) {
    got <- tryCatch(install_test2inf(), error = function(e) {
      message("Automatic download failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(got) && file.exists(got)) return(got)
  }

  stop("Could not find the test2inf binary.\n\n",
       "BadgerInfeR does not ship one: the bundle is ~23 MB per platform and ",
       "JuliaC has no cross-compilation, so each platform's bundle must be ",
       "built on that platform.\n\n",
       "This machine is ", .Platform$OS.type, "/", R.version$arch,
       ", so it needs a bundle built there.\n\n",
       "To get one:\n",
       "  * build it:   test2inf-cli/build.sh   (needs Julia 1.12+), then set\n",
       "                TEST2INF_BIN to dist/bin/", exe_name, "\n",
       "  * or install a prebuilt bundle into inst/bin/ of this package.\n\n",
       "Ship the whole dist/ directory, not just the executable: the runtime ",
       "libraries beside it are required.")
}

#' Report what platform a binary would be needed for
#'
#' `BadgerInfeR` is a thin wrapper around a compiled executable, so it is only
#' as portable as that executable. There is no cross-compilation: a bundle has
#' to be built on the platform it will run on.
#'
#' @return Invisibly, a list describing the current platform and whether a
#'   usable binary was found.
#' @export
badger_platform <- function() {
  exe_name <- if (.Platform$OS.type == "windows") "test2inf.exe" else "test2inf"
  found <- tryCatch(badger_binary(), error = function(e) NA_character_)
  info <- list(os = .Platform$OS.type,
               sysname = unname(Sys.info()[["sysname"]]),
               arch = R.version$arch,
               expected_binary = exe_name,
               binary = found,
               available = !is.na(found))
  cat("BadgerInfeR platform check\n")
  cat("  OS:      ", info$sysname, " (", info$os, ")\n", sep = "")
  cat("  arch:    ", info$arch, "\n", sep = "")
  cat("  binary:  ", if (info$available) info$binary else "NOT FOUND", "\n", sep = "")
  if (!info$available) {
    cat("\n  Build one with test2inf-cli/build.sh on this platform, then set\n")
    cat("  TEST2INF_BIN, or install a prebuilt bundle into inst/bin/.\n")
  }
  invisible(info)
}

#' Check that the binary runs
#'
#' Invokes `--help` and returns invisibly. Useful as a first call in a script,
#' so a missing or wrong-architecture binary fails immediately rather than
#' partway through a long fit.
#'
#' @return `TRUE`, invisibly.
#' @export
check_badger_binary <- function() {
  out <- suppressWarnings(
    system2(badger_binary(), "--help", stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("The test2inf binary did not run:\n", paste(out, collapse = "\n"))
  }
  invisible(TRUE)
}

# Run the binary, returning its stdout. Errors carry the binary's own output,
# which is where the model diagnostics (divergences, step size) are reported.
run_test2inf <- function(args) {
  out <- suppressWarnings(
    system2(badger_binary(), as.character(args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("test2inf failed (exit ", status, "):\n", paste(out, collapse = "\n"))
  }
  out
}
