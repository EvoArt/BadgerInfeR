# Getting a platform-appropriate binary onto the user's machine.
#
# The bundle is ~115 MB per platform (the whole Julia runtime; trimming it
# safely needs the transitive closure of the import graph, which nothing here
# computes) and JuliaC cannot cross-compile, so one bundle per platform is
# built by CI and attached to a GitHub release. The package itself ships no
# binary; it fetches the right one on first use into a user cache directory.

TEST2INF_RELEASE_TAG <- "latest"
TEST2INF_REPO <- "EvoArt/test2inf-cli"

#' The platform tag this machine needs
#'
#' Matches the `label` values in `.github/workflows/build-cli.yml`.
#' @return A string such as `"windows-x86_64"`, or `NA` if unsupported.
#' @export
badger_platform_tag <- function() {
  sysname <- unname(Sys.info()[["sysname"]])
  arch <- R.version$arch
  arm <- grepl("aarch64|arm64", arch, ignore.case = TRUE)
  if (sysname == "Windows") return("windows-x86_64")
  if (sysname == "Linux")   return(if (arm) NA_character_ else "linux-x86_64")
  # Intel macOS has no prebuilt bundle: GitHub's macos-13 runner never
  # started, so CI cannot produce one. Report it as unsupported rather than
  # hand back an ARM tag that would download a bundle this machine cannot run.
  if (sysname == "Darwin")  return(if (arm) "macos-arm64" else NA_character_)
  NA_character_
}

#' Where downloaded bundles are kept
#' @return Path to the cache directory (not created until needed).
#' @export
badger_cache_dir <- function() {
  root <- Sys.getenv("BADGERINFER_CACHE", "")
  if (nzchar(root)) return(root)
  tools::R_user_dir("BadgerInfeR", which = "cache")
}

#' Download and install the test2inf binary
#'
#' Fetches the bundle built for this platform from the project's GitHub
#' releases, verifies its checksum, and unpacks it into [badger_cache_dir()].
#'
#' Called automatically by [badger_hmm()] on first use unless
#' `BADGERINFER_NO_DOWNLOAD` is set, in which case it must be run explicitly.
#'
#' @param tag Release tag to install from; `"latest"` by default.
#' @param force Reinstall even if a binary is already present.
#' @param quiet Suppress progress messages.
#' @return Path to the installed executable, invisibly.
#' @export
install_test2inf <- function(tag = TEST2INF_RELEASE_TAG, force = FALSE,
                             quiet = FALSE) {
  plat <- badger_platform_tag()
  if (is.na(plat)) {
    stop("No prebuilt bundle exists for ", Sys.info()[["sysname"]], "/",
         R.version$arch, ".\n",
         "Build one with test2inf-cli/build.sh (needs Julia 1.12+) and set ",
         "TEST2INF_BIN to dist/bin/test2inf.")
  }

  exe_name <- if (.Platform$OS.type == "windows") "test2inf.exe" else "test2inf"
  # Cache per RELEASE TAG, not just per platform. Keying on platform alone
  # meant an existing binary was reused whatever release it came from, so
  # upgrading silently kept the old one -- which hid a fix for a binary that
  # would not start on some CPUs. "latest" gets its own slot and is always
  # re-fetched, since its contents change without the key changing.
  dest <- file.path(badger_cache_dir(), paste0(plat, "-", tag))
  exe <- file.path(dest, "bin", exe_name)
  stale <- identical(tag, "latest")
  if (!force && !stale && file.exists(exe)) {
    return(invisible(normalizePath(exe)))
  }

  base <- if (identical(tag, "latest")) {
    sprintf("https://github.com/%s/releases/latest/download", TEST2INF_REPO)
  } else {
    sprintf("https://github.com/%s/releases/download/%s", TEST2INF_REPO, tag)
  }
  asset <- sprintf("test2inf-%s.tar.gz", plat)
  url <- file.path(base, asset)

  tmp <- tempfile(fileext = ".tar.gz")
  on.exit(unlink(tmp), add = TRUE)
  if (!quiet) message("Downloading ", asset, " ...")
  ok <- tryCatch({
    utils::download.file(url, tmp, mode = "wb", quiet = quiet)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(tmp) || file.size(tmp) == 0) {
    stop("Could not download ", url, "\n",
         "If the repository is private, download the asset manually and set ",
         "TEST2INF_BIN to the extracted bin/", exe_name, ".")
  }

  # Verify against the published checksum. A truncated or substituted download
  # would otherwise produce silently wrong inference rather than an error.
  verify_checksum(url, tmp, quiet = quiet)

  unlink(dest, recursive = TRUE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tmp, exdir = dest)

  if (!file.exists(exe)) {
    stop("The downloaded bundle does not contain bin/", exe_name, ".")
  }
  if (.Platform$OS.type != "windows") {
    Sys.chmod(exe, "0755")
  }
  if (!quiet) message("Installed to ", dest)
  invisible(normalizePath(exe))
}

verify_checksum <- function(url, path, quiet = FALSE) {
  sum_url <- paste0(url, ".sha256")
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  got <- tryCatch({
    utils::download.file(sum_url, tmp, quiet = TRUE)
    trimws(sub("\\s.*$", "", readLines(tmp, warn = FALSE)[1]))
  }, error = function(e) NA_character_, warning = function(w) NA_character_)

  if (is.na(got) || !nzchar(got)) {
    warning("No published checksum for this asset; the download was not ",
            "verified.", call. = FALSE)
    return(invisible(FALSE))
  }
  # base R has md5sum but no sha256, so this needs `digest`. Treat its absence
  # as "unverified" rather than an error: refusing to install over a missing
  # optional dependency would be worse than warning.
  have <- if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(path, algo = "sha256", file = TRUE)
  } else {
    NA_character_
  }
  if (is.na(have)) {
    warning("Package 'digest' is not installed, so the download's checksum ",
            "could not be verified. Install it for a verified download.",
            call. = FALSE)
    return(invisible(FALSE))
  }
  if (!identical(tolower(have), tolower(got))) {
    stop("Checksum mismatch for the downloaded bundle.\n",
         "  expected: ", got, "\n  got:      ", have, "\n",
         "Refusing to install: a corrupted binary would produce wrong results.")
  }
  if (!quiet) message("Checksum verified.")
  invisible(TRUE)
}
