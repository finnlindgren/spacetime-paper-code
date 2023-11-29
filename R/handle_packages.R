#' Warn about packages without the required versions
#'
#' @param needed a vector with named character entries indicating package name
#' and version pairs. Warnings are given for packages that are not installed or
#' have older versions, and installation instructions are shown.
#' @param attach logical; if `TRUE`, run `library(packagename)` for each needed
#' package and give an error if some package isn't installed.
#'
#' @return logical; TRUE if all the packages are installed and have the
#' required versions
handle_packages <- function(needed, attach = FALSE) {
  installed <- names(needed) %in% installed.packages()
  names(installed) <- names(needed)
  version <- vapply(
    names(needed),
    function(x) {
      if (installed[x]) {
        getNamespaceVersion(x)
      } else {
        NA_character_
      }
    },
    ""
  )
  names(version) <- names(needed)
  version_ok <- vapply(
    names(needed),
    function(x) {
      compareVersion(version[x], needed[x]) >= 0
    },
    TRUE
  )
  names(version_ok) <- names(needed)

  special_installation <-
    c(
      "INLA" = 'install.packages("INLA", repos = c(options("repos", INLA="https://inla.r-inla-download.org/R/testing")))'
      # "INLAspacetime" = 'remotes::install_github("eliaskrainski/INLAspacetime")'
    )

  for (pkg in names(needed)[!version_ok]) {
    inst <- if (pkg %in% names(special_installation)) {
      special_installation[pkg]
    } else {
      paste0('install.packages("', pkg, '")')
    }
    if (!installed[pkg]) {
      warning(
        paste0(
          "The ", pkg, " package is not installed. Install with\n",
          "  ", inst
        ),
        call. = FALSE
      )
    } else {
      warning(
        paste0(
          pkg,
          " >= ",
          needed[pkg],
          " needed, but ",
          version[pkg],
          " is installed. Upgrade with\n",
          "  ",
          inst
        ),
        call. = FALSE
      )
    }
  }

  if (attach) {
    if (!all(version_ok)) {
      stop("Packages not fully installed.")
    }
    pkgs_ <- lapply(names(needed), library, character.only = TRUE)
  }

  # Return TRUE if all is ok
  all(version_ok)
}
