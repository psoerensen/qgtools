# ------------------------------------------------------------------
# Kernel validators
# ------------------------------------------------------------------

#' Validate a kernel object based on inferred type
#' @keywords internal
validate_kernel <- function(kernel, kernel_id) {
  if (!is.list(kernel)) {
    return(paste0(kernel_id, ": kernel must be a list-like object"))
  }

  # --- 1. IID kernels FIRST ---
  if (inherits(kernel, "IIDlist") || (!is.null(kernel$type) && kernel$type == "iid")) {
    return(validate_IIDlist(kernel, kernel_id))
  }

  issues <- character(0)

  # --- 2. Structured kernels ---
  if (!is.null(kernel$ldfiles) || !is.null(kernel$fnLD)) {
    issues <- c(issues, validate_LDlist(kernel, kernel_id))

  } else if (!is.null(kernel$pedigree) || !is.null(kernel$fnPED) || !is.null(kernel$idsP)) {
    issues <- c(issues, validate_PEDlist(kernel, kernel_id))

  } else if (!is.null(kernel$fnG) || !is.null(kernel$idsG)) {
    issues <- c(issues, validate_GRMlist(kernel, kernel_id))

  } else if (!is.null(kernel$bedfiles) || !is.null(kernel$fnBED) || !is.null(kernel$markers)) {
    issues <- c(issues, validate_Glist(kernel, kernel_id))

  } else if (!is.null(kernel$M) || !is.null(kernel$features)) {
    issues <- c(issues, validate_Mlist(kernel, kernel_id))

  } else {
    issues <- c(issues, paste0(kernel_id, ": unable to infer kernel type"))
  }

  unique(issues)
}

# ------------------------------------------------------------------
# GRMlist
# ------------------------------------------------------------------

validate_GRMlist <- function(GRMlist, kernel_id = "GRM") {
  issues <- character(0)

  if (is.null(GRMlist$idsG)) {
    issues <- c(issues, paste0(kernel_id, ": GRMlist$idsG is required"))
  }

  if (is.null(GRMlist$fnG)) {
    issues <- c(issues, paste0(kernel_id, ": GRMlist$fnG is required (disk-backed GRM)"))
  } else if (!is.character(GRMlist$fnG)) {
    issues <- c(issues, paste0(kernel_id, ": GRMlist$fnG must be character path(s)"))
  }

  issues
}

# ------------------------------------------------------------------
# Glist
# ------------------------------------------------------------------

validate_Glist <- function(Glist, kernel_id = "G") {
  issues <- character(0)

  if (is.null(Glist$ids)) {
    issues <- c(issues, paste0(kernel_id, ": Glist$ids is required"))
  }

  if (is.null(Glist$rsids) && is.null(Glist$markers)) {
    issues <- c(issues, paste0(kernel_id, ": Glist$rsids or Glist$markers is required"))
  }

  has_plink   <- !is.null(Glist$bedfiles) || !is.null(Glist$fnBED)
  has_generic <- !is.null(Glist$fn) && !is.null(Glist$n) && !is.null(Glist$m)

  if (!has_plink && !has_generic) {
    issues <- c(
      issues,
      paste0(kernel_id, ": Glist must be PLINK-backed or generic matrix-backed")
    )
  }

  if (has_generic) {
    if (!is.character(Glist$fn)) {
      issues <- c(issues, paste0(kernel_id, ": Glist$fn must be a path string"))
    }
    if (!is.numeric(Glist$n) || !is.numeric(Glist$m)) {
      issues <- c(issues, paste0(kernel_id, ": Glist$n and Glist$m must be numeric"))
    }
  }

  issues
}

# ------------------------------------------------------------------
# PEDlist
# ------------------------------------------------------------------

validate_PEDlist <- function(PEDlist, kernel_id = "PED") {
  issues <- character(0)

  if (is.null(PEDlist$ids) && is.null(PEDlist$idsP)) {
    issues <- c(issues, paste0(kernel_id, ": PEDlist$ids or idsP is required"))
  }

  if (is.null(PEDlist$fnPED) &&
      is.null(PEDlist$pedfile) &&
      is.null(PEDlist$pedigree)) {
    issues <- c(
      issues,
      paste0(kernel_id, ": provide fnPED/pedfile (disk) or pedigree (data.frame)")
    )
  }

  issues
}

# ------------------------------------------------------------------
# Mlist
# ------------------------------------------------------------------

validate_Mlist <- function(Mlist, kernel_id = "M") {
  issues <- character(0)

  if (is.null(Mlist$ids)) {
    issues <- c(issues, paste0(kernel_id, ": Mlist$ids is required"))
  }

  if (is.null(Mlist$features) && is.null(Mlist$markers)) {
    issues <- c(issues, paste0(kernel_id, ": Mlist$features or Mlist$markers is required"))
  }

  if (is.null(Mlist$fn) && is.null(Mlist$M)) {
    issues <- c(issues, paste0(kernel_id, ": provide disk-backed fn or in-memory matrix M"))
  }

  issues
}

# ------------------------------------------------------------------
# LDlist
# ------------------------------------------------------------------

validate_LDlist <- function(LDlist, kernel_id = "LD") {
  issues <- character(0)

  if (is.null(LDlist$ldfiles) && is.null(LDlist$fnLD)) {
    issues <- c(issues, paste0(kernel_id, ": LDlist needs ldfiles or fnLD"))
  }

  if (is.null(LDlist$msize)) {
    issues <- c(issues, paste0(kernel_id, ": LDlist$msize is required"))
  }

  issues
}

#' Validate IID kernel
#' @keywords internal
validate_IIDlist <- function(IIDlist, kernel_id = "IID") {
  issues <- character(0)

  if (is.null(IIDlist$ids)) {
    issues <- c(issues, paste0(kernel_id, ": ids is required for IID kernel"))
  } else if (!(is.character(IIDlist$ids) || is.numeric(IIDlist$ids))) {
    issues <- c(issues, paste0(kernel_id, ": ids must be character or numeric"))
  }

  issues
}

# ------------------------------------------------------------------
# Model-level validation
# ------------------------------------------------------------------

#' Validate qgtools model specification
#' @keywords internal
validate_qg_model <- function(models, vcov = NULL, Klist = NULL) {
  issues <- character(0)

  issues <- c(issues, validate_vcov_schema(vcov))
  issues <- c(issues, validate_vcov_structure(vcov))
  issues <- c(issues, validate_vcov_trait_structure(vcov))
  issues <- c(issues, validate_vcov_sets(vcov, Klist))

  # No vcov => fixed-effects-only model
  if (is.null(vcov)) return(issues)

  # Collect random-effect indices from models
  model_indices <- unique(unlist(lapply(models, function(m) {
    if (length(m$random) == 0) return(character(0))
    vapply(m$random, `[[`, character(1), "index")
  })))

  for (nm in names(vcov)) {
    v <- vcov[[nm]]

    # Structural checks
    for (fld in c("index", "kernel", "traits")) {
      if (is.null(v[[fld]])) {
        issues <- c(issues, paste0("vcov$", nm, ": missing field '", fld, "'"))
        next
      }
    }

    if (!v$index %in% model_indices) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": index '", v$index, "' not found in model random effects")
      )
    }

    if (is.null(Klist) || is.null(Klist[[v$kernel]])) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": kernel '", v$kernel, "' not found in Klist")
      )
      next
    }

    issues <- c(
      issues,
      validate_kernel(Klist[[v$kernel]], v$kernel)
    )

    model_traits <- names(models)
    if (!all(v$traits %in% model_traits)) {
      issues <- c(
        issues,
        paste0(
          "vcov$", nm, ": unknown trait(s): ",
          paste(setdiff(v$traits, model_traits), collapse = ", ")
        )
      )
    }
  }

  # Check ID consistency across kernels sharing index
  by_index <- split(vcov, vapply(vcov, `[[`, character(1), "index"))

  for (idx in names(by_index)) {
    kernels <- vapply(by_index[[idx]], `[[`, character(1), "kernel")
    ids <- lapply(kernels, function(k) {
      K <- Klist[[k]]
      if (!is.null(K$ids))  return(K$ids)
      if (!is.null(K$idsG)) return(K$idsG)
      if (!is.null(K$idsP)) return(K$idsP)
      NULL
    })

    ids <- ids[!vapply(ids, is.null, logical(1))]
    if (length(ids) > 1 && !all(vapply(ids[-1], identical, logical(1), ids[[1]]))) {
      issues <- c(
        issues,
        paste0("Index '", idx, "': kernels have incompatible IDs")
      )
    }
  }

  # 4. vcov structure validation
  issues <- c(
    issues,
    validate_vcov_trait_structure(vcov),
    validate_vcov_sets_kernel(vcov, Klist)
  )


  unique(issues)
}

# ------------------------------------------------------------
# vcov schema validation
# ------------------------------------------------------------

validate_vcov_schema <- function(vcov) {
  issues <- character(0)

  if (!is.list(vcov)) {
    return("vcov must be a list")
  }

  required_fields <- c("index", "traits", "kernel", "structure")

  for (nm in names(vcov)) {
    v <- vcov[[nm]]

    if (!is.list(v)) {
      issues <- c(issues, paste0("vcov$", nm, " must be a list"))
      next
    }

    missing <- setdiff(required_fields, names(v))
    if (length(missing) > 0) {
      issues <- c(
        issues,
        paste0("vcov$", nm, " missing field(s): ", paste(missing, collapse = ", "))
      )
    }

    if (!is.character(v$index) || length(v$index) != 1) {
      issues <- c(issues, paste0("vcov$", nm, "$index must be a single character"))
    }

    if (!is.character(v$traits) || length(v$traits) < 1) {
      issues <- c(issues, paste0("vcov$", nm, "$traits must be character vector"))
    }

    if (!is.character(v$kernel) || length(v$kernel) != 1) {
      issues <- c(issues, paste0("vcov$", nm, "$kernel must be a single character"))
    }

    if (!is.character(v$structure) || length(v$structure) != 1) {
      issues <- c(issues, paste0("vcov$", nm, "$structure must be a single character"))
    }
  }

  issues
}

.allowed_vcov_structures <- c(
  "diagonal",
  "unstructured",
  "identity",
  "shared"
)

validate_vcov_structure <- function(vcov) {
  issues <- character(0)

  for (nm in names(vcov)) {
    s <- vcov[[nm]]$structure
    if (!s %in% .allowed_vcov_structures) {
      issues <- c(
        issues,
        paste0(
          "vcov$", nm, "$structure='", s,
          "' not in allowed: ",
          paste(.allowed_vcov_structures, collapse = ", ")
        )
      )
    }
  }

  issues
}

validate_vcov_trait_structure <- function(vcov) {
  issues <- character(0)

  for (nm in names(vcov)) {
    v <- vcov[[nm]]
    ntr <- length(v$traits)

    if (v$structure == "unstructured" && ntr < 2) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": unstructured requires >= 2 traits")
      )
    }
  }

  issues
}

validate_vcov_sets <- function(vcov, Klist) {
  issues <- character(0)

  for (nm in names(vcov)) {
    v <- vcov[[nm]]
    if (is.null(v$sets)) next

    kernel <- Klist[[v$kernel]]
    if (is.null(kernel)) next

    # kernel type inference
    is_marker_kernel <- !is.null(kernel$markers) || !is.null(kernel$rsids)

    if (!is_marker_kernel) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": sets only allowed for marker-based kernels (Glist/Mlist)")
      )
      next
    }

    if (!is.list(v$sets) || is.null(names(v$sets))) {
      issues <- c(
        issues,
        paste0("vcov$", nm, "$sets must be a *named* list")
      )
      next
    }

    markers <- kernel$markers %||% kernel$rsids

    for (s in names(v$sets)) {
      bad <- setdiff(v$sets[[s]], markers)
      if (length(bad) > 0) {
        issues <- c(
          issues,
          paste0("vcov$", nm, "$sets$", s,
                 " contains unknown markers: ",
                 paste(head(bad, 5), collapse = ", "))
        )
      }
    }
  }

  issues
}

#' Validate vcov trait-structure compatibility
#' @keywords internal
validate_vcov_trait_structure <- function(vcov) {
  issues <- character(0)

  for (nm in names(vcov)) {
    v <- vcov[[nm]]

    if (is.null(v$structure)) {
      issues <- c(issues, paste0("vcov$", nm, ": structure is required"))
      next
    }

    structure <- v$structure
    traits <- v$traits

    allowed <- c("diagonal", "unstructured", "identity", "scalar")
    if (!structure %in% allowed) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": invalid structure '", structure,
               "' (allowed: ", paste(allowed, collapse = ", "), ")")
      )
      next
    }

    ntraits <- length(traits)

    if (structure == "unstructured" && ntraits < 2) {
      issues <- c(
        issues,
        paste0("vcov$", nm, ": unstructured requires >= 2 traits")
      )
    }
  }

  unique(issues)
}

#' Validate vcov sets vs kernel compatibility
#' @keywords internal
validate_vcov_sets_kernel <- function(vcov, Klist) {
  issues <- character(0)

  for (nm in names(vcov)) {
    v <- vcov[[nm]]
    if (is.null(v$sets)) next

    kernel <- v$kernel
    K <- Klist[[kernel]]

    is_marker_kernel <-
      !is.null(K$markers) ||
      !is.null(K$rsids) ||
      !is.null(K$bedfiles) ||
      !is.null(K$fnBED)

    if (!is_marker_kernel) {
      issues <- c(
        issues,
        paste0(
          "vcov$", nm, ": sets specified but kernel '", kernel,
          "' is not marker-based"
        )
      )
    }
  }

  unique(issues)
}


