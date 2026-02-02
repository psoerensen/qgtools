test_that("sets not allowed for non-marker kernels", {
  vcov <- list(
    bad = list(
      index = "id",
      traits = c("bmi", "height"),
      kernel = "GRM",
      structure = "diagonal",
      sets = list(A = 1:10)
    )
  )

  Klist <- list(GRM = list(idsG = 1:10, fnG = "x"))

  expect_true(
    any(grepl("not marker-based", validate_vcov_sets_kernel(vcov, Klist)))
  )
})


test_that("invalid vcov structure is caught", {
  vcov <- list(
    bad = list(
      index = "id",
      traits = "bmi",
      kernel = "G",
      structure = "weird"
    )
  )

  expect_true(
    any(grepl("invalid structure", validate_vcov_trait_structure(vcov)))
  )
})
