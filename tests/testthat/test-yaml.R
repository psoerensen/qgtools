test_that("YAML write/read roundtrip preserves model structure", {

  formulas <- list(
    bmi = as.formula(bmi ~ age + (1 | id) + (time | id), env = emptyenv()),
    height = as.formula(height ~ age + (1 | id), env = emptyenv())
  )

  vcov <- list(
    polygenic = list(
      index = "id",
      traits = c("bmi", "height"),
      kernel = "GRM",
      structure = "unstructured"
    )
  )

  Klist <- list(GRM = list(idsG = 1:10, fnG = "dummy"))

  tmp <- tempfile(fileext = ".yaml")

  obj_written <- write_qg_yaml(
    formulas = formulas,
    vcov = vcov,
    Klist = Klist,
    file = tmp
  )

  obj_read <- read_qg_yaml(tmp)

  expect_true(is.list(obj_read))
  expect_true("models" %in% names(obj_read))
  expect_true("vcov" %in% names(obj_read))
  expect_true("kernels" %in% names(obj_read))
})

test_that("Random regression terms survive YAML", {

  formulas <- list(
    y = as.formula(y ~ age + (1 | id) + (time | id), env = emptyenv())
  )

  tmp <- tempfile(fileext = ".yaml")

  write_qg_yaml(formulas, file = tmp)
  obj <- read_qg_yaml(tmp)

  rnd <- obj$models$y$random

  expect_equal(length(rnd), 2)

  expect_true(any(vapply(rnd, function(x) x$intercept, logical(1))))
  expect_true(any(vapply(rnd, function(x) "time" %in% x$slopes, logical(1))))
})

test_that("Klist kernels are written correctly", {

  formulas <- list(y = as.formula(y ~ 1 + (1 | id), env = emptyenv()))

  Klist <- list(
    GRM = list(idsG = 1:10, fnG = "dummy"),
    PED = list(idsP = 1:10, fnPED = "dummy")
  )

  tmp <- tempfile(fileext = ".yaml")

  write_qg_yaml(formulas, Klist = Klist, file = tmp)
  obj <- read_qg_yaml(tmp)

  expect_equal(sort(obj$kernels), c("GRM", "PED"))
})
