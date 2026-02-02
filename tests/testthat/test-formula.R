test_that("formula parser handles complex random-effect structures", {

  f <- function(x) as.formula(x, env = emptyenv())

  fms <- list(
    y1 = f("y1 ~ 1 + (1 | id)"),
    y2 = f("y2 ~ age + (time | id)"),
    y3 = f("y3 ~ age + (1 + time | id)"),
    y4 = f("y4 ~ age + (0 + time | id)"),
    y5 = f("y5 ~ age + (time1 + time2 | id)"),
    y6 = f("y6 ~ age + (poly(time,2) | id)"),
    y7 = f("y7 ~ age + (1 | id) + (time | id)"),
    y8 = f("y8 ~ age + (1 | id) + (time | herd)"),
    y9 = f("y9 ~ age + (1 + time | id) + (1 | herd)")
  )

  parsed <- qg_parse_formulas(fms)

  ## Basic structure
  expect_true(is.list(parsed))
  expect_equal(length(parsed), 9)

  ## Intercept logic
  expect_true(parsed$y2$random[[1]]$intercept)
  expect_false(parsed$y4$random[[1]]$intercept)

  ## Slopes
  expect_true("time" %in% parsed$y3$random[[1]]$slopes)
  expect_equal(parsed$y6$random[[1]]$slopes, "poly(time, 2)")

  ## Multiple random terms
  expect_equal(length(parsed$y7$random), 2)

  ## Multiple indices
  expect_equal(parsed$y8$random[[2]]$index, "herd")
  expect_equal(parsed$y9$random[[2]]$index, "herd")
})


test_that("vcov validation catches bad structure", {
  vcov <- list(
    bad = list(
      index = "id",
      traits = "bmi",
      kernel = "G",
      structure = "unstructured"
    )
  )

  expect_true(
    any(grepl("unstructured requires", validate_vcov_trait_structure(vcov)))
  )
})
