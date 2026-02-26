test_that("decision2S works for lower sided", {
  decLower <- decision2S(
    qc = c(0, 0.8),
    pc = c(0.05, 0.5),
    lower.tail = TRUE
  )

  expect_true(has_lower(decLower))
  expect_false(has_upper(decLower))

  # Complementary attribute is still present for consistency:
  expect_function(attr(decLower, "upper"))

  expect_class(decLower, c("decision2S", "decision2S_1sided", "function"))
  expect_class(lower(decLower), c("decision2S_atomic", "function"))

  expect_snapshot(print(decLower))

  flat_prior1 <- mixnorm(c(1, 0, 100), sigma = 10)
  flat_prior2 <- mixnorm(c(1, 1, 20), sigma = 5)
  expect_equal(decLower(flat_prior1, flat_prior2), 1)

  lower_part <- lower(decLower)
  expect_equal(lower_part(flat_prior1, flat_prior2), 1)

  dist <- decLower(flat_prior1, flat_prior2, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")

  priorT <- mixnorm(c(1, Inf, 10), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -Inf, 10), sigma = 88, param = "mn")
  successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
  expect_equal(lower(successCrit)(priorT, priorP), 1)
  expect_equal(upper(successCrit)(priorT, priorP), 1)
  expect_equal(lower(successCrit)(priorP, priorT), 1)
  expect_equal(upper(successCrit)(priorP, priorT), 0)
})

test_that("decision2S works for upper sided", {
  decUpper <- decision2S(
    qc = c(1.2, 2),
    pc = c(0.05, 0.5),
    lower.tail = FALSE
  )

  expect_false(has_lower(decUpper))
  expect_true(has_upper(decUpper))

  # Complementary attribute is still present for consistency:
  expect_function(attr(decUpper, "lower"))

  expect_class(decUpper, c("decision2S", "decision2S_1sided", "function"))
  expect_class(upper(decUpper), c("decision2S_atomic", "function"))

  expect_snapshot(print(decUpper))

  flat_prior1 <- mixnorm(c(1, 0, 100), sigma = 10)
  flat_prior2 <- mixnorm(c(1, 1, 20), sigma = 5)
  expect_equal(decUpper(flat_prior1, flat_prior2), 0)

  upper_part <- upper(decUpper)
  expect_equal(upper_part(flat_prior1, flat_prior2), 0)

  dist <- decUpper(flat_prior1, flat_prior2, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")

  priorT <- mixnorm(c(1, Inf, 10), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -Inf, 10), sigma = 88, param = "mn")
  successCrit <- decision2S(c(0.95, 0.5), c(0, 50), TRUE)
  expect_equal(lower(successCrit)(priorT, priorP), 0)
  expect_equal(upper(successCrit)(priorT, priorP), 1)
  expect_equal(lower(successCrit)(priorP, priorT), 1)
  expect_equal(upper(successCrit)(priorP, priorT), 1)
})

test_that("decision2S works for two sided", {
  decMixed <- decision2S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = c(TRUE, FALSE, TRUE)
  )

  expect_true(has_lower(decMixed))
  expect_true(has_upper(decMixed))
  expect_class(decMixed, c("decision2S", "decision2S_2sided", "function"))
  expect_class(lower(decMixed), c("decision2S_atomic", "function"))
  expect_class(upper(decMixed), c("decision2S_atomic", "function"))

  expect_snapshot(print(decMixed))

  flat_prior1 <- mixnorm(c(1, 0, 100), sigma = 10)
  flat_prior2 <- mixnorm(c(1, 1, 20), sigma = 5)
  expect_equal(decMixed(flat_prior1, flat_prior2), 0)

  lower_part <- lower(decMixed)
  upper_part <- upper(decMixed)
  expect_equal(lower_part(flat_prior1, flat_prior2), 1)
  expect_equal(upper_part(flat_prior1, flat_prior2), 0)

  dist <- decMixed(flat_prior1, flat_prior2, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")
})
