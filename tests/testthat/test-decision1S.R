test_that("decision1S works for lower sided", {
  dec <- decision1S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = TRUE
  )

  expect_true(has_lower(dec))
  expect_false(has_upper(dec))

  expect_class(dec, c("decision1S", "decision1S_1sided", "function"))
  expect_class(lower(dec), c("decision1S_atomic", "function"))

  expect_snapshot(print(dec))

  flat_prior <- mixnorm(c(1, 0, 100), sigma = 10)
  expect_equal(dec(flat_prior), 1)

  dist <- dec(flat_prior, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")
})

test_that("decision1S works for upper sided", {
  dec <- decision1S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = FALSE
  )

  expect_false(has_lower(dec))
  expect_true(has_upper(dec))

  expect_class(dec, c("decision1S", "decision1S_1sided", "function"))
  expect_class(upper(dec), c("decision1S_atomic", "function"))

  expect_snapshot(print(dec))

  flat_prior <- mixnorm(c(1, 0, 100), sigma = 10)
  expect_equal(dec(flat_prior), 0)

  dist <- dec(flat_prior, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")
})

test_that("decision1S works for two sided", {
  decMixed <- decision1S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = c(TRUE, FALSE, TRUE)
  )

  expect_true(has_lower(decMixed))
  expect_true(has_upper(decMixed))
  expect_class(decMixed, c("decision1S", "decision1S_2sided", "function"))
  expect_class(lower(decMixed), c("decision1S_atomic", "function"))
  expect_class(upper(decMixed), c("decision1S_atomic", "function"))

  expect_snapshot(print(decMixed))

  flat_prior <- mixnorm(c(1, 0, 100), sigma = 10)
  expect_equal(decMixed(flat_prior), 0)

  lower_part <- lower(decMixed)
  upper_part <- upper(decMixed)
  expect_equal(lower_part(flat_prior), 1)
  expect_equal(upper_part(flat_prior), 0)

  dist <- decMixed(flat_prior, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")
})
