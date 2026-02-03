test_that("decision1S works for lower sided", {
  dec <- decision1S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = TRUE
  )

  expect_true(is_lower(dec))
  expect_false(is_upper(dec))

  expect_class(dec, c("decision1S", "decision1S_1sided", "function"))
  expect_class(lower(dec), c("decision1S_atomic", "function"))

  expect_snapshot(print(dec))
})

test_that("decision1S works for upper sided", {
  dec <- decision1S(
    qc = c(0, 0.8, 1.2),
    pc = c(0.05, 0.5, 0.05),
    lower.tail = FALSE
  )

  expect_false(is_lower(dec))
  expect_true(is_upper(dec))

  expect_class(dec, c("decision1S", "decision1S_1sided", "function"))
  expect_class(upper(dec), c("decision1S_atomic", "function"))

  expect_snapshot(print(dec))
})

## Mixed lower.tail usage also works
## This is relevant to calculate the probability of "interim" outcomes.
## Example:
## P(theta <= 0) > 0.05
## P(theta > theta_target) > 0.5

decMixed <- decision1S(
  qc = c(0, 0.8, 1.2),
  pc = c(0.05, 0.5, 0.05),
  lower.tail = c(TRUE, FALSE, TRUE)
)

test_that("print method works for mixed lower.tail", {
  expect_snapshot(print(decMixed))
})

test_that("two-sided decision function works", {
  flat_prior <- mixnorm(c(1, 0, 100), sigma = 10)
  expect_equal(decMixed(flat_prior), 0)

  lower_part <- lower(decMixed)
  upper_part <- upper(decMixed)
  expect_equal(lower_part(flat_prior), 1)
  expect_equal(upper_part(flat_prior), 0)

  dist <- decMixed(flat_prior, dist = TRUE)
  expect_snapshot_value(dist, style = "deparse")
})
