## Mixed lower.tail usage also works
## This is relevant to calculate the probability of "interim" outcomes.
## Example:
## P(theta1 - theta2 <= 0) > 0.05
## P(theta1 - theta2 > theta_target) > 0.5

decMixed <- decision2S(
  qc = c(0, 0.8, 1.2),
  pc = c(0.05, 0.5, 0.05),
  lower.tail = c(TRUE, FALSE, TRUE)
)

test_that("print method works for mixed lower.tail", {
  expect_snapshot(print(decMixed))
})

test_that("two-sided decision function works", {
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
