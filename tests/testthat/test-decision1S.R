## Mixed lower.tail usage also works
## This is relevant to calculate the probability of "interim" outcomes.
## Example:
## P(theta <= 0) > 0.05
## P(theta > theta_target) > 0.5

theta_target <- 0.8
decMixed <- decision1S(
  qc = c(0, theta_target),
  pc = c(0.05, 0.5),
  lower.tail = c(TRUE, FALSE)
)

test_that("print method works for mixed lower.tail", {
  expect_snapshot(print(decMixed))
})

test_that("length method works", {
  expect_equal(length(decMixed), 2L)
})

test_that("subsetting works", {
  dec1 <- decMixed[1]
  expect_equal(length(dec1), 1L)
  expect_identical(attr(dec1, "lower.tail"), TRUE)
  
  dec2 <- decMixed[2]
  expect_equal(length(dec2), 1L)
  expect_identical(attr(dec2, "lower.tail"), FALSE)
})
