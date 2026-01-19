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

test_that("length method works", {
  expect_equal(length(decMixed), 3L)
})

test_that("subsetting works", {
  dec1 <- decMixed[1]
  expect_equal(length(dec1), 1L)
  expect_identical(attr(dec1, "lower.tail"), TRUE)

  dec2 <- decMixed[2]
  expect_equal(length(dec2), 1L)
  expect_identical(attr(dec2, "lower.tail"), FALSE)

  dec12 <- decMixed[1:2]
  expect_equal(length(dec12), 2L)
  expect_identical(attr(dec12, "lower.tail"), c(TRUE, FALSE))

  # This case is important: both criterions have same lower.tail,
  # therefore it should be simplified to a single logical.
  dec13 <- decMixed[c(1, 3)]
  expect_equal(length(dec13), 2L)
  expect_identical(attr(dec13, "lower.tail"), TRUE)
})
