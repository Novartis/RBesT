test_that("decision1S_boundary works for normal outcome", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  dec <- decision1S(pc = 0.95, qc = 0.6, lower.tail = TRUE)
  result <- decision1S_boundary(prior, n = 50, decision = dec)
  expect_equal(result, -0.7359727, tolerance = 1e-5)

  dec <- decision1S(pc = 0.8, qc = 1.5, lower.tail = FALSE)
  result <- decision1S_boundary(prior, n = 50, decision = dec)
  expect_equal(result, 2.105779, tolerance = 1e-5)
})

test_that("Mixed lower.tail usage works for decision boundary calculation", {
  prior <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  dec_lower <- decision1S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  result_lower <- decision1S_boundary(
    prior,
    n = 50,
    decision = dec_lower
  )

  dec_upper <- decision1S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  result_upper <- decision1S_boundary(
    prior,
    n = 50,
    decision = dec_upper
  )
  assert_true(result_lower != result_upper)

  decMixed <- decision1S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  result <- decision1S_boundary(prior, 50, decMixed)
  expected <- c(lower_than = result_lower, higher_than = result_upper)
  expect_equal(result, expected)
})
