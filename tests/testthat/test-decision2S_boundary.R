test_that("decision2S_boundary works for normal outcome", {
  priorT <- mixnorm(c(1, 0, 0.001), sigma = 88, param = "mn")
  priorP <- mixnorm(c(1, -49, 20), sigma = 88, param = "mn")
 
  successCrit <- decision2S(c(0.95, 0.5), c(0, 50), FALSE)
  futilityCrit <- decision2S(c(0.90), c(40), TRUE)

  successBoundary <- decision2S_boundary(priorP, priorT, 10, 20, successCrit)
  futilityBoundary <- decision2S_boundary(priorP, priorT, 10, 20, futilityCrit)

  gridVals <- -25:25 - 49

  successBounds <- successBoundary(gridVals)
  futilityBounds <- futilityBoundary(gridVals)
  expect_snapshot_value(successBounds, style = "deparse")
  expect_snapshot_value(futilityBounds, style = "deparse")

  # Now we define the criterion for the gray zone using mixed lower.tail.
  grayzoneCrit <- decision2S(
    c(0.95, 0.5, 0.9), 
    c(0, 50, 40), 
    c(TRUE, TRUE, FALSE)
  ) 
  grayzoneBoundary <- decision2S_boundary(priorP, priorT, 10, 20, grayzoneCrit)
  grayzoneBoundsLower <- grayzoneBoundary$lower_than(gridVals)
  grayzoneBoundsHigher <- grayzoneBoundary$higher_than(gridVals)

  expect_snapshot_value(grayzoneBoundsLower, style = "deparse")
  expect_snapshot_value(grayzoneBoundsHigher, style = "deparse")
  
  # In this case there is no gray zone:
  expect_false(any(grayzoneBoundsHigher < grayzoneBoundsLower))
})

test_that("Mixed lower.tail usage works for decision boundary calculation", {
  skip_on_cran()

  prior1 <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  prior2 <- mixnorm(rob = c(0.2, 1, 2), inf = c(0.8, 3, 2), sigma = 5)

  dec_lower <- decision2S(pc = 0.5, qc = 1.5, lower.tail = TRUE)
  boundary_fn_lower <- decision2S_boundary(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_lower
  )

  gridVals <- -25:25 - 49
  result_lower <- boundary_fn_lower(gridVals)

  dec_upper <- decision2S(pc = 0.6, qc = 0.5, lower.tail = FALSE)
  boundary_fn_upper <- decision2S_boundary(
    prior1,
    prior2,
    n1 = 50,
    n2 = 50,
    decision = dec_upper
  )
  result_upper <- boundary_fn_upper(gridVals)
  
  decMixed <- decision2S(
    qc = c(1.5, 0.5),
    pc = c(0.5, 0.6),
    lower.tail = c(TRUE, FALSE)
  )
  boundary_fn_mixed <- decision2S_boundary(prior1, prior2, 50, 50, decMixed)
  
  result_mixed_lower <- boundary_fn_mixed$lower_than(gridVals)
  result_mixed_upper <- boundary_fn_mixed$higher_than(gridVals)
  
  # For 2-sample, mixed lower.tail returns a list
  expect_equal(result_mixed_lower, result_lower)
  expect_equal(result_mixed_upper, result_upper)
})
