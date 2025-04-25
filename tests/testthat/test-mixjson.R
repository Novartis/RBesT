test_that("write_mix_json and read_mix_json work correctly", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  write_mix_json(nm, tempfile, pretty = TRUE, digits = 1)

  # Check that the file was created
  expect_true(file.exists(tempfile))

  # Deserialize the JSON file back into a mixture object
  mix <- read_mix_json(tempfile)

  # Check that the deserialized object matches the original
  expect_true(all(nm == mix))
})

test_that("write_mix_json handles EM objects correctly", {
  # Create a mixture object with EM attributes
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  class(nm) <- c("EM", class(nm)) # Simulate an EM object

  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  expect_message(
    write_mix_json(nm, tempfile, pretty = TRUE, digits = 1),
    "Dropping EM information from mixture object before serialization."
  )

  # Deserialize the JSON file back into a mixture object
  mix <- read_mix_json(tempfile)

  # Check that the deserialized object matches the original (without EM attributes)
  expect_true(all(unclass(nm) == mix))
})

test_that("read_mix_json rescaling works correctly", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)
  nm[1,1] <- 0.15
  
  # Use withr to create a temporary file
  tempfile <- withr::local_tempfile(fileext = ".json")

  # Serialize the mixture object to JSON
  write_mix_json(nm, tempfile, pretty = TRUE, digits = 2)

  # Deserialize the JSON file with rescaling
  mixrescaled <- read_mix_json(tempfile, rescale = TRUE)

  # Check that the weights sum to 1 after rescaling
  expect_equal(sum(mixrescaled[1, ]), 1)

  # Deserialize the JSON file without rescaling
  mix_no_rescale <- read_mix_json(tempfile, rescale = FALSE)

  # Check that the weights do not sum to 1 without rescaling
  expect_false(sum(mix_no_rescale[1, ]) == 1)
})

test_that("write_mix_json warns about missing digits argument", {
  # Create a mixture object
  nm <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

  # Use withr to create a temporary file
  temp_file <- withr::local_tempfile(fileext = ".json")

  # Expect a warning when digits argument is not provided
  expect_warning(
    write_mix_json(nm, temp_file, pretty = TRUE),
    "JSON serialization by default restricts number of digits"
  )
})
