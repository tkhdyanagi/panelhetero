y <- rbind(1:5, 6:10, 11:15)

# errors for data
expect_error(hpjecdf(
  data = NULL,
  acov_order = 0,
  acor_order = 1
))

expect_error(hpjecdf(
  data = 1,
  acov_order = 0,
  acor_order = 1
))

expect_error(hpjecdf(
  data = 1:10,
  acov_order = 0,
  acor_order = 1
))

expect_error(hpjecdf(
  data = "x",
  acov_order = 0,
  acor_order = 1
))

# errors for acov_order
expect_error(hpjecdf(
  data = y,
  acov_order = -1,
  acor_order = 1
))

expect_error(hpjecdf(
  data = y,
  acov_order = 1:10,
  acor_order = 1
))

expect_error(hpjecdf(
  data = y,
  acov_order = y,
  acor_order = 1
))

expect_error(hpjecdf(
  data = y,
  acov_order = "x",
  acor_order = 1
))

# errors for acor_order
expect_error(hpjecdf(
  data = y,
  acov_order = 0,
  acor_order = -1
))

expect_error(hpjecdf(
  data = y,
  acov_order = 0,
  acor_order = 1:10
))

expect_error(hpjecdf(
  data = y,
  acov_order = 0,
  acor_order = y
))

expect_error(hpjecdf(
  data = y,
  acov_order = 0,
  acor_order = "x"
))
