y <- rbind(1:5, 6:10, 11:15)


# errors for data
expect_error(tojmoment(
  data = NULL,
  acov_order = 0,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = 1,
  acov_order = 0,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = 1:10,
  acov_order = 0,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = "x",
  acov_order = 0,
  acor_order = 1,
  R = 1000
))

# errors for acov_order
expect_error(tojmoment(
  data = y,
  acov_order = -1,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = 1:10,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = y,
  acor_order = 1,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = "x",
  acor_order = 1,
  R = 1000
))

# errors for acor_order
expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = -1,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = 1:10,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = y,
  R = 1000
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = "x",
  R = 1000
))

# errors for R
expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = 1,
  R = -1
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = 1,
  R = 1:10
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = 1,
  R = y
))

expect_error(tojmoment(
  data = y,
  acov_order = 0,
  acor_order = 1,
  R = "x"
))