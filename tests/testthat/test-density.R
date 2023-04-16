y <- rbind(1:5, 6:10, 11:15)

# errors for data
expect_error(nekd(
  data = NULL,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = 1,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = 1:10,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = "x",
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

# errors for acov_order
expect_error(nekd(
  data = y,
  acov_order = -1,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 1:10,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = y,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = "x",
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

# errors for acor_order
expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = -1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1:10,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = y,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = "x",
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

# errors for mean_bw
expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = -1,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 1:10,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = y,
  acov_bw = 0.1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = "x",
  acov_bw = 0.1,
  acor_bw = 0.1
))


# errors for acov_bw
expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = -1,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 1:10,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = y,
  acor_bw = 0.1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = "x",
  acor_bw = 0.1
))

# errors for acor_bw
expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = -1
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = 1:10
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = y
))

expect_error(nekd(
  data = y,
  acov_order = 0,
  acor_order = 1,
  mean_bw = 0.1,
  acov_bw = 0.1,
  acor_bw = "x"
))
