# a matrix
y <- rbind(1:5, 6:10, 11:15)

#-------------------------------------------------------------------------------
# testing errors for `nemoment`
#-------------------------------------------------------------------------------

# errors for data
expect_error(nemoment(data = NULL, acov_order = 0, acor_order = 1, R = 1000))

expect_error(nemoment(data = 1, acov_order = 0, acor_order = 1, R = 1000))

expect_error(nemoment(data = 1:10, acov_order = 0, acor_order = 1, R = 1000))

expect_error(nemoment(data = "x", acov_order = 0, acor_order = 1, R = 1000))

# errors for acov_order
expect_error(nemoment(data = y, acov_order = -1, acor_order = 1, R = 1000))

expect_error(nemoment(data = y, acov_order = 1:10, acor_order = 1, R = 1000))

expect_error(nemoment(data = y, acov_order = y, acor_order = 1, R = 1000))

expect_error(nemoment(data = y, acov_order = "x", acor_order = 1, R = 1000))

# errors for acor_order
expect_error(nemoment(data = y, acov_order = 0, acor_order = -1, R = 1000))

expect_error(nemoment(data = y, acov_order = 0, acor_order = 1:10, R = 1000))

expect_error(nemoment(data = y, acov_order = 0, acor_order = y, R = 1000))

expect_error(nemoment(data = y, acov_order = 0, acor_order = "x", R = 1000))

# errors for R
expect_error(nemoment(data = y, acov_order = 0, acor_order = 1, R = -1))

expect_error(nemoment(data = y, acov_order = 0, acor_order = 1, R = 1:10))

expect_error(nemoment(data = y, acov_order = 0, acor_order = 1, R = y))

expect_error(nemoment(data = y, acov_order = 0, acor_order = 1, R = "x"))



#-------------------------------------------------------------------------------
# testing errors for `hpjmoment`
#-------------------------------------------------------------------------------

# errors for data
expect_error(hpjmoment(data = NULL, acov_order = 0, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = 1, acov_order = 0, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = 1:10, acov_order = 0, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = "x", acov_order = 0, acor_order = 1, R = 1000))

# errors for acov_order
expect_error(hpjmoment(data = y, acov_order = -1, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = y, acov_order = 1:10, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = y, acov_order = y, acor_order = 1, R = 1000))

expect_error(hpjmoment(data = y, acov_order = "x", acor_order = 1, R = 1000))

# errors for acor_order
expect_error(hpjmoment(data = y, acov_order = 0, acor_order = -1, R = 1000))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = 1:10, R = 1000))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = y, R = 1000))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = "x", R = 1000))

# errors for R
expect_error(hpjmoment(data = y, acov_order = 0, acor_order = 1, R = -1))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = 1, R = 1:10))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = 1, R = y))

expect_error(hpjmoment(data = y, acov_order = 0, acor_order = 1, R = "x"))


#-------------------------------------------------------------------------------
# testing errors for `neecdf`
#-------------------------------------------------------------------------------

# errors for data
expect_error(neecdf(data = NULL, acov_order = 0, acor_order = 1))

expect_error(neecdf(data = 1, acov_order = 0, acor_order = 1))

expect_error(neecdf(data = 1:10, acov_order = 0, acor_order = 1))

expect_error(neecdf(data = "x", acov_order = 0, acor_order = 1))

# errors for acov_order
expect_error(neecdf(data = y, acov_order = -1, acor_order = 1))

expect_error(neecdf(data = y, acov_order = 1:10, acor_order = 1))

expect_error(neecdf(data = y, acov_order = y, acor_order = 1))

expect_error(neecdf(data = y, acov_order = "x", acor_order = 1))

# errors for acor_order
expect_error(neecdf(data = y, acov_order = 0, acor_order = -1))

expect_error(neecdf(data = y, acov_order = 0, acor_order = 1:10))

expect_error(neecdf(data = y, acov_order = 0, acor_order = y))

expect_error(neecdf(data = y, acov_order = 0, acor_order = "x"))



#-------------------------------------------------------------------------------
# testing errors for `hpjecdf`
#-------------------------------------------------------------------------------

# errors for data
expect_error(hpjecdf(data = NULL, acov_order = 0, acor_order = 1))

expect_error(hpjecdf(data = 1, acov_order = 0, acor_order = 1))

expect_error(hpjecdf(data = 1:10, acov_order = 0, acor_order = 1))

expect_error(hpjecdf(data = "x", acov_order = 0, acor_order = 1))

# errors for acov_order
expect_error(hpjecdf(data = y, acov_order = -1, acor_order = 1))

expect_error(hpjecdf(data = y, acov_order = 1:10, acor_order = 1))

expect_error(hpjecdf(data = y, acov_order = y, acor_order = 1))

expect_error(hpjecdf(data = y, acov_order = "x", acor_order = 1))

# errors for acor_order
expect_error(hpjecdf(data = y, acov_order = 0, acor_order = -1))

expect_error(hpjecdf(data = y, acov_order = 0, acor_order = 1:10))

expect_error(hpjecdf(data = y, acov_order = 0, acor_order = y))

expect_error(hpjecdf(data = y, acov_order = 0, acor_order = "x"))



#-------------------------------------------------------------------------------
# testing errors for `nekd`
#-------------------------------------------------------------------------------

# errors for data
expect_error(nekd(data = NULL, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = 1, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = 1:10, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = "x", acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for acov_order
expect_error(nekd(data = y, acov_order = -1, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 1:10, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = y, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = "x", acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for acor_order
expect_error(nekd(data = y, acov_order = 0, acor_order = -1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1:10, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = y, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = "x", mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for mean_bw
expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = -1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 1:10, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = y, acov_bw = 0.1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = "x", acov_bw = 0.1, acor_bw = 0.1))


# errors for acov_bw
expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = -1, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 1:10, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = y, acor_bw = 0.1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = "x", acor_bw = 0.1))

# errors for acor_bw
expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = -1))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 1:10))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = y))

expect_error(nekd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = "x"))



#-------------------------------------------------------------------------------
# testing errors for `hpjkd`
#-------------------------------------------------------------------------------

# errors for data
expect_error(hpjkd(data = NULL, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = 1, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = 1:10, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = "x", acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for acov_order
expect_error(hpjkd(data = y, acov_order = -1, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 1:10, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = y, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = "x", acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for acor_order
expect_error(hpjkd(data = y, acov_order = 0, acor_order = -1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1:10, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = y, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = "x", mean_bw = 0.1, acov_bw = 0.1, acor_bw = 0.1))

# errors for mean_bw
expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = -1, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 1:10, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = y, acov_bw = 0.1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = "x", acov_bw = 0.1, acor_bw = 0.1))


# errors for acov_bw
expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = -1, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 1:10, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = y, acor_bw = 0.1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = "x", acor_bw = 0.1))

# errors for acor_bw
expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = -1))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = 1:10))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = y))

expect_error(hpjkd(data = y, acov_order = 0, acor_order = 1, mean_bw = 0.1, acov_bw = 0.1, acor_bw = "x"))
