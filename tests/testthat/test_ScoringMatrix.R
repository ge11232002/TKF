test_that("test_PAMn", {
    PAM250R <- PAMn(GONNET, 250, method = "R")
    PAM250C <- PAMn(GONNET, 250, method = "C")
    expect_equal(PAM250C, PAM250R)
})

