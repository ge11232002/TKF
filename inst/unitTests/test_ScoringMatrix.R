test_PAMn <- function(){
  data(GONNET)
  PAM250R <- PAMn(GONNET, 250, method="R")
  PAM250C <- PAMn(GONNET, 250, method="C")
  checkEqualsNumeric(PAM250R, PAM250C)
}

