
dfbeta <- function(formula, data){

  mph <- coxph(formula, data, model = TRUE)
  ret <- as.matrix(resid(mph, type = 'dfbeta'))
  attr(ret, 'score') <- resid(mph, type = 'score')
  imat <- coxph.detail(mph)$imat
  if('numeric' %in% class(imat)){
    attr(ret, 'info') <- sum(imat)
  }else{
    attr(ret, 'info') <- rowSums(imat, dims = 2)
  }

  attr(ret, 'SE') <- sqrt(diag(vcov(mph)))
  attr(ret, 'logHR') <- coef(mph)
  ret

}
