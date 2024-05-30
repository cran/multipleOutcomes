
FisherInformation <- function(score1, score2){

  shared_ids <- intersect(rownames(score1), rownames(score2))
  if(length(shared_ids) == 0){
    ret <- matrix(0, nrow = ncol(score1), ncol = ncol(score2))
  }else{
    ret <- cov(score1[shared_ids, , drop = FALSE],
               score2[shared_ids, , drop = FALSE])
  }

  attr(ret, 'n') <- length(shared_ids)

  invisible(ret)
}
