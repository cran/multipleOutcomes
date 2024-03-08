#' Object Summaries
#'
#' @param object an object returned by `multipleOutcomes()`.
#' @param model_index `NULL` if displaying summary of all fitted models; otherwise,
#' an integer indicating the fitted model.
#' @param ... for debugging only
#'
#' @return a list
#' @export
#'
summary.multipleOutcomes <- function(object, model_index = NULL, ...){

  se <- sqrt(diag(object$mcov))
  tval <- coef(object) / se

  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = pchisq(tval^2, df = 1, lower.tail = FALSE))
  colnames(TAB) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')

  id_map <- object$id_map
  if(!is.null(model_index)){
    if(model_index > length(id_map)){
      stop('Invalid model index.')
    }
    TAB <- TAB[id_map[[model_index]], , drop = FALSE]
    rownames(TAB) <- names(id_map[[model_index]])
  }else{
    rownames(TAB) <-
      lapply(
        1:length(id_map),
        function(mid){
          paste0('model_', mid, '_', names(id_map[[mid]]))
        }
      ) %>% unlist()
  }

  res <- list(call = object$call,
              coefficients = TAB)

  class(res) <- 'summary.multipleOutcomes'
  res

}
