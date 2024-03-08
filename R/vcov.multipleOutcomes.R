#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' @param object an object returned by `multipleOutcomes()`.
#' @param model_index `NULL` if displaying covariance matrix of all fitted models;
#' otherwise, an integer indicating the fitted model.
#' @param ... for debugging only
#'
#' @return a matrix of covariance of all estimates
#' @export
#'
vcov.multipleOutcomes <- function(object, model_index = NULL, ...){

  ret <- object$mcov

  if(!is.null(model_index)){
    id_map <- object$id_map
    if(model_index > length(id_map)){
      stop('Invalid model index.')
    }

    ret <- ret[id_map[[model_index]], id_map[[model_index]]]
    rownames(ret) <- names(id_map[[model_index]])
    colnames(ret) <- names(id_map[[model_index]])
  }

  ret

}
