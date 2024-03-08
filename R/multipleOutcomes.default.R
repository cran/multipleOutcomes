#' @export
multipleOutcomes.default <- function(..., family, data, data_index = NULL, score_epsilon = 1e-6){

  formulas <- list(...)
  formulas[sapply(formulas, function(x){!inherits(x, 'formula')})] <- NULL

  stopifnot(all(family %in% c('gaussian', 'binomial')))

  if(length(family) == 1){
    family <- rep(family, length(formulas))
  }

  if(is.data.frame(data)){
    data <- list(data)
  }

  if(is.null(data_index)){
    if(length(data) == 1){
      data_index <- rep(1, length(formulas))
    }else{
      stop('data_index should be a integer vector of length ', length(formula))
    }
  }
  stopifnot(all(sapply(data, is.data.frame)))

  stopifnot(all(sapply(data_index, function(x){abs(x - round(x)) < .Machine$double.eps^0.5})))
  stopifnot(all(data_index <= length(data) & data_index >= 1))
  stopifnot(length(data_index) == length(formulas))

  inverseHessianMatrix <- function(family, design_matrix, model = NULL){
    n_sam <- nrow(design_matrix)
    n_par <- ncol(design_matrix)
    hess <- matrix(0, n_par, n_par)
    if(family == 'gaussian'){
      hess <- -t(design_matrix / n_sam) %*% design_matrix
    }else{
      hess <- t(design_matrix / n_sam) %*% (design_matrix * (-fitted(model) * (1 - fitted(model))))
    }

    ret <- solve(hess)
    attr(ret, 'n') <- n_sam

    invisible(ret)
  }

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

  IDMapping <- function(df, i, var_name = NULL){
    if(i == 1){
      ids <- 1:df[1]
    }else{
      ids <- (sum(df[1:(i - 1)]) + 1):sum(df[1:i])
    }

    if(!is.null(var_name) && length(var_name) == length(ids)){
      names(ids) <- var_name
    }

    ids
  }

  n <- sapply(data, nrow)

  n_models <- length(formulas)
  models <- vector('list', n_models)
  inv_hess <- vector('list', n_models)
  df <- rep(NA, n_models)
  bet <- NULL
  vars <- vector('list', n_models)

  for(i in 1:n_models){
    models[[i]] <- glm(formulas[[i]], family[i], data[[data_index[i]]])

    df[i] <- coef(models[[i]]) %>% length
    bet <- c(bet, coef(models[[i]]) %>% unname)
    vars[[i]] <- coef(models[[i]]) %>% names
  }

  n_par <- length(bet)
  mcov <- matrix(0, n_par, n_par)
  n_shared_sample_sizes <- matrix(0, n_models, n_models)
  id_map <- NULL

  for(i in 1:n_models){
    n_shared_sample_sizes[i, i] <- nrow(data[[data_index[i]]])
    dm1 <- model.matrix(formulas[[i]], model.frame(formulas[[i]], data[[data_index[i]]]))

    score1 <- resid(models[[i]], type = 'response') * dm1
    stopifnot((colSums(score1) / nrow(score1)) %>% abs %>% max < score_epsilon)

    inv_hess[[i]] <- inverseHessianMatrix(family[i], dm1, models[[i]])

    id1 <- IDMapping(df, i, vars[[i]])
    id_map[[i]] <- id1
    for(j in i:length(formulas)){
      n_shared_sample_sizes[i, j] <- intersect(rownames(data[[data_index[[i]]]]), rownames(data[[data_index[[j]]]])) %>% length()
      n_shared_sample_sizes[j, i] <- n_shared_sample_sizes[i, j]

      dm2 <- model.matrix(formulas[[j]], model.frame(formulas[[j]], data[[data_index[j]]]))

      score2 <- resid(models[[j]], type = 'response') * dm2
      stopifnot((colSums(score2) / nrow(score2)) %>% abs %>% max < score_epsilon)
      if(is.null(inv_hess[[j]])){
        inv_hess[[j]] <- inverseHessianMatrix(family[j], dm2, models[[j]])
      }

      id2 <- IDMapping(df, j)
      info <- FisherInformation(score1, score2)
      mcov[id1, id2] <-
        inv_hess[[i]] %*%
        info %*%
        inv_hess[[j]] *
        (attr(info, 'n') / nrow(score1) / nrow(score2))
      mcov[id2, id1] <- t(mcov[id1, id2])

      rm(dm2, score2)
    }

    models[i] <- list(NULL)
    inv_hess[i] <- list(NULL)

  }

  fit <-
    list(
      coefficients = bet,
      mcov = mcov,
      id_map = id_map,
      n_shared_sample_sizes = n_shared_sample_sizes
    )

  fit$call <- match.call()
  class(fit) <- 'multipleOutcomes'
  invisible(fit)

}
