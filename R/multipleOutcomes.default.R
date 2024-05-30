#' @export
multipleOutcomes.default <- function(..., family, data, data_index = NULL, score_epsilon = 1e-6){

  input <- checkInput(..., family = family, data = data, data_index = data_index)
  formulas <- input$formulas
  family <- input$family
  data <- input$data
  data_index <- input$data_index
  
  inverseHessianMatrixGlm <- function(family, design_matrix, model = NULL){
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
  
  inverseJacobianMatrixGmm <- function(func, u, data){
    func2 <- function(eval_at, data){
      apply(func(eval_at, data), 2, sum)
    }
    
    n_sam <- nrow(data)
    jac <- numDeriv::jacobian(func2, u, data = data) / n_sam
    ret <- solve(jac)
    attr(ret, 'n') <- n_sam
    
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
    
    if(family[i] %in% c('binomial', 'gaussian')){
      models[[i]] <- glm(formulas[[i]], family[i], data[[data_index[i]]])
    }else{
      if(family[i] %in% 'coxph'){
        models[[i]] <- coxph(formulas[[i]], data[[data_index[i]]], model = TRUE)
      }else{ # gmm
        models[[i]] <- gmm4(formulas[[i]], data[[data_index[i]]], theta0 = eval(formals(formulas[[i]])$start))
        attr(models[[i]], 'score') <- formulas[[i]](models[[i]]@theta, data[[data_index[i]]])
      }
    }
    
    if(family[i] %in% 'gmm'){
      cf <- models[[i]]@theta
    }else{
      cf <- coef(models[[i]])
    }

    df[i] <- length(cf)
    bet <- c(bet, unname(cf))
    vars[[i]] <- names(cf)
  }

  n_par <- length(bet)
  mcov <- matrix(0, n_par, n_par)
  n_shared_sample_sizes <- matrix(0, n_models, n_models)
  id_map <- NULL

  for(i in 1:n_models){
    n_shared_sample_sizes[i, i] <- nrow(data[[data_index[i]]])

    if(family[i] %in% c('binomial', 'gaussian')){
      dm1 <- model.matrix(formulas[[i]], model.frame(formulas[[i]], data[[data_index[i]]]))
      score1 <- resid(models[[i]], type = 'response') * dm1
      inv_hess[[i]] <- inverseHessianMatrixGlm(family[i], dm1, models[[i]])
      rm(dm1)
    }else{
      if(family[i] %in% 'coxph'){
        score1 <- -as.matrix(resid(models[[i]], type = 'score'))
        imat <- coxph.detail(models[[i]])$imat
        if('numeric' %in% class(imat)){
          inv_hess[[i]] <- solve(sum(imat) / nrow(score1))
        }else{
          inv_hess[[i]] <- solve(rowSums(imat, dims = 2) / nrow(score1))
        }
      }else{ # gmm
        score1 <- attr(models[[i]], 'score')
        inv_hess[[i]] <- inverseJacobianMatrixGmm(formulas[[i]], models[[i]]@theta, data[[data_index[i]]])
      }
    }

    stopifnot((colSums(score1) / nrow(score1)) %>% abs %>% max < score_epsilon)

    id1 <- IDMapping(df, i, vars[[i]])
    id_map[[i]] <- id1
    for(j in i:length(formulas)){
      n_shared_sample_sizes[i, j] <- intersect(rownames(data[[data_index[[i]]]]), rownames(data[[data_index[[j]]]])) %>% length()
      n_shared_sample_sizes[j, i] <- n_shared_sample_sizes[i, j]

      if(family[j] %in% c('binomial', 'gaussian')){
        dm2 <- model.matrix(formulas[[j]], model.frame(formulas[[j]], data[[data_index[j]]]))
        score2 <- resid(models[[j]], type = 'response') * dm2
        if(is.null(inv_hess[[j]])){
          inv_hess[[j]] <- inverseHessianMatrixGlm(family[j], dm2, models[[j]])
        }
        rm(dm2)
      }else{
        if(family[i] %in% 'coxph'){
          score2 <- -as.matrix(resid(models[[j]], type = 'score'))
          if(is.null(inv_hess[[j]])){
            imat <- coxph.detail(models[[j]])$imat
            if('numeric' %in% class(imat)){
              inv_hess[[j]] <- solve(sum(imat) / nrow(score2))
            }else{
              inv_hess[[j]] <- solve(rowSums(imat, dims = 2) / nrow(score2))
            }
          }
        }else{ # gmm
          score2 <- attr(models[[j]], 'score')
          inv_hess[[j]] <- inverseJacobianMatrixGmm(formulas[[j]], models[[j]]@theta, data[[data_index[j]]])
        }
      }

      stopifnot((colSums(score2) / nrow(score2)) %>% abs %>% max < score_epsilon)

      id2 <- IDMapping(df, j)
      info <- FisherInformation(score1, score2)
      mcov[id1, id2] <-
        inv_hess[[i]] %*%
        info %*%
        t(inv_hess[[j]]) *
        (attr(info, 'n') / nrow(score1) / nrow(score2))
      mcov[id2, id1] <- t(mcov[id1, id2])

      rm(score2)
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
