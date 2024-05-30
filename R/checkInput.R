
checkInput <- function(..., family, data, data_index){
  
  formulas <- list(...)
  for(i in seq_along(formulas)){
    if(!inherits(formulas[[i]], 'function') && !inherits(formulas[[i]], 'formula')){
      stop(str_glue('{formulas[[i]]} is neither a function nor formula.'))
    }
    
    if(inherits(formulas[[i]], 'function')){
      if(is.null(formals(formulas[[i]])$start)){
        stop('Any function specified in ... should have an argument "start" with a default value.')
      }
    }
  }
  
  if(!all(family %in% c('gaussian', 'binomial', 'coxph', 'gmm'))){
    stop('Only gaussian, binomial, coxph, and gmm are supported families.')
  }
  
  if(length(family) == 1){
    family <- rep(family, length(formulas))
  }
  
  if(length(family) != length(formulas)){
    stop('family should be a vector of length ', length(formulas))
  }
  
  if(!setequal(which(sapply(formulas, function(x) inherits(x, 'function'))), which(family %in% 'gmm'))){
    stop('family should be gmm if a function is specified. ')
  }
  
  if(is.data.frame(data)){
    data <- list(data)
  }
  
  if(is.null(data_index)){
    if(length(data) == 1){
      data_index <- rep(1, length(formulas))
    }
  }
  
  if(length(data_index) != length(formulas)){
    stop('data_index should be a integer vector of length ', length(formula))
  }
  
  stopifnot(all(sapply(data, is.data.frame)))
  
  stopifnot(all(sapply(data_index, function(x){abs(x - round(x)) < .Machine$double.eps^0.5})))
  stopifnot(all(data_index <= length(data) & data_index >= 1))
  
  list(formulas = formulas, family = family, data = data, data_index = data_index)
  
}
