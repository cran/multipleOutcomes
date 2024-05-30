#' Prognostic Variables Assisted Treatment Effect Detection
#' @description
#' `pated` is a wrapper function of `multipleOutcomes` for testing treatment effect 
#' in randomized clinical trials. It assumes that prognostic variables are fully 
#' randomized. This assumption can help enhancing statistical power of conventional 
#' approaches in detecting the treatment effect. Specifically, the sensitivity 
#' of the conventional models specified in `...` are improved by `pated`. 
#' 
#' @param ... formulas of models to be fitted, or moment functions for gmm. 
#' @param family a character vector of families to be used in the models.
#' Currently only `gaussian`, `binomial`, `coxph` and `gmm` are supported.
#' `long` for longitudinal data may be supported in the future.
#' `family` can be of length 1 if all models are fitted in thesame family; 
#' otherwise family should be specified for each of the models in `...`.
#' @param data a data frame if all models are fitted on the same dataset;
#' otherwise a list of data frames for fitting models in `...`. Note that a
#' dataset can be used to fit multiple models, thus, `length(data)` is unnecessary
#' to be equal to the number of models in `...`. The row names in a data frame
#' are treated as sample IDs. Consequently, for any two records in different
#' data frames that correspond to the same sample, their row names should be
#' consistent.
#'
#' @return a data frame of testing results.
#' @export
#'
#' @examples
#' # see vignettes
pated <- function(..., family, data){
  fit <- multipleOutcomes(..., family = family, data = data, data_index = NULL)

  parseFormula <- function(...){
    formulas <- list(...)
    formulas[sapply(formulas, function(x){!inherits(x, 'formula')})] <- NULL
    outcome <- NULL
    arm <- NULL
    for(fstr in formulas){
      outcome <- c(outcome, trimws(sub("\\s*~.*", "", deparse1(fstr))))
      arm <- c(arm, trimws(sub("^.*?~", "", deparse1(fstr))))
    }
    data.frame(outcome, arm)
  }

  fml <- parseFormula(...)

  id <- sapply(1:length(fit$id_map), function(i) fit$id_map[[i]][fml$arm[i]])
  id1 <- id[1]
  id2 <- id[-1]
  Delta <- coef(fit)[id1]
  delta <- coef(fit)[id2]
  mcov <- vcov(fit)
  opt_c <- as.vector(- solve(mcov[id2, id2, drop = FALSE]) %*% mcov[id2, id1, drop = FALSE])

  estimate <- Delta + sum(delta * opt_c)
  stderr <- sqrt(
    mcov[id1, id1] +
      2 * as.vector(t(opt_c) %*% mcov[id2, id1, drop = FALSE]) +
      as.vector(t(opt_c) %*% mcov[id2, id2, drop = FALSE] %*% opt_c)
  )

  pvalue <- pchisq((estimate / stderr)^2, df = 1, lower.tail = FALSE)
  
  treatment <-
    data.frame(
      term = fml$outcome[1],
      estimate = estimate,
      stderr = stderr,
      pvalue = pvalue,
      method = 'PATED',
      corr = NA
    )

  nonconfounder <-
    data.frame(
      term = fml$outcome,
      estimate = coef(fit)[id],
      stderr = sqrt(diag(vcov(fit))[id]),
      pvalue = pchisq((coef(fit)[id] / sqrt(diag(vcov(fit)))[id])^2, df = 1, lower.tail = FALSE),
      method = c('Standard', rep('Prognostic', nrow(fml) - 1)),
      corr = cov2cor(mcov)[id, id1]
    )

  ret <- rbind(treatment, nonconfounder)
  ret$method <- paste0(ret$method, ifelse(ret$method == 'Prognostic' & ret$pvalue < .05, '*', ''))
  attr(ret, 'Rel. Eff.') <- mcov[id1, id1] / stderr^2
  ret

}
