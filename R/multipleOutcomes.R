#' Fitting Regression Models for Multiple Outcomes and Returning the Matrix of
#' Covariance
#' @description
#' `multipleOutcomes` can fit different types of models for multiple outcomes 
#' simultaneously and return model parameters and variance-covariance matrix 
#' for further analysis. 
#' 
#' @importFrom stats coef cov fitted formula glm model.frame model.matrix pchisq
#' printCoefmat resid vcov cov2cor
#' @importFrom momentfit gmm4
#' @importFrom stringr str_glue
#' @importFrom dplyr %>% mutate
#' @importFrom survival coxph coxph.detail
#'
#' @param ... formulas of models to be fitted, or moment functions for gmm. 
#' @param family a character vector of families to be used in the models.
#' Currently only `gaussian`, `binomial`, `coxph` and `gmm` are supported.
#' `long` for longitudinal data may be supported in the future.
#' `family` can be of length 1 if all models are fitted in the
#' same family; otherwise family should be specified for each of the models in `...`.
#' @param data a data frame if all models are fitted on the same dataset;
#' otherwise a list of data frames for fitting models in `...`. Note that a
#' dataset can be used to fit multiple models, thus, `length(data)` is unnecessary
#' to be equal to the number of models in `...`. The row names in a data frame
#' are treated as sample IDs. Consequently, for any two records in different
#' data frames that correspond to the same sample, their row names should be
#' consistent.
#' @param data_index `NULL` if `data` is a data frame; otherwise, a vector in
#' integer specifying mapping a model in `...` to a data frame in `data` (a list).
#' @param score_epsilon whatever.
#'
#' @return It returns an object of class "multipleOutcomes", which is a list
#' containing the following components:
#' \tabular{ll}{
#' \code{coefficients} \tab an unnamed vector of coefficients of all fitted models.
#' Use `id_map` for variable mapping. \cr
#' \tab \cr
#' \code{mcov} \tab a unnamed matrix of covariance of `coefficients`. Use `id_map`
#' for variable mapping. \cr
#' \tab \cr
#' \code{id_map} \tab a list mapping the elements in `coefficients` and `mcov` to
#' variable names. \cr
#' \tab \cr
#' \code{n_shared_sample_sizes} \tab a matrix of shared sample sizes between datasets
#' being used to fit the models. \cr
#' \tab \cr
#' \code{call} \tab the matched call.\cr
#' }
#' @export
#'
#' @examples
#' ## More examples can be found in the vignettes.
#' library(mvtnorm)
#' genData <- function(seed = NULL){
#'
#'   set.seed(seed)
#'   n <- 400
#'   sigma <- matrix(c(1, .6, .6, 1), 2)
#'   x <- rmvnorm(n, sigma = sigma)
#'   gam <- c(.1, -.2)
#'   z <- rbinom(n, 1, plogis(1-1/(1+exp(-.5+x%*%gam+.1*rnorm(n)))))
#'
#'   bet <- c(-.2,.2)
#'   #y <- rbinom(n, 1, plogis(1-1/(1+exp(-.5+x%*%bet + .2*z-.3*rnorm(n)))))
#'   y <- -.5+x%*%bet + .2*z-.3*rnorm(n)
#'
#'   data.frame(y = y, z = z, x1 = x[, 1], x2 = x[, 2])
#'
#' }
#'
#' dat <- genData(123456)
#' dat1 <- head(dat,200)
#' dat2 <- tail(dat,200)
#' ## fitting four models simultaneously.
#' fit <-
#'   multipleOutcomes(
#'     y ~ z + x1 - 1,
#'     z ~ x1 + x2,
#'     z ~ x1 - 1,
#'     y ~ x2,
#'     ## z can be fitted with a linear or logistic regression
#'     family = c('gaussian', 'binomial', 'gaussian','gaussian'),
#'     data = list(dat1, dat2),
#'     ## each dataset is used to fit two models
#'     data_index = c(1, 1, 2, 2)
#'   )
#'
#'   ## unnamed coefficients of all model parameters
#'   coef(fit)
#'
#'   ## named coefficients of a specific model
#'   coef(fit, 2)
#'
#'   ## unnamed covariance matrix of all model parameters
#'   vcov(fit)
#'
#'   ## named covariance matrix of a specific model
#'   vcov(fit, 1)
#'
#'   ## summary of all parameter estimates
#'   summary(fit)
#'
#'   ## summary of parameters in a specific model
#'   summary(fit, 4)
#'
multipleOutcomes <- function(..., family, data, data_index = NULL, score_epsilon = 1e-6){

  UseMethod('multipleOutcomes')

}
