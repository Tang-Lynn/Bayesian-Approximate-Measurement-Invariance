# CJ: Please read the tidyverse style guide for R code;
# CJ: for example, there should always be a space between comma's and mathematical operators,
# CJ: variable names should be snake_case not camelCase, and function names should be verbs when possible
# CJ: https://style.tidyverse.org/syntax.html
# testing
library(lavaan)  # CFA model
library(bain) # example from this package
library(psych)    # Descriptives
library(mvtnorm) # compute Bayes factors
library(Matrix) # matrix
# Arguments
# 1. x, the output of the two-group CFA model with
#    the factor variance = 1 and the factor mean = 0.
# CJ: I would suggest calling this argument x,
# CJ: because even though at this moment the function only accepts lavaan
# CJ: objects, that will likely change in the future.
# CJ: Use S3 method dispatch to call the correct function (see my suggestion below).
# 2. type, a string that specifies the type of the test.
#    If type = 'metric', metric invariance is tested;
#    If type = 'scalar', scalar invariance is tested.
# 3. data is the data that the two-group CFA model fits.
# CJ: I believe this is inside the lavaan object, so don't need this
# 4. a, a numeric value that specifies the standardized difference tolerance.
# CJ: Just call this tolerance then!
# 5. minBF, a numeric value that specifies the minimum value of Bayes factors.
#   if Bayes factor >= minBF, full invariance is supported. Otherwise,
#   each pair of loadings or intercepts and partial invariance will be tested.
# 6. times, an integer number that specifies the times of the minimal training
#    sample size J in the fraction b. For example, if times = 1, the
#    fraction uses J; if times = 2, the fraction uses 2J.
# CJ: Is the times argument analogous to the fraction argument in bain? If so,
# CJ: I would recommend also having a fraction argument here.
# 7. seed, an integer number that creates a repeatable random number sequence.
# note the example was attached after this function
# CJ: The random seed should always be set outside the function, so it is explicit
# CJ: So please remove the seed argument

#' @title Bayesian Measurement Invariance Test
#' @description FUNCTION_DESCRIPTION
#' @param x An object for which a method exists.
#' @param type Character, indicating the type of measurement invariance to test.
#' Choose from `c('metric', 'scalar'). Default: `'metric'`.
#' @param tolerance Numeric, specifies the standardized difference tolerance.
#' Default: 0.2
#' @param minBF Numeric, specifies the minimum value of Bayes factors.
#' If Bayes factor >= minBF, full invariance is supported. Otherwise,
#' each pair of loadings or intercepts and partial invariance will be tested.
#' Default: 3
#' @param fraction A number representing the fraction of information in the data
#' used to construct the prior distribution. The default value 1 denotes the
#' minimal fraction, 2 denotes twice the minimal fraction, etc.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#' # Fit a multiple group latent regression model
#' model <- 'A  =~ Ab + Al + Af + An + Ar + Ac'
#' res <- cfa(model, data = sesamesim, std.lv = TRUE, group = "sex")
#' # Set a random seed
#' set.seed(2020)
#' # Calculate approximate metric invariance
#' result1 <- bmi(res)
#'  }
#' }
#' @rdname bmi
#' @export
bmi <-
  function(x,
           type = "metric",
           tolerance = .2,
           minBF = 3,
           fraction = 1) {
    UseMethod("bmi", x)
  }

#' @method bmi lavaan
#' @export
bmi.lavaan <-
  function(x,
           type = "metric",
           tolerance = .2,
           minBF = 3,
           fraction = 1) {
    # Capture function call
    cl <- match.call()
    # Add arguments
    #Part1: prepare for recomputing the model
    # Replace function call with internal mi function
    cl[[1L]] <- str2lang(paste0("bmi_", type))
    # Evaluate MI
    out <- eval.parent(cl)
    class(out) <- c("bayesian_invariance", class(out))
    return(out)
  }# end function

# CJ: You still have to pass the correct arguments here
bmi_metric <- function(...) {
  cl <- match.call()
  dots <- list(...)
  suppressMessages(attach(dots))
  the_items <- lavaan::lavNames(x, type = "ov")
  if (!length(the_items) > 0) {
    warning("Could not establish partial invariance; all items are variant.")
    out <- dots[c("BF_vec", "BF_item_mat")]
    out$warnings <-
      "Could not establish partial invariance; all items are variant."
    class(out) <- c("bmi_object", class(out))
    return(out)
  }
  if (is.null(dots[["BF_vec"]]))
    BF_vec <- NULL
  if (is.null(dots[["BF_item_mat"]]))
    BF_item_mat <-
    matrix(
      nrow = 0,
      ncol = length(the_items),
      dimnames = list(NULL, the_items)
    )
  partab <- lavaan::partable(x)
  indicatorNum <- length(lavaan::lavNames(x, type = "ov"))
  n_groups <- length(lavInspect(x, what = "group.label"))
  # fraction b
  fraction_b <-
    fraction * (indicatorNum + 1) / (2 * lavInspect(x, what = "nobs"))
  #Part2: recreate the model in lavaan to align the loadings
  # compute the factor's variance of group 1 and 2 that the product of loadings per group = 1
  # run lavaan again
  BF_Args <- prepare_bf_args(x, items = the_items, fraction_b = fraction_b, tolerance = tolerance)
  # BF <-

  BF <- do.call(calc_bf_mi, BF_Args)

  BF_item <- sapply(seq(indicatorNum), function(i) {
      calc_bf_mi(
        lower = BF_Args$lower[i],
        upper = BF_Args$upper[i],
        mprior = 0,
        mpost = BF_Args$mpost[i],
        sprior = BF_Args$spost[i, i, drop = FALSE],
        spost = BF_Args$sprior[i, i, drop = FALSE]
      )
    })
  BF_drop <- NULL
  if (BF < minBF) {
    drop_in_order <- the_items[order(BF_item, decreasing = FALSE)]
    BF_drop <- sapply(1:(length(drop_in_order)-1), function(i) {
      # Recalculate the Bayes factor, dropping an increasing number of items
      update_items <- the_items[!the_items %in% drop_in_order[1:i]]
      BF_Args <- prepare_bf_args(x, items = update_items, fraction_b = fraction_b, tolerance = tolerance)
      do.call(calc_bf_mi, BF_Args)
    })
  }
  out <- list(BF = BF, BF_item = BF_item, BF_drop = BF_drop)
  out[["cl"]] <- cl
  class(out) <- c("bmi_object", class(out))
  return(out)
}

prepare_bf_args <- function(x, items, fraction_b, tolerance) {
  loadings_original <- get_loadings_by_group(x)
  var_by_group <- get_var_by_group(x, items = items)
  lavAdjusted <- lav_fix_var(x, var_by_group)
  #Part3: get the adjusted loadings, that is, comparable loadings
  # get the adjusted loadings
  loadings_adjusted <- get_loadings_by_group(lavAdjusted, items = items)
  indicatorNum <- nrow(loadings_adjusted)
  # get the posterior and prior means
  # CJ: The mean of the posterior is the difference between the two group loadings
  meanPosterior <- loadings_adjusted[, 1] - loadings_adjusted[, 2]
  meanPrior <- rep(0, indicatorNum)
  #Part4: get upper and lower bounds
  bounds <- get_bounds(lavAdjusted, items = items, var_by_group, tolerance)

  #Part5: get posterior and prior covariance
  # get covariance of group 1 and 2
  # get posterior covariance
  covariance <- lav_cov_matrix(lavAdjusted, items = items)
  # CJ: A helps us get the covariance of the difference between indicators
  A <- cbind(diag(indicatorNum),-1 * diag(indicatorNum))
  covPosterior <- A %*% covariance %*% t(A)

  # get the prior covariance
  # Make block diagonal matrix of fractions to divide by
  divby <-
    as.matrix(do.call(
      bdiag,
      lapply(fraction_b, matrix, nrow = indicatorNum, ncol = indicatorNum)
    ))
  # Divide covariance by fractions
  covPrior <- covariance / divby
  # Set NaNs to zero (divided by zero)
  covPrior[is.nan(covPrior)] <- 0
  covPrior <- A %*% covPrior %*% t(A)
  return(list(lower = bounds[["lowerBound"]],
              upper = bounds[["upperBound"]],
              mprior = meanPrior,
              mpost = meanPosterior,
              sprior = covPrior,
              spost = covPosterior))
}

#
# recursive_partial_bf <- function(x, bmi_object){
#   BF <- bmi_object$BF
#   BF_item <- bmi_object$BF_item
#   cl <- bmi_object$call
#   browser()
#   #Part6: compute BF for each pair of the loadings
#   the_items <- lavaan::lavNames(x, type = "ov")
#   if(!length(the_items) > 0){
#     stop("Could not establish partial invariance; all items are variant.")
#   }
#   drop_item <- the_items[which.min(BF_item)]
#   x <- lav_drop_item(x, drop_item)
#   cl[[x]] <- x
#   cl$
#   if(partBF < minBF){
#     return("something")
#   }
#
#     # Part7: compute partial BF. If patial BF always < minBF,
#     #       the loop will stop with length(PartBF_item) == 1.
#     # CJ: WHile loops are risky. This should probably be a recursive function.
#
#     loadings_partial <- loadings_partial[-which.min(BF_item), ]
#
#     # the number of the partial loadings
#       partNum <- length(partLoadG2)
#       # recomputer factor's variance of group 1 and 2 that the product of partial loadings per group = 1
#       partVarG1 <- (prod(partLoadG1)**(1/partNum))**2
#       partVarG2 <- (prod(partLoadG2)**(1/partNum))**2
#       # adjust factor's variance
#       # CJ: Instead of specifying model syntax as text, I think it's better to just provide the parameter table?
#       adjustPartFactor <- paste0('F', ' ~~ ' ,'c(',partVarG1,',',
#                                  partVarG2,')* ','F')
#       # set model for lavaan
#       partAdjustedModel <- c(Equation,adjustPartFactor)
#       # run lavaan again
#       lavPart <- cfa(partAdjustedModel, dat, std.lv = TRUE, group = groupVariable)
#       ## get the adjusted loadings
#       # get the adjusted estimates
#       estiPart <- parameterEstimates(lavPart)
#       # get the adjusted loadings
#       partLoad <- estiPart[estiPart$op == "=~", "est"]
#       # get the adjusted loadings of group 1 and 2
#       AdjustedPartLoadG1 <- partLoad[1:indicatorNum]
#       AdjustedPartLoadG2 <- partLoad[(indicatorNum + 1):(2*indicatorNum)]
#       ## get the upper bound
#       # get the covariance of the observed variables
#       partCovx <- lavInspect(lavPart, what = "cov.ov")
#       # The max value of the loading in group 1 and 2
#       partLoadMaxG1 <- sqrt(diag(matrix(unlist(partCovx[1]),nrow=indicatorNum))/partVarG1)
#       partLoadMaxG2 <- sqrt(diag(matrix(unlist(partCovx[2]),nrow=indicatorNum))/partVarG2)
#       # get the bound
#       partLoadMax <- 0
#       for(i in 1:indicatorNum){
#         if(partLoadMaxG1[i]>partLoadMaxG2[i]){
#           partLoadMax[i] <- tolerance *partLoadMaxG2[i]
#
#         }else{
#           partLoadMax[i] <- tolerance *partLoadMaxG1[i]
#         }
#       }
#       # get the upper bound
#       partUpperBound <- unlist(partLoadMax)
#       ## get the posterior covariance
#       # get the covariance of group 1 and 2
#       partCov1 <- lavInspect(lavPart, "vcov")[varName1, varName1]
#       partCov2 <- lavInspect(lavPart, "vcov")[varName2, varName2]
#       # get the posterior covariance
#       partCov <- as.matrix(bdiag(partCov1,partCov2))
#       partCovPosterior<- A%*%partCov%*%t(A)
#       # get the prior covariance
#       partCovPrior1 <- partCov1/fraction_b[1]
#       partCovPrior2 <- partCov2/fraction_b[2]
#       partCovPrior <- as.matrix(bdiag(partCovPrior1,partCovPrior2))
#       partCovPrior <- A%*%partCovPrior%*%t(A)
#       ## get partial loadings, upper bound, and prior and posterior covarance
#       PartBF_item2 <- unlist(BF_item)
#       # the number of the deleted loadings
#       DeletedNum <- indicatorNum - partNum
#       for (j in 1:DeletedNum){
#         # the partial loadings of group 1 and 2
#         AdjustedPartLoadG1 <- AdjustedPartLoadG1[- which.min(PartBF_item2)]
#         AdjustedPartLoadG2 <- AdjustedPartLoadG2[- which.min(PartBF_item2)]
#         # the partial upper bound
#         partUpperBound <- partUpperBound[- which.min(PartBF_item2)]
#         # the partial posterior covariance
#         partCovPosterior <- partCovPosterior[- which.min(PartBF_item2),- which.min(PartBF_item2)]
#         # the partial prior covariance
#         partCovPrior <- partCovPrior[- which.min(PartBF_item2),- which.min(PartBF_item2)]
#         # control this loop process
#         PartBF_item2 <- PartBF_item2[- which.min(PartBF_item2)]
#       }
#       # get the partial posterior and prior means
#       partMeanPosterior <- AdjustedPartLoadG1 - AdjustedPartLoadG2
#       partMeanPrior <- rep(0, partNum)
#       # get the partial lower bound
#       partLowerBound <- -partUpperBound
#       ## compute partial BF
#       # CJ: Here you set the same seed every time. The result is: no more random sampling
#       # set.seed(seed)
#       partBF <- (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]))/
#         (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]))
#       # control the loop process
#       PartBF_item <- PartBF_item[- which.min(PartBF_item)]
#     }#end while
#     ## arrange partial BF
#     partBF <- as.data.frame(partBF)
#     # CJ: WHy are you rounding this? THat loses information
#     #partBF <- round(partBF,3)
#     colnames(partBF) <- 'partial'
#     # get the deleted indicator
#     deletedIndicator <- setdiff(names(BF_item), names(PartBF_item))
#     ## print BF of all loadings, BF of each loadings, BF of partial loadings
#     # judge partial BF >= minBF or not
#     # CJ: Do not do this inside the function to compute measurement invariance.
#     # CJ: Make a separate function that reports the results, and which does this
#     if(partBF >= minBF) {
#       cat("\n With the minimum Bayes factor >=",minBF,",","partial approximate metric invariance is supported.
#         \n The following indicators are free :\n",paste0(deletedIndicator,collapse = " , "),"\n\n")
#     } else {
#       # if partial BF < minBF, partial BF do not exist.
#       partBF <- 'null'
#       cat("\n With the minimum Bayes factor >=", minBF, ",loadings are not
#           \n about equal across groups, respectively.Please check your model or minBF! \n\n")
#     }
#     ## arrange a table
#     BF_item <- as.data.frame(BF_item)
#     # CJ: DO NOT ROUND! It loses information
#     #BF_item <- round(BF_item,3)
#     out <- cbind(BF,partBF,BF_item)
#     # CJ: Do not print output from the function that calculates invariance.
#     # CJ: Make a separate print function
#     return(out)
#   }# end if (BF < minBF)
# } # end testing metric invariance

# CJ: You still have to pass the correct arguments here
bmi_scalar <- function(...) {
  #Part2: recreate model in lavaan to align intercepts
  # get the original intercepts
  Nint <- estiOut[estiOut$op == "~1", "est"]
  # get the intercepts of group 1 and 2
  intG1 <- Nint[1:indicatorNum]
  intG2 <- Nint[(indicatorNum + 2):(2 * indicatorNum + 1)]
  # compute the factor's mean of group1 and 2 that the sum of intercepts per group = 0
  meanG1 <- sum(intG1) / sum(loadG1)
  meanG2 <- sum(intG2) / sum(loadG2)
  # adjust factor's mean
  adjustFactorM <- paste0('F', ' ~ ' , 'c(', meanG1, ',',
                          meanG2, ')* ', '1')
  # set model for lavaan
  Model <- c(Equation, adjustFactor, adjustFactorM)
  # run lavaan again
  lavAdjusted <-
    cfa(Model, dat, std.lv = TRUE, group = groupVariable)
  #Part3: get the adjusted intercepts, that is, comparable intercepts
  # get the adjusted estimates
  estiAdjusted <- parameterEstimates(lavAdjusted)
  # get the adjusted intercepts
  AdjustedInt <- estiAdjusted[estiAdjusted$op == "~1", "est"]
  # get the adjusted intercepts of group 1 and 2
  AdjustedInt1 <- AdjustedInt[1:indicatorNum]
  AdjustedInt2 <-
    AdjustedInt[(indicatorNum + 2):(2 * indicatorNum + 1)]
  # get the posterior and prior means
  meanPosterior <- AdjustedInt1 - AdjustedInt2
  meanPrior <- rep(0, indicatorNum)
  #Part4: get upper and lower bounds
  # data from group 1 and 2
  da <- lavInspect(lavAdjusted, what = "data")
  dataG1 <- as.data.frame(da[1])
  dataG2 <- as.data.frame(da[2])
  # robust max and min intercept in group 1
  # CJ: Don't use an extra package just to compute the quantile
  difG1 <-
    apply(sapply(dataG1, quantile, probs = c(.2, .8)), 2, diff)

  # robust max and min intercept in group 2
  difG2 <-
    apply(sapply(dataG2, quantile, probs = c(.2, .8)), 2, diff)
  # get the bound
  # CJ: Instead of a for loop, use vectorization and boolean indexing
  dif <- tolerance * difG1[i]
  dif[difG1 > difG2] <- (tolerance * difG2)[difG1 > difG2]

  # the upper and lower bounds
  upperBound <- unlist(dif)
  lowerBound <- -upperBound
  #Part5: get posterior and prior covariance
  # get covariance of group 1 and 2
  varName1 <- paste0(indicatorName, '~1')
  varName2 <- paste0(indicatorName, '~1.g2')
  cov1 <- lavInspect(lavAdjusted, "vcov")[varName1, varName1]
  cov2 <- lavInspect(lavAdjusted, "vcov")[varName2, varName2]
  # get posterior covariance
  covariance <- as.matrix(bdiag(cov1, cov2))
  A1 <- diag(indicatorNum)
  A2 <- diag(indicatorNum) * (-1)
  A <- cbind(A1, A2)
  covPosterior <- A %*% covariance %*% t(A)

  # get prior covariance
  covPrior1 <- cov1 / fraction_b[1]
  covPrior2 <- cov2 / fraction_b[2]
  covPrior <- as.matrix(bdiag(covPrior1, covPrior2))
  covPrior <- A %*% covPrior %*% t(A)
  ## compute Bayes factor
  set.seed(seed)
  BF <-
    (
      pmvnorm(
        lower = lowerBound,
        upper = upperBound,
        mean = meanPosterior,
        sigma = covPosterior
      )[1] * (
        1 - pmvnorm(
          lower = lowerBound,
          upper = upperBound,
          mean = meanPrior,
          sigma = covPrior
        )[1]
      )
    ) /
    (
      pmvnorm(
        lower = lowerBound,
        upper = upperBound,
        mean = meanPrior,
        sigma = covPrior
      )[1] * (
        1 - pmvnorm(
          lower = lowerBound,
          upper = upperBound,
          mean = meanPosterior,
          sigma = covPosterior
        )[1]
      )
    )
  BF <- as.data.frame(BF)
  colnames(BF) <- 'scalar'
  rownames(BF) <- 'Bayes factor'
  #BF <- round(BF,3)
  BF1 <- unlist(BF)
  if (is.nan(BF1)) {
    cat('Approximate scalar invariance : Bayes factor cannot be computed. \n\n')
    print(BF)
  }
  else if (BF >= minBF) {
    cat(
      'When the minimum Bayes factor should  >=',
      minBF,
      ',',
      'approximate scalar invariance is supported. \n\n'
    )
    print(BF)
  } else {
    #Part6: compute BF for each pair of the intercept
    BF_item <- list()
    for (i in 1:indicatorNum) {
      # get each posterior mean
      eachMeanPosterior <- meanPosterior[i]
      # get each upper and lower bound
      eachUpperBound <- upperBound[i]
      eachLowerBound <- -eachUpperBound
      # get each posterior covariance
      eachCovPosterior <- covPosterior[i, ]
      eachCovPosterior <- eachCovPosterior[i]
      eachCovPosterior <- matrix(eachCovPosterior, 1, 1)
      # get each prior covariance
      eachCovPrior <- covPrior[i, ]
      eachCovPrior <- eachCovPrior[i]
      eachCovPrior <- matrix(eachCovPrior, 1, 1)
      # Bayes Factor for each loading
      set.seed(seed)
      BF_item[i] <-
        (
          pmvnorm(
            lower = eachLowerBound,
            upper = eachUpperBound,
            mean = eachMeanPosterior,
            sigma = eachCovPosterior
          )[1] * (
            1 - pmvnorm(
              lower = eachLowerBound,
              upper = eachUpperBound,
              mean = 0,
              sigma = eachCovPrior
            )[1]
          )
        ) /
        (
          pmvnorm(
            lower = eachLowerBound,
            upper = eachUpperBound,
            mean = 0,
            sigma = eachCovPrior
          )[1] * (
            1 - pmvnorm(
              lower = eachLowerBound,
              upper = eachUpperBound,
              mean = eachMeanPosterior,
              sigma = eachCovPosterior
            )[1]
          ) / (
            1 - pmvnorm(
              lower = eachLowerBound,
              upper = eachUpperBound,
              mean = 0,
              sigma = eachCovPrior
            )[1]
          )
        )
      #BF_item[[i]] <- round(BF_item[[i]],3)
      # name BF with the indicator's name
      names(BF_item)[i] <- paste(indicatorName[i],  sep = "")
    } # end for(i in 1:indicatorNum)
    partBF <- BF
    PartBF_item <- unlist(BF_item)
    partIntG1 <- intG1
    partIntG2 <- intG2
    #Part7 : compute partial BF. If patial BF always < minBF,
    #        the loop will exit with length(PartBF_item) == 1.
    while (partBF < minBF & length(PartBF_item) != 1) {
      ## recompute model
      # get partial intercepts
      partIntG1 <- partIntG1[-which.min(PartBF_item)]
      partIntG2 <- partIntG2[-which.min(PartBF_item)]
      # the number of the partial intercepts
      partNum <- length(partIntG2)
      # recompute factor's mean of group 1 and 2 that the sum of the partial intercepts per group = 0
      partMeanG1 <- sum(partIntG1) / sum(loadG1)
      partMeanG2 <- sum(partIntG2) / sum(loadG2)
      # adjust factor's mean
      adjustPartFactorM <- paste0('F', ' ~ ' , 'c(', partMeanG1, ',',
                                  partMeanG2, ')* ', '1')
      # set model for lavaan
      adjustedModel <- c(Equation, adjustFactor, adjustPartFactorM)
      # run lavaan again
      lavPart <-
        cfa(adjustedModel, dat, std.lv = TRUE, group = groupVariable)
      ## get the adjusted intercepts
      # get the adjusted estimates
      estiPart <- parameterEstimates(lavPart)
      # get the adjusted intercepts
      partInt <- estiPart[estiPart$op == "~1", "est"]
      # get the adjusted intercepts of group 1 and 2
      AdjustedPartIntG1 <- partInt[1:indicatorNum]
      AdjustedPartIntG2 <-
        partInt[(indicatorNum + 2):(2 * indicatorNum + 1)]
      ## get the upper bound
      partUpperBound <- upperBound
      ## get posterior and prior covariance
      # get covariance of group 1 and 2
      partCov1 <- lavInspect(lavPart, "vcov")[varName1, varName1]
      partCov2 <- lavInspect(lavPart, "vcov")[varName2, varName2]
      # posterior covariance
      partCov <- as.matrix(bdiag(partCov1, partCov2))
      partCovPosterior <- A %*% partCov %*% t(A)
      # prior covariance
      partCovPrior1 <- partCov1 / fraction_b[1]
      partCovPrior2 <- partCov2 / fraction_b[2]
      partCovPrior <-
        as.matrix(bdiag(partCovPrior1, partCovPrior2))
      partCovPrior <- A %*% partCovPrior %*% t(A)
      ## get partial intercepts, upper bound, and prior and posterior covarance
      PartBF_item2 <- unlist(BF_item)
      # the number of the deleted loadings
      DeletedNum <- indicatorNum - partNum
      for (j in 1:DeletedNum) {
        # the partial loadings of group 1 and 2
        AdjustedPartIntG1 <-
          AdjustedPartIntG1[-which.min(PartBF_item2)]
        AdjustedPartIntG2 <-
          AdjustedPartIntG2[-which.min(PartBF_item2)]
        # the partial upper bound
        partUpperBound <-
          partUpperBound[-which.min(PartBF_item2)]
        # the partial posterior covariance
        partCovPosterior <-
          partCovPosterior[-which.min(PartBF_item2), -which.min(PartBF_item2)]
        # the partial prior covariance
        partCovPrior <-
          partCovPrior[-which.min(PartBF_item2), -which.min(PartBF_item2)]
        # control this loop process
        PartBF_item2 <- PartBF_item2[-which.min(PartBF_item2)]
      }
      # get partial posterior and prior means
      partMeanPosterior <- AdjustedPartIntG1 - AdjustedPartIntG2
      partMeanPrior <- rep(0, partNum)
      #  get the partial lower bound
      partLowerBound <- -partUpperBound
      ## compute partial BF
      set.seed(seed)
      partBF <-
        (
          pmvnorm(
            lower = partLowerBound,
            upper = partUpperBound,
            mean = partMeanPosterior,
            sigma = partCovPosterior
          )[1] * (
            1 - pmvnorm(
              lower = partLowerBound,
              upper = partUpperBound,
              mean = partMeanPrior,
              sigma = partCovPrior
            )[1]
          )
        ) /
        (
          pmvnorm(
            lower = partLowerBound,
            upper = partUpperBound,
            mean = partMeanPrior,
            sigma = partCovPrior
          )[1] * (
            1 - pmvnorm(
              lower = partLowerBound,
              upper = partUpperBound,
              mean = partMeanPosterior,
              sigma = partCovPosterior
            )[1]
          )
        )
      # control the loop process
      PartBF_item <- PartBF_item[-which.min(PartBF_item)]
    }#end while
    ## arrange partial BF
    partBF <- as.data.frame(partBF)
    #partBF <- round(partBF,3)
    colnames(partBF) <- 'partial'
    # get deleted indicator
    deletedIndicator <-
      setdiff(names(BF_item), names(PartBF_item))
    ## print BF of all intercepts, BF of each intercepts, BF of partial intercepts
    ## arrange a table
    BF_item <- as.data.frame(BF_item)
    #BF_item <- round(BF_item,3)
    out <- cbind(BF, partBF, BF_item)
    # CJ: Do not print in this function
    return(out)
  }# end if (BF < minBF)
}

#' @method print bayesian_invariance
#' @export
print.bayesian_invariance <- function(x, ...) {
  # CJ: DO all of the printing here

  # judge partial BF > minBF or not
  if (partBF >= minBF) {
    cat(
      "\n With the minimum Bayes factor >=",
      minBF,
      ",",
      "partial approximate scalar invariance is supproted.
        \n The following indicators are free :\n",
      paste0(deletedIndicator, collapse = " , "),
      "\n\n"
    )
  } else {
    # if partial BF < minBF, partial BF do not exist.
    partBF <- 'null'
    cat(
      "\n With the minimum Bayes factor >=",
      minBF,
      ", intercepts are not
          \n about equal across group, respectively.Please check your model or minBF! \n\n"
    )
  }
}

get_bounds <- function(x, items = NULL, var_by_group, tolerance) {
  # get the covariance of the observed variables
  covx <- lavInspect(x, what = "cov.ov")
  if(!is.null(items)) covx <- lapply(covx, function(thiscov){thiscov[items, items, drop = FALSE]})
  # The max value of the loadings in group 1 and 2
  loadmax_grp <- do.call(cbind, lapply(seq_along(covx), function(i) {
    sqrt(diag(covx[[i]]) / var_by_group[i])
  }))
  colnames(loadmax_grp) <- names(covx)
  # the upper and lower bounds
  upperBound <- loadmax_grp[, 1]
  upperBound[loadmax_grp[, 1] > loadmax_grp[, 2]] <-
    loadmax_grp[, 2][loadmax_grp[, 1] > loadmax_grp[, 2]]
  upperBound <- tolerance * upperBound
  lowerBound <- -1 * upperBound
  list(lowerBound = lowerBound, upperBound = upperBound)
}

calc_bf_mi <- function(lower, upper, mprior, mpost, sprior, spost) {
  pmvpost <-
    pmvnorm(
      lower = lower,
      upper = upper,
      mean = mpost,
      sigma = spost
    )
  pmvprior <-
    pmvnorm(
      lower = lower,
      upper = upper,
      mean = mprior,
      sigma = sprior
    )
  BF <- (pmvpost * (1 - pmvprior)) / (pmvprior * (1 - pmvpost))
  if (is.nan(BF)) {
    # SHould these be error messages instead?
    message('Approximate metric invariance : Bayes factor cannot be computed. \n\n')
  }
  # Run a test: Make the BF NaN, and see if the rest of the function works
  return(BF)
}

get_lav_data <- function(x, ...) {
  dat <- lavInspect(object = x, what = "data")
  if (inherits(dat, "list")) {
    if (length(dat) > 1) {
      grp <- lavInspect(object = x, what = "group")
      grp_levels <- lavInspect(object = x, what = "group.label")
      dat <-
        do.call(rbind, lapply(seq_along(grp_levels), function(i) {
          out <- as.data.frame(dat[[i]])
          out[[grp]] <- grp_levels[i]
          out
        }))
    } else {
      dat <- dat[[1]]
    }
  }
  dat
}

get_loadings_by_group <- function(x, items = NULL) {
  ests <- lavaan::inspect(x, what = "est")
  if (!is.null(ests[["lambda"]])) {
    out <- matrix(ests[["lambda"]], ncol = 1)
  } else {
    groupnames <- names(ests)
    out <- do.call(cbind, lapply(ests, `[[`, "lambda"))
    colnames(out) <- groupnames
  }
  if(!is.null(items)) out <- out[items, , drop = FALSE]
  return(out)
}

get_var_by_group <- function(x, items = NULL) {
  loadings <- get_loadings_by_group(x)
  if (is.null(items))
    items <- rownames(loadings)
  loadings <-
    loadings[rownames(loadings) %in% items, , drop = FALSE]
  apply(loadings, 2, function(i) {
    (prod(i) ** (1 / nrow(loadings))) ** 2
  })
}

lav_drop_item <- function(x, item) {
  thepars <- partable(x)
  thepars <-
    thepars[!(thepars$lhs == item |
                thepars$rhs == item), c("lhs", "op", "rhs", "group", "free", "ustart")]
  lavaan::update(x, model = thepars)
}

lav_cov_matrix <- function(x, items = NULL){
  vcv = x@vcov$vcov
  rownames(vcv) <- colnames(vcv) <- x@ParTable$rhs[x@ParTable$free > 0]
  indx <- (x@ParTable$op == "=~")[x@ParTable$free > 0]
  if(is.null(items)){
    return(vcv[indx, indx, drop = FALSE])
  } else {
    out <- vcv[indx, indx, drop = FALSE]
    indx <- which(rownames(out) %in% items)
    return(out[indx, indx, drop = FALSE])
  }
}

lav_fix_var <- function(x, var_by_group) {
  partab <- lavaan::partable(x)
  lv_nam <- lavaan::lavNames(partab, type = "lv")
  add_this <- data.frame(
    lhs = lv_nam,
    op = "~~",
    rhs = lv_nam,
    group = 1:length(var_by_group),
    free = 0,
    ustart = var_by_group
  )
  new_partab <- suppressWarnings(lav_partable_merge(
    partab,
    add_this,
    remove.duplicated = TRUE,
    fromLast = TRUE
  )[, c("lhs", "op", "rhs", "group", "free",
        "ustart")])
  lavaan::update(x, model = new_partab)
}

# section 1: example
# CFA model
model22 <- 'A  =~ Ab + Al + Af + An + Ar + Ac'
# Fit the multiple group latent regression model
lavout <-
  cfa(model22,
      data = sesamesim,
      std.lv = TRUE,
      group = "sex")
# approximate metric invariance
set.seed(2020)
result1 <-
  bmi(
    lavout,
    type = 'metric',
    minBF = 120,
    tolerance = 0.2,
    fraction = 2
  )
# partial approximate metric invariance
# result2 <- BMI1(lavout,type = 'metric',data = sesamesim, tolerance = 0.2,minBF = 65,times=2)
# # approximate scalar invariance
# result3 <- BMI1(lavout,type = 'scalar',data = sesamesim, tolerance = 0.2,minBF = 3,times=2)
# # partial approximate scalar invariance
# result4 <- BMI1(lavout,type = 'scalar',data = sesamesim, tolerance = 0.2,minBF = 10,times=2)

library(microbenchmark)
microbenchmark(
  lavInspect(lavout, "vcov"),
  lav_cov_matrix(lavout)
)
