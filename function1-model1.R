# CJ: Please read the tidyverse style guide for R code;
# CJ: for example, there should always be a space between comma's and mathematical operators,
# CJ: variable names should be snake_case not camelCase, and function names should be verbs when possible
# CJ: https://style.tidyverse.org/syntax.html
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
bmi <- function(x, type = "metric", tolerance = .2, minBF = 3, fraction = 1){
  UseMethod("bmi", x)
}

#' @method bmi lavaan
#' @export
bmi.lavaan <- function(x, type = "metric", tolerance = .2, minBF = 3, fraction = 1){
  dat <- as.data.frame(lavInspect(object = x, what = "data"))
  names(dat) <- lavNames(x)
  #Part1: prepare for recomputing the model
  # get the estimates
  estiOut <- parameterEstimates(x)
  # get the indicaors' name
  indicatorName <- unique(estiOut[estiOut$op == "=~", "rhs"])
  # get the number of indicator
  indicatorNum <- length(indicatorName)
  # add the indicators one by one
  indicatorSum <- paste0(indicatorName,collapse = " + ")
  # create the CFA equation
  Equation <- paste0('F', ' =~ ' ,indicatorSum)
  # adjust latent variable
  # adjustFactor <- paste0('F', ' ~~ ' ,'c(',varG1,',',
  #                      varG2,')* ','F')
  # get the group variable's name
  groupVariable <- lavInspect(x, what = "group")
  # testing metric invariance
  # CJ: Make a separate function for metric invariance
  out <- switch(type,
                "metric" = bmi_metric(...), # CJ: You still have to pass the correct arguments here
                "scalar" = bmi_scalar(...)) # CJ: You still have to pass the correct arguments here
  class(out) <- c("bayesian_invariance", class(out))
  return(out)
}# end function

# CJ: You still have to pass the correct arguments here
bmi_metric <- function(...){
  #Part2: recreate the model in lavvan to align the loadings
  # get the original loadings
  Nload <- estiOut[ estiOut$op == "=~", "est"]
  # get the loadings of group 1 and 2
  loadG1 <- Nload[1:indicatorNum]
  loadG2 <- Nload[(indicatorNum + 1):(2*indicatorNum)]
  # compute the factor's variance of group 1 and 2 that the product of loadings per group = 1
  varG1 <- (prod(loadG1)**(1/indicatorNum))**2
  varG2 <- (prod(loadG2)**(1/indicatorNum))**2
  # adjust the factor's variance of group 1 and 2
  adjustFactor <- paste0('F', ' ~~ ' ,'c(',varG1,',',
                         varG2,')* ','F')
  # set model for lavaan
  Model <- c(Equation,adjustFactor)
  # run lavvan again
  lavAdjusted <- cfa(Model, dat, std.lv = TRUE, group = groupVariable)
  #Part3: get the adjusted loadings, that is, comparable loadings
  # get the adjusted estimates
  estiAdjusted <- parameterEstimates(lavAdjusted)
  # get the adjusted loadings
  AdjustedLoad <- estiAdjusted[estiAdjusted$op == "=~", "est"]
  # get the adjusted loadings of group 1 and 2
  AdjustedLoadG1 <- AdjustedLoad[1:indicatorNum]
  AdjustedLoadG2 <- AdjustedLoad[(indicatorNum + 1):(2*indicatorNum)]
  # get the posterior and prior means
  meanPosterior <- AdjustedLoadG1 - AdjustedLoadG2
  meanPrior <- rep(0, indicatorNum)
  #Part4: get upper and lower bounds
  # get the covariance of the observed variables
  covx <- lavInspect(lavAdjusted, what = "cov.ov")
  # The max value of the loadings in group 1 and 2
  loadMaxG1 <- sqrt(diag(matrix(unlist(covx[1]),nrow=indicatorNum))/varG1)
  loadMaxG2 <- sqrt(diag(matrix(unlist(covx[2]),nrow=indicatorNum))/varG2)
  # get the bound
  loadMax <- 0
  for(i in 1:indicatorNum){
    if(loadMaxG1[i]>loadMaxG2[i]){
      loadMax[i] <- a*loadMaxG2[i]

    }else{
      loadMax[i] <- a*loadMaxG1[i]
    }
  }
  # the upper and lower bounds
  upperBound <- unlist(loadMax)
  lowerBound <- -upperBound
  #Part5: get posterior and prior covariance
  # get covariance of group 1 and 2
  varName1 <- paste0('F', '=~' ,indicatorName)
  varName2 <- paste0('F', '=~' ,indicatorName,'.g2')
  cov1 <- lavInspect(lavAdjusted, "vcov")[varName1, varName1]
  cov2 <- lavInspect(lavAdjusted, "vcov")[varName2, varName2]
  # get posterior covariance
  covariance <- as.matrix(bdiag(cov1,cov2))
  A1<- diag(indicatorNum)
  A2<- diag(indicatorNum)*(-1)
  A <- cbind(A1,A2)
  covPosterior<- A%*%covariance%*%t(A)
  # get sample size
  sampsizes <- lavInspect(lavAdjusted, what = "nobs")
  # fraction b
  b <- times*(indicatorNum+1)/(2*sampsizes)
  # get the prior covariance
  covPrior1 <- cov1/b[1]
  covPrior2 <- cov2/b[2]
  covPrior <- as.matrix(bdiag(covPrior1,covPrior2))
  covPrior <- A%*%covPrior%*%t(A)
  ## compute Bayes factor
  set.seed(seed)
  BF <- (pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPosterior, sigma=covPosterior)[1]*(1- pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPrior, sigma=covPrior)[1]))/
    (pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPrior, sigma=covPrior)[1]*(1- pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPosterior, sigma=covPosterior)[1]))
  BF <- as.data.frame(BF)
  colnames(BF) <- 'metric'
  rownames(BF) <- 'Bayes factor'
  #BF <- round(BF,3)
  BF1 <- unlist(BF)
  if(is.nan(BF1)){
    cat('Approximate metric invariance : Bayes factor cannot be computed. \n\n')
    print(BF)
  }
  else if(BF >= minBF){
    cat('When the minimum Bayes factor should  >=',minBF,',','approximate metric invariance is supported. \n\n')
    print(BF)
  } else {
    #Part6: compute BF for each pair of the loadings
    eachBF <- list()
    for(i in 1:indicatorNum){
      # get each posterior mean
      eachMeanPosterior <- meanPosterior[i]
      # get each upper and lower bound
      eachUpperBound <- upperBound[i]
      eachLowerBound <- -eachUpperBound
      # get each posterior covariance
      eachCovPosterior <- covPosterior[i,]
      eachCovPosterior <- eachCovPosterior[i]
      eachCovPosterior <- matrix(eachCovPosterior,1,1)
      # get each prior covariance
      eachCovPrior <- covPrior[i,]
      eachCovPrior <- eachCovPrior[i]
      eachCovPrior <- matrix(eachCovPrior,1,1)
      # Bayes Factor for each loading
      set.seed(seed)
      eachBF[i] <- (pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=eachMeanPosterior, sigma=eachCovPosterior)[1]*(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]))/
        (pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]*(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=eachMeanPosterior, sigma=eachCovPosterior)[1])/(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]))
      #eachBF[[i]] <- round(eachBF[[i]],3)
      # name each BF with the indicator's name
      names(eachBF)[i] <- paste(indicatorName[i],  sep = "")
    } # end for(i in 1:indicatorNum)
    partBF <- BF
    PartEachBF <- unlist(eachBF)
    partLoadG1 <- loadG1
    partLoadG2 <- loadG2
    # Part7: compute partial BF. If patial BF always < minBF,
    #       the loop will stop with length(PartEachBF) == 1.
    # CJ: WHile loops are risky. Make sure to set a maximum number of iterations
    while(partBF < minBF & length(PartEachBF) != 1){
      ## recompute model
      # get the partial loadings
      partLoadG1 <- partLoadG1[- which.min(PartEachBF)]
      partLoadG2 <- partLoadG2[- which.min(PartEachBF)]
      # the number of the partial loadings
      partNum <- length(partLoadG2)
      # recomputer factor's variance of group 1 and 2 that the product of partial loadings per group = 1
      partVarG1 <- (prod(partLoadG1)**(1/partNum))**2
      partVarG2 <- (prod(partLoadG2)**(1/partNum))**2
      # adjust factor's variance
      # CJ: Instead of specifying model syntax as text, I think it's better to just provide the parameter table?
      adjustPartFactor <- paste0('F', ' ~~ ' ,'c(',partVarG1,',',
                                 partVarG2,')* ','F')
      # set model for lavaan
      partAdjustedModel <- c(Equation,adjustPartFactor)
      # run lavvan again
      lavPart <- cfa(partAdjustedModel, dat, std.lv = TRUE, group = groupVariable)
      ## get the adjusted loadings
      # get the adjusted estimates
      estiPart <- parameterEstimates(lavPart)
      # get the adjusted loadings
      partLoad <- estiPart[estiPart$op == "=~", "est"]
      # get the adjusted loadings of group 1 and 2
      AdjustedPartLoadG1 <- partLoad[1:indicatorNum]
      AdjustedPartLoadG2 <- partLoad[(indicatorNum + 1):(2*indicatorNum)]
      ## get the upper bound
      # get the covariance of the observed variables
      partCovx <- lavInspect(lavPart, what = "cov.ov")
      # The max value of the loading in group 1 and 2
      partLoadMaxG1 <- sqrt(diag(matrix(unlist(partCovx[1]),nrow=indicatorNum))/partVarG1)
      partLoadMaxG2 <- sqrt(diag(matrix(unlist(partCovx[2]),nrow=indicatorNum))/partVarG2)
      # get the bound
      partLoadMax <- 0
      for(i in 1:indicatorNum){
        if(partLoadMaxG1[i]>partLoadMaxG2[i]){
          partLoadMax[i] <- tolerance *partLoadMaxG2[i]

        }else{
          partLoadMax[i] <- tolerance *partLoadMaxG1[i]
        }
      }
      # get the upper bound
      partUpperBound <- unlist(partLoadMax)
      ## get the posterior covariance
      # get the covariance of group 1 and 2
      partCov1 <- lavInspect(lavPart, "vcov")[varName1, varName1]
      partCov2 <- lavInspect(lavPart, "vcov")[varName2, varName2]
      # get the posterior covariance
      partCov <- as.matrix(bdiag(partCov1,partCov2))
      partCovPosterior<- A%*%partCov%*%t(A)
      # get the prior covariance
      partCovPrior1 <- partCov1/b[1]
      partCovPrior2 <- partCov2/b[2]
      partCovPrior <- as.matrix(bdiag(partCovPrior1,partCovPrior2))
      partCovPrior <- A%*%partCovPrior%*%t(A)
      ## get partial loadings, upper bound, and prior and posterior covarance
      PartEachBF2 <- unlist(eachBF)
      # the number of the deleted loadings
      DeletedNum <- indicatorNum - partNum
      for (j in 1:DeletedNum){
        # the partial loadings of group 1 and 2
        AdjustedPartLoadG1 <- AdjustedPartLoadG1[- which.min(PartEachBF2)]
        AdjustedPartLoadG2 <- AdjustedPartLoadG2[- which.min(PartEachBF2)]
        # the partial upper bound
        partUpperBound <- partUpperBound[- which.min(PartEachBF2)]
        # the partial posterior covariance
        partCovPosterior <- partCovPosterior[- which.min(PartEachBF2),- which.min(PartEachBF2)]
        # the partial prior covariance
        partCovPrior <- partCovPrior[- which.min(PartEachBF2),- which.min(PartEachBF2)]
        # control this loop process
        PartEachBF2 <- PartEachBF2[- which.min(PartEachBF2)]
      }
      # get the partial posterior and prior means
      partMeanPosterior <- AdjustedPartLoadG1 - AdjustedPartLoadG2
      partMeanPrior <- rep(0, partNum)
      # get the partial lower bound
      partLowerBound <- -partUpperBound
      ## compute partial BF
      set.seed(seed)
      partBF <- (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]))/
        (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]))
      # control the loop process
      PartEachBF <- PartEachBF[- which.min(PartEachBF)]
    }#end while
    ## arrange partial BF
    partBF <- as.data.frame(partBF)
    # CJ: WHy are you rounding this? THat loses information
    #partBF <- round(partBF,3)
    colnames(partBF) <- 'partial'
    # get the deleted indicator
    deletedIndicator <- setdiff(names(eachBF), names(PartEachBF))
    ## print BF of all loadings, BF of each loadings, BF of partial loadings
    # judge partial BF >= minBF or not
    # CJ: Do not do this inside the function to compute measurement invariance.
    # CJ: Make a separate function that reports the results, and which does this
    if(partBF >= minBF) {
      cat("\n With the minimum Bayes factor >=",minBF,",","partial approximate metric invariance is supported.
        \n The following indicators are free :\n",paste0(deletedIndicator,collapse = " , "),"\n\n")
    } else {
      # if partial BF < minBF, partial BF do not exist.
      partBF <- 'null'
      cat("\n With the minimum Bayes factor >=", minBF, ",loadings are not
          \n about equal across groups, respectively.Please check your model or minBF! \n\n")
    }
    ## arrange a table
    eachBF <- as.data.frame(eachBF)
    # CJ: DO NOT ROUND! It loses information
    #eachBF <- round(eachBF,3)
    out <- cbind(BF,partBF,eachBF)
    # CJ: Do not print output from the function that calculates invariance.
    # CJ: Make a separate print function
    return(out)
  }# end if (BF < minBF)
} # end testing metric invariance

# CJ: You still have to pass the correct arguments here
bmi_scalar <- function(...){
    #Part2: recreate model in lavvan to align intercepts
    # get the original intercepts
    Nint <- estiOut[estiOut$op == "~1", "est"]
    # get the intercepts of group 1 and 2
    intG1 <- Nint[1:indicatorNum]
    intG2 <- Nint[(indicatorNum + 2):(2*indicatorNum + 1)]
    # compute the factor's mean of group1 and 2 that the sum of intercepts per group = 0
    meanG1 <- sum(intG1)/sum(loadG1)
    meanG2 <- sum(intG2)/sum(loadG2)
    # adjust factor's mean
    adjustFactorM <- paste0('F', ' ~ ' ,'c(',meanG1,',',
                            meanG2,')* ','1')
    # set model for lavaan
    Model <- c(Equation,adjustFactor, adjustFactorM)
    # run lavvan again
    lavAdjusted <- cfa(Model, dat, std.lv = TRUE, group = groupVariable)
    #Part3: get the adjusted intercepts, that is, comparable intercepts
    # get the adjusted estimates
    estiAdjusted <- parameterEstimates(lavAdjusted)
    # get the adjusted intercepts
    AdjustedInt <- estiAdjusted[estiAdjusted$op == "~1", "est"]
    # get the adjusted intercepts of group 1 and 2
    AdjustedInt1 <- AdjustedInt[1:indicatorNum]
    AdjustedInt2 <- AdjustedInt[(indicatorNum + 2):(2*indicatorNum + 1)]
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
    difG1 <- apply(sapply(dataG1, quantile, probs = c(.2, .8)), 2, diff)

    # robust max and min intercept in group 2
    difG2 <- apply(sapply(dataG2, quantile, probs = c(.2, .8)), 2, diff)
    # get the bound
    # CJ: Instead of a for loop, use vectorization and boolean indexing
    dif <- tolerance*difG1[i]
    dif[difG1 > difG2] <- (tolerance*difG2)[difG1 > difG2]

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
    covariance <- as.matrix(bdiag(cov1,cov2))
    A1<- diag(indicatorNum)
    A2<- diag(indicatorNum)*(-1)
    A <- cbind(A1,A2)
    covPosterior<- A%*%covariance%*%t(A)
    # get sample size
    sampsizes <- lavInspect(lavAdjusted, what = "nobs")
    # fraction b
    b <- times*(indicatorNum+1)/(2*sampsizes)
    # get prior covariance
    covPrior1 <- cov1/b[1]
    covPrior2 <- cov2/b[2]
    covPrior <- as.matrix(bdiag(covPrior1,covPrior2))
    covPrior <- A%*%covPrior%*%t(A)
    ## compute Bayes factor
    set.seed(seed)
    BF <- (pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPosterior, sigma=covPosterior)[1]*(1- pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPrior, sigma=covPrior)[1]))/
      (pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPrior, sigma=covPrior)[1]*(1- pmvnorm(lower=lowerBound, upper = upperBound, mean=meanPosterior, sigma=covPosterior)[1]))
    BF <- as.data.frame(BF)
    colnames(BF) <- 'scalar'
    rownames(BF) <- 'Bayes factor'
    #BF <- round(BF,3)
    BF1 <- unlist(BF)
    if(is.nan(BF1)){
      cat('Approximate scalar invariance : Bayes factor cannot be computed. \n\n')
      print(BF)
    }
    else if(BF >= minBF){
      cat('When the minimum Bayes factor should  >=',minBF,',','approximate scalar invariance is supported. \n\n')
      print(BF)
    } else {
      #Part6: compute BF for each pair of the intercept
      eachBF <- list()
      for(i in 1:indicatorNum){
        # get each posterior mean
        eachMeanPosterior <- meanPosterior[i]
        # get each upper and lower bound
        eachUpperBound <- upperBound[i]
        eachLowerBound <- -eachUpperBound
        # get each posterior covariance
        eachCovPosterior <- covPosterior[i,]
        eachCovPosterior <- eachCovPosterior[i]
        eachCovPosterior <- matrix(eachCovPosterior,1,1)
        # get each prior covariance
        eachCovPrior <- covPrior[i,]
        eachCovPrior <- eachCovPrior[i]
        eachCovPrior <- matrix(eachCovPrior,1,1)
        # Bayes Factor for each loading
        set.seed(seed)
        eachBF[i] <- (pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=eachMeanPosterior, sigma=eachCovPosterior)[1]*(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]))/
          (pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]*(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=eachMeanPosterior, sigma=eachCovPosterior)[1])/(1-pmvnorm(lower=eachLowerBound, upper = eachUpperBound, mean=0, sigma=eachCovPrior)[1]))
        #eachBF[[i]] <- round(eachBF[[i]],3)
        # name BF with the indicator's name
        names(eachBF)[i] <- paste(indicatorName[i],  sep = "")
      } # end for(i in 1:indicatorNum)
      partBF <- BF
      PartEachBF <- unlist(eachBF)
      partIntG1 <- intG1
      partIntG2 <- intG2
      #Part7 : compute partial BF. If patial BF always < minBF,
      #        the loop will exit with length(PartEachBF) == 1.
      while(partBF < minBF & length(PartEachBF) != 1){
        ## recompute model
        # get partial intercepts
        partIntG1 <- partIntG1[- which.min(PartEachBF)]
        partIntG2 <- partIntG2[- which.min(PartEachBF)]
        # the number of the partial intercepts
        partNum <- length(partIntG2)
        # recompute factor's mean of group 1 and 2 that the sum of the partial intercepts per group = 0
        partMeanG1 <- sum(partIntG1)/sum(loadG1)
        partMeanG2 <- sum(partIntG2)/sum(loadG2)
        # adjust factor's mean
        adjustPartFactorM <- paste0('F', ' ~ ' ,'c(',partMeanG1,',',
                                    partMeanG2,')* ','1')
        # set model for lavaan
        adjustedModel <- c(Equation,adjustFactor, adjustPartFactorM)
        # run lavvan again
        lavPart <- cfa(adjustedModel, dat, std.lv = TRUE, group = groupVariable)
        ## get the adjusted intercepts
        # get the adjusted estimates
        estiPart <- parameterEstimates(lavPart)
        # get the adjusted intercepts
        partInt <- estiPart[estiPart$op == "~1", "est"]
        # get the adjusted intercepts of group 1 and 2
        AdjustedPartIntG1 <- partInt[1:indicatorNum]
        AdjustedPartIntG2 <- partInt[(indicatorNum + 2):(2*indicatorNum + 1)]
        ## get the upper bound
        partUpperBound <- upperBound
        ## get posterior and prior covariance
        # get covariance of group 1 and 2
        partCov1 <- lavInspect(lavPart, "vcov")[varName1, varName1]
        partCov2 <- lavInspect(lavPart, "vcov")[varName2, varName2]
        # posterior covariance
        partCov <- as.matrix(bdiag(partCov1,partCov2))
        partCovPosterior<- A%*%partCov%*%t(A)
        # prior covariance
        partCovPrior1 <- partCov1/b[1]
        partCovPrior2 <- partCov2/b[2]
        partCovPrior <- as.matrix(bdiag(partCovPrior1,partCovPrior2))
        partCovPrior <- A%*%partCovPrior%*%t(A)
        ## get partial intercepts, upper bound, and prior and posterior covarance
        PartEachBF2 <- unlist(eachBF)
        # the number of the deleted loadings
        DeletedNum <- indicatorNum - partNum
        for (j in 1:DeletedNum){
          # the partial loadings of group 1 and 2
          AdjustedPartIntG1 <- AdjustedPartIntG1[- which.min(PartEachBF2)]
          AdjustedPartIntG2 <- AdjustedPartIntG2[- which.min(PartEachBF2)]
          # the partial upper bound
          partUpperBound <- partUpperBound[- which.min(PartEachBF2)]
          # the partial posterior covariance
          partCovPosterior <- partCovPosterior[- which.min(PartEachBF2),- which.min(PartEachBF2)]
          # the partial prior covariance
          partCovPrior <- partCovPrior[- which.min(PartEachBF2),- which.min(PartEachBF2)]
          # control this loop process
          PartEachBF2 <- PartEachBF2[- which.min(PartEachBF2)]
        }
        # get partial posterior and prior means
        partMeanPosterior <- AdjustedPartIntG1 - AdjustedPartIntG2
        partMeanPrior <- rep(0, partNum)
        #  get the partial lower bound
        partLowerBound <- -partUpperBound
        ## compute partial BF
        set.seed(seed)
        partBF <- (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]))/
          (pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPrior, sigma=partCovPrior)[1]*(1-pmvnorm(lower=partLowerBound, upper = partUpperBound, mean=partMeanPosterior, sigma=partCovPosterior)[1]))
        # control the loop process
        PartEachBF <- PartEachBF[- which.min(PartEachBF)]
      }#end while
      ## arrange partial BF
      partBF <- as.data.frame(partBF)
      #partBF <- round(partBF,3)
      colnames(partBF) <- 'partial'
      # get deleted indicator
      deletedIndicator <- setdiff(names(eachBF), names(PartEachBF))
      ## print BF of all intercepts, BF of each intercepts, BF of partial intercepts
      ## arrange a table
      eachBF <- as.data.frame(eachBF)
      #eachBF <- round(eachBF,3)
      out <- cbind(BF,partBF,eachBF)
      # CJ: Do not print in this function
      return(out)
    }# end if (BF < minBF)
}

#' @method print bayesian_invariance
#' @export
print.bayesian_invariance <- function(x, ...){
  # CJ: DO all of the printing here

  # judge partial BF > minBF or not
  if(partBF >= minBF) {
    cat("\n With the minimum Bayes factor >=",minBF,",","partial approximate scalar invariance is supproted.
        \n The following indicators are free :\n",paste0(deletedIndicator,collapse = " , "),"\n\n")
  } else {
    # if partial BF < minBF, partial BF do not exist.
    partBF <- 'null'
    cat("\n With the minimum Bayes factor >=", minBF, ", intercepts are not
          \n about equal across group, respectively.Please check your model or minBF! \n\n")
  }
}

# section 1: example
# CFA model
model22 <- 'A  =~ Ab + Al + Af + An + Ar + Ac'
# Fit the multiple group latent regression model
lavout <- cfa(model22, data = sesamesim, std.lv = TRUE, group = "sex")
# approximate metric invariance
set.seed(2020)
result1 <- BMI1(lavout,type = 'metric',data = sesamesim, tolerance = 0.2,minBF = 10,times=2)
# partial approximate metric invariance
result2 <- BMI1(lavout,type = 'metric',data = sesamesim, tolerance = 0.2,minBF = 65,times=2)
# approximate scalar invariance
result3 <- BMI1(lavout,type = 'scalar',data = sesamesim, tolerance = 0.2,minBF = 3,times=2)
# partial approximate scalar invariance
result4 <- BMI1(lavout,type = 'scalar',data = sesamesim, tolerance = 0.2,minBF = 10,times=2)
