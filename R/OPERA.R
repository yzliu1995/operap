#' Calculating the row sums
#'
#' This function calculates the row sums of a given input, which can be a number, a vector, or a matrix. Before performing the calculation, the function checks whether the input is a vector or a matrix. If it's a vector, the function treats it as a single column matrix and does nothing. If it's a matrix, the function calculates the row sums.
#'
#' @param x The input, which can be a number, a vector, or a matrix.
#' @param ... see [base::rowSums()] for details
#' @return A vector of sums.
#'
#' @examples
#' RowSums(array(1:60, c(3,4,5)), dims = 1L)
#' # Output: 590 610 630
#'
#' @export
#'
RowSums <- function(x, ...){
  # checks if it is a matrix
  # if it is not a matrix, the function converts it to a matrix
  if (!is.matrix(x) && !is.array(x))
    x <- as.matrix(x)
  # checks if it is a column vector
  # if yes, no further operations are needed
  if (dim(x)[2] == 1L && length(dim(x))==2)
    return(x)
  return(rowSums(x, ...))
}


#' Calculating the mode
#'
#' This function takes a vector as input and calculates the mode, which is the most frequently occurring value in the vector. The function achieves this by utilizing two auxiliary functions: [base::tabulate()], which calculates the frequency of each unique value in the vector, and [base::match()], which returns the position of each unique value when matched.
#'
#' @param x A vector.
#'
#' @return The mode of the input vector.
#'
#' @examples
#' Mode(c(1, 2, 3, 2, 4, 2))
#' # Output: The mode is 2
#'
#' @seealso
#'
#' [base::tabulate()]: Calculates the frequency for each unique value in a vector.
#'
#' [base::match()]: Returns the position of each unique value when matched in a vector.
#'
#' Other functions for finding statistical measures: [base::mean()], [base::min()], [base::max()], [stats::median()]
#'
#' @export
#'
Mode <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

#' Getting all the edges and the subsequent vertices of the current vertex
#'
#' This function takes a starting vertex, which represents the current levels for all ordinal variables or risk factors, and a vector specifying the number of categories in each ordinal variable or risk factor. It then generates the edges for the current vertex, which represents the ordered relationships between vertices. The function returns a list that contains a vector of all edges of the current vertex and a matrix of the next vertices of the current vertex.
#'
#' @param v The starting vertex representing the current levels for all ordinal variables or risk factors.
#' @param ncat A vector specifying the number of categories in each variable or risk factor.
#'
#' @return A list containing a vector of all the edges of the current vertex and a matrix of the next vertices of the current vertex.
#'
#' @examples
#' nextEdges(v = c(1, 2, 1), ncat = c(2, 3, 2))
#'
#' @seealso
#'
#' \code{\link{edgesHasse}}
#'
#' @export
#'
nextEdges <- function(v, ncat){
  # variable names
  l <- letters[1:length(ncat)]
  currentVertex <- paste0(l, v, collapse='')
  # edges
  e <- base::rep(NA, sum(v < ncat)*2)
  # next vertices
  vs <- matrix(base::rep(NA, sum(v < ncat)*length(ncat)), nrow = sum(v < ncat), ncol = length(ncat))
  j = 1
  # iterates over each variable
  for(i in 1:length(ncat)){
    v_new <- v
    # if the current variable of the current vertex has not reached the max level
    if(v[i] < ncat[i]){
      # moves one level up for the current variable
      v_new[i] <-  v[i] + 1
      e[2*j-1] <- currentVertex
      nextVertex <- paste0(l, v_new, collapse='')
      vs[j,] <- v_new
      e[2*j] <- nextVertex
      j = j + 1
    }
  }
  return(list(e, vs))
}

#' Generating the edges
#'
#' The function generates edges in the Hasse diagram based on the given inputs. It then returns a vector of edges representing the ordered relationships between vertices in the Hasse diagram.
#'
#' @param ncat A vector specifying the number of categories in each ordinal variable or risk factor with a total ordering.
#' @param e The current edges in the Hasse diagram, starting from an empty vector.
#' @param vs The current vertex.  When it is null, it starts from (1, 1, ..., 1) depending on the number of variables.
#'
#' @return A vector of edges representing the connections between vertices in the Hasse diagram.
#'
#' @seealso
#'
#' \code{\link{nextEdges}}
#'
#' @examples
#' edgesHasse(ncat = c(2, 3, 2), e = vector())
#'
#' @export
#'
edgesHasse <- function(ncat, e, vs = NULL){
  # The current vertex, starting from (1, 1, ..., 1) depending on the number of variables
  if(is.null(vs)){
    vs <- matrix(base::rep(1, length(ncat)), 1, length(ncat))
  }
  temp <- c()
  for(i in 1:length(ncat)){
    # the total number of edges
    temp[i]  <- (ncat[i] - 1) * prod(ncat[-i])
  }
  # whether all edges have been found
  if(length(e) == sum(temp) * 2){
    return(e)
  }else{
    # saves new vertices
    vs_new <- c()
    # goes through all the vertices
    for(i in 1:dim(vs)[1]){
      if(sum(vs[i,]) < prod(ncat)){
        # finds the next edges for the current vertex
        r <- nextEdges(v = vs[i,], ncat = ncat)
        # returns the next vertices and all edges of the current vertex
        e  <- append(e, r[[1]], after = length(e))
        # updates the vertices
        vs_new <- rbind(vs_new, r[[2]])
      }
    }
    # removes any duplicated vertices
    vs_new <- vs_new[!duplicated(vs_new),]
    if(length(vs_new) > 0){
      # recursively searches the edges for the next vertices
      edgesHasse(ncat = ncat, e = e, vs = vs_new)
    }
  }
}

#' Generating the partial ordering matrix
#'
#' The function performs a transformation on the Hasse diagram to generate a matrix that expresses the partial ordering constraints on coefficients. The resulting matrix represents the constraints imposed by the partial ordering relations.
#'
#' @param Hasse The partial ordering relations expressed by a directed graph.
#' @param cellName The labels of the nodes in the Hasse diagram.
#'
#' @return A matrix expressing the partial ordering constraints on coefficients.
#'
#' @importFrom igraph as_edgelist ecount
#' @export
#'
as_inequalities <- function(Hasse, cellName){

  PO <- matrix(0, nrow = ecount(Hasse), ncol = length(cellName))
  colnames(PO) <- cellName
  edgelist <- as_edgelist(Hasse)
  for (i in 1:ecount(Hasse)){
    PO[i,edgelist[i,1]] <- -1
    PO[i,edgelist[i,2]] <- 1
  }
  return(PO)
}

#' Simulating a dataset
#'
#' The function generates a data frame representing one dataset based on the provided inputs.
#'
#' @param beta.mat The betas or specified coefficients in matrix form for all nodes.
#' @param p.mat The probabilities of a patient being sampled from nodes in matrix form.
#' @param n The total number of patients.
#' @param seed The seed used for reproducibility.
#' @param cen.prob The proportion of censored patients.
#' @param label The labels of nodes.
#' @param ordered Whether the dataset should be sorted by survival time in ascending order.
#' @param type The outcome type, which can be either survival outcome or binary outcome.
#' @param firstCon Whether the first variable is continuous.
#'
#'
#' @return A data frame representing one dataset.
#'
#' @importFrom reticulate array_reshape
#' @importFrom stats rexp rbinom runif pchisq
#' @import dplyr
#'
#'
#' @export
#'
dataset.gen <- function(beta.mat, p.mat,
                        n = 800, seed = 2021, cen.prob = 0.8,
                        label, ordered = TRUE, type = "surv", firstCon = FALSE){
  # sets seed for reproducibility
  set.seed(seed = seed)
  # the number of levels for each ordinal variable
  p <- dim(p.mat)
  # the probabilities of a patient being sampled from nodes
  pN <- c()
  for(k in 1:length(label)){
    pN <- c(pN, p.mat[t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  label[k]), split = "")[[1]])))])
  }

  idx <- sample(1:(prod(p)), n, replace = TRUE, prob = pN)
  samples <- label[idx]
  vars <- letters[1:length(p)]

  # gets all the unique coefficients from the smallest to the largest
  uniqueBetas <- sort(unique(as.numeric(array_reshape(beta.mat, dim = c(1, prod(dim(beta.mat)))))))
  # gets all the levels of all the variables for each observation
  values <- as.matrix(do.call(rbind.data.frame,lapply(samples, function(s){t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  s), split = "")[[1]])))})))
  # gets the true stages for all observations
  stages <- sapply(samples, function(s){which(unique(uniqueBetas)  == beta.mat[t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  s), split = "")[[1]])))])-1})
  # gets all the variables specifying all the levels for the risk factors for all observations
  variables <- apply(values, 1, function(v){paste0(vars, v, collapse = "")})
  # survival outcome
  if(type == "surv"){
    # failure times that follow exponential distribution with rate \lambda (i.e., mean 1/\lambda)(x) = \lambda exp{-\lambda x}
    y <- rexp(n, rate = exp(beta.mat[values]))

    eta <-  seq(from = min(beta.mat)/2 , to = max(beta.mat)*2 , length.out = 200)

    targetdelta <- base::rep(0, n)
    targety <- y

    for (j in 1:length(eta)){
      # exponential censoring time
      # C <- rexp(n, rate = exp(eta[j]))
      C <- rexp(n, rate = exp(eta[j])*(1 - 1/(stages+2)))
      tilda_y <- pmin(C,y)

      # y > C indicates censoring
      delta <- as.numeric(y<=C)
      # checks censoring rate
      # 1 - mean(delta) indicates current censoring rate
      # 1 - mean(targetdelta) indicates the closest censoring rate to the target censoring rate so far
      if (abs(1 - mean(delta) - cen.prob) < abs(1 - mean(targetdelta) - cen.prob)) {
        targety <- tilda_y
        targetdelta <- delta
      }
    }

    dat <- cbind.data.frame(values, targety, targetdelta, stages, variables)
    # adds column names to values
    colnames(dat) <- c(vars, "time", "status", "tStage", "variable")

    finalD <- data.frame(dat)

    # orders the dataset by time
    if(ordered){
      finalD <- finalD[order(finalD$time),]
    }

  }

  # binary outcome
  if(type == "bin"){
    # gets the outcomes following the binomial distribution with p = 1/(1+exp(-\beta))
    y <- rbinom(n, 1, prob = 1/(1 + exp(-beta.mat[values])))
    # gets the failure rate in each node
    falRate <- as.data.frame.matrix(table(variables, y))
    falRate <- falRate %>% mutate(prop  = .$`1`/(.$`1`  + .$`0`))
    # assures the extreme failure rates satisfying 0.1 <= p_i <= 0.9
    while(sum(falRate$prop < 0.1) + sum(falRate$prop > 0.9) > 0){
      y <- rbinom(n, 1, prob = 1/(1 + exp(-beta.mat[values])))
      falRate <- as.data.frame.matrix(table(variables, y))
      falRate <- falRate %>% mutate(prop  = .$`1`/(.$`1`  + .$`0`))
    }
    dat <- cbind.data.frame(values, y, stages, variables)
    # adds column names to values
    colnames(dat) <- c(vars, "outcome", "tStage", "variable")

    finalD <- data.frame(dat)

    # orders the dataset by outcome, which is not actually necessary
    if(ordered){
      finalD <- finalD[order(finalD$outcome),]
    }
  }

  if(firstCon){
    finalD <- finalD %>% mutate(a = ifelse(replicate(nrow(finalD), runif(1)) <= 0.5, 2*a - 1, 2*a))
    finalD <- finalD %>% mutate(variable = apply(finalD[,1:length(vars), drop = F], 1, function(v){paste0(vars, v, collapse = "")}))
  }

  return(finalD)
}

#' Generating the design matrix for risk categories
#'
#' The function generates a matrix `Z` that indicates the risk categories or nodes to which each observation in the dataset belongs. Each row of the matrix corresponds to an observation, and each column represents a node. The value in each cell of the matrix indicates whether the observation belongs to the corresponding node (1 if yes, 0 if no).
#'
#' @param dat The dataset that includes a `variable` column specifying each observation's node.
#' @param colNames The node names.
#'
#' @return A matrix `Z` representing the node to which each observation resides.
#'
#' @export
#'
Zmatrix <- function(dat, colNames){
  n <- nrow(dat)
  Z <- matrix(0, nrow = n, ncol = length(colNames))
  colnames(Z) <- colNames
  Z[cbind(1:nrow(Z),  match(dat$variable, colNames))] <- 1
  return(Z)
}

#' Generating a list regarding color assignments and true stages
#'
#' The function generates a list that contains the following elements: colors for nodes, unique coefficients, the color palette, and the true stage for nodes.
#'
#' @param beta.mat The betas in matrix form.
#' @param label The labels of nodes.
#' @param c The basic color palette.
#'
#'
#' @return A list containing colors for nodes, unique coefficients, the color palette, and the true stage for nodes.
#'
#' @importFrom colorspace sequential_hcl
#'
#' @export
#'
generatePalette <- function(beta.mat, label, c = "Viridis"){

  uniqueBetas <- unique(as.numeric(array_reshape(beta.mat, dim = c(1, prod(dim(beta.mat))))))
  pal <- sequential_hcl(palette = c, n = length(uniqueBetas), rev = TRUE)
  uniqueS <- 1:length(uniqueBetas)
  value <- as.matrix(do.call(rbind.data.frame,lapply(label, function(l){t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  l), split = "")[[1]])))})))
  nodePal <- sapply(label, function(l){pal[which(uniqueBetas == beta.mat[t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  l), split = "")[[1]])))])]})
  trueS <- sapply(label, function(l){uniqueS[which(uniqueBetas == beta.mat[t(cbind(as.numeric(strsplit(gsub("[^0-9.]", "",  l), split = "")[[1]])))])]})
  names(trueS) <- label

  return(list(nodePal, uniqueBetas, pal, trueS))
}

#' Generating a dummy matrix for a categorical variable
#'
#' The function generates a dummy matrix where each column represents a category or level of the categorical variable. The matrix contains binary values (0 or 1) indicating whether each observation falls into a particular category. If a reference level is specified, the dummy matrix will include one less column compared to the total number of levels.
#'
#' @param X The categorical variable provided as a vector.
#' @param ref Whether a reference level is needed. If set to `TRUE`, the first level is chosen as the reference level.
#'
#' @return A dummy matrix representing the categorical variable.
#'
#' @export
#'
dM <- function(X, ref = TRUE){
  if(length(unique(X)) == 1){
    stop("Only one category!")
  }
  M <- matrix(0, nrow = length(X), ncol = length(levels(as.factor(X))))
  M[cbind(1:nrow(M), match(as.factor(X), levels(as.factor(X))))] <- 1
  if(!ref){
    return(M)
  }else{
    # chooses the first level as the reference level
    m <- M[,2:ncol(M)]
    if(is.null(dim(m))){
      return(t(t(m)))
    }else{
      return(m)
    }
  }
}

#' Generating a list of dummy variables
#'
#' The function generates a list where each element corresponds to a categorical variable in the data matrix. Each element of the list is a dummy variable, represented as a binary vector or column, indicating the presence or absence of each category or level of the categorical variable.
#'
#' @param dat The data matrix containing categorical variables.
#'
#' @return A list of dummy variables, where each element represents a categorical variable and is a binary vector indicating the presence or absence of each category or level.
#'
#' @export
#'
listDummies <- function(dat){
  l <- split(dat, col(dat))
  names(l) <- paste0("x", 1:ncol(dat))
  return(l)
}

#' Calculating the edge misclassification rate
#'
#' The function computes the edge misclassification rate, which measures the accuracy of the estimated edges compared to the true edges. It quantifies the proportion of misclassified edges in relation to the total number of edges.
#'
#' @param edgeMatrix A vector of partial ordering relationship node variables.
#' @param trueS The true stages.
#' @param estimateS The estimated stages.
#'
#' @return The edge misclassification rate.
#'
#' @export
#'
edgeMisclassfication <- function(edgeMatrix, trueS, estimateS){
  tM <- edgeMatrix %*% t(t(trueS))
  eM <- edgeMatrix %*% t(t(estimateS))
  return(1 - round(sum(tM - eM == 0)/length(eM), 3))
}

#' Calculating the Kaplan-Meier survival probabilities
#'
#' The function computes the K-M survival probabilities for the specified `times`. The K-M survival probabilities estimate the probability of survival at each time point based on the observed survival data.
#'
#' @param obj A `survfit` object representing the survival data.
#' @param times A vector of time points at which the survival probabilities are to be computed.
#'
#' @return The K-M survival probabilities for the specified `times`.
#'
summaryC <- function(obj, times){
  # gets the survival table
  survTable <- summary(obj, times = times)
  survProbs <- survTable$surv
  if(length(survProbs) != length(times)){
    survProbs <- c(survProbs, rep(survProbs[length(survProbs)], length(times) - length(survProbs)))
  }
  return(survProbs)
}

#' Calculating the integrated Brier score
#'
#' The function computes the integrated Brier score, which is a measure of the accuracy of survival predictions. It assesses the agreement between the observed survival outcomes and the predicted survival probabilities over a range of time points. The integrated Brier score provides an overall evaluation of the predictive performance.
#'
#' @param obj An object of class `Surv` representing the survival data.
#' @param pred A list of `survfit` objects representing the predicted survival probabilities.
#'
#' @return The integrated Brier score and the starting and ending time points.
#'
#' @importFrom prodlim prodlim
#' @importFrom survival Surv
#' @importFrom stats predict
#'
#' @references
#' Graf, E., Schmoor, C., Sauerbrei, W., & Schumacher, M. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statistics in Medicine, 18(17-18), 2529-2545. Retrieved from \href{https://rdrr.io/cran/ipred/src/R/sbrier.R}{here}
#'
#' @export
#'
sbrierC <- function(obj, pred){

  btime = range(obj[,1])
  # gets # of times
  N <- nrow(obj)
  # gets the times and censoring of the data, order them w.r.t. time
  time <- obj[,1]
  ot <- order(time)
  cens <- obj[ot,2]
  time <- time[ot]

  # gets the times to compute the integrated Brier score over
  btime <- time[time >= btime[1] & time <= btime[2]]
  pred <- pred[ot]
  survs <-  matrix(unlist(lapply(pred, summaryC, times = btime)),
                   nrow=length(btime), ncol=N)

  # estimates censoring distribution
  # deals with ties
  hatcdist <- prodlim::prodlim(Surv(time, cens) ~ 1, reverse = TRUE)
  csurv <- predict(hatcdist, times = time, type = "surv")
  csurv[csurv == 0] <- Inf

  # conditional survival for new time points
  csurv_btime <- predict(hatcdist, times = btime, type = "surv")
  csurv_btime[is.na(csurv_btime)] <- min(csurv_btime, na.rm = TRUE)
  csurv_btime[csurv_btime == 0] <- Inf

  bsc <- rep(0, length(btime))

  # computes the integrated Brier score
  for (j in 1:length(btime)) {
    help1 <- as.integer(time <= btime[j] & cens == 1)
    help2 <- as.integer(time > btime[j])
    bsc[j] <-  mean((0 - survs[j,])^2*help1*(1/csurv) +
                      (1-survs[j,])^2*help2*(1/csurv_btime[j]))
  }

  # applies trapezoid rule
  idx <- 2:length(btime)
  RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx]) / 2)
  RET <- RET / diff(range(btime))

  names(RET) <- "integrated Brier score"
  attr(RET, "time") <- range(btime)
  return(RET)
}

#' Getting the constrained estimates for cancer stages
#'
#' Finds the estimates for cancer stages given the partial ordering constraints.
#'
#' @param coeff0 A vector of initial coefficients for nodes.
#' @param theta0 A vector of initial coefficients for covariates that need adjusting.
#' @param cen A vector of censoring statuses, or the binary outcome.
#' @param y A vector of outcomes.
#' @param stage A vector of stages after initial classification.
#' @param Cov A matrix of covariates.
#' @param withCov Whether any covariates need adjusting.
#' @param maxiter The maximum iteration in each approximation step.
#' @param eps A numeric value denoting the accuracy level.
#' @param getLikelihood Whether it outputs the value of the negative likelihood.
#' @param fullModel Whether it is a full model, representing the model before pruning to the current model.
#' @param fullZ A matrix of stages for the full model.
#' @param type Outcome type, either survival outcome or binary outcome.
#' @param seed Seed used for generating the simulation dataset.
#'
#' @return The estimates for cancer stages given the partial ordering constraints.
#'
#' @import Rcpp quadprog
#'
#' @export
#'
cffStg <- function(coeff0, theta0, cen, y, stage, Cov,
                   withCov = FALSE, maxiter = 20, eps = 10^-4,
                   getLikelihood = FALSE, fullModel = FALSE,
                   fullZ = NULL, type = "surv", seed = 0){

  nIter <- 0

  if(!fullModel){

    numStages <- length(unique(stage))

    if(numStages > 1){
      fullCnstrn <- cbind(0, diag(numStages-1)) - diag(numStages)[1:(numStages-1),]
    }else{
      fullCnstrn <- diag(1)
    }

    coeff <- coeff0
    Z <- as.matrix(data.frame(listDummies(dM(stage, ref = F))))


  }else{

    numStages <- length(unique(stage))

    if(numStages > 1){
      fullCnstrn <- cbind(0, diag(numStages-1)) - diag(numStages)[1:(numStages-1),]
      fullCnstrn <- cbind(fullCnstrn, 0)
    }else{
      fullCnstrn <- diag(1)
      fullCnstrn <- cbind(fullCnstrn, 0)
    }
    coeff <- coeff0
    Z <- fullZ

  }

  if(withCov){
    Lambda <- rep(0, length(coeff0)+length(theta0))
  }else{
    Lambda <- rep(0, length(coeff0))
  }

  r <- IRLSA(coeff0, theta0, cen, y, Z,
             Cov, fullCnstrn, Lambda, eps, maxiter, seed,
             withCov, type)

  k <- 1

  while(wInf(r[[1]]) & (k <= 30)){

    k <- k + 1
    r <- IRLSA(coeff0/k + eps, theta0/k + eps, cen, y, Z,
               Cov, fullCnstrn, Lambda, eps, maxiter, seed,
               withCov, type)

  }

  coeff <- r[[1]]
  names(coeff) <- colnames(Z)

  if(withCov){
    theta <- r[[2]]
    names(theta) <- colnames(Cov)
  }else{
    theta <- 0
  }

  if(!getLikelihood){
    if(withCov){
      return(list(coeff, theta))
    }else{
      return(coeff)
    }
  }else{
    if(withCov){
      return(list(logLK(Z, cen, coeff, y, Cov, theta, type, withCov), coeff, theta))
    }else{
      return(list(logLK(Z, cen, coeff, y, Cov, theta, type, withCov), coeff))
    }
  }
}

#' Fitting OPERA without pruning
#'
#' This function implements OPERA without pruning
#'
#' @param Data The dataset
#' @param TimeN The variable name for survival times
#' @param cenN The variable name for censoring times
#' @param yN The variable name for binary outcomes
#' @param covN The variable names for covariates
#' @param Z A matrix representing the node to which each observation resides.
#' @param Cnstrn A matrix of constraints
#' @param withCov Whether any covariates need adjusting
#' @param maxiter Maximum number of iterations in each approximation step
#' @param eps A numeric accuracy level
#' @param AIC Whether AIC is used
#' @param sPoint A numeric value for an initial coefficient
#' @param type A string of outcome type, either survival ("surv") or binary ("bin") outcome
#' @param seed Seed used for generating the simulation dataset
#' @param GIC The penalty term used in the criterion to identify the optimal stage
#' @param N The size of the search space for tuning parameter
#' @param preS Whether the coefficients from the regression model are used as the initial coefficients
#' @param perc The minimum proportion of patients in each stage
#' @param minObs The minimum number of patients in each stage
#'
#' @return The initial staging result from OPERA without pruning
#'
#' @importFrom survival coxph strata
#' @importFrom stats as.formula glm xtabs binomial
#'
#' @importFrom Matrix bdiag
#'
#' @export
#'
operai <- function(Data, TimeN, cenN, yN, covN, Z, Cnstrn,
                   withCov = F, maxiter = 10, eps = 10^-4,
                   AIC = TRUE, sPoint = 0.001, type = "surv",
                   seed = 0, GIC = 2, N = NULL, preS = T, perc = 0.1, minObs = 30, ...){

  if(AIC){
    GIC <- 2
  }

  Zp <- Z
  Cnstrnp <- Cnstrn

  if(withCov){
    Cov  <- as.matrix(Data[, covN, drop = F])
  }else{
    Cov <- matrix(0, 1, 1)
  }


  if(type == "surv"){
    Times <- Data[, TimeN]
    cen <- Data[, cenN]
    y <- Data[, TimeN]

  }

  if(type == "bin"){
    cen <- Data[, yN]
    y <- Data[, yN]

  }

  if(withCov){
    if(type == "surv"){
      suppressWarnings(betas <- coxph(as.formula(paste0("Surv(time, status) ~ ", paste(colnames(Cov), collapse = "+"), "+variable")), data = Data)$coef)
    }

    if(type == "bin"){
      suppressWarnings(betas <- glm(as.formula(paste0("outcome ~", paste(colnames(Cov),collapse = "+"), "+variable-1")), data = Data, family = binomial)$coef)
    }

    betas[is.na(betas)] <- sPoint

    preTheta <- betas[1:(ncol(Cov))]
    preTheta[is.na(preTheta)] <- sPoint
    names(preTheta) <- colnames(Cov)

    coeffs <- betas[(ncol(Cov)+1):length(betas)]

    names(coeffs) <- gsub("variable", "", names(coeffs))

    preCoeff <- c(sPoint, sapply(colnames(Z), function(x) ifelse(x %in% names(coeffs), coeffs[x], sPoint)))

    names(preCoeff)[1] <- "mu"
  }else{

    if(type == "surv"){
      suppressWarnings(coeffs <- coxph(Surv(time, status) ~ variable, data = Data)$coef)
    }

    if(type == "bin"){
      suppressWarnings(coeffs <- glm(outcome ~ variable - 1, family = binomial, data = Data)$coef)
    }

    coeffs[is.na(coeffs)] <- sPoint

    names(coeffs) <- gsub("variable", "", names(coeffs))

    preCoeff <- c(sPoint, sapply(colnames(Z), function(x) ifelse(x %in% names(coeffs), coeffs[x], sPoint)))

    names(preCoeff)[1] <- "mu"
    preTheta <- 0
  }

  coeff_p <- preCoeff
  theta_p <- preTheta

  if(withCov){
    covName <- colnames(Cov)
    if(preS){
      theta <- theta_p
    }else{
      theta <- base::rep(sPoint, length(covName))
    }
    names(theta) <- covName
  }else{
    theta <- 0
  }

  # # of nodes
  numCoeff <- dim(Z)[2]
  # nodes
  stagingV <- colnames(Z)
  # saves final classification results
  result <- base::rep(NA, numCoeff)
  names(result) <- stagingV
  # adds mu to Z
  Z <- cbind(1,Z)
  colnames(Z)[1] <- "mu"
  # nodes plus mu
  allnames <- colnames(Z)

  S <- 1
  # initializes coefficients
  if(preS){
    coeff <- coeff_p
  }else{
    coeff <- c(sPoint, base::rep(sPoint, numCoeff))
  }

  names(coeff) <- allnames

  # adds the column corresponding to mu into the matrix representing the partial ordering constraints
  nOrd <- cbind(0, Cnstrn)
  colnames(nOrd) <- allnames
  # assures all node coefficients non-negative
  fullCnstrn <- rbind(nOrd, cbind(0, diag(1, numCoeff)))

  # saves the difference in coefficients between the last and the second last iteration
  allCs <- c()

  while(sum(is.na(result)) > 0){

    best_coeff <- coeff
    best_theta <- theta

    best_abic <- NULL

    # finds the upper bound for the tuning parameter lambda

    # starts from 1
    lambda_now <- 1
    lambda_max <- lambda_now

    coeff_tmp <- best_coeff

    if(withCov){
      coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
    }else{
      coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
    }

    # IRWLS
    resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                        Z, Cov, fullCnstrn, coeffLambda,
                        eps, maxiter, seed, withCov, type)

    if(is.infinite(resultIRLS[[length(resultIRLS)]])){
      coeff_tmp <- best_coeff
      break
    }

    coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)

    k1 <- 0
    # as long as at least one coefficient for a node is not equal to 0
    # finds the lambda_max that shrinks all the coefficients for all the nodes to zero
    while((sum(coeff_tmp[2:(numCoeff+1)]!=0) > 0) & (k1 <= 30)){

      if(k1 == 0){
        coeff_tmp <- best_coeff
      }

      # updates lambda_now
      lambda_now <- lambda_now*2
      # places lambda's only on the coefficients of nodes/risk factors
      if(withCov){
        coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
      }else{
        coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
      }

      # IRWLS
      resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                          Z, Cov, fullCnstrn, coeffLambda,
                          eps, maxiter, seed, withCov, type)

      if(is.infinite(resultIRLS[[length(resultIRLS)]])){
        coeff_tmp <- best_coeff
        break
      }

      # rounding is necessary since most values are not exactly zeros but close to zero
      coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)
      k1 <- k1 + 1
    }

    # if 1 is big enough to shrink all node coefficients to zero
    if(k1 >= 1){
      lambda_max <- lambda_now
    }else{

      # finds lambda_max by taking a half of lambda each time
      lambda_now <- 1
      lambda_max <- lambda_now

      k2 <- 0
      # as long as all node coefficients equal 0
      while((sum(coeff_tmp[2:(numCoeff+1)]!=0) == 0) & (k2 <= 30)){

        if(k2 == 0){
          coeff_tmp <- best_coeff
        }



        # updates lambda_now
        lambda_now <- lambda_now/2

        if(withCov){
          coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
        }else{
          coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
        }

        # IRWLS
        resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                            Z, Cov, fullCnstrn, coeffLambda,
                            eps, maxiter, seed, withCov, type)

        if(is.infinite(resultIRLS[[length(resultIRLS)]])){
          coeff_tmp <- best_coeff
          break
        }

        coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)


        k2 <- k2 + 1
      }

      if(k2 >= 1){
        lambda_max <- lambda_now*2
      }
    }

    coeff_tmp <- best_coeff

    if(withCov){
      coeffLambda <- c(0, base::rep(0, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
    }else{
      coeffLambda <- c(0, base::rep(0, numCoeff), base::rep(0, S-1))
    }

    resultIRLS <-  IRLSA(coeff_tmp, theta, cen, y,
          Z, Cov, fullCnstrn, coeffLambda,
          eps, maxiter, seed, withCov, type)

    coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)

    numberZeros <- sum(coeff_tmp[2:(numCoeff+1)] == 0)

    # finds lambda_min by doubling

    # starts from 1
    lambda_now <- 1
    lambda_min <- lambda_now

    if(withCov){
      coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
    }else{
      coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
    }

    # IRLS
    resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                        Z, Cov, fullCnstrn, coeffLambda,
                        eps, maxiter, seed, withCov, type)

    if(is.infinite(resultIRLS[[length(resultIRLS)]])){
      coeff_tmp <- best_coeff
      break
    }

    coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)

    k3 <- 0
    # as long as no node coefficients equals zero
    while((sum(coeff_tmp[2:(numCoeff+1)] == 0) <= numberZeros) & (k3 <= 30)){

      if(k3 == 0){
        coeff_tmp <- best_coeff
      }

      # updates lambda_now
      lambda_now <- lambda_now*2

      if(withCov){
        coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
      }else{
        coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
      }

      # IRLS
      resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                          Z, Cov, fullCnstrn, coeffLambda,
                          eps, maxiter, seed, withCov, type)

      if(is.infinite(resultIRLS[[length(resultIRLS)]])){
        coeff_tmp <- best_coeff
        break
      }

      coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)

      k3 <- k3 + 1

    }

    if(k3 >= 1){
      lambda_min <- lambda_now/2
    }else{

      lambda_now <- 1
      lambda_min <- lambda_now

      k4 <- 0
      # as long as at least one node coefficient equals 0
      # as long as # of zero coefficients is more than when lambda equals zero
      while((sum(coeff_tmp[2:(numCoeff+1)] == 0) > numberZeros) & (k4 <= 30)){
        if(k4 == 0){
          coeff_tmp <- best_coeff
        }

        # updates lambda_now
        lambda_now <- lambda_now/2

        if(withCov){
          coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1), base::rep(0,length(covName)))
        }else{
          coeffLambda <- c(0, base::rep(lambda_now, numCoeff), base::rep(0, S-1))
        }

        # IRLS
        resultIRLS <- IRLSA(coeff_tmp, theta, cen, y,
                            Z, Cov, fullCnstrn, coeffLambda,
                            eps, maxiter, seed, withCov, type)

        if(is.infinite(resultIRLS[[length(resultIRLS)]])){
          coeff_tmp <- best_coeff
          break
        }

        coeff_tmp <- round(resultIRLS[[1]], -log10(eps)-1)

        k4 <- k4 + 1

      }
      if(k4 >= 1){
        lambda_min <- lambda_now
      }
    }

    if(is.null(N)){
      Np <- 2*length(coeff_tmp)
    }else{
      Np <- N
    }

    # finds the lambda with the highest AIC or BIC
    if(lambda_min < lambda_max){
      lambda_val <- c(0, gs(lambda_min, lambda_max, Np))
    }else{
      lambda_val <- c(0, lambda_min)
    }

    coeff <- best_coeff

    for (i in 1:length(lambda_val)){

      if(withCov){
        lambdaCoeff <- c(0, base::rep(lambda_val[i],numCoeff), base::rep(0,S-1), base::rep(0,length(covName)))
      }else{
        lambdaCoeff <- c(0, base::rep(lambda_val[i],numCoeff), base::rep(0,S-1))
      }


      # IRLS
      resultIRLS <- IRLSA(coeff, theta, cen, y,
                          Z, Cov, fullCnstrn, lambdaCoeff,
                          eps, maxiter, seed, withCov, type)

      # updates the best coefficients
      if(!is.infinite(resultIRLS[[length(resultIRLS)]])){
        coeff <- resultIRLS[[1]]
        names(coeff) <- colnames(Z)
        coeff <- coeff - min(coeff[2:(2+numCoeff-1)])
        coeff <- round(coeff, -log10(eps)-1)

        if(withCov){
          theta <- resultIRLS[[2]]
          theta <- round(theta, -log10(eps)-1)
          names(theta) <- colnames(Cov)
          abic_now <- 2*logLK(Z = Z, cen = cen, beta = coeff, y = y, withCov = withCov, cov = Cov, theta = theta, type = type) + GIC*sum(coeff[2:(numCoeff+1)] != 0) + GIC * length(theta)

        }else{
          abic_now <- 2*logLK(Z = Z, cen = cen, beta = coeff, y = y, withCov = withCov, cov = Cov, theta = theta, type = type) + GIC*sum(coeff[2:(numCoeff+1)] != 0)
        }

        if(is.null(best_abic)){
          best_abic <- abic_now
          best_coeff <- coeff
          if(withCov){
            best_theta <- theta
            best_C <-  resultIRLS[[3]]
          }else{
            best_C <-  resultIRLS[[2]]
          }
          best_lambda <- lambda_val[i]
        }else{
          if (abic_now < best_abic){
            best_abic <- abic_now
            best_coeff <- coeff
            if(withCov){
              best_theta <- theta
              best_C <-  resultIRLS[[3]]
            }else{
              best_C <-  resultIRLS[[2]]
            }
            best_lambda <- lambda_val[i]
          }
        }


      }

    }

    allCs <- c(allCs, best_C)
    print(paste0("Find the stage ", S))
    coeff  <- best_coeff
    print(coeff)
    print(paste0('The best lambda = ', round(best_lambda, -log10(eps))))

    if(withCov){
      theta <- best_theta
      print(theta)
    }

    # updates constraints
    result[names(which(coeff[2:(numCoeff+1)]==0))] <- S
    newnOrd <- cbind(nOrd[, c('mu',names(which(coeff[2:(numCoeff+1+S-1)]!=0)))],
                     RowSums(nOrd[,names(which(coeff[2:(numCoeff+1)]==0))]))

    colnames(newnOrd)[dim(newnOrd)[2]] = paste0('S',S)
    newnOrd <- as.matrix(dplyr::distinct(data.frame(newnOrd)))
    newnOrd <- newnOrd[RowSums(newnOrd != 0) != 0, ]
    nOrd <- newnOrd

    newZ <- cbind(Z[, c('mu',names(which(coeff[2:(numCoeff+1+S-1)]!=0)))],
                  RowSums(Z[,names(which(coeff[2:(numCoeff+1)]==0))]))
    colnames(newZ)[dim(newZ)[2]] = paste0('S',S)
    Z <- newZ

    SOrd <- matrix(0,S,S)

    for (i in 1:S-1){
      SOrd[i,i:(i+1)] <- c(-1,1)
    }
    SOrd[S,S] <- -1
    colnames(SOrd) <- sapply(1:S , function(x) paste0('S',x))

    numCoeff <- sum(coeff[2:(numCoeff+1)]!=0)

    fullCnstrn <- rbind(nOrd, cbind(0,as.matrix(Matrix::bdiag(diag(1,numCoeff),SOrd))))
    fullCnstrn <- as.matrix(dplyr::distinct(data.frame(fullCnstrn)))

    allnames <- colnames(Z)

    if(!preS){
      coeff <- base::rep(sPoint, 1 + numCoeff + S)
      names(coeff) <- allnames
    }else{
      if(withCov){
        if(type  == "surv"){
          suppressWarnings(Model <- coxph(Surv(Times, cen) ~ ., data = data.frame(cbind(Z, Cov, Times, cen))))
        }
        if(type  == "bin"){
          suppressWarnings(Model <- glm(y ~ . - 1, data = data.frame(cbind(Z, Cov, y)), family = binomial))
        }
      }else{
        if(type  == "surv"){
          suppressWarnings(Model <- coxph(Surv(Times, cen) ~ ., data = data.frame(cbind(Z, Times, cen))))
        }
        if(type  == "bin"){
          suppressWarnings(Model <- glm(y ~ . - 1, data = data.frame(cbind(Z, y)), family = binomial))
        }

      }

      coeff <- Model$coefficients
      coeff[is.na(coeff)] <- -sPoint

      if(withCov){
        theta <- coeff[colnames(Cov)]
      }else{
        theta <- 0
      }
      coeff <- coeff[allnames]
    }

    S = S + 1
  }

  result <- result[stagingV]

  totalS <- max(result)
  print("The initial staging result:")
  print(result)

  if( (sum(allCs<=eps)/(length(allCs)))  < 1){
    warning(paste0("The algorithm does not converge at the stage ", paste(which(allCs > eps), collapse = ",")))
  }

  if(max(result) > 1){

    Data[,"cStage"] <- result[as.character(Data[,"variable"])] - 1

    counts <- xtabs(~cStage, data = Data)

    if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))){

      result_f <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                       Data = Data, withCov = withCov, maxiter = maxiter,
                       eps = eps, GIC = GIC, initial_result = result,
                       useIbs = F, useLRT = F, useAIC = T,
                       coarse_pruning = F, fine_pruning = F,
                       fine_pruning_quad = F, prefix_stage = NULL,
                       threshold = NULL, type = type, perc = perc, minObs = minObs)[[1]]
    }else{
      result_f <- result
    }

  }else{
    result_f <- result
  }

  return(list("The result" = result, "After excluding small stages" = result_f))
}


#' Calculating the estimated coefficients for stages and the corresponding value for the chosen criterion
#'
#' This function calculates the estimated coefficients for stages and the corresponding value for the chosen criterion. It takes a dataset, censoring statuses, outcomes, current staging result, covariates, and other parameters as inputs.
#'
#' @param Data The given dataset.
#' @param cen A vector of censoring statuses, or the binary outcome.
#' @param y A vector of outcomes.
#' @param curR The current staging result.
#' @param Cov A matrix of covariates.
#' @param withCov A boolean indicating whether any covariates need adjusting.
#' @param useAIC A boolean indicating whether AIC (Akaike Information Criterion) is used
#' @param useIbs A boolean indicating whether the (integrated) Brier score is applied.
#' @param type The type of outcome, either survival outcome or binary outcome.
#' @param seed The seed used for generating the simulation dataset.
#' @param fullModel A boolean indicating whether it is a full model, representing the model before pruning to the current model.
#' @param fullZ A matrix of stages for the full model.
#' @param GIC The penalty term used in the criterion. The default value is 2 when AIC is used.
#'
#' @return The estimated coefficients for stages and the corresponding value for the chosen criterion.
#'
#' @importFrom survival survfit
#'
#' @export
#'
getsEst <- function(Data, cen, y, curR, Cov, withCov, type, useIbs = F, useAIC = F, seed = 0,
                    fullModel = FALSE, fullZ = NULL, GIC = 2){


  Data[,"cStage"] <- curR[as.character(Data[,"variable"])] - 1

  if(withCov){
    covName <- colnames(Cov)
  }

  r <- max(curR)

  if(useIbs){
    if(type == "surv"){

      if(!withCov){

        pred <- survfit(Surv(Data[,"time"], Data[,"status"]) ~ Data[,"cStage"])
        ibs <- sbrierC(Surv(Data[,"time"], Data[,"status"]),
                       as.vector(lapply(Data[,"cStage"],
                                        function(x) pred[as.numeric(x) +  1])) )

        bestIBS  <- ibs[[1]]

      }else{


        dList <- list()
        dList$time <- Data[, "time"]
        dList$status <- Data[, "status"]
        for(k in 1:length(covName)){
          dList[[covName[k]]] <- Data[, covName[k]]
        }
        suppressWarnings(dfit <- coxph(Surv(time, status) ~ . + strata(Data[,"cStage"]), data = dList))
        pred <- survfit(dfit)

        ibs <- sbrierC(Surv(Data[,"time"], Data[,"status"]),
                       as.vector(lapply(Data[,"cStage"],
                                        function(x) pred[as.numeric(x) +  1])) )

        bestIBS  <- ibs[[1]]
      }
    }

    if(type == "bin"){
      if(!withCov){


        dList <- listDummies(dM(Data[,"cStage"], ref = F))
        dList$outcome <- Data[, "outcome"]

        suppressWarnings(dfit <- glm(outcome ~ . - 1, data = dList))

        dfit <- cffStg(coeff0 = dfit$coef, theta0 = 0,
                       cen = cen, y = y,
                       stage = Data$cStage,
                       Cov = matrix(0,1,1), withCov = withCov, maxiter = 20, type = type,
                       seed = seed, getLikelihood = F)

        bestIBS <-  mean((1/(1+exp(- dM(Data[,"cStage"], ref = F) %*% dfit)) - dList$outcome)^2)

      }else{

        dList <- listDummies(dM(Data[,"cStage"], ref = F))
        dList$outcome <- Data[, "outcome"]

        for(k in 1:length(covName)){
          dList[[covName[k]]] <- Data[, covName[k]]
        }

        suppressWarnings(dfit <- glm(outcome ~ . - 1, data = dList))

        dfit <- cffStg(coeff0 = dfit$coef[-c((length(dfit$coef)-ncol(Cov)+1):length(dfit$coef))],
                       theta0 = dfit$coef[(length(dfit$coef)-ncol(Cov)+1):length(dfit$coef)],
                       cen = cen, y = y,
                       stage = Data$cStage,
                       Cov = Cov, withCov = withCov, maxiter = 20, type = type,
                       seed = seed, getLikelihood = F)
        bestIBS <- mean((1/(1+exp(- (dM(Data[,"cStage"], ref = F) %*% dfit[[1]] + Cov %*% dfit[[2]]) )) - dList$outcome)^2)
      }
    }

    return(bestIBS)
  }

  if(useAIC){


    if(withCov){

      if(type == "surv"){

        if(fullModel){
          dList <- fullZ

        }else{
          dList <- listDummies(dM(Data[,"cStage"]))
        }

        dList$time <- Data[, "time"]
        dList$status <- Data[, "status"]
        for(k in 1:length(covName)){
          dList[[covName[k]]] <- Data[, covName[k]]
        }
        suppressWarnings(dfit <- coxph(Surv(time, status) ~ ., data = dList))

        Coeff0 <- dfit$coef[1:(length(dfit$coef) - length(covName))]

        if(fullModel){
          Coeff0[is.na(Coeff0)] <- 0
        }else{
          Coeff0 <- c(0, Coeff0)
        }

        Theta0 <- dfit$coef[(length(dfit$coef) - length(covName)+1):length(dfit$coef)]
      }

      if(type == "bin"){

        if(fullModel){
          dList <- fullZ
        }else{
          dList <- listDummies(dM(Data[,"cStage"], ref = F))
        }

        dList$outcome <- Data[, "outcome"]

        for(k in 1:length(covName)){
          dList[[covName[k]]] <- Data[, covName[k]]
        }

        suppressWarnings(dfit <- glm(outcome ~ . - 1, data = dList))

        Coeff0 <- dfit$coef[1:(length(dfit$coef) - length(covName))]
        Theta0 <- dfit$coef[(length(dfit$coef) - length(covName)+1):length(dfit$coef)]
      }

      dfit <- cffStg(coeff0 = Coeff0, theta0 = Theta0,
                     cen = cen, y = y,
                     stage = Data$cStage,
                     Cov = Cov, withCov = withCov, maxiter = 20, type = type,
                     seed = seed, getLikelihood = T, fullModel = fullModel, fullZ = as.matrix(data.frame(fullZ)))

      best_beta <- dfit[[2]]

      best_theta <- dfit[[3]]

      best_L <-  dfit[[1]]

      best_aic <- 2*dfit[[1]] + GIC * (length(best_beta) + length(best_theta))

      return(list(best_beta, best_L, best_aic, best_theta))

    }else{

      if(type == "surv"){
        if(fullModel){
          dList <- fullZ
        }else{
          dList <- listDummies(dM(Data[,"cStage"]))
        }
        dList$time <- Data[, "time"]
        dList$status <- Data[, "status"]
        suppressWarnings(dfit <- coxph(Surv(time, status) ~ ., data = dList))
        Coeff0 <- dfit$coef

        if(fullModel){
          Coeff0[is.na(Coeff0)] <- 0
        }else{
          Coeff0 <- c(0, Coeff0)
        }
      }

      if(type == "bin"){

        if(fullModel){
          dList <- fullZ
        }else{
          dList <- listDummies(dM(Data[,"cStage"], ref = F))
        }

        dList$outcome <- Data[, "outcome"]

        suppressWarnings(dfit <- glm(outcome ~ . - 1, data = dList))

        Coeff0 <- dfit$coef
      }

      dfit <- cffStg(coeff0 = Coeff0, theta0 = 0,
                     cen = cen, y = y,
                     stage = Data$cStage,
                     Cov = matrix(0,1,1),
                     withCov = FALSE, maxiter = 20,
                     type = type,
                     seed = seed, getLikelihood = T, fullModel = fullModel, fullZ = as.matrix(data.frame(fullZ)))

      best_beta <- dfit[[2]]

      best_aic <- 2*dfit[[1]] + GIC * length(best_beta)

      best_L <-  dfit[[1]]

      return(list(best_beta, best_L, best_aic))
    }


  }

}


#' Fitting the pruning step for the OPERA algorithm
#'
#' This function performs the pruning step for the OPERA algorithm.
#'
#' @param cen A vector specifying the censoring statuses, or the binary outcome.
#' @param y A vector representing the outcome variable.
#' @param Z A matrix of nodes.
#' @param Cnstrn A matrix of constraints.
#' @param Cov A matrix of covariates.
#' @param Data A given dataset.
#' @param withCov A boolean indicating whether any covariates need adjusting.
#' @param maxiter The maximum number of iterations in each approximation step.
#' @param eps The desired accuracy level.
#' @param GIC The penalty term used in the criterion. The default value is 2 when AIC is used.
#' @param useAIC A boolean indicating whether AIC (Akaike Information Criterion) is used.
#' @param useIbs A boolean indicating whether the (integrated) Brier score is applied.
#' @param useLRT A boolean indicating whether the likelihood ratio test is applied.
#' @param initial_result The initial result from OPERA.
#' @param coarse_pruning A boolean indicating whether coarse pruning is applied.
#' @param fine_pruning A boolean indicating whether fine pruning using exhaustive search is applied.
#' @param fine_pruning_quad A boolean indicating whether fine pruning using quadratic programming constraint is applied.
#' @param prefix_stage The number of stages needed.
#' @param threshold The type-I error rate for the likelihood ratio test.
#' @param type The outcome type, either survival outcome or binary outcome.
#' @param seed The seed used for generating the simulation dataset.
#' @param getsAll A boolean indicating whether pruning is performed until only two stages are left with the values for all the criteria
#' @param checkSampleSize A boolean indicating whether the number of patients in each stage is checked.
#' @param perc The minimum proportion of patients required in each stage
#' @param minObs The minimum number of patients required in each stage
#'
#' @return The cancer staging result after pruning.
#'
#'
#' @importFrom stats pchisq
#'
#' @export
#'
operap <- function(cen, y, Z, Cnstrn, Cov, Data,
                   withCov = F, maxiter = 10, eps = 10^-4, GIC = NULL,
                   initial_result, useIbs  = FALSE, useLRT  = TRUE, useAIC = FALSE, coarse_pruning = FALSE,
                   fine_pruning = TRUE, fine_pruning_quad = FALSE,
                   prefix_stage = 5, threshold = 0.01, type = "surv", seed = 0, getsAll = F, perc = 0.1, minObs = 30, checkSampleSize = T, reScale = 1, ...){


  if(getsAll & (fine_pruning | fine_pruning_quad | coarse_pruning)){
    useIbs  <- TRUE
    useLRT  <- TRUE
    useAIC <- TRUE
  }
  # saves a copy of Z
  Zp <- Z
  Cnstrnp  <- Cnstrn

  Z <- Zp

  result <- initial_result
  # gets the number of stages
  totalS <- max(result)

  if(!is.null(prefix_stage)){
    if(prefix_stage <= 1 | prefix_stage > totalS){
      warning(paste0("The number of stages has to be larger than 1 and no larger than ", totalS, "!"))
      return(initial_result)
    }
  }


  if(is.null(GIC)){
    GIC <- 2
  }

  # gets the number of coefficients based on the number of nodes
  numCoeff <- dim(Z)[2]

  stagingV <- colnames(Z)

  if(withCov){
    covName <- colnames(Cov)
  }else{
    Cov <- matrix(0, 1, 1)
  }

  if(checkSampleSize & (sum(fine_pruning + fine_pruning_quad + coarse_pruning)  == 0)){

    allEsts <- getsEst(Data = Data, cen = cen, y = y, curR = initial_result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
            seed = seed, GIC = GIC)

    best_beta <- allEsts[[1]]

    best_aic  <- allEsts[[3]]

    if(withCov){
      best_theta <- allEsts[[4]]
    }

    best_result <- initial_result

    print(paste("The model criterion equals: ", round(best_aic, -log10(eps)-1), " before performing the coarse pruning."))

    print("Checking the number of patients in each stage:")

    Data[,"cStage"] <- best_result[as.character(Data[,"variable"])]
    counts <- xtabs(~cStage, data = Data)



    if( ((sum(counts < as.integer(perc*nrow(Data))) == 0) | (sum(counts < minObs) == 0))){

      result <- best_result
      print("Each stage has sufficient patients.")

    }else{
      print(paste0("There are no sufficient patients in the following stages: ", paste0(names(which((counts <  as.integer(perc*nrow(Data))) &  (counts <  minObs) )),  collapse = ", ")))
      print("Merging has begun:")

      r <- unname(best_result)

      toBeM <- as.integer(names(which((counts <  as.integer(perc*nrow(Data))) &  (counts <  minObs) )))

      if(1 %in% toBeM){
        r <- ifelse(r %in% c(1, 2), 1.5, r)
      }

      if(max(r) %in% toBeM){
        r <- ifelse(r %in% c(max(r) - 1, max(r)), max(r) - 0.5, r)
      }

      toBeMN <- toBeM[!toBeM  %in% c(1, max(r))]

      bestAIC <- NULL

      if(length(toBeMN) > 0){
        howToMerge <- as.matrix(expand.grid(replicate(length(toBeMN), c(-0.5, 0.5),  simplify = F)))
        for(kk in 1:nrow(howToMerge)){
          for(cc in 1:ncol(howToMerge)){
            if(howToMerge[kk, cc] > 0){
              r <- ifelse(r %in% c(toBeMN[cc], toBeMN[cc]+1), toBeMN[cc]+0.5, r)
            }else{
              r <- ifelse(r %in% c(toBeMN[cc], toBeMN[cc]-1), toBeMN[cc]-0.5, r)
            }
          }

          r <- as.numeric(as.factor(r))

          names(r) <- stagingV

          print("The current stages: ")
          print(r)

          eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = r, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                           seed = seed, GIC = GIC)

          curbeta <- eEsts[[1]]

          aic <- eEsts[[3]]

          if(withCov){
            theta <- eEsts[[4]]
          }

          print(paste("The current model criterion equals: ", round(aic, -log10(eps)-1)))

          if(is.null(bestAIC)){
            bestAIC <- aic
            curBest <- r
          }else{
            if(aic < bestAIC){
              bestAIC <- aic
              curBest <- r
            }
          }

        }
      }else{
        r <- as.numeric(as.factor(r))

        names(r) <- stagingV

        print("The current stages: ")
        print(r)
        curBest <- r
      }

      Data[,"cStage"] <- curBest[as.character(Data[,"variable"])]

      counts <- xtabs(~cStage, data = Data)

      if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))  & ( length(unique(curBest)) > 2 ) ){

        curBest <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                       Data = Data, withCov = withCov, maxiter = maxiter,
                       eps = eps, GIC = GIC, initial_result = curBest,
                       useIbs = F, useLRT = F, useAIC = F,
                       coarse_pruning = F, fine_pruning = F,
                       fine_pruning_quad = F, prefix_stage = NULL,
                       threshold = NULL, type = type, perc = perc, minObs = minObs, checkSampleSize = T)[[1]]
      }

    }

    print("After merging, the best result: ")
    print(curBest)
    finalResult <- list(curBest)
  }

  if(max(totalS) > 2 & coarse_pruning){

    if(is.null(prefix_stage)){
      reducedS <- totalS - 2
    }else{
      reducedS <- totalS - prefix_stage
    }

    if(reducedS > 0){

      currentReducedS <- 0

      best_result <- initial_result

      if(useIbs){

        bestIBS <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov,
                           type = type, useAIC = F, useIbs = T,
                           seed = seed, GIC = GIC)

        print(paste0("The initial ibs equals: ", round(bestIBS, -log10(eps)), " before pruning."))
      }

      allEsts <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                         seed = seed, GIC = GIC)

      best_beta <- allEsts[[1]]

      best_aic  <- allEsts[[3]]

      if(withCov){
        best_theta <- allEsts[[4]]
      }

      print("The coarse pruning has started.")

      if(useIbs){
        ibss <- c(bestIBS)
      }

      if(useLRT){
        pvs <- c(0)
      }

      if(useAIC){
        print(paste0("The initial AIC equals: ", round(best_aic, -log10(eps)), " before pruning."))
        aicss <- c(best_aic)
      }


      Data[,"cStage"] <- best_result[as.character(Data[,"variable"])] - 1

      counts <- xtabs(~cStage, data = Data)


      if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))  & ( length(unique(best_result)) > 2 ) ){

        best_result_f <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                                Data = Data, withCov = withCov, maxiter = maxiter,
                                eps = eps, GIC = GIC, initial_result = best_result,
                                useIbs = F, useLRT = F, useAIC = F,
                                coarse_pruning = F, fine_pruning = F,
                                fine_pruning_quad = F, prefix_stage = NULL,
                                threshold = NULL, type = type, perc = perc, minObs = minObs, checkSampleSize = T)[[1]]
      }else{
        best_result_f <- best_result
      }

      rs <- list(best_result_f)

      while(currentReducedS < reducedS){

        result <- best_result

        # updates the total number of stages
        totalS <- length(unique(result))

        eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov,
                           type = type, useAIC = T, useIbs = F,
                           seed = seed, GIC = GIC)

        oldL <- eEsts[[2]]

        best_L <- NULL
        # merges each stage to the previous stage or the next stage
        # starts from the second stage and ends at the last second stage
        for(s1 in 2:(totalS-1)){
          for(s2 in c(-1,1)){

            r <- ifelse(result == s1, s1+s2, result)

            r <- as.numeric(as.factor(r))
            names(r) <- stagingV

            eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = r, Cov = Cov, withCov = withCov,
                             type = type, useAIC = T, useIbs = F,
                             seed = seed, GIC = GIC)
            curbeta <- eEsts[[1]]
            L <- eEsts[[2]]

            if(withCov){
              theta <- eEsts[[4]]
            }


            if(sum(is.infinite(curbeta))==0){
              if(is.null(best_L)){
                best_L <- L
                best_beta <- curbeta
                if(withCov){
                  best_theta <- theta
                }
                best_result <- r
              }else{
                if (L <= best_L){
                  best_L <- L
                  best_beta <- curbeta
                  if(withCov){
                    best_theta <- theta
                  }
                  best_result <- r
                }
              }
            }

          }
        }

        print(paste0("The staging result after pruning down one stage with the largest likelihood is"))
        print(best_result)

        currentReducedS <- currentReducedS + 1


        Data[,"cStage"] <- best_result[as.character(Data[,"variable"])] - 1


        counts <- xtabs(~cStage, data = Data)


        if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))  & ( length(unique(best_result)) > 2 ) ){

          best_result_f <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                           Data = Data, withCov = withCov, maxiter = maxiter,
                           eps = eps, GIC = GIC, initial_result = best_result,
                           useIbs = F, useLRT = F, useAIC = T,
                           coarse_pruning = F, fine_pruning = F,
                           fine_pruning_quad = F, prefix_stage = NULL,
                           threshold = NULL, type = type, perc = perc, minObs = minObs)[[1]]
        }else{
          best_result_f <- best_result
        }

        if(useLRT){

          pvalue <- 1 - pchisq(2*abs(oldL - best_L), df = 1)

          print(paste("The p-value equals", round(pvalue, -log10(eps))))

          pvs <- c(pvs, pvalue)

          finalResult <- list("The result" = rs[[length(rs)]], "The p-values" = pvs)

          if(!getsAll){
            if(!is.null(threshold)){
              if(pvalue <= threshold){
                break
              }
            }
          }
        }


        rs[[length(rs)+1]] <- best_result_f

        if(!getsAll){
          if(useLRT  & (!is.null(threshold))){
            finalResult <- list("The result" = rs[[length(rs)]], "The p-values" = pvs)
          }
        }

        if(useLRT  & (is.null(threshold))){
          finalResult <- list("The results" = rs, "The p-values" = pvs)
        }

        if(useIbs){

          curIBS <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov,
                             type = type, useAIC = F, useIbs = T,
                             seed = seed, GIC = GIC)

          ibss <- c(ibss, curIBS)

          print(paste("The (integrated) brier score equals", round(curIBS, -log10(eps))))

          finalResult <- list("The result" = rs[[which.min(ibss)]], "The brier scores" = ibss)
        }

        if(useAIC){

          allEsts <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                             seed = seed, GIC = GIC)

          curaic  <- allEsts[[3]]

          aicss <- c(aicss, curaic)
          print(paste("The AIC equals", round(curaic, -log10(eps))))

          finalResult <- list("The result" = rs[[which.min(aicss)]], "The AICs or equivalances" = aicss)
        }

        if(!is.null(prefix_stage)){
          finalResult <- list("The result" = rs[[length(rs)]])
        }

        if(getsAll){
          finalResult <- list("The results" = rs, "The p-values" = pvs, "The brier scores" = ibss, "The AICs or equivalances" = aicss)
        }

      }
    }
  }

  if(max(totalS) > 2 & (fine_pruning_quad | fine_pruning)){

    if(fine_pruning_quad){
      allCs <- c()
      whichCs <- c()
    }


    if(is.null(prefix_stage)){
      reducedS <- totalS - 2
    }else{
      reducedS <- totalS - prefix_stage
    }

    if(reducedS > 0){

      currentReducedS <- 0

      best_result <- initial_result


      if(useIbs){

        bestIBS <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov,
                           type = type, useAIC = F, useIbs = T,
                           seed = seed, GIC = GIC)

        print(paste0("The initial (integrated) brier score equals: ", round(bestIBS, -log10(eps)), " before pruning."))
      }

      if(fine_pruning_quad){
        print("The fine pruning has started using the quadratic programming constraint.")
      }

      if(fine_pruning){
        print("The fine pruning has started using the exhaustive enumeration.")
      }


      allEsts <- getsEst(Data = Data, cen = cen, y = y, curR = best_result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                         seed = seed, GIC = GIC)

      best_beta <- allEsts[[1]]

      best_aic  <- allEsts[[3]]

      if(withCov){
        best_theta <- allEsts[[4]]
      }


      PO <- Cnstrn

      if(useIbs){
        ibss <- c(bestIBS)
      }

      if(useLRT){
        pvs <- c(0)
      }

      if(useAIC){
        print(paste0("The initial AIC equals: ", round(best_aic, -log10(eps)), " before pruning."))
        aicss <- c(best_aic)
      }


      Data[,"cStage"] <- best_result[as.character(Data[,"variable"])] - 1

      counts <- xtabs(~cStage, data = Data)


      if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))  & ( length(unique(best_result)) > 2 ) ){

        best_result_f <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                              Data = Data, withCov = withCov, maxiter = maxiter,
                              eps = eps, GIC = GIC, initial_result = best_result,
                              useIbs = F, useLRT = F, useAIC = T,
                              coarse_pruning = F, fine_pruning = F,
                              fine_pruning_quad = F, prefix_stage = NULL,
                              threshold = NULL, type = type, perc = perc, minObs = minObs)[[1]]
      }else{
        best_result_f <- best_result
      }

      rs <- list(best_result_f)

      result <- best_result

      while(currentReducedS < reducedS){

        # update the total number of stages
        totalS <- length(unique(result))

        best_L <- NULL

        if(fine_pruning_quad){
        # splits each stage into two partitions and merges one to the previous stage and the other to the next stage
        # starts from the second stage and ends at the last second stage
          for(s in 2:(totalS-1)){
          # gather all nodes at current stage
          nds <- names(result[result == s])
          # constraints for nodes
          consts <- PO[, nds, drop = FALSE]
          consts <- consts[!duplicated(consts), ,drop = FALSE]

          if(is.null(dim(consts))){
            consts <- base::matrix(consts)
            colnames(consts) <- nds
          }

          allS <- 1:totalS
          otherS <- allS[! allS %in% c(s)]
          otherConsts <- matrix(0, dim(consts)[1], length(otherS))
          colnames(otherConsts) <- paste0(c("S"), otherS)
          consts <- cbind(consts, otherConsts)
          consts <- consts[RowSums(consts) == 0,,drop = FALSE]
          if(is.null(dim(consts))){
            consts <- t(as.matrix(consts))
          }
          consts <- consts[apply(consts, 1, function(x) sum(x == 0)) != dim(consts)[2], ,drop = FALSE]
          if(is.null(dim(consts))){
            consts <- t(as.matrix(consts))
          }
          # gets the constraints
          for(o in 1:(length(otherS)-1)){
            for(q in 1:(length(otherS)-o)){
              newConst <- t(as.matrix(base::rep(0, dim(consts)[2])))
              colnames(newConst) <- colnames(consts)
              newConst[,colnames(otherConsts)[o]] = -1
              newConst[,colnames(otherConsts)[o+q]] = 1
              consts <- rbind(consts, newConst)
            }
          }

          NSconsts1 <- cbind(diag(dim(consts)[2] - length(otherS)), matrix(0, nrow = (dim(consts)[2] - length(otherS)), ncol = length(otherS)))
          colnames(NSconsts1) <- colnames(consts)
          NSconsts1[,paste0("S", s-1)] <- base::rep(-1, dim(NSconsts1)[1])

          NSconsts2 <- cbind(-diag(dim(consts)[2] - length(otherS)), matrix(0, nrow = dim(consts)[2] - length(otherS), ncol = length(otherS)))
          colnames(NSconsts2) <- colnames(consts)
          NSconsts2[,paste0("S", s+1)] <- base::rep(1, dim(NSconsts2)[1])

          allConsts <- rbind(consts, NSconsts1, NSconsts2)

          # gets the quadratic constraint
          quadConsts <- rbind(NSconsts1, NSconsts2)

          interConsts <- t(as.matrix(base::rep(0, dim(quadConsts)[2])))
          colnames(interConsts) <- colnames(quadConsts)

          interConsts[,paste0("S", s-1)] <- -1
          interConsts[,paste0("S", s+1)] <- 1

          qConsts1 <- t(quadConsts) %*% quadConsts
          qConsts2 <- (t(interConsts) %*% interConsts) * (dim(interConsts)[2] - length(otherS))

          qConsts <- qConsts2 - qConsts1

          if(withCov){
            qConsts <- as.matrix(Matrix::bdiag(qConsts, diag(0,length(covName))))
          }

          # starts pruning

          # initializes Z
          Z <- Zp[,nds,drop =FALSE]
          for(eachS in 1:length(otherS)){
            nextZ <- as.matrix(RowSums(Zp[,names(result[result == otherS[eachS]]),drop =FALSE]))
            colnames(nextZ) <- paste0("S", otherS[eachS], collapse = "")
            Z <- cbind(Z, nextZ)
          }

          n.iter <- 0
          diff.beta <- Inf

          eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov,
                           type = type, useAIC = T, useIbs = F,
                           seed = seed, GIC = GIC)

          coeffs <- eEsts[[1]]
          curbeta <- base::rep(0, dim(allConsts)[2])
          curbeta[(length(curbeta) - max(result) + 2):length(curbeta)] <- coeffs[-s]
          curbeta[1:(length(curbeta) - max(result) + 1)] <- coeffs[s]
          names(curbeta) <- colnames(allConsts)

          if(withCov){
            diff.theta <- Inf
            theta <- eEsts[[4]]
            names(theta) <- covName
          }

          # adds more constraints
          oneMoreConst <- rep(0, ncol(allConsts))
          names(oneMoreConst)  <- colnames(allConsts)
          oneMoreConst[length(oneMoreConst)] <- -1

          secondMoreConst <- rep(0, ncol(allConsts))
          names(secondMoreConst)  <- colnames(allConsts)
          secondMoreConst[length(secondMoreConst) - (max(result)-2)] <- 1

          oldConsts <- allConsts
          oneConsts <- rbind(allConsts,  oneMoreConst)
          twoConsts <-  rbind(allConsts,  oneMoreConst, secondMoreConst)

          repeatConvs <- 0
          listOfBetas <- list()

          if(withCov){
            listOfthetas  <- list()
            thetaInitial <- theta
          }

          vectorOfConvs  <-  c()
          betaInitial <- curbeta
          vectorOfL <- c()

          while(repeatConvs < 3){

            curbeta <- betaInitial
            if(withCov){
              theta <- thetaInitial
            }else{
              Cov <- matrix(0,1,1)
              theta <- 0
            }


            # IRWLS
            if(repeatConvs ==  0){
              irwls_pruning <- IRLSPAT(coeff = curbeta/reScale, theta = theta, cen = cen, y = y, Z = Z,
                                       cov = Cov, LCnstrn = oldConsts, qCnstrn = qConsts,
                                       mCnstrn = repeatConvs,
                                       minCoeff = min(curbeta), maxCoeff = max(curbeta),
                                       eps = eps, maxiter = maxiter, withCov = withCov,
                                       seed = seed, type = type)
            }

            if(repeatConvs ==  1){
              irwls_pruning <- IRLSPAT(coeff = curbeta/reScale, theta = theta, cen = cen, y = y, Z = Z,
                                       cov = Cov, LCnstrn = oneConsts, qCnstrn = qConsts,
                                       mCnstrn = repeatConvs,
                                       minCoeff = min(curbeta), maxCoeff = max(curbeta),
                                       eps = eps, maxiter = maxiter, withCov = withCov,
                                       seed = seed, type = type)

            }

            if(repeatConvs ==  2){

              irwls_pruning <- IRLSPAT(coeff = curbeta/reScale, theta = theta, cen = cen, y = y, Z = Z,
                                       cov = Cov, LCnstrn = twoConsts, qCnstrn = qConsts,
                                       mCnstrn = repeatConvs,
                                       minCoeff = min(curbeta), maxCoeff = max(curbeta),
                                       eps = eps, maxiter = maxiter, withCov = withCov,
                                       seed = seed, type = type)

            }


            diff.beta <- irwls_pruning[[2]]

            repeatConvs <- repeatConvs + 1

            if(withCov){
              listOfBetas[[repeatConvs]] <- irwls_pruning[[1]][1:length(curbeta)]
              listOfthetas[[repeatConvs]] <- irwls_pruning[[1]][-c(1:length(curbeta))]
            }else{
              listOfBetas[[repeatConvs]] <- irwls_pruning[[1]]
            }
            vectorOfConvs  <-  c(vectorOfConvs, diff.beta)
            vectorOfL <- c(vectorOfL, irwls_pruning[[3]])
          }
          # chooses the set of coefficients with the largest likelihood
          curbeta <- listOfBetas[[which.min(vectorOfL)]]
          eachC <- vectorOfConvs[which.min(vectorOfL)]
          whichEachC <- which.min(vectorOfL)

          names(curbeta) <- colnames(Z)

          if(sum(is.infinite(curbeta))==0){

            result_p  <- result

            compareWithS <- result_p[names(curbeta[1])]

            needToBeMerged_p <- names(curbeta[1:(length(curbeta) - totalS + 1)])

            for(nt in 1:length(needToBeMerged_p)){
              if(curbeta[needToBeMerged_p[nt]] == curbeta[paste0("S", compareWithS - 1)]){
                result_p[needToBeMerged_p[nt]] <- result_p[needToBeMerged_p[nt]] - 1
              }else{
                result_p[needToBeMerged_p[nt]] <- result_p[needToBeMerged_p[nt]] + 1
              }
            }

            # changes the result
            result_p[result_p > compareWithS] <- result_p[result_p > compareWithS] - 1

            if(withCov){
              theta <- listOfthetas[[which.min(vectorOfL)]]

              names(theta) <- colnames(Cov)
            }

            L_now <- getsEst(Data = Data, cen = cen, y = y, curR = result_p, Cov = Cov,
                             withCov = withCov, type = type, useAIC = T,
                             seed = seed, GIC = GIC)[[2]]

            # updates the best aic or bic
            if (is.null(best_L)){
              best_L <- L_now
              best_beta <- curbeta
              if(withCov){
                best_theta <- theta
              }
              needToBeMerged <- needToBeMerged_p
              best_result <- result_p
              best_C <- eachC
              best_which <- whichEachC
            }else{
              if(L_now < best_L){
                best_L <- L_now
                best_beta <- curbeta
                needToBeMerged <- needToBeMerged_p
                best_result <- result_p
                if(withCov){
                  best_theta <- theta
                }
                best_C <- eachC
                best_which <- whichEachC
              }
            }

          }

        }


          # stops pruning
          if(is.null(best_L)){
            break
          }

          allCs <- c(allCs, best_C)
          whichCs <- c(whichCs, best_which)

          currentReducedS <- currentReducedS + 1

          result <- best_result

          print(paste0("The staging result after pruning down one stage with the largest likelihood is"))
          print(result)

        }

        if(fine_pruning){
          # splits each stage into two partitions and merges one to the previous stage and the other to the next stage
          # starts from the second stage and ends at the last second stage
          for(s in 2:(totalS-1)){
            # gathers all nodes at current stage
            nds <- names(result[result == s])
            # gets the constraints for nodes
            consts <- PO[, nds, drop = FALSE]
            consts <- consts[!duplicated(consts), ,drop = FALSE]
            # fixes the constraints by setting drop = F
            if(is.null(dim(consts))){
              consts <- base::matrix(0, nrow = 1, ncol = length(nds))
              colnames(consts) <- nds
            }

            consts <- consts[RowSums(consts) == 0, ,drop = FALSE]

            if(is.null(dim(consts))){
              consts <- base::matrix(0, nrow = 1, ncol = length(nds))
              colnames(consts) <- nds
            }

            consts <- consts[apply(consts, 1, function(x) sum(x == 0)) != dim(consts)[2], ,drop = FALSE]

            if(is.null(dim(consts))){
              consts <- base::matrix(0, nrow = 1, ncol = length(nds))
              colnames(consts) <- nds
            }

            allCmbns <- as.matrix(expand.grid(replicate(length(nds), c(s-1,s+1), simplify = F)))
            colnames(allCmbns) <- nds

            allCmbns <- allCmbns[!apply(allCmbns %*% t(consts), 1, function(x) sum(x<0)>0), ,drop = FALSE]

            #  starts the pruning
            for(i in 1:nrow(allCmbns)){
              r <- result

              r[as.character(nds)] <- allCmbns[i,]

              r <- as.numeric(as.factor(r))
              names(r) <- stagingV

              eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = r, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                                 seed = seed, GIC = GIC)

              curbeta <- eEsts[[1]]
              L <- eEsts[[2]]

              if(withCov){
                theta <- eEsts[[4]]
              }

              if(sum(is.infinite(curbeta))==0){
                if(is.null(best_L)){
                  best_L <- L
                  best_beta <- curbeta
                  if(withCov){
                    best_theta <- theta
                  }

                  needToBeMerged <- nds

                  best_result <- r
                }else{
                  if (L <= best_L){
                    best_L <- L
                    best_beta <- curbeta
                    if(withCov){
                      best_theta <- theta
                    }
                    needToBeMerged <- nds
                    best_result <- r

                  }
                }
              }
            }


          }


          # stops pruning
          if(is.null(best_L)){
            break
          }

          # changes the result
          result <- best_result

          currentReducedS <- currentReducedS + 1


          Data[,"cStage"] <- result[as.character(Data[,"variable"])] - 1

          print(paste0("The staging result after pruning down one stage with the largest likelihood is"))
          print(result)
        }


        if(useLRT){

          eEsts <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                           seed = seed, GIC = GIC)

          Data[,"cStage"] <- result[as.character(Data[,"variable"])] - 1
          secondList <- listDummies(dM(Data[,"cStage"], ref = F))
          secondList[[paste0("x", max(result)+1)]] <- ifelse(as.character(Data[,"variable"])  %in% needToBeMerged, 1, 0)

          firstL <- eEsts[[2]]

          eEsts2  <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                             seed = seed, fullModel = T, fullZ = secondList, GIC = GIC)


          secondL <- eEsts2[[2]]

          pvalue <- 1 - pchisq(2*abs(secondL - firstL), df = 1)

          print(paste("The p-value equals", round(pvalue, -log10(eps))))

          pvs <- c(pvs, pvalue)

          finalResult <- list("The result" = rs[[length(rs)]], "The p-value" = pvs[length(pvs)])

          if(!getsAll){
            if(!is.null(threshold)){
              if(pvalue <= threshold){
                break
              }
            }
          }
        }

        Data[,"cStage"] <- result[as.character(Data[,"variable"])] - 1

        counts <- xtabs(~cStage, data = Data)


        if( ((sum(counts < as.integer(perc*nrow(Data))) > 0) & (sum(counts < minObs) > 0))  & ( length(unique(result)) > 2 ) ){

          result_f <- operap(cen = cen, y = y, Z = Zp, Cnstrn = Cnstrnp, Cov = Cov,
                                Data = Data, withCov = withCov, maxiter = maxiter,
                                eps = eps, GIC = GIC, initial_result = result,
                                useIbs = F, useLRT = F, useAIC = T,
                                coarse_pruning = F, fine_pruning = F,
                                fine_pruning_quad = F, prefix_stage = NULL,
                                threshold = NULL, type = type, perc = perc, minObs = minObs)[[1]]
        }else{
          result_f <- result
        }


        rs[[length(rs)+1]] <- result_f

        if(!getsAll){
          if(useLRT  & (!is.null(threshold))){
            finalResult <- list("The result" = rs[[length(rs)]], "The p-values" = pvs)
          }
        }

        if(useLRT  & (is.null(threshold))){
          finalResult <- list("The results" = rs, "The p-values" = pvs)
        }

        if(useIbs){

          curIBS <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov,
                             type = type, useAIC = F, useIbs = T,
                             seed = seed, GIC = GIC)

          ibss <- c(ibss, curIBS)

          print(paste("The (integrated) brier score equals", round(curIBS, -log10(eps))))

          finalResult <- list("The result" = rs[[which.min(ibss)]], "The brier scores" = ibss)
        }

        if(useAIC){

          allEsts <- getsEst(Data = Data, cen = cen, y = y, curR = result, Cov = Cov, withCov = withCov, type = type, useAIC = T,
                             seed = seed, GIC = GIC)

          curaic  <- allEsts[[3]]

          aicss <- c(aicss, curaic)
          print(paste("The AIC equals", round(curaic, -log10(eps))))

          finalResult <- list("The result" = rs[[which.min(aicss)]], "The AICs or equivalances" = aicss)
        }

        if(!is.null(prefix_stage)){
          finalResult <- list("The result"  = rs[[length(rs)]])
        }



      }

      if(getsAll){
        finalResult <- list("The results" = rs, "The p-values" = pvs, "The brier scores" = ibss, "The AICs or equivalances" = aicss)
      }

      if(fine_pruning_quad){
        finalResult$`Convergence Rate` <- sum(allCs  <=  eps)/length(allCs)
      }

    }
  }

  if((max(totalS) > 2) & (is.null(prefix_stage) | ((!is.null(prefix_stage))& (sum(totalS > prefix_stage) == 1) ))){
    return(finalResult)
  }else{
    return(initial_result)
  }
}


#' Constructing an object having the coefficients
#'
#' This function constructs an object called "param" used for better storage, printing, and visualization of the coefficients.
#'
#' @param beta.vec The values of the coefficients. If there is a single tuning parameter, it should be a vector. If there are multiple tuning parameters, it should be a matrix, with each row associated with one tuning parameter.
#' @param lambda The tuning parameter associated with this beta. It can be a scalar or a vector.
#'
#' @return The object that stores the coefficients.
#'
#' @export
#'
as.param <- function(beta.vec, lambda = NULL){
  beta.param <- beta.vec
  attr(beta.param, "lambda") <- lambda
  class(beta.param) <- "param"
  return(beta.param)
}

#' Fitting the lasso tree method
#'
#' This function is the core function for lasso tree method.
#'
#' @param ncat A vector specifying the number of levels for each risk factor.
#' @param cen An indicator for censoring, or the binary outcome.
#' @param y The time to failure, or the binary outcome.
#' @param lambda The tuning parameter of the fused group lasso penalty. It can be a single value or a vector.
#' @param dat The dataset.
#' @param label Combinations of risk factors.
#' @param Apo Partial order constraints.
#' @param inibeta A vector of initial guesses for the beta coefficients.
#' @param trace A boolean variable indicating whether to show the process of calculation.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance level.
#' @param type The outcome type, either survival outcome or binary outcome.
#' @param giveNB A boolean indicating whether the number of boundaries is given.
#' @param fullA A boolean indicating whether the second derivatives are fully derived.
#'
#'
#' @return The constrained estimates using lasso tree method
#'
#'
#' @seealso
#'
#' \code{\link{as.param}}
#'
#' @export
#'
lasso.tree <- function(ncat, cen, y, lambda, dat, label, Apo, inibeta = NULL,
                       trace=F, maxiter=200, eps=1e-4,
                       type = "surv", giveNB = F, fullA = F) {
  n <- nrow(dat)
  # number of boundaries
  if(giveNB){
    nb <- nrow(Apo)
  }else{
    nb <- 0
    for(j in 1:length(ncat)){
      nb <- nb + prod(ncat[-j])*(ncat[j]-1)
    }
  }

  # creates the design matrix x
  X <- Zmatrix(dat = dat, colNames = label)

  if(type == "surv"){
    # calculates risk set
    indfail <- seq(1, n)[cen == 1]
    indset <- matrix(F, nrow = n, ncol = length(indfail))
    for (j in 1:length(indfail))
      indset[,j] <- (y >= y[indfail[j]])
  }

  if(type == "bin"){
    indset <- NULL
  }

  # fits cox model/logistic regression to get w
  if(type == "surv"){
    suppressWarnings(model <- coxph(Surv(y,cen)~X))
  }

  if(type ==  "bin"){
    suppressWarnings(model <- glm(y ~ X-1, family = binomial))
  }

  if (is.null(inibeta)){
    inibeta <- model$coefficients
    inibeta[is.na(inibeta)] <- eps
  }

  coeffs <- model$coefficients

  w <- NULL
  for (i in 1:(nb))
    w <- c(w , Apo[i,which(Apo[i,]!=0)] %*% as.vector(coeffs)[which(Apo[i,]!=0)])
  w[is.na(w)] <- 0
  w <- pmax(w,0)
  w <- sqrt(rowSums(matrix(w^2,nb,1))+0.0001)^(-1)

  if (length(lambda) == 1){
    if(type == "surv"){
      beta.out <- lasso_tree_single(X = X, indset = indset, Apo = Apo, nb = nb, w = w,
                                    cen = cen, y = y, lambda = lambda,
                                    inibeta = inibeta, maxiter = maxiter,
                                    eps = eps, trace = trace, type = type, fullA = fullA)
    }else{
      beta.out <- lasso_tree_single(X = X, indset = matrix(0,1,1), Apo = Apo, nb = nb, w = w,
                                    cen = cen, y = y, lambda = lambda,
                                    inibeta = inibeta, maxiter = maxiter,
                                    eps = eps, trace = trace, type = type, fullA = fullA)
    }


  }else{
    if(type == "surv"){
      beta.out <- lasso_tree_multi(X = X, indset = indset, Apo = Apo, nb = nb, cen = cen, w = w,
                                 y = y, lambda = lambda,inibeta = inibeta,
                                 maxiter = maxiter,
                                 eps = eps, trace = trace, type = type, fullA = fullA)
    }else{
      beta.out <- lasso_tree_multi(X = X, indset = matrix(0,1,1), Apo = Apo, nb = nb, cen = cen, w = w,
                                   y = y, lambda = lambda,inibeta = inibeta,
                                   maxiter = maxiter,
                                   eps = eps, trace = trace, type = type, fullA = fullA)
    }

  }
  beta.out <- as.param(beta.vec = beta.out, lambda = lambda)
  return(beta.out)
}

#' Calculating AIC or BIC value for a given beta
#'
#' This function calculates the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC) value for a given beta.
#'
#' @param beta An object of class "param".
#' @param ncat A vector specifying the number of levels for each risk factor.
#' @param dat The dataset.
#' @param label The combinations of risk factors.
#' @param y A vector representing the time to event or binary outcome.
#' @param cen A vector indicating the censoring status, or the binary outcome.
#' @param type The outcome type, either survival outcome or binary outcome.
#' @param nlambda The size of the search space for the tuning parameter.
#' @param useBIC A boolean indicating whether BIC is used.
#' @param eps The tolerance level.
#'
#' @return The AIC or BIC value. If beta contains a single lambda value, the function returns a single AIC or BIC value. If beta contains a vector of lambda values, the function returns a vector of AIC or BIC values.
#'
#' @export
#'
a_bic <- function(beta, ncat, dat, label, y, cen, type, nlambda, eps, useBIC = F){
  n <- nrow(dat)

  # creates the design matrix x
  X <- Zmatrix(dat = dat, colNames = label)

  if(nlambda <= 1){
  	if(useBIC){
    	ans <- 2*logLK(Z = X, cen = cen, beta = beta, y = y, withCov = F, cov = matrix(0,1,1), theta = 0, type = type) + log(nrow(X))*length(unique(beta))
    }else{
    	ans <- 2*logLK(Z = X, cen = cen, beta = beta, y = y, withCov = F, cov = matrix(0,1,1), theta = 0, type = type) + 2*length(unique(beta))
    }
  }else{
    ans <- NULL
    for (j in 1:nlambda){
      b <- beta[j, ]

      if(useBIC){
      	ans <- c(ans, 2*logLK(Z = X, cen = cen, beta = b, y = y, withCov = F, cov = matrix(0,1,1), theta = 0, type = type) + log(nrow(X))*length(unique(b)))
      }else{
      	ans <- c(ans, 2*logLK(Z = X, cen = cen, beta = b, y = y, withCov = F, cov = matrix(0,1,1), theta = 0, type = type) + 2*length(unique(b)))
      }


      print(paste0("The total number of all the parameters equals ",length(unique(b))))
      print(paste0("AIC/BIC equals ", round(ans[length(ans)], -log10(eps) - 1)))
      names(b) <- label
      print(round(b, -log10(eps) - 1))
    }
  }
  return(ans)
}

#' Fitting the Lasso Tree Method with Parameter Tuning
#'
#' This function is the core function for fitting the Lasso Tree method with parameter tuning.
#'
#' @param ncat A vector specifying the number of levels for each risk factor.
#' @param y The time to failure, or the binary outcome.
#' @param cen An indicator for censoring, or the binary outcome.
#' @param dat The dataset.
#' @param label The combinations of risk factors.
#' @param Apo Theee partial order constraints.
#' @param type The outcome type, either survival outcome or binary outcome.
#' @param inibeta A vector of initial values for the beta coefficients.
#' @param trace A boolean variable indicating whether to show the calculation process.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance level.
#' @param giveNB A boolean indicating whether the number of boundaries is given.
#' @param fullA A boolean indicating whether the second derivatives are fully derived.
#' @param useBIC A boolean indicating whether BIC is used.
#' @param N_lasso The size of the search space for tuning parameter
#'
#' @return A list that includes the best solution of beta according to the AIC or BIC, a vector of AIC or BIC values corresponding to the input lambda, and a `param` object called `beta.seq` that includes all the coefficients.
#'
#' @export
#'
lasso.tree.bic <- function(ncat, y, cen, Apo, label, dat,
                           inibeta = NULL, trace=F, maxiter=200,
                           eps=1e-4, type = "surv", giveNB = F, fullA = F,  useBIC = F,
                           N_lasso = NULL,...){

  lambda0Beta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = 0,
                            dat = dat, label = label,
                            Apo = Apo, inibeta = inibeta, trace = trace,
                            maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)

  print(paste0("The number of the unique parameters equals ", length(unique(lambda0Beta))))
  lam <- 1

  initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                         dat = dat, label = label,
                         Apo = Apo, inibeta = inibeta, trace = trace,
                         maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)

  if(length(unique(initBeta)) > 1){
    k <- 0
    while(length(unique(initBeta)) > 1){
      initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                             dat = dat, label = label,
                             Apo = Apo, inibeta = inibeta, trace = trace,
                             maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)
      lam <- lam * 2
      k <- k + 1

      if(k > 10){
        break
      }
    }
    maxlam <- lam/2
  }else{
    k <- 0
    while(length(unique(initBeta)) == 1){
      initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                             dat = dat, label = label,
                             Apo = Apo, inibeta = inibeta, trace = trace,
                             maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)
      lam <- lam / 2
      k <- k + 1
      if(k > 10){
        break
      }
    }
    maxlam <- lam*4
  }
  print("The maximum lambda equals: ")

  print(maxlam)

  lam <- 10^-2
  initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                         dat = dat, label = label,
                         Apo = Apo, inibeta = inibeta, trace = trace,
                         maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)

  if(length(unique(initBeta)) <  length(unique(lambda0Beta))){
    k <- 0
    while(length(unique(initBeta)) <  length(unique(lambda0Beta))){
      initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                             dat = dat, label = label,
                             Apo = Apo, inibeta = inibeta, trace = trace,
                             maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)

      lam <- lam / 2
      k <- k + 1
      if(k >= 10){
        break
      }
    }
    minlam <- lam*2
  }else{
    k <- 0
    while(length(unique(initBeta)) ==  length(unique(lambda0Beta))){
      initBeta <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lam,
                             dat = dat, label = label,
                             Apo = Apo, inibeta = inibeta, trace = trace,
                             maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)
      lam <- lam * 2
      k <- k + 1
      if(k >= 10){
        break
      }
    }
    minlam <- lam/4
  }
  print("The minimum lambda equals: ")
  print(minlam)

  if(is.null(N_lasso)){
    Np <- ncol(Apo)
  }else{
    Np <- N_lasso
  }

  lambda <- c(0, exp(seq(log(minlam),log(maxlam),length.out = Np)))

  nlambda <- length(lambda)

  beta.seq <- lasso.tree(ncat = ncat, cen = cen, y = y, lambda = lambda,
                         dat = dat, label = label,
                         Apo = Apo, inibeta = inibeta, trace = trace,
                         maxiter = maxiter, eps = eps, type = type, giveNB = giveNB, fullA = fullA)


  bics <- a_bic(beta = beta.seq, ncat = ncat, dat = dat, label = label,
              y = y, cen = cen, type = type, nlambda = nlambda, eps = eps, useBIC = useBIC)

  ind.opt <- which.min(bics)
  beta.opt <- as.param(beta.vec = beta.seq[ind.opt, ], lambda = lambda[ind.opt])

  stage <- rep(NA, length(beta.opt))

  for (j in 1:length(sort(unique(beta.opt)))){
    stage[which(beta.opt == sort(unique(beta.opt))[j])] <- j
  }

  names(stage) <-  colnames(Apo)

  return(stage)
}


#' Generating Risk Categories as Separate Variables
#'
#' This function constructs a dataset with separate variables indicating different ordinal risk factors based on a given dataset with a variable named `variable` that represents risk categories.
#'
#' @param dat The dataset that contains a variable named `variable` with values such as `a1b2` indicating the first level for the risk factor `a` and the second level for the risk factor `b`.
#' @param numRiskFs The total number of risk factors.
#'
#' @return The dataset with separate variables indicating different ordinal risk factors.
#'
#' @keywords dataset, risk factors, ordinal variables
#' @export
#'
createRiskGroups <- function(dat, numRiskFs){
  for(k in 1:numRiskFs){
    dat[, letters[k]] <- sapply(dat$variable, function(x) {as.numeric(strsplit(x,"")[[1]][2*k])})
  }
  return(dat)
}

#' Generating Node Labels for `visNetwork` without Stages Specified
#'
#' This function generates node labels for `visNetwork` based on a given risk category, the names of all risk factors, and the names for all levels of each risk factor.
#'
#' @param x The dataset that contains a variable named `variable` with values such as `a1b2` indicating the first level for the risk factor `a` and the second level for the risk factor `b`.
#' @param riskNames The names of all risk factors.
#' @param riskFactors The list of names for all levels of all risk factors.
#'
#' @return The node labels with the name for the specified level for each risk factor.
#'
#' @export
#'
NTitle <- function(x, riskNames, riskFactors){
  titles <- c("<p>")
  for(k in 1:length(riskNames)){
    titles <- paste0(titles, paste0(riskNames[k], " - "),
                     riskFactors[[k]][as.numeric(strsplit(x,"[A-Za-z]")[[1]][k+1])])
    if(k != length(riskNames)){
      titles <- paste0(titles, "<br>")
    }else{
      titles <- paste0(titles, "</p>")
    }
  }
  return(titles)
}

#' Generating Node Labels for `visNetwork` with Stages Specified
#'
#' This function generates node labels for `visNetwork` based on a given risk category, the names of all risk factors, and the names for all levels of each risk factor. It also incorporates stages specified by the user.
#'
#' @param x The dataset that contains a variable named `variable` with values such as `a1b2` indicating the first level for the risk factor `a` and the second level for the risk factor `b`.
#' @param riskNames The names of all risk factors.
#' @param riskFactors The list of names for all levels of all risk factors.
#' @param S The vector indicating the corresponding stages.
#'
#' @return The node labels with the name for the specified level for each risk factor and the stage for each risk factor.
#'
#' @export
#'
NTitleS <- function(x, riskNames, riskFactors, S){
  titles <- c("<p>")
  for(k in 1:length(riskNames)){
    titles <- paste0(titles, paste0(riskNames[k], " - "),
                     riskFactors[[k]][as.numeric(strsplit(x,"[A-Za-z]")[[1]][k+1])])
    titles <- paste0(titles, "<br>")

  }
  titles <- paste0(titles, "Stage: ", S[x], "</p>")
  return(titles)
}


#' Performing cancer staging
#'
#' This function performs cancer staging
#'
#' @param ncat A vector specifying the number of categories in each variable or risk factor.
#' @param dat The dataset
#' @param TimeN The variable name for survival times
#' @param yN The variable name for a binary outcome
#' @param cenN The variable name for censoring times or a binary outcome
#' @param covN The variable names for covariates
#' @param withCov Whether any covariates need adjusting
#' @param filepath The file path to save all the results
#' @param riskVariables The variable names for all risk factors
#' @param riskFactors The list of names for all levels of all risk factors.
#' @param riskNames The names of all risk factors.
#' @param useLassoT Whether lasso tree method is used
#' @param useOPERA Whether OPERA is used
#' @param usePruning Whether pruning is used
#' @param edges edges needed when a risk factor lacks total orderign
#' @param ifE whether edges are given
#' @param type The outcome type, either survival outcome or binary outcome.
#' @param getsAll If all criteria in pruning are used
#' @param useAIC If AIC is used as a criterion in pruning
#' @param useIbs If IBS is used as a criterion in pruning
#' @param useLRT If likelihood ratio test is used as a criterion in pruning
#' @param coarse_pruning A boolean indicating whether coarse pruning is applied.
#' @param fine_pruning A boolean indicating whether fine pruning using exhaustive search is applied.
#' @param fine_pruning_quad A boolean indicating whether fine pruning using quadratic programming constraint is applied.
#' @param prefix_stage The number of stages needed.
#' @param threshold The type-I error rate for the likelihood ratio test.
#' @param eps The tolerance level except for lasso tree method
#' @param eps_lasso The tolerance level for lasso tree method
#' @param maxiter The maximum number of iterations in each approximation step except for lasso tree method
#' @param maxiter_lasso The maximum number of iterations in each approximation step for lasso tree method
#' @param GIC The penalty term used in the criterion except for lasso tree method.
#' @param useBIC Whether BIC is used for lasso tree method. If false, AIC is used.
#' @param perc The minimum proportion of patients required in each stage
#' @param minObs The minimum number of patients required in each stage
#' @param plt Whether results are needed to be plotted
#' @param pic_width The width of the figure of the network of risk categories
#' @param pic_height The height of the figure of the network of risk categories
#' @param legend_size The size of the legend in the figure of the network of risk categories
#' @param x_pos The horizontal position of the legend in the figure of the network of risk categories
#' @param yratio The ratio proportional to the number of levels in each factor to control the vertical position of the legend in the figure of the network of risk categories
#' @param yaxis_min The minimum value of the y-axis for the Kaplan-Meier curves
#' @param xpos_bs The horizontal position of the label for the Brier score
#' @param x_label The label of the x-axis for the Kaplan-Meier curves
#'
#' @return The staging result
#'
#' @importFrom stringr str_replace_all
#' @importFrom igraph simplify plot.igraph graph layout_as_tree
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics title legend
#' @importFrom visNetwork visNetwork visNodes visOptions visEdges visLayout visLegend visSave
#' @importFrom survminer ggsurvplot ggsurvtable
#' @importFrom ggplot2 annotate theme_bw theme element_text
#' @importFrom tidyr drop_na unite
#' @importFrom utils as.roman
#' @export

runOpera <- function(ncat, dat, TimeN, yN, cenN, covN = NULL, withCov = F,
                     filepath = "./operawp_results/",
                     riskVariables,
                     riskFactors,
                     riskNames,
                     useLassoT = F, useOPERA = F, usePruning = F,
                     edges = NULL, ifE = F, type = "surv",
                     getsAll = F, useAIC = F, useIbs = F, useLRT = T,
                     prefix_stage = NULL,  threshold = NULL,
                     coarse_pruning = F, fine_pruning = F,
                     fine_pruning_quad = F,
                     eps = 10^-4, eps_lasso = 10^-2, maxiter = 10, maxiter_lasso = 200,
                     GIC = NULL, useBIC = F, perc = 0.1, minObs = 30,
                     plt = T, pic_width = 8, pic_height = 8, legend_size = 0.7,
                     x_pos = -1, yratio = 0.1,  yaxis_min = 0,
                     xpos_bs = 0.02, x_label = "Time (Days)", ...){
  # @import Rcpp grDevices graphics stats utils colorspace cowplot devtools doParallel plyr dplyr foreach ggplot2 htmlwidgets igraph Matrix networkD3 prodlim quadprog readr reshape2 reticulate visNetwork survival ipred R.utils tidyr survminer stringr

  print(paste0("The input dataset has a total number of ",  nrow(dat), " observations."))

  if(type == "surv"){
    dat  <- dat[,  c(TimeN, cenN, riskVariables, covN), drop = F]  %>% drop_na() %>% arrange(TimeN)
    dat  <- dat[order(dat[, TimeN]),]
  }

  if(type == "bin"){
    dat  <- dat[,  c(yN, riskVariables, covN), drop = F] %>% drop_na()
  }
  print(paste0("The input dataset has a total number of ",  nrow(dat), " observations after removing missing values."))

  for(v in 1:length(riskVariables)){
    dat[, LETTERS[v]] <- paste0(letters[v],
                               as.numeric(factor(dat[, riskVariables[v]], ordered = TRUE, levels = riskFactors[[v]])))
  }

  dat <- dat[, LETTERS[1:length(riskVariables)]] %>% unite("variable", sep = "", remove = TRUE) %>% cbind(dat)

  colnames(dat) <- str_replace_all(colnames(dat), "[^[:alnum:]]", "_")

  if(type == "surv"){
    TimeN <- str_replace_all(TimeN, "[^[:alnum:]]", "_")
    cenN <- str_replace_all(cenN, "[^[:alnum:]]", "_")

  }

  if(type == "bin"){
    yN <- str_replace_all(yN, "[^[:alnum:]]", "_")
  }


  if(withCov){
    covN <-  str_replace_all(covN, "[^[:alnum:]]", "_")
  }


  suppressWarnings(dir.create(filepath, recursive = T))

  if(is.null(GIC)){
    GIC <- 2
  }


  if(withCov){
    Cov  <- as.matrix(dat[, covN, drop = F])
  }else{
    Cov <- matrix(0, 1, 1)
  }

  input_cov <- Cov

  if(type == "surv"){
    Times <- dat[, TimeN]
    cen <- dat[, cenN]
    y <- dat[, TimeN]
    dat$time <- dat[, TimeN]
    dat$status <- dat[, cenN]
  }

  if(type == "bin"){
    cen <- dat[, yN]
    y <- dat[, yN]
    dat$outcome <- dat[, yN]
  }


  if(!ifE){
    edges <- edgesHasse(ncat = ncat, e = c(), vs = matrix(base::rep(1, length(ncat)), 1, length(ncat)))
  }

  N <- data.frame(id = unique(edges), label = unique(edges))

  Hasse_ABC = igraph::simplify(graph(
    edges = edges, directed = T))

  if(plt){
    pdf(file = paste0(filepath, "tree_like_structure.pdf"),
        width = pic_width,
        height = pic_height)
  }

  plot.igraph(Hasse_ABC, layout = layout_as_tree, ...)
  suppressWarnings(title("The Tree-like Structure of Risk Factors", ...))

  ypos <- 1

  for(k in 1:length(ncat)){
    legend(x = x_pos, y = ypos, title = riskNames[k], legend = paste0(1:ncat[k], ": ", riskFactors[[k]]), horiz = FALSE, cex =  legend_size, box.lty = 0)
    ypos <- ypos - ncat[k]*yratio
  }

  if(plt){
    dev.off()
  }

  N$title <- sapply(N$id, function(x) NTitle(x = x,
                                             riskNames = riskNames,
                                             riskFactors = riskFactors))

  E <- data.frame(
    from = edges[seq_along(edges) %% 2 > 0],
    to = edges[seq_along(edges) %% 2 == 0]
  )

  E <- E %>% dplyr::distinct(from, to, .keep_all = T)

  network <- visNetwork(N, E, width = "100%", height = "1000px") %>%
    visNodes(
      color = list(background = "lightblue",
                   border = "darkblue",
                   highlight = "yellow")) %>%
    visOptions(highlightNearest = TRUE) %>%
    visEdges(
      arrows = "middle",
      color = list(color = "blue", highlight = "red")) %>%
    visLayout(randomSeed = 10)

  if(plt){
    visSave(network, paste0(filepath, "tree-like_structure.html"))
  }

  dat <- createRiskGroups(dat = dat, numRiskFs = length(ncat))
  input_PO <- as_inequalities(Hasse = Hasse_ABC, cellName = N$id)
  input_Z <- Zmatrix(dat,  N$id)

  if(useLassoT){

    resultNoP <- list(lasso.tree.bic(ncat = ncat, cen = cen, y = y, Apo = input_PO, label = N$id, dat = dat,
                                     inibeta = NULL, maxiter = maxiter_lasso, eps = eps_lasso, type = type,
                                     giveNB = ifE, useBIC = useBIC,  ...))

    dat[,"cStage"] <- resultNoP[[1]][as.character(dat[,"variable"])] - 1


    counts <- xtabs(~cStage, data = dat)

    if( ((sum(counts < as.integer(perc*nrow(dat))) > 0) & (sum(counts < minObs) > 0))){

      resultNoP[[2]] <- operap(cen = cen, y = y, Z = input_Z, Cnstrn = input_PO, Cov = Cov,
                      Data = dat, withCov = withCov, maxiter = maxiter,
                      eps = eps, GIC = GIC, initial_result = resultNoP[[1]],
                      useIbs = F, useLRT = F, useAIC = T,
                      coarse_pruning = F, fine_pruning = F,
                      fine_pruning_quad = F, prefix_stage = NULL,
                      threshold = NULL, type = type, perc = perc, minObs = minObs)[[1]]
    }else{
      resultNoP[[2]] <- resultNoP[[1]]
    }

  }

  if(useOPERA){
    resultNoP <- operai(Data = dat, TimeN = TimeN, cenN = cenN, yN= yN, covN = covN, Z = input_Z, Cnstrn = input_PO,
                  withCov = withCov, maxiter = maxiter, eps = eps,
                  type = type, GIC = GIC, perc = perc, minObs = minObs, ...)

  }

  if(usePruning & max(resultNoP[[2]]) > 2){
    resultP <- operap(cen = cen, y = y, Z = input_Z, Cnstrn = input_PO, Cov = Cov, Data = dat,
                                  withCov = withCov, maxiter = maxiter, eps = eps, GIC = GIC,
                                  initial_result = resultNoP[[1]], type = type, perc = perc, minObs = minObs,
                                  getsAll = getsAll, useAIC = useAIC,
                                  useIbs = useIbs, useLRT = useLRT,
                                  prefix_stage = prefix_stage,
                                  threshold = threshold, coarse_pruning = coarse_pruning,
                                  fine_pruning = fine_pruning, fine_pruning_quad = fine_pruning_quad, ...)


    finalR <- list("Before Pruning" = resultNoP[[2]], "After Pruning" = resultP)

    if(useOPERA){
      filenames <- c("No_Pruning_with_OPERA")
    }

    if(useLassoT){
      if(useBIC){
        filenames <- c("No_Pruning_with_Lasso_Tree_with_BIC")
      }else{
        filenames <- c("No_Pruning_with_Lasso_Tree_with_AIC")
      }
    }

    results <- list(finalR[["Before Pruning"]])



    if(coarse_pruning){
      if(useAIC){
        filenames <- c(filenames, "Coarse_Pruning_Using_AIC_or_Its_Equilvalence")
      }
      if(useIbs){
        filenames <- c(filenames, "Coarse_Pruning_Using_the_Brier_Score")
      }

      if(useLRT & !is.null(threshold)){
        filenames <- c(filenames, paste0("Coarse_Pruning_Using_the_LRT_with_Alpha_", threshold))
      }

      if(!is.null(prefix_stage)){
        filenames <- c(filenames, paste0("Coarse_Pruning_Using_the_Prespecied_Number_of_Stages_", prefix_stage))
      }
    }

    if(fine_pruning){
      if(useAIC){
        filenames <- c(filenames, "Fine_Pruning_with_Exhaustive_Search_Using_AIC_or_Its_Equivalence")
      }
      if(useIbs){
        filenames <- c(filenames, "Fine_Pruning_with_Exhaustive_Search_Using_the_Brier_Score")
      }
      if(useLRT & !is.null(threshold)){
        filenames <- c(filenames, paste0("Fine_Pruning_with_Exhaustive_Search_Using_the_LRT_with_Alpha_", threshold))
      }
      if(!is.null(prefix_stage)){
        filenames <- c(filenames, paste0("Fine_Pruning_with_Exhaustive_Search_Using_the_Prespecied_Number_of_Stages_", prefix_stage))
      }
    }

    if(fine_pruning_quad){
      if(useAIC){
        filenames <- c(filenames, "Fine_Pruning_with_the_Quadratic_Constraint_Using_AIC_or_Its_Equivalence")
      }
      if(useIbs){
        filenames <- c(filenames, "Fine_Pruning_with_the_Quadratic_Constraint_Using_the_Brier_Score")
      }
      if(useLRT & !is.null(threshold)){
        filenames <- c(filenames, paste0("Fine_Pruning_with_the_Quadratic_Constraint_Using_the_LRT_with_Alpha_", threshold))
      }
      if(!is.null(prefix_stage)){
        filenames <- c(filenames, paste0("Fine_Pruning_with_the_Quadratic_Constraint_Using_the_Prespecied_Number_of_Stages_", prefix_stage))
      }
    }



    if(getsAll & !is.null(threshold)){
      if(coarse_pruning){
        results[[length(results)+1]] <- finalR[["After Pruning"]][["The result using AIC"]]
      }else{
        results[[length(results)+1]] <- finalR[["After Pruning"]][[1]][[which.min(finalR[["After Pruning"]][["The AICs or equivalances"]])]]

      }
      results[[length(results)+1]] <- finalR[["After Pruning"]][[1]][[which.min(finalR[["After Pruning"]][["The brier scores" ]])]]
      results[[length(results)+1]] <- finalR[["After Pruning"]][[1]][[suppressWarnings(min(min(which(finalR[["After Pruning"]][["The p-values"]][-c(1)] <= threshold)+1)-1, length(finalR[["After Pruning"]][["The p-values"]])))]]
    }else{
      results[[length(results)+1]] <- finalR[["After Pruning"]][[1]]
    }

  }else{
    finalR <- list("Before Pruning" = resultNoP[[2]])

    if(useOPERA){
      filenames <- c("No_Pruning_with_OPERA")
    }

    if(useLassoT){
      if(useBIC){
        filenames <- c("No_Pruning_with_Lasso_Tree_with_BIC")
      }else{
        filenames <- c("No_Pruning_with_Lasso_Tree_with_AIC")
      }
    }

      results <- list(finalR[["Before Pruning"]])
  }

  for(j in 1:length(filenames)){

    resultClassification <- results[[j]]
    dat <- dat %>% mutate(stage = resultClassification[variable])

    if(plt){
      pdf(file = paste0(filepath, filenames[j], "_tree_like_structure.pdf"),   # The directory you want to save the file in
          width = pic_width,
          height =  pic_height)
    }

    plot.igraph(Hasse_ABC, layout = layout_as_tree,
                vertex.color =  sapply(resultClassification, function(x){sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE)[x]}), ...)

    suppressWarnings(title(str_replace_all(filenames[j], "_", " "), ...))

    ypos <- 1

    for(k in 1:length(ncat)){
      legend(x = x_pos, y = ypos, title = riskNames[k], legend = paste0(1:ncat[k], ": ",  riskFactors[[k]]), horiz=FALSE, cex=legend_size, box.lty = 0)
      ypos <- ypos - ncat[k]*yratio
    }

    legend("bottomright", legend = sort(unique(resultClassification)),
           col = sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE),
           horiz = FALSE, pch = base::rep(15,length(unique(resultClassification))), cex = legend_size,  bty = "n", title="The Classified Stages")

    if(plt){
      dev.off()
    }

    N$title <- sapply(N$id, function(x) NTitleS(x = x, riskNames = riskNames, riskFactors = riskFactors, S = resultClassification))

    N$color.background <- sapply(N$id, function(x){sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE)[resultClassification[x]]})
    N$color.highlight.background <- sapply(N$id, function(x){sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE)[resultClassification[x]]})

    E <- data.frame(
      from = edges[seq_along(edges) %% 2 > 0],
      to = edges[seq_along(edges) %% 2 == 0]
    )

    LegendN <- data.frame(label = as.character(as.roman(1:length(unique(resultClassification)))), shape = "circle",
                          color.highlight.background = sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE),
                          color.background = sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE),
                          size = 50)

    network <- visNetwork(N, E, width = "100%", height = "1000px") %>%
      visOptions(highlightNearest = TRUE) %>%
      visEdges(
        arrows = "middle",
        color = list(color = "blue", highlight = "red")) %>%
      visLayout(randomSeed = 10) %>% visLegend(useGroups = FALSE, addNodes = LegendN, main = "Stage", width = 0.1) %>% visNodes(font = list(size = 30))

    if(plt){
      visSave(network, paste0(filepath, filenames[j], "_tree_like_structure.html"))
    }


    if(plt){

      if(type == "surv"){

        dList <- list()
        dList$time <- dat[, TimeN]
        dList$status <- dat[, cenN]
        dList$stage <- dat[, "stage"]

        if(withCov){
          covName <- colnames(input_cov)
          for(k in 1:length(covName)){
            dList[[covName[k]]] <- dat[, covName[k]]
          }
          dfit <- coxph(Surv(time, status) ~ . + strata(stage), data = dList)
        }else{
          dfit <- coxph(Surv(time, status) ~ strata(stage), data = dList)
        }

        pred <- survfit(dfit)

        ibs <- sbrierC(Surv(dat[,TimeN], dat[,cenN]),
                       as.vector(lapply(dat[,"stage"],
                                        function(x) pred[as.numeric(x)])) )

        pdf(file = paste0(filepath, filenames[j], "_K-M.pdf"),
            width = 8,
            height = 8)
        survp <- ggsurvplot(
          fit = survfit(Surv(time, status) ~ stage,
                        data = dat, conf.type = "log-log"),
          data = dat,
          surv.median.line = "hv",
          font.x = c(15, "bold.italic", "black"),
          font.y = c(15, "bold.italic", "black"),
          xlab = x_label,
          ylab = "The Estimated Survival Probabilities",
          ylim = c(yaxis_min,1),
          censor = T,
          size = 1,                 # change line size
          palette = sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE),  # custom color palettes
          conf.int = T,          # Add confidence interval
          pval = TRUE,              # Add p-value
          pval.size = 5,            # p-value text size
          pval.method =  TRUE,      # test name
          risk.table = T,        # Add risk table
          legend = "none",
          title = str_replace_all(filenames[j], "_", " "),
          font.tickslab = c(15),
          ggtheme = theme_bw() + theme(
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size =  12)
          )
        )

        if(max(dat$stage) > 1){
          survp$table <- ggsurvtable(survfit(Surv(time, status) ~ stage,
                                             data = dat, conf.type = "log-log"),
                                     data = dat,
                                     color = "stage",
                                     survtable = "risk.table",
                                     palette = sequential_hcl(palette = "Hawaii", n = length(unique(resultClassification)), rev = TRUE),
                                     fontsize = 3,
                                     risk.table.type = "abs_pct",
                                     xlab = x_label,
                                     ylab = "",
                                     legend.title = "",
                                     font.tickslab = c(3, "bold"))
        }else{
          survp$table <- ggsurvtable(survfit(Surv(time, status) ~ stage,
                                             data = dat, conf.type = "log-log"),
                                     data = dat,
                                     survtable = "risk.table",
                                     risk.table.type = "abs_pct",
                                     fontsize = 3,
                                     xlab = x_label,
                                     ylab = "",
                                     legend.title = "",
                                     font.tickslab = c(3, "bold"))
        }

        survp$plot <- survp$plot+
          ggplot2::annotate("text",
                            x = xpos_bs, y = 0.05, # x and y coordinates of the text
                            label = paste0("The Brier Score = ", round(ibs, -log10(eps))), size = 5)
        print(survp, newpage = FALSE)
        dev.off()
      }
    }

  }

  return(finalR)
}
