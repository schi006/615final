batchCoxprepare <- function(formula, data) {
  # Parameter identification as in  `survival::coxph()`.
  Call <- match.call()                                  #save input
  indx <- match(c("formula", "data"),
                names(Call), nomatch = 0)               #check if input is Ok
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")

  mf <- eval(temp, parent.frame())                      #eval check if input is Ok, and mf save data
  Y <- model.extract(mf, "response")
  #the following check if Y is Ok
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival formula")
  type <- attr(Y, "type")

  if (type != "right" && type != "counting")
    stop(paste("Cox model doesn't support \"", type, "\" survival data",
               sep = ""))

  dataprepare<-cbindC(as.matrix(mf[,-1]),Y)                          #cbind X and Y;fist column of mf is Y
  colnames(dataprepare)<-c("event","times",names(mf[,-1]))
  return(dataprepare)
}

batchCoxbatch <- function(formula, data, learning.rate, beta){       #core function

  batchData <- batchCoxprepare(formula = formula, data = data)
  # calculate the log-likelihood for this batch sample

  batchData <- batchData[order(-batchData[,2]), ] #sort time

  # scores occure in nominator and denominator
  scores <- elementE(matMultC(batchData[, -c(1, 2)], beta) )
  nominator <- nominatorC(batchData[, -c(1, 2)],scores)
  denominator <- cumsum(scores)
  # sum over non-censored observations
  partial_sum <- (batchData[, -c(1, 2)] - nominator/denominator)*batchData[,1]
  # each column indicates one explanatory variable
  U_batch <- colSums(partial_sum)
  names(U_batch)=colnames(batchData)[c(-1,-2)]
  return(beta + learning.rate * U_batch)
}

batchCoxcheck <- function(formula, data, rate,
                          initial, epsilon) {

  stopifnot(is.list(data) & length(data) > 0)
  stopifnot(length(unique(unlist(lapply(data, ncol)))) == 1)
  # check names and types for every variables
  stopifnot(is.function(rate))
  stopifnot(is.numeric(epsilon))
  stopifnot(is.numeric(initial))

  # check length of the start parameter
  if (length(initial) == 1) { #if input one value for beta，then beta is a vector of input beta
    initial <-
      rep(initial,
          length(
            unlist(
              strsplit(
                as.character(
                  formula
                )[3],
                split = "\\+")
            )
          )
      )

  }

  return(initial)
}



batchCox <- function(formula, data, rate = function(x){1/x},
                     initial = 0, epsilon = 1e-5, maxit = 500) {
  # check arguments
  beta_start <-
    batchCoxcheck(
      formula,
      data,
      rate,
      initial,
      epsilon
    )
  n <- length(data)
  diff <- epsilon + 1 #initial value of diff can be any value larger than epsilon
  i <- 1
  beta_new <- list()     # steps are saved in a list so that they can
  beta_old <- beta_start # be traced in the future
  # estimate
  while(i <= maxit & diff > epsilon) {
    beta_new[[i]] <-
      unlist(
        batchCoxbatch(
          formula = formula,
          beta = beta_old,
          learning.rate = rate(i),
          data = data[[ceiling(runif(1,0,n))]] #two disadvantages: 1.data are input together
        )                                      #here; 2.prepare function run everytime
      )

    diff <- sqrt(sum((beta_new[[i]] - beta_old)^2))
    if(is.na(diff))  break
    beta_old <- beta_new[[i]]
    i <- i + 1
  }
  # return results
  list(
    Call = match.call(),
    epsilon = epsilon,
    rate = rate,
    steps = i,
    coefficients = c(list(beta_start), beta_new)
  )
}
CoxData <- function(n, lambda, rho, x, beta, cens.rate){

  # real Weibull times
  u <- runif(n)
  Treal <- (- log(u) / (lambda * elementE(matMultC(x,beta))))^(1 / rho)

  # censoring times
  Censoring <- rexp(n, cens.rate)

  # follow-up times and event indicators
  time <- pmin(Treal, Censoring)
  status <- as.numeric(Treal <= Censoring)

  # data set
  data.frame(id = 1:n, time = time, status = status, x = x)
}
predict.batchCox <- function(object,newdata){
  # Parameter identification as in  `survival::coxph()`.
  # prepare x
  Call <- object$Call                                   #save input
  indx <- match(c("formula", "data"),
                names(Call), nomatch = 0)               #check if input is Ok
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  newdata$time<-0
  newdata$status<-0
  temp[[indx[2]]]<-newdata
  mf <- eval(temp, parent.frame())                      #eval check if input is Ok, and mf save data
  X<-as.matrix(mf[,-1])
  #prepare beta
  beta<-object$coefficients[[length(object$coefficients)]]
  matMultC(X,beta)
}
