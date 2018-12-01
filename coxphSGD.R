####################################################################################################
#cbind
cbindCpp <- 'NumericMatrix cbindC(NumericMatrix mf, NumericMatrix Y) {
  int acoln = mf.ncol();
  NumericMatrix out = no_init_matrix(Y.nrow(), acoln + 2);
  for (int j = 0; j < acoln + 2; j++) {
    if (j < 2) {
      out(_, j) = Y(_, 1-j);
    } else {
      out(_, j) = mf(_, j - 2);
    }
  }
  return out;
}'
cppFunction(cbindCpp)
#matrix multiplication
#install.packages("inline")
#library("inline")
matMultC<-cxxfunction(signature(a="NumericMatrix",
                                b="NumericMatrix"),
                         plugin="RcppEigen",
                         body="
NumericMatrix aa(a);
NumericMatrix bb(b);
const Eigen::Map<Eigen::MatrixXd> am(as<Eigen::Map<Eigen::MatrixXd> >(aa));
const Eigen::Map<Eigen::MatrixXd> bm(as<Eigen::Map<Eigen::MatrixXd> >(bb));
Eigen::MatrixXd c=am*bm;
return(wrap(c));
")
####################################################################################################
coxphSGDprepare <- function(formula, data) {
  # Parameter identification as in  `survival::coxph()`.
  Call <- match.call()                                  #保存调用函数时的输入
  indx <- match(c("formula", "data"),
                names(Call), nomatch = 0)               #检查输入是否符合函数的要求
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  
  mf <- eval(temp, parent.frame())                      #eval判断输入是否合法，mf储存数据
  Y <- model.extract(mf, "response")
  #这一部分判断Y是否合法
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  
  if (type != "right" && type != "counting")
    stop(paste("Cox model doesn't support \"", type, "\" survival data",
               sep = ""))
  
  # collect times, status, variables and reorder samples
  # to make the algorithm more clear to read and track
  cbindC(as.matrix(mf[,-1]),Y) -> data2return                        #mf的第一列应该是Y，之后是X
  #colnames(data2return)[1:2]<-c("event","times")
}

coxphSGDbatch <- function(formula, data, learning.rate, beta){
  
  # collect times, status, variables and reorder samples
  # to make the algorithm more clear to read and track
  batchData <- coxphSGDprepare(formula = formula, data = data) # sorts times lol
  # calculate the log-likelihood for this batch sample

########################################################################################  
  batchData <- batchData[order(-batchData[,2]), ] # dplyr::arrange(-times) / sorts time again but with different order

  # scores occure in nominator and denominator
  scores <- apply(batchData[, -c(1, 2)], 1,
                  function(element) exp(element %*% beta) )
  nominator <- apply(batchData[, -c(1, 2)], 2,
                     function(element) cumsum(scores*element) )
  denominator <- cumsum(scores)
  # sum over non-censored observations
  partial_sum <- (batchData[, -c(1, 2)] - nominator/denominator)*batchData[,1]
  # each column indicates one explanatory variable
  U_batch <- colSums(partial_sum)
  return(beta + learning.rate * U_batch)
}

coxphSGDcheck <- function(formula, data, learn.rates,
                          beta.zero, epsilon) {
  
  stopifnot(is.list(data) & length(data) > 0)
  stopifnot(length(unique(unlist(lapply(data, ncol)))) == 1)
  # + check names and types for every variables
  stopifnot(is.function(learn.rates))
  stopifnot(is.numeric(epsilon))
  stopifnot(is.numeric(beta.zero))
  
  # check length of the start parameter
  if (length(beta.zero) == 1) {
    beta.zero <-
      rep(beta.zero,
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
  
  return(beta.zero)
}



coxphSGD <- function(formula, data, learn.rates = function(x){1/x},
                     beta.zero = 0, epsilon = 1e-5, max.iter = 500,
                     verbose = FALSE) {
  # check arguments
  beta_start <-
    coxphSGDcheck(
      formula,
      data,
      learn.rates,
      beta.zero,
      epsilon
    )
  n <- length(data)
  diff <- epsilon + 1
  i <- 1
  beta_new <- list()     # steps are saved in a list so that they can
  beta_old <- beta_start # be traced in the future
  # estimate
  while(i <= max.iter & diff > epsilon) {
    beta_new[[i]] <-
      unlist(
        coxphSGDbatch(
          formula = formula,
          beta = beta_old,
          learning.rate = learn.rates(i),
          data = data[[ifelse(i%%n==0, n, i%%n)]]
        )
      )
    
    diff <- sqrt(sum((beta_new[[i]] - beta_old)^2))
    beta_old <- beta_new[[i]]
    i <- i + 1
    if (verbose) {
      cat("\r iteration: ", i, "\r")
    }
  }
  # return results
  list(
    Call = match.call(),
    epsilon = epsilon,
    learn.rates = learn.rates,
    steps = i,
    coefficients = c(list(beta_start), beta_new)
  )
}
