setwd("/Users/sunyichi/Documents/GitHub/615/")

# simulation generate time to event data
libraray(survival)
dataCox <- function(n, lambda, rho, x, beta, cens.rate){
  
  # real Weibull times
  u <- runif(n)
  Treal <- (- log(u) / (lambda * exp(x %*% beta)))^(1 / rho)
  
  # censoring times
  Censoring <- rexp(n, cens.rate)
  
  # follow-up times and event indicators
  time <- pmin(Treal, Censoring)
  status <- as.numeric(Treal <= Censoring)
  
  # data set
  data.frame(id = 1:n, time = time, status = status, x = x)
}



coxphSGDprepare <- function(formula, data) {
  # Parameter identification as in  `survival::coxph()`.
  Call <- match.call()
  indx <- match(c("formula", "data"), names(Call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  
  mf <- eval(temp, parent.frame())
  Y <- model.extract(mf, "response")
  
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  
  if (type != "right" && type != "counting")
    stop(paste("Cox model doesn't support \"", type, "\" survival data",
               sep = ""))
  # collect times, status, variables and reorder samples
  # to make the algorithm more clear to read and track
  cbind(event = unclass(Y)[,2], # 1 indicates event, 0 indicates cens   //cbind=>c++: https://stackoverflow.com/questions/31913437/r-fast-cbind-matrix-using-rcpp
        times = unclass(Y)[,1],
        mf[, -1]) -> data2return
  data2return[order(data2return$times), ] # dplyr::arrange(times)   //order=>sort algorithm
}



coxphSGDbatch <- function(formula, data, learning.rate, beta){
  # collect times, status, variables and reorder samples
  # to make the algorithm more clear to read and track
  batchData <- coxphSGDprepare(formula = formula, data = data) # sorts times lol
  # calculate the log-likelihood for this batch sample
  batchData <- batchData[order(-batchData$times), ] # dplyr::arrange(-times) / sorts time again but with different order
  # scores occure in nominator and denominator
  
  scores <- apply(batchData[, -c(1, 2)], 1,
                  function(element) exp(element %*% beta)) 
  nominator <- apply(batchData[, -c(1, 2)], 2,
                     function(element) cumsum(scores*element) )   #apply=>c++ matrix calculation
  denominator <- cumsum(scores)  
  # sum over non-censored observations
  partial_sum <- (batchData[, -c(1, 2)] - nominator/denominator)*batchData[, "event"]  #c++ matrix
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


coxphSGD <- function(formula, data, learn.rates = function(x){1/x},  #learn.rates function 使learning rate随interation递减，除了1/x还有其他形式
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
          data = data[[ifelse(i%%n==0, n, i%%n)]]   # 
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

#############################simulation
set.seed(456)
x <- matrix(sample(0:1, size = 20000, replace = TRUE), ncol = 2)
dCox <- dataCox(10^4, lambda = 3, rho = 2, x,
                beta = c(2,2), cens.rate = 5)
batch_id <- sample(1:90, size = 10^4, replace = TRUE)
dCox_split <- split(dCox, batch_id)
results <-
  coxphSGD(formula     = Surv(time, status) ~ x.1+x.2,
           data        = dCox_split,
           epsilon     = 1e-5,
           learn.rates = function(x){1/(100*sqrt(x))},
           beta.zero   = c(0,0),
           max.iter    = 10*90)
coeff_by_iteration <-
  as.data.frame(
    do.call(
      rbind,
      results$coefficients
    )
  )




library(reshape2)
set.seed(456)
x <- matrix(sample(0:1, size = 20000, replace = TRUE), ncol = 2)
head(x)
dCox <- dataCox(10^4, lambda = 3, rho = 2, x,
                beta = c(2,2), cens.rate = 5)
head(dCox)

batch_id <- sample(1:90, size = 10^4, replace = TRUE)
dCox_split <- split(dCox, batch_id)
results <-
  coxphSGD(formula     = Surv(time, status) ~ x.1+x.2,
           data        = dCox_split,
           epsilon     = 1e-5,
           learn.rates = function(x){1/(100*sqrt(x))},
           beta.zero   = c(0,0),
           max.iter    = 10*90)

coeff_by_iteration <-
  as.data.frame(
    do.call(
      rbind,
      results$coefficients
    )
  )
head(coeff_by_iteration)


coxph_loglik <- function(beta, formula, data) {
  coxph(formula, init=beta, control=list('iter.max'=0), data =data)$loglik[2]
}
coxph_loglik <- function(beta, formula, data) {
  coxph(formula, init=beta, control=list('iter.max'=0), data =data)$loglik[2]
}
calculate_outer_cox_3 <- function(dCox){
  ## contours
  outer_res <- outer(seq(0,4, length = 25),
                     seq(0,4, length = 25),
                     Vectorize( function(beta1,beta2){
                       coxph_loglik(beta=c(beta1,beta2), Surv(time, status)~x.1+x.2-1, dCox)
                     } )
  )
  outer_res_melted <- melt(outer_res)
  outer_res_melted$Var1 <- as.factor(outer_res_melted$Var1)
  levels(outer_res_melted$Var1) <- as.character(seq(0,4, length = 25))
  outer_res_melted$Var2 <- as.factor(outer_res_melted$Var2)
  levels(outer_res_melted$Var2) <- as.character(seq(0,4, length = 25))
  outer_res_melted$Var1 <- as.numeric(as.character(outer_res_melted$Var1))
  outer_res_melted$Var2 <- as.numeric(as.character(outer_res_melted$Var2))
  return(outer_res_melted)
}
calculate_outer_cox_3(dCox) -> outerCox
save(outerCox, file = 'dev/outerCox.rda')
#d2ggplot <- coeff_by_iteration
beta.zero <- c(0,0)
solution <- c(2,2)
library(ggplot2)
ggplot() +
  stat_contour(aes(x=outerCox$Var1,
                   y=outerCox$Var2,
                   z=outerCox$value),
               bins = 40, alpha = 0.25) +
  geom_path(aes(coeff_by_iteration[['x.1']],
                coeff_by_iteration[['x.2']]),
            #group = d2ggplot$version,
            #colour = d2ggplot$version),
            size = 1) +
  theme_bw(base_size = 20) +
  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        legend.position = "top") +
  scale_colour_brewer(palette="Dark2",
                      name = 'Algorithm \n & Steps') +
  geom_point(aes(x = beta.zero[1], y = beta.zero[2]),
             col = "black",
             size = 4, shape = 17) +
  geom_point(aes(x = solution[1], y = solution[2]),
             col = "black", size = 4, shape = 15) +
  geom_point(aes(x = summary(coxph(Surv(time, status) ~ x.1+x.2, data = dCox))$coeff[1,1],
                 y = summary(coxph(Surv(time, status) ~ x.1+x.2, data = dCox))$coeff[2,1]),
             col = "black", size = 4, shape = 13) +
  xlab("X1") +
  ylab("X2") -> p


###################################real data application
testCox <- readRDS("testCox.rds")
trainCox<- readRDS("trianCox.rds")
class(trainCox)

length(trainCox)

head(
  trainCox[[1]][ # train data is a list of data frames that corresponds to batches
    c(7, 10),    # pick few rows
    c(210,302,356,898,911,1092:1093) # pick few columns, and survival outcomes
    ]
)

formulaSGD <- archivist::aread('MarcinKosinski/MasterThesis/064277e1c2a1fbea36d7d0ac518b9c8d')

# Surv(times, patient.vital_status ~ all posible genes
# how many genes are there in the whole formula?
length(
  strsplit(
    as.character(formulaSGD)[3],
    split = "+",
    fixed = T
  )[[1]]
)

coxphSGD(
  formula = formulaSGD,
  data = trainCox,
  learn.rates = function(x){1/x},
  max.iter = 490 # 98*5 / 5 full epoches
) -> model

library(survminer)
fit_braf <- survfit(
  Surv(times, patient.vital_status) ~ BRAF,
  data = do.call(rbind, trainCox)
)
ggsurvplot(
  fit_braf,
  palette = c("#E7B800", "#2E9FDF"),
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =             # change legend labels.
    c("No mutation in BRAF", "Mutation in BRAF"),    
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,2500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
)

fit_egfr <- survfit(
  Surv(times, patient.vital_status) ~ EGFR,
  data = do.call(rbind, trainCox)
)
ggsurvplot(
  fit_egfr,
  palette = c("#E7B800", "#2E9FDF"),
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =             # change legend labels.
    c("No mutation in EGFR", "Mutation in EGFR"),
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,2500),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 500,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)
















