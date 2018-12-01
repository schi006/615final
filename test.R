setwd("C:/Users/10636/Documents/GitHub/615final")
a=c(3,2)
b=c(5,1)
coxphSGDprepare <- function(formula, data) {
  # Parameter identification as in  `survival::coxph()`.
  Call <- match.call()
  indx <- match(c("formula", "data"),
                names(Call), nomatch = 0)
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  
  mf <- eval(temp, parent.frame())
  #print(mf)
  mfm<-as.matrix(mf[,-1])
  print(is.matrix(mf))
  Y <- model.extract(mf, "response")
  print(is.matrix(Y))
  #print(Y)
  #print(unclass(Y))
}
coxphSGDprepare(formula     = Surv(time, status) ~ x.1+x.2,
                data        = dCox)


library(Rcpp)
cppFunction('void sort1129(NumericVector A){
            std::sort(A.begin(),A.end());
          }')


library(Rcpp)
?cppFunction
cppFunction('NumericVector sort1129(NumericVector v){
 NumericVector idx(v.size());
 itoa(idx.begin(),idx.end(),0);
 std::sort(idx.begin(),idx.end(),
 [&v](size_t i1,size_t i2){return v[i1]<v[i2];});
 return idx;
}')



cppFunction('int add(int x, int y, int z){
           int sum=x+y+z;
           return sum;
           }')

set.seed(42)
a <- matrix(rnorm(2186), 1093)
b <- matrix(runif(66673), 1093)
#cbind(a,b)
#cbindC(b,a)
library(microbenchmark)
microbenchmark(cbind(a,b),cbindC(a,b))

sourceCpp("orderC.cpp")
r<-runif(1093)
microbenchmark(orderC(r),order(r))

coxphSGDbatch(formula     = Surv(time, status) ~ x.1+x.2,
                     data        = dCox,
                     learning.rate = function(x){1/(100*sqrt(x))},
                     beta   = c(0,0))


a <- matrix(rnorm(1093), 1093)
b <- matrix(runif(1093), 1)+10
c<-matMultC(a,b)
d<-
microbenchmark(apply(b, 1,function(element) exp(element %*% a) ),matMultC(b,a),b%*%a)
identical(apply(b, 1,function(element) exp(element %*% a) ),matMultC(b,a))

M1 <- matrix(sample(1e3),ncol=50)
M2 <- matrix(sample(1e3),nrow=50)

identical(matMultC(M1,M2), M1 %*% M2)
microbenchmark(
  +   matMultC(M1,M2),
  +   M1 %*% M2,
  +   times=10000L)