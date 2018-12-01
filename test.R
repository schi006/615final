setwd("C:/1.1/bios615/project")
a=c(3,2)
b=c(5,1)
coxphSGDprepare <- function(formula, data) {
  # Parameter identification as in  `survival::coxph()`.
  Call <- match.call()
  indx <- match(c("formula", "data"),
                names(Call), nomatch = 0)
  print(Call)
  print('\n')
  print('\n')
  print(indx)
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  
  mf <- eval(temp, parent.frame())
  #print(mf)
  Y <- model.extract(mf, "response")
  print(Y)
  print(unclass(Y))
}
coxphSGDprepare(data=2,fomula=a~b,data=1,data=8)


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
a <- matrix(rnorm(50), 25)
b <- matrix(runif(100), 25)
cbind(a,b)
cbindC(b,a)
