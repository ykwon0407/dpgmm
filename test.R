set.seed(100)
Rcpp::sourceCpp('~/project/bayesian/generator.cpp')
list_data = generator_cpp(5, A = 0.6, tol=0.0001)
theta=list_data[[1]]
data=list_data[[2]]

# posterior 계산해보고
# algorithm별로 차이점 생각한뒤
# code 짜기
# posteriror distribution 구해서 histogram 그리기