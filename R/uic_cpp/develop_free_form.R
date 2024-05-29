# clear R working global environment
rm(list = ls())

# C++ testing
## Compilation dependencies
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

## 1. Element wise multiplication 
## 1.1. using Rcpp:
code1 <- "Rcpp::NumericVector element_multip1(const Rcpp::NumericVector& v_one,
                                              const Rcpp::NumericVector& v_two) {
    if (v_one.size() == 1 && v_two.size() != 1) {
      double tmp = v_one[0];
      return tmp * v_two;
    } else if (v_one.size() != 1 && v_two.size() == 1) {
      double tmp = v_two[0];
      return tmp * v_one;
    } else {
      return v_one * v_two;
    }
  }
"
### compile:
Rcpp::cppFunction(
  code = code1
)
## 1.2. using arma:
code2 <- "arma::vec element_multip2(const arma::vec& v_one, 
                                    const arma::vec& v_two) {
  if (v_one.size() == 1 && v_two.size() != 1) {
    double tmp = v_one[0];
    return tmp * v_two;  // Scalar multiplication
  } else if (v_one.size() != 1 && v_two.size() == 1) {
    double tmp = v_two[0];
    return tmp * v_one;  // Scalar multiplication
  } else {
    return v_one % v_two;  // Element-wise multiplication
  }
}"
### compile:
Rcpp::cppFunction(
  code = code2,
  depends = "RcppArmadillo"
)
## 1.3. using arma returning a vector
code3 <- "Rcpp::NumericVector element_multip3(const arma::vec& v_one, 
                                              const arma::vec& v_two) {
  arma::vec result;
  if (v_one.size() == 1 && v_two.size() != 1) {
    double tmp = v_one[0];
    result = tmp * v_two;  // Scalar multiplication
  } else if (v_one.size() != 1 && v_two.size() == 1) {
    double tmp = v_two[0];
    result = tmp * v_one;  // Scalar multiplication
  } else {
    result = v_one % v_two;  // Element-wise multiplication
  }
  return Rcpp::wrap(result);  // Convert arma::vec to Rcpp::NumericVector
}"
### compile:
Rcpp::cppFunction(
  code = code3,
  depends = "RcppArmadillo"
)
### test:
element_multip1(v_one = rep(2, 10), v_two = 1:10)
element_multip1(v_one = 2, v_two = 1:10)
element_multip2(v_one = rep(2, 10), v_two = 1:10)[,1]
element_multip2(v_one = 2, v_two = 1:10)[,1]
microbenchmark::microbenchmark(
  "rcpp11" = element_multip1(v_one = rep(2, 10), v_two = 1:10),
  "rcpp12" = element_multip1(v_one = 2, v_two = 1:10),
  "arma11" = element_multip2(v_one = rep(2, 10), v_two = 1:10)[,1],
  "arma12" = element_multip2(v_one = 2, v_two = 1:10)[,1],
  "arma21" = element_multip2(v_one = rep(2, 10), v_two = 1:10),
  "arma22" = element_multip2(v_one = 2, v_two = 1:10)
  )

## 2. Matrix Multiplication using Rcpp:
code <- 'NumericMatrix matrix_mult1(const NumericMatrix A, const NumericMatrix B) {
  // Get the number of rows and columns for the resulting matrix
  int nrow = A.nrow();
  int ncol = B.ncol();
  
  // Create an empty NumericMatrix to store the result
  NumericMatrix C(nrow, ncol);
  
  // Perform matrix multiplication
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double sum = 0;
      for (int k = 0; k < A.ncol(); k++) {
        sum += A(i, k) * B(k, j);
      }
      C(i, j) = sum;
    }
  }
  
  return C; // Return the result matrix
}'
### compile:
Rcpp::cppFunction(
  code = code,
  includes = 'using namespace Rcpp;'
)
## 2. Matrix Multiplication using arma:
code <- 'arma::mat matrix_mult2(const arma::mat& A, const arma::mat& B) {
  // Perform matrix multiplication using Armadillo
  arma::mat C = A * B;
  return C;
}'
### compile:
Rcpp::cppFunction(
  code = code,
  depends = "RcppArmadillo"
)
### test:
m_one <- matrix(1,   nrow = 4, ncol = 4)
m_two <- matrix(1:16, nrow = 4, byrow = T)
matrix_mult1(A = m_one, B = m_two)
matrix_mult2(A = m_one, B = m_two)
m_one %*% m_two
crossprod(m_one, m_two) # surely works with multiplication of square matrices
### compare:
m_one <- matrix(1,   nrow = 4, ncol = 4)
m_two <- matrix(1:16, nrow = 4, byrow = T)
mat_mult_RvC <- microbenchmark::microbenchmark(
  "C1" = matrix_mult1(A = m_one, B = m_two),
  "C2" = matrix_mult2(A = m_one, B = m_two),
  "R1" = m_one %*% m_two
)
m_one <- matrix(1,   nrow = 1e5, ncol = 4)
m_two <- matrix(1:16, nrow = 4, byrow = T)
mat_mult_RvC2 <- microbenchmark::microbenchmark(
  "C1" = matrix_mult1(A = m_one, B = m_two),
  "C2" = matrix_mult2(A = m_one, B = m_two),
  "R1" = m_one %*% m_two
)
m_one <- matrix(1,   nrow = 1e5, ncol = 10)
m_two <- matrix(1:16, nrow = 10, byrow = T)
mat_mult_RvC3 <- microbenchmark::microbenchmark(
  "C1" = matrix_mult1(A = m_one, B = m_two),
  "C2" = matrix_mult2(A = m_one, B = m_two),
  "R1" = m_one %*% m_two
)
