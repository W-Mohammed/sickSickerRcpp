# Using arma::cube:----

## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.Trt")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.Trt")

# Transition matrix:
m_t_p <- matrix(
  data = c(1 - t_p["p.HS1"] - t_p["p.HD"], t_p["p.HS1"], 0, t_p["p.HD"],
           t_p["p.S1H"], 1 - t_p["p.S1H"] - t_p["p.S1S2"] - t_p["p.S1D"], 
           t_p["p.S1S2"],  t_p["p.S1D"], 0, 0, 1 - t_p["p.S2D"], t_p["p.S2D"], 0,
           0, 0, 1),
  nrow = n.s,
  ncol = n.s,
  byrow = TRUE
)
## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)
## C++ code:----
code <- 
  'arma::cube testV_Cpp( arma::rowvec v_S_t,
                        arma::mat& t_p,
                        int& n_I,
                        int& n_S) {
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube test(n_I, n_S, n_S) ;
  
  // assign probabilities from t_p to the slices
  for(arma::uword i = 0; i < test.n_slices; i++) {
  
    // Extract the transition probabilities from health state i as a row vector:
    arma::rowvec t_p_i = t_p.row(i) ;
    
    // Replicate t_p_i n_I times vertically, and 1 time horizontally:
    test.slice(i) = arma::repmat(t_p_i, n_I, 1) ;
  }

  return(test) ;
}'

## Compile C++ function:----
Rcpp::cppFunction(
  code = code,
  depends = depends,
  #plugins = plugins,
  includes = includes
)

## Test Rcpp function:----
test = testV_Cpp(
  v_S_t =  v_M_1,
  t_p = matrix(
    data = rep(x = 1:4, each = n.s), 
    nrow = n.s, 
    ncol = n.s, 
    byrow = TRUE),
  n_I = n.i,
  n_S = n.s
)
# Using arma::cumsum:----
## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.Trt")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.Trt")

# Transition matrix:
m_t_p <- matrix(
  data = c(1 - t_p["p.HS1"] - t_p["p.HD"], t_p["p.HS1"], 0, t_p["p.HD"],
           t_p["p.S1H"], 1 - t_p["p.S1H"] - t_p["p.S1S2"] - t_p["p.S1D"], 
           t_p["p.S1S2"],  t_p["p.S1D"], 0, 0, 1 - t_p["p.S2D"], t_p["p.S2D"], 0,
           0, 0, 1),
  nrow = n.s,
  ncol = n.s,
  byrow = TRUE
)
## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)
## C++ code:----
code2 <- 
  'arma::mat test_Cpp( arma::mat m_P_t) {
    arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
    return(m_CumProb_t) ; 
}'

## Compile C++ function:----
Rcpp::cppFunction(
  code = code2,
  depends = depends,
  includes = includes
)

# Using custom findInterval:----
## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c('#include <algorithm>')
## C++ code:----
code3 <- '
Rcpp::IntegerVector findInterval2(Rcpp::NumericVector x, Rcpp::NumericVector breaks) {
  
  Rcpp::IntegerVector out(x.size());

  Rcpp::NumericVector::iterator it, pos;
  Rcpp::IntegerVector::iterator out_it;

  for(it = x.begin(), out_it = out.begin(); it != x.end(); 
      ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }

  return out;
}
'
code4 <- '
Rcpp::IntegerVector findInterval3(Rcpp::NumericVector x, Rcpp::NumericVector breaks) {
  
  Rcpp::IntegerVector out(x.size());

  Rcpp::NumericVector::iterator it, pos;
  Rcpp::IntegerVector::iterator out_it;
  
  // Store the iterator for breaks.begin() to avoid multiple calls
  Rcpp::NumericVector::iterator breaks_begin = breaks.begin();
  Rcpp::NumericVector::iterator breaks_end = breaks.end();

  for(it = x.begin(), out_it = out.begin(); it != x.end(); 
      ++it, ++out_it) {
    pos = std::upper_bound(breaks_begin, breaks_end, *it);
    *out_it = std::distance(breaks_begin, pos);
  }

  return out;
}
'
code5 <- '
Rcpp::IntegerVector findInterval4(Rcpp::NumericVector x, Rcpp::NumericVector breaks) {
  
  Rcpp::IntegerVector out(x.size());

  Rcpp::NumericVector::iterator it, pos;
  Rcpp::IntegerVector::iterator out_it;
  
  // Store the iterator for breaks.begin() to avoid multiple calls
  Rcpp::NumericVector::iterator breaks_begin = breaks.begin();
  Rcpp::NumericVector::iterator breaks_end = breaks.end();
  Rcpp::IntegerVector::iterator out_begin = out.begin();
  Rcpp::NumericVector::iterator x_begin = x.begin();

  for(it = x_begin, out_it = out_begin; it != x.end(); 
      ++it, ++out_it) {
    pos = std::upper_bound(breaks_begin, breaks_end, *it);
    *out_it = std::distance(breaks_begin, pos);
  }

  return out;
}
'
## Compile C++ function:----
cpp_functions_defs <- list(
  code3, code4, code5
)
for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    #depends = depends,
    includes = includes
  )
}
## Test function:----
x <- 2:18
v <- c(5, 10, 15) # create two bins [5,10) and [10,15)
test3R <- findInterval(x = x, vec = v)
test3Rcpp1 <- findInterval2(x = x, breaks = v)
test3Rcpp2 <- findInterval3(x = x, breaks = v)
test3Rcpp3 <- findInterval4(x = x, breaks = v)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "findInterval" = findInterval(x = x, vec = v),
  "findInterval2" = findInterval2(x = x, breaks = v), # fastest
  "findInterval3" = findInterval3(x = x, breaks = v),
  "findInterval4" = findInterval4(x = x, breaks = v)
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)


# Using arma::diagmat:----
## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.Trt")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.Trt")

# Transition matrix:
m_t_p <- matrix(
  data = c(1 - t_p["p.HS1"] - t_p["p.HD"], t_p["p.HS1"], 0, t_p["p.HD"],
           t_p["p.S1H"], 1 - t_p["p.S1H"] - t_p["p.S1S2"] - t_p["p.S1D"], 
           t_p["p.S1S2"],  t_p["p.S1D"], 0, 0, 1 - t_p["p.S2D"], t_p["p.S2D"], 0,
           0, 0, 1),
  nrow = n.s,
  ncol = n.s,
  byrow = TRUE
)
## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)
## C++ code:----
code6 <- 
  'arma::mat test_Cpp( arma::colvec v_S_t) {
    arma::mat m_costs_diag = arma::diagmat(v_S_t) ;
    return(m_costs_diag) ; 
}'

## Compile C++ function:----
Rcpp::cppFunction(
  code = code6,
  depends = depends#,
  #includes = includes
)
## Test function:----
v <- c(5, 10, 15)
test4R <- diag(v)
test4Rcpp <- test_Cpp(v)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "diagR" = diag(v),
  "diagCpp" = test_Cpp(v)
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)


# Sub-setting values:----
## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.Trt")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.Trt")

# Transition matrix:
m_t_p <- matrix(
  data = c(1 - t_p["p.HS1"] - t_p["p.HD"], t_p["p.HS1"], 0, t_p["p.HD"],
           t_p["p.S1H"], 1 - t_p["p.S1H"] - t_p["p.S1S2"] - t_p["p.S1D"], 
           t_p["p.S1S2"],  t_p["p.S1D"], 0, 0, 1 - t_p["p.S2D"], t_p["p.S2D"], 0,
           0, 0, 1),
  nrow = n.s,
  ncol = n.s,
  byrow = TRUE
)
## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)
## C++ code:----
code7 <- 
  'arma::colvec test7_Cpp( arma::vec& c_vec, arma::uvec indices) {
    arma::colvec v_costs = c_vec.elem(indices - 1) ;
    return(v_costs) ; 
}'
code8 <- 
  'arma::colvec test8_Cpp( arma::vec& c_vec, arma::vec& v_indices) {
    arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_indices);
    arma::colvec v_costs = c_vec.elem(uv_indices - 1) ;
    return(v_costs) ; 
}'
code9 <- 
  'arma::colvec test9_Cpp( arma::vec& c_vec, arma::vec& v_indices) {
    arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_indices - 1);
    arma::colvec v_costs = c_vec.elem(uv_indices) ;
    return(v_costs) ; 
}'
## Compile C++ function:----
cpp_functions_defs <- list(
  code7, code8, code9
)
for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    depends = depends,
    #includes = includes
  )
}
## Test function:----
v_indices <- c(1, 2, 3, 3, 1, 1)
test5R <- c_vec[v_indices]
test5Rcpp_1 <- test7_Cpp(c_vec, v_indices)
test5Rcpp_2 <- test8_Cpp(c_vec, v_indices)
test5Rcpp_3 <- test9_Cpp(c_vec, v_indices)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "sub-setting_R" = c_vec[v_indices],
  "sub-setting_RCpp1" = test7_Cpp(c_vec, v_indices),
  "sub-setting_RCpp2" = test8_Cpp(c_vec, v_indices),
  "sub-setting_RCpp3" = test9_Cpp(c_vec, v_indices)
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)

