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

# log and exp:----
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
v_M_1 <- sample(
  x = 1:n.s, 
  prob = rep(0.25, 4), 
  replace = TRUE, 
  size = n.i
)

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
code_time_in_state <- 
  'void time_in_state( arma::mat& m_t_states,
                        const arma::colvec& v_S_t,
                        const std::vector<int>& v_tracked_states) {
  // m_t_states: matrix storing times spent in each health state.
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_tracked_states: vector of health states in which to track time.

  // Track time if in tracked states:
  for(int s : v_tracked_states) {
    // Rescale states numerical identifier (s) to C++ index (i):
    int i = s - 1 ;
    
    // Create a logical vector where state matches s:
    arma::uvec other_states = arma::find(v_S_t != s) ;

    // Copy column data to a temporary vector:
    arma::vec state_col = m_t_states.col(i) ;
    
    // Increment the time for all individuals in the column.
    state_col += 1 ;
    
    // Ensure other_states has valid indices within the range of tmp_col
    if (!other_states.empty()) {
      // Reset time in other states to 0 using .elem():
      state_col.elem(other_states).fill(0) ;
    }
    
    // Record time back to states time matrix:
    m_t_states.col(i) = state_col ;
  }
}'
Rcpp::cppFunction(
  code = code_time_in_state,
  depends = depends,
  #includes = includes
)

## Code 10 - arma::log:----
code10 <- 
  'arma::mat test_log( arma::mat& m_t_states,
                      arma::mat& m_probs) {
  
  // m_t_states: duration spent in each health state.
  // m_probs: matrix containing transition probabilities.

  arma::mat tmp_m_rates = - arma::log(1 - m_probs) ;
  
  return(tmp_m_rates) ;

}'
## Code 11 - arma::exp:----
code11 <- 
  'arma::mat test_exp( arma::mat& m_t_states,
                      arma::mat& m_t_rates) {
                      
  // m_t_states: duration spent in each health state.
  // m_t_rates: matrix containing transition rates.
  
  arma::mat tmp_pob_mat = 1 - (arma::exp(- m_t_rates)) ;
  
  return(tmp_pob_mat) ;
}'
## Code 12 - arma::log and arma::exp:----
code12 <- 
  'arma::mat test_log_exp( arma::mat& m_t_states,
                      arma::mat& m_probs) {
  
  // m_t_states: duration spent in each health state.
  // m_probs: matrix containing transition probabilities.

  arma::mat tmp_pob_mat = - arma::log(1 - m_probs) ;
  tmp_pob_mat = 1 - (arma::exp(- tmp_pob_mat)) ;
  
  return(tmp_pob_mat) ;

}'

## Compile C++ function:----
cpp_functions_defs <- list(
  code10, code11, code12
)
for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    depends = depends,
    #includes = includes
  )
}
## Test function:----
m_t_r <- -log(1 - m_t_p)
m_t_states <- matrix(
  data = 0,
  nrow = n.i,
  ncol = n.s
)
time_in_state(
  m_t_states = m_t_states,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)

rates <- test_log(
  m_t_states = matrix(),
  m_probs = m_t_p
)
probs <- test_exp(
  m_t_states = matrix(),
  m_t_rates = m_t_r
)
rates_probs <- test_log_exp(
  m_t_states = matrix(),
  m_probs = m_t_p
)

## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "R_log" = -log(1 - m_t_p),
  "Rcpp_log" =  test_log(
    m_t_states = matrix(),
    m_probs = m_t_p
  ),
  "R_exp" = exp(m_t_r),
  "Rcpp_exp" = test_exp(
    m_t_states = matrix(),
    m_t_rates = m_t_r
  ),
  "R_log_exp" = 1 - exp(- -log(1 - m_t_p)),
  "Rcpp_log_exp" = test_log_exp(
    m_t_states = matrix(),
    m_probs = m_t_p
  )
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)

# discounting_weights:----
## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%

## Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)
## C++ code:----
## Code 13 - discounting_weights:----
code13 <- 
  'std::vector<double> estimate_discounting_weights ( 
                                const double discount_rate,
                                const double time_horizon) {
                                
  // Prepare the discount rate for the function:
  std::vector<double> base = as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::rep((1 + discount_rate), (time_horizon + 1)))} ) ;
  std::vector<double> expo = as< std::vector<double> >(Rcpp::NumericVector {wrap(Rcpp::seq(0, time_horizon))} ) ;

  std::vector<double> results(expo.size()) ;
  std::transform(base.begin(), base.end(), expo.begin(), results.begin(),
                 [&](double lhs, double rhs) -> double {
                   return (1/ std::pow(lhs, rhs)) ;
                 }) ;
  return(results);
}'
## Code 14 - discounting_weights:----
code14 <- 
  'std::vector<double> estimate_discounting_weights2(
                                                    const double discount_rate, 
                                                    const double time_horizon) {
    // Prepare object of size time_horizon + 1:                                                
    std::vector<double> weights(static_cast<size_t>(time_horizon) + 1);
    
    // Fill the vector with time periods 0 to time_horizon:
    std::iota(weights.begin(), weights.end(), 0);

    // Calculate discounting weights
    std::transform(weights.begin(), weights.end(), weights.begin(),
                   [discount_rate](double time) -> double {
                       return 1.0 / std::pow(1.0 + discount_rate, time);
                   });
    
    return weights;
}'
## Code 15 - discounting_weights:----
code15 <- 
  'arma::vec estimate_discounting_weights3(const double discount_rate, 
                                       const double time_horizon) {
    // discount_rate: the discount rate. 
    // time_horizon: the model time horizon.
     
    // Prepare object of size time_horizon + 1:
    arma::vec weights = arma::linspace(0, time_horizon, static_cast<size_t>(time_horizon) + 1);
    
    // Calculate discounting weights using element-wise operation
    weights.transform([discount_rate](double time) -> double {
        return 1.0 / std::pow(1.0 + discount_rate, time);
    });
    
    return(weights) ;
}'
## Compile C++ function:----
cpp_functions_defs <- list(
  code13, code14, code15
)
for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    depends = depends,
    #includes = includes
  )
}
## Test function:----
R_d_w <- 1 / (1 + d.e) ^ (0:n.t) 
Cpp_d_w <- estimate_discounting_weights(
  discount_rate = 0.03, 
  time_horizon = n.t
)
Cpp_d_w2 <- estimate_discounting_weights2(
  discount_rate = 0.03, 
  time_horizon = n.t
)
Cpp_d_w3 <- estimate_discounting_weights3(
  discount_rate = 0.03, 
  time_horizon = n.t
)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "R_d_w" =  1 / (1 + d.e) ^ (0:n.t) ,
  "Rcpp_d_w" = estimate_discounting_weights(
    discount_rate = 0.03, 
    time_horizon = n.t
  ),
  "Rcpp_d_w2" =  estimate_discounting_weights2(
    discount_rate = d.e, 
    time_horizon = n.t
  ),
  "Rcpp_d_w3" =  estimate_discounting_weights3(
    discount_rate = 0.03, 
    time_horizon = n.t
  )
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)

# MicroSim:----
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
rp.S1S2 <- 0.2                 # increase of the mortality rate with every additional year being sick

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
c_vec <- c(c.H, c.S1, c.S2, 0, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.D", "c.Trt")
c_T_vec <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, 0)
names(c_T_vec) <- c("c.H", "c.S1", "c.S2", "c.D")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, 0, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.D", "u.Trt")
u_T_vec <- c(u.H, u.Trt, u.S2, 0)
names(u_T_vec) <- c("u.H", "u.S1", "u.S2", "u.D")

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

## Compile C++ function:----
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_1.0.cpp"
  )
)
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_2.0.cpp"
  )
)
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_3.0.cpp"
  )
)
## Test function:----
results_1 <- MicroSimV_Cpp_1(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  b_Trt = FALSE,
  n_Seed = 1
)
results_1_T <- MicroSimV_Cpp_1(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  b_Trt = TRUE,
  n_Seed = 1
)
results_2 <- MicroSimV_Cpp_2(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_2_T <- MicroSimV_Cpp_2(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_T <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_2 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_T_2 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_T <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_3_LYs <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0,
  d_dE = 0,
  n_Seed = 1
)
results_3_4_dLYs <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "MicroSim_Cpp_1" = MicroSimV_Cpp_1( # 908ms median 955ms mean
    v_S_t = v_M_1,
    t_P = t_p,
    v_C = c_vec,
    v_U = u_vec,
    n_I = n.i,
    n_S = n.s,
    n_T = n.t,
    n_Cl = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    b_Trt = FALSE,
    n_Seed = 1
  ),
  "MicroSim_Cpp_2" = MicroSimV_Cpp_2( # 542ms median 577ms mean
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  ),
  "MicroSimV_Cpp_3_1" = MicroSimV_Cpp_3(
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    v_tracked_states = vector(),
    v_states_from = vector(),
    v_states_to = vector(),
    v_states_comp = vector(),
    v_increase_rate = vector(),
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  ),
  "MicroSimV_Cpp_3_2" = MicroSimV_Cpp_3(
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    v_tracked_states = c(2, 3),
    v_states_from = c(2, 3),
    v_states_to = c(4, 4),
    v_states_comp = c(2, 3),
    v_increase_rate = c(rp.S1S2, rp.S1S2),
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  )
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)