# Parameters:----
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

# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# Using arma::cube:----

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
