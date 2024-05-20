# Parameters:----
n.i   <- 1e6                   # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names
seed_no <- 1234

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

# R function:----
samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}
# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
code_ProbsV <- 
  'arma::mat ProbsV_Cpp( arma::rowvec& v_S_t,
                       int& n_I,
                       int& n_S,
                       arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  
  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
    
    // Get the transition probabilities from health state i as a row vector:
    // Replicate t_P.row(i) n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P.row(i), n_I, 1) ;
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Identify individuals occupying state i at time t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    arma::uvec v_O_t = arma::find(v_S_t == (i + 1)) ;
    
    // Grab transition probabilities from state i
    // Assign probabilities for state i occupancy according to v_O_t:
    m_P_t.rows(v_O_t) = c_P_t.slice(i).rows(v_O_t) ;
  }
  
  // check if vector of probabilities sum to 1
  // need to round up to the 1e-6, otherwise it breaks
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / 
  std::pow(10, n) ;
  bool notSumToOne = arma::any(sum_row != 1.000000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else {
    return(m_P_t) ;
  }
}'
Rcpp::cppFunction(
  code = code_ProbsV,
  depends = depends,
  includes = includes
)
m_P_t <- ProbsV_Cpp(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
## Code 1 - Sick-Sicker SampleV:----
code1 <- 
  'arma::mat SampleV_Cpp( arma::mat m_P_t,
                       int& n_I,
                       int& n_S,
                       int m = 1) {
  // m_P_t: numeric matrix containing the transition probabilities for individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // m: number of health states to sample.
  
  // declare some assistive functions from R:
  Function printR("print") ;
  Function roundR("round") ;
  
  // create a matrix for sampled health states (m_M_t):
  arma::mat m_M_t(n_I, m, arma::fill::ones) ;
  
  // create a matrix m_CumProb_t:
  arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
  
  // recheck the probabilities:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec v_CumProb_t = as< Rcpp::NumericVector >(Rcpp::wrap(
    roundR(m_CumProb_t.col(m_CumProb_t.n_cols - 1), 6))) ;
  // the cumulative probability in the last column is expected to be equal to 1
  if(any(v_CumProb_t != 1.0000)) {
    stop("error in multinom: probabilities do not sum to 1") ;
  }
  
  // sample from a uniform U~(0, 1) distribution:
  arma::mat m_U(n_I, n_S, fill::ones) ; // a matrix to save values sampled from the U~(0, 1)
  for(int i = 0; i < m; i++) {
    // in each row, sample one random value for n_I individuals and repeat that value n_S times
    m_U = arma::repmat( randu<colvec>(n_I), 1, n_S ) ;
    // identify the first individual-specific health-state with a cumulative probability higher than the their corresponding sampled value
    // using a logical (true/false or 1/0), matrix get the value to be added to 1 (the starting health-state)
    // one plus the sum of each row of the logical matrix gives the health-state for the corresponding individuals at time t + 1
    m_M_t.col(i) = m_M_t.col(i) + arma::sum( (m_U > m_CumProb_t), 1 ) ;
  }
  
  return(m_M_t) ;
}'
## Code 2 - SampleV:----
code2 <- 
  'arma::mat SampleV_Cpp2( arma::mat& m_P_t,
                       int& n_I,
                       int& n_S,
                       int m = 1) {
  // m_P_t: numeric matrix containing the transition probabilities for individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // m: number of health states to sample.
  
  // create m_CumProb_t matrix with row-wise cumlative transition probabilities:
  arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
  
  // create a matrix for sampled health states (m_M_t):
  arma::mat m_M_t(n_I, m, arma::fill::ones) ;
  
  // recheck the probabilities:
  // need to round up to the 1e-6, otherwise it breaks
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec v_CumProb_t = arma::round(m_CumProb_t.col(m_CumProb_t.n_cols - 1) * 
    std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = arma::any(v_CumProb_t != 1.000000) ;
  
  if(notSumToOne) {
    stop("error in multinom: probabilities do not sum to 1") ;
  }
  
  // Initialise a matrix to save values sampled from the U~(0, 1)  
  arma::mat m_U(n_I, n_S, fill::ones) ;
  for(int i = 0; i < m; i++) {
    // in each row, sample one random value for n_I individuals and repeat that value n_S times
    m_U = arma::repmat( arma::randu<colvec>(n_I), 1, n_S ) ;
    
    // identify the first individual-specific health-state with a cumulative probability higher than the their corresponding sampled value
    // using a logical (true/false or 1/0), matrix get the value to be added to 1 (the starting health-state)
    // one plus the sum of each row of the logical matrix gives the health-state for the corresponding individuals at time t + 1
    m_M_t.col(i) = m_M_t.col(i) + arma::sum( (m_U > m_CumProb_t), 1 ) ;
  }
  
  return(m_M_t) ;
}'
## Code 3 - SampleV:----
code3 <- 
  'arma::colvec SampleV_Cpp3( arma::mat& m_P_t,
                       int& n_I,
                       int& n_S) {
  // m_P_t: numeric matrix containing the transition probabilities for individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  
  // create m_CumProb_t matrix with row-wise cumlative transition probabilities:
  arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
  
  // create a column vector for sampled health states (v_s_t):
  arma::colvec v_s_t = arma::ones<arma::colvec>(n_I) ;
  
  // recheck the probabilities:
  // need to round up to the 1e-6, otherwise it breaks
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec v_CumProb_t = arma::round(m_CumProb_t.col(m_CumProb_t.n_cols - 1) * 
    std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = arma::any(v_CumProb_t != 1.000000) ;
  
  if(notSumToOne) {
    stop("error in multinom: probabilities do not sum to 1") ;
  }
  
  // Initialise a matrix to save values sampled from the U~(0, 1)  
  arma::mat m_U(n_I, n_S, fill::ones) ;
  // in each row, sample one random value for n_I individuals and repeat that value n_S times
  m_U = arma::repmat( arma::randu<colvec>(n_I), 1, n_S ) ;
    
  // identify the first individual-specific health-state with a cumulative probability higher than the their corresponding sampled value
  // using a logical (true/false or 1/0), matrix get the value to be added to 1 (the starting health-state)
  // one plus the sum of each row of the logical matrix gives the health-state for the corresponding individuals at time t + 1
  v_s_t = v_s_t + arma::sum( (m_U > m_CumProb_t), 1 ) ;
  
  return(v_s_t) ;
}'

# Compile C++ code:----
cpp_functions_defs <- list(
  code1, code2, code3
)

for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    depends = depends,
    includes = includes#,
    #plugins = plugins
  )
}

# Test Rcpp functions:----
set.seed(seed_no)
test <- samplev(
  m_P_t,
  1
)
set.seed(seed_no)
test1 <- SampleV_Cpp(
  m_P_t = m_P_t,
  n_I = n.i,
  n_S = n.s,
  m = 1
)
set.seed(seed_no)
test2 <- SampleV_Cpp2(
  m_P_t = m_P_t,
  n_I = n.i,
  n_S = n.s,
  m = 1
)
set.seed(seed_no)
test3 <- SampleV_Cpp3(
  m_P_t = m_P_t,
  n_I = n.i,
  n_S = n.s
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 100,
  "SampleV_R" = samplev(
    m_P_t,
    1
  ),
  "SampleV_Cpp" = SampleV_Cpp(
    m_P_t = m_P_t,
    n_I = n.i,
    n_S = n.s,
    m = 1
  ),
  "SampleV_Cpp2" = SampleV_Cpp2(
    m_P_t = m_P_t,
    n_I = n.i,
    n_S = n.s,
    m = 1
  ),
  "SampleV_Cpp3" = SampleV_Cpp3(
    m_P_t = m_P_t,
    n_I = n.i,
    n_S = n.s
  )
)

# Print the results
print(results)

# Unit: milliseconds
#         expr      min        lq      mean    median        uq      max neval
#    SampleV_R 135.3824 152.93265 161.83898 158.11995 165.97400 240.1300   100
#  SampleV_Cpp  94.6879  98.69145 105.03883 100.61445 106.54750 162.2282   100
# SampleV_Cpp2  53.0879  55.41015  58.92750  56.98580  58.93545 113.3415   100
# SampleV_Cpp3  51.4090  53.49470  56.64896  54.29895  56.71630 110.1587   100

# For a visual comparison
boxplot(results)