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
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
## Code 1 - Sick-Sicker ProbsV:----
code1 <- 
  'arma::mat ProbsV_Cpp( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      NumericVector& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: vector containing transition probabilities.
  
  // declare other R functions:
  Function cat("cat") ;
  Function printR("print") ;
  Function roundR("round") ;
  
  // create a matrix for the state transition probabilities (m_P_t):
  arma::mat m_P_t(n_I, n_S) ;
  
  // define probabilities:
  arma::mat m_H  = arma::repmat( arma::rowvec {1 - t_P["p.HS1"] - t_P["p.HD"], t_P["p.HS1"], 0, t_P["p.HD"]}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {t_P["p.S1H"], 1 - t_P["p.S1H"] - t_P["p.S1S2"] - t_P["p.S1D"], t_P["p.S1S2"], t_P["p.S1D"]}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {0, 0, 1 - t_P["p.S2D"], t_P["p.S2D"]}, n_I, 1 ) ;
  arma::mat m_D  = arma::repmat( arma::rowvec {0, 0, 0, 1}, n_I, 1 ) ;
  
  // update v.p.it with the appropriate probabilities:
  m_P_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 4 ) ) = m_D.rows  ( arma::find( v_S_t == 4 ) ) ;
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 2 - ProbsV:----
code2 <- 
  'arma::mat ProbsV_Cpp1_1( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      NumericVector& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: vector containing transition probabilities.
  
  // declare other R functions:
  Function cat("cat") ;
  Function printR("print") ;
  Function roundR("round") ;
  
  // create a matrix for the state transition probabilities (m_P_t):
  arma::mat m_P_t(n_I, n_S) ;
  
  // define probabilities:
  arma::mat m_H  = arma::repmat( arma::rowvec {1 - t_P["p.HS1"] - t_P["p.HD"], t_P["p.HS1"], 0, t_P["p.HD"]}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {t_P["p.S1H"], 1 - t_P["p.S1H"] - t_P["p.S1S2"] - t_P["p.S1D"], t_P["p.S1S2"], t_P["p.S1D"]}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {0, 0, 1 - t_P["p.S2D"], t_P["p.S2D"]}, n_I, 1 ) ;
  arma::mat m_D  = arma::repmat( arma::rowvec {0, 0, 0, 1}, n_I, 1 ) ;
  
  // update v.p.it with the appropriate probabilities:
  m_P_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 4 ) ) = m_D.rows  ( arma::find( v_S_t == 4 ) ) ;
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 3 - ProbsV:----
code3 <- 
  'arma::mat ProbsV_Cpp1_2( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      NumericVector& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: vector containing transition probabilities.
  
  // declare other R functions:
  Function roundR("round") ;

  // create a matrix for the state transition probabilities (m_P_t):
  arma::mat m_P_t(n_I, n_S) ;
  
  // define probabilities:
  arma::mat m_H  = arma::repmat( arma::rowvec {1 - t_P["p.HS1"] - t_P["p.HD"], t_P["p.HS1"], 0, t_P["p.HD"]}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {t_P["p.S1H"], 1 - t_P["p.S1H"] - t_P["p.S1S2"] - t_P["p.S1D"], t_P["p.S1S2"], t_P["p.S1D"]}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {0, 0, 1 - t_P["p.S2D"], t_P["p.S2D"]}, n_I, 1 ) ;
  arma::mat m_D  = arma::repmat( arma::rowvec {0, 0, 0, 1}, n_I, 1 ) ;
  
  // update v.p.it with the appropriate probabilities:
  m_P_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  m_P_t.rows ( arma::find( v_S_t == 4 ) ) = m_D.rows  ( arma::find( v_S_t == 4 ) ) ;
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 4 - ProbsV:----
code4 <- 
  'arma::mat ProbsV_Cpp2( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  
  // declare other R functions:
  Function roundR("round") ;

  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
  
    // Extract the transition probabilities from health state i as a row vector:
    arma::rowvec t_P_i = t_P.row(i) ;
    
    // Replicate t_P_i n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P_i, n_I, 1) ;
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Grab transition probabilities from state i
    arma::mat t_P_i = c_P_t.slice(i) ;
    
    // Assign probabilities by current states occupancy as in v_S_t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = t_P_i.rows(arma::find(v_S_t == (i + 1))) ;

  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'

## Code 5 - ProbsV:----
code5 <- 
  'arma::mat ProbsV_Cpp2_1( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.

  // declare other R functions:
  Function roundR("round") ;

  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
  
    // Extract the transition probabilities from health state i as a row vector:
    arma::rowvec t_P_i = t_P.row(i) ;
    
    // Replicate t_P_i n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P_i, n_I, 1) ;
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Grab transition probabilities from state i
    arma::mat t_P_i = c_P_t.slice(i) ;
    
    // Assign probabilities by current states occupancy as in v_S_t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = t_P_i.rows(arma::find(v_S_t == (i + 1))) ;

  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 6 - ProbsV:----
code6 <- 
  'arma::mat ProbsV_Cpp3( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  
  // declare other R functions:
  Function roundR("round") ;

  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
  
    // Extract the transition probabilities from health state i as a row vector:
    arma::rowvec t_P_i = t_P.row(i) ;
    
    // Replicate t_P_i n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P_i, n_I, 1) ;
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Grab transition probabilities from state i
    arma::mat t_P_i = c_P_t.slice(i) ;
    
    // Identify individuals occupying state i at time t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    arma::uvec v_O_t = arma::find(v_S_t == (i + 1)) ;
    
    // Assign probabilities for state i occupancy according to v_O_t:
    m_P_t.rows(v_O_t) = t_P_i.rows(v_O_t) ;

  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 7 - ProbsV:----
code7 <- 
  'arma::mat ProbsV_Cpp3_1( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.

  // declare other R functions:
  Function roundR("round") ;

  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
  
    // Extract the transition probabilities from health state i as a row vector:
    arma::rowvec t_P_i = t_P.row(i) ;
    
    // Replicate t_P_i n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P_i, n_I, 1) ;
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Grab transition probabilities from state i
    arma::mat t_P_i = c_P_t.slice(i) ;
    
    // Identify individuals occupying state i at time t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    arma::uvec v_O_t = arma::find(v_S_t == (i + 1)) ;
    
    // Assign probabilities for state i occupancy according to v_O_t:
    m_P_t.rows(v_O_t) = t_P_i.rows(v_O_t) ;

  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 8 - ProbsV:----
code8 <- 
  'arma::mat ProbsV_Cpp4( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.

  // declare other R functions:
  Function roundR("round") ;
  
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
    // Grab transition probabilities from state i c_P_t.slice(i):
    // Identify individuals occupying state i at time t arma::find(v_S_t == (i + 1)):
    // Adjust states occupancy to allign with C++ indexing, which starts at 0 v_S_t == (i + 1):
    // Assign probabilities for state i occupancy:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = c_P_t.slice(i).rows(arma::find(v_S_t == (i + 1))) ;

  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 9 - ProbsV:----
code9 <- 
  'arma::mat ProbsV_Cpp4_1( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.

  // declare other R functions:
  Function roundR("round") ;
  
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
    // Grab transition probabilities from state i c_P_t.slice(i):
    // Identify individuals occupying state i at time t arma::find(v_S_t == (i + 1)):
    // Adjust states occupancy to allign with C++ indexing, which starts at 0 v_S_t == (i + 1):
    // Assign probabilities for state i occupancy:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = c_P_t.slice(i).rows(arma::find(v_S_t == (i + 1))) ;
  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  colvec sums = as< Rcpp::NumericVector >(Rcpp::wrap(roundR(arma::sum(m_P_t, 1), 6))) ;
  bool notSumToOne = any(sums != 1.0000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else
    return(m_P_t) ;
}'
## Code 10 - ProbsV:----
code10 <- 
  'arma::mat ProbsV_Cpp5( arma::rowvec& v_S_t,
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
    // Grab transition probabilities from state i c_P_t.slice(i):
    // Identify individuals occupying state i at time t arma::find(v_S_t == (i + 1)):
    // Adjust states occupancy to allign with C++ indexing, which starts at 0 v_S_t == (i + 1):
    // Assign probabilities for state i occupancy:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = c_P_t.slice(i).rows(arma::find(v_S_t == (i + 1))) ;
  }
  
  return(m_P_t) ;
}'
## Code 11 - ProbsV:----
code11 <- 
  'arma::mat ProbsV_Cpp6( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P ) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  
  // declare other R functions:
  Function printR("print") ;
  Function roundR("round") ;
  
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
    // Grab transition probabilities from state i c_P_t.slice(i):
    // Identify individuals occupying state i at time t arma::find(v_S_t == (i + 1)):
    // Adjust states occupancy to allign with C++ indexing, which starts at 0 v_S_t == (i + 1):
    // Assign probabilities for state i occupancy:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = c_P_t.slice(i).rows(arma::find(v_S_t == (i + 1))) ;
  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  int n = 6; // Rounding to the 6th decimal place
  arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = any(sum_row != 1.000000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else {
    return(m_P_t) ;
  }
}'

## Code 12 - ProbsV:----
code12 <- 
  'arma::mat ProbsV_Cpp7( arma::rowvec& v_S_t,
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
    // Grab transition probabilities from state i c_P_t.slice(i):
    // Identify individuals occupying state i at time t arma::find(v_S_t == (i + 1)):
    // Adjust states occupancy to allign with C++ indexing, which starts at 0 v_S_t == (i + 1):
    // Assign probabilities for state i occupancy:
    m_P_t.rows(arma::find(v_S_t == (i + 1))) = c_P_t.slice(i).rows(arma::find(v_S_t == (i + 1))) ;
  }
  
  // check if vector of probabilities sum to 1:
  // need to round up to the 1e-6, otherwise it breaks, e.g, once it was 1.00000001
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = any(sum_row != 1.000000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else {
    return(m_P_t) ;
  }
}'
## Code 13 - ProbsV:----
code13 <- 
  'arma::mat ProbsV_Cpp8( arma::rowvec& v_S_t,
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
  arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = arma::any(sum_row != 1.000000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else {
    return(m_P_t) ;
  }
}'

# Compile C++ code:----
cpp_functions_defs <- list(
  code1, code2, code3, code4, code5, code6, code7, code8, code9, code10, code11,
  code12, code13
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
test1 <- ProbsV_Cpp(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = t_p
)
test1_1 <- ProbsV_Cpp1_1(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = t_p
)
test1_2 <- ProbsV_Cpp1_2(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = t_p
)
test2 <- ProbsV_Cpp2(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test2_1 <- ProbsV_Cpp2_1(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test3 <- ProbsV_Cpp3(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test3_1 <- ProbsV_Cpp3_1(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test4 <- ProbsV_Cpp4(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test4_1 <- ProbsV_Cpp4_1(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test5 <- ProbsV_Cpp5(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test6 <- ProbsV_Cpp6(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test7 <- ProbsV_Cpp7(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
test8 <- ProbsV_Cpp8(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "ProbsV_Cpp" = ProbsV_Cpp(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = t_p
  ),
  "ProbsV_Cpp1_1" = ProbsV_Cpp1_1(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = t_p
  ),
  "ProbsV_Cpp1_2" = ProbsV_Cpp1_2(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = t_p
  ),
  "ProbsV_Cpp2" = ProbsV_Cpp2(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp2_1" = ProbsV_Cpp2_1(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp3" = ProbsV_Cpp3(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp3_1" = ProbsV_Cpp3_1(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp4" = ProbsV_Cpp4(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp4_1" = ProbsV_Cpp4_1(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp5" = ProbsV_Cpp5(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp6" = ProbsV_Cpp6(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp7" = ProbsV_Cpp7(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  ),
  "ProbsV_Cpp8" = ProbsV_Cpp8(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)