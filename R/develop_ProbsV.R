# Parameters:----
n.i   <- 1e6                # number of simulated individuals
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
Probs <- function(M_it) { 
  # M_it:    health state occupied by individual i at cycle t (character variable)
  
  m.p.it <- matrix(NA, n.s, n.i)     # create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # assign names to the vector
  
  # update the v.p with the appropriate probabilities   
  m.p.it[,M_it == "H"]  <- c(1 - p.HS1 - p.HD, p.HS1, 0, p.HD)                  # transition probabilities when healthy
  m.p.it[,M_it == "S1"] <- c(p.S1H, 1- p.S1H - p.S1S2 - p.S1D, p.S1S2, p.S1D)   # transition probabilities when sick
  m.p.it[,M_it == "S2"] <- c(0, 0, 1 - p.S2D, p.S2D)                            # transition probabilities when sicker
  m.p.it[,M_it == "D"]  <- c(0, 0, 0, 1)                                        # transition probabilities when dead   
  ifelse(colSums(m.p.it) == 1, return(t(m.p.it)), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
}       

# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
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
  includes = includes#,
  #plugins = plugins
)

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
## Code 14 - ProbsV:----
code14 <- 
  'arma::mat ProbsV_Cpp9( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P,
                      arma::mat& m_t_states, 
                      std::vector<int>& v_states_from, 
                      std::vector<int>& v_states_to,
                      std::vector<int>& v_states_comp,
                      std::vector<double>& v_increase_rate) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  // m_t_states: matrix containing time spent in tracked states.
  // v_states_from: vector containing the states from which transitions occur
  // v_states_to: vector containing the states to which transitions occur
  // v_states_comp: vector containing the states complementing the transition probabilities
  // v_increase_rate: change in transition rates with every additional cycle in tracked states

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
    
    // Get the health state to which the transition occur, re-scale to C++:
    int state_from = i + 1 ;

    // Use std::find to check if i is in v_states_from
    auto it = std::find(v_states_from.begin(), v_states_from.end(), state_from) ;
    if (it != v_states_from.end()) {
      // Found state_from in v_states_from, get the index:
      auto index = std::distance(v_states_from.begin(), it) ;

      // Check if state_from appears more than once in the v_states_from:
      for(size_t y = index; y <  v_states_from.size(); y++) {
        if(state_from == v_states_from[y]) {
        
          // Get the health state to which the transition occur, re-scale to C++:
          int state_to = v_states_to[y] - 1 ;

          // Convert the probabilities in the target column to rates:
          arma::vec v_R_t = - arma::log(1 - c_P_t.slice(i).col(state_to)) ;

          // Get cycle change in rates:
          double delta_rate = v_increase_rate[y] ;

          // Get accumulated change over time in i state:
          arma::vec v_dR_t = 1 + m_t_states.col(i) * delta_rate ;

          // Convert the rates back to probabilities after adjusting them:
          arma::colvec v_P_t = 1 - arma::exp(- v_R_t % v_dR_t) ;

          c_P_t.slice(i).col(state_to) =  v_P_t ; 
        }
      }
      
      // Calculate complement and update the target column
      int target = v_states_comp[index] - 1 ;
      c_P_t.slice(i).col(target) = 1 - (arma::sum(c_P_t.slice(i), 1) - c_P_t.slice(i).col(target)) ;
    }
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
## Code 15 - ProbsV (ProbsV_Cpp9 with printing):----
code15 <- 
  'arma::mat ProbsV_Cpp10( arma::rowvec& v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P,
                      arma::mat& m_t_states, 
                      std::vector<int>& v_states_from, 
                      std::vector<int>& v_states_to,
                      std::vector<int>& v_states_comp,
                      std::vector<double>& v_increase_rate) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  // m_t_states: matrix containing time spent in tracked states.
  // v_states_from: vector containing the states from which transitions occur
  // v_states_to: vector containing the states to which transitions occur
  // v_states_comp: vector containing the states complementing the transition probabilities
  // v_increase_rate: change in transition rates with every additional cycle in tracked states

  Function printR("print") ;
  
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
    
    // Get the health state to which the transition occur, re-scale to C++:
    int state_from = i + 1 ;
    printR("from") ;
    printR(state_from) ;

    // Use std::find to check if i is in v_states_from
    auto it = std::find(v_states_from.begin(), v_states_from.end(), state_from) ;
    if (it != v_states_from.end()) {
      // Found state_from in v_states_from, get the index:
      auto index = std::distance(v_states_from.begin(), it) ;
      printR("slice_before") ;
      printR(c_P_t.slice(i)) ;
      
      // Check if state_from appears more than once in the v_states_from:
      for(size_t y = index; y <  v_states_from.size(); y++) {
        if(state_from == v_states_from[y]) {
        
          // Get the health state to which the transition occur, re-scale to C++:
          int state_to = v_states_to[y] - 1 ;
          printR("to index") ;
          printR(state_to) ;
          
          // Convert the probabilities in the target column to rates:
          printR("transition matrix probabilities:") ;
          printR(c_P_t.slice(i).col(state_to)) ;
          arma::vec v_R_t = - arma::log(1 - c_P_t.slice(i).col(state_to)) ;
          printR("rates") ;
          printR(v_R_t) ;
          
          // Get cycle change in rates:
          double delta_rate = v_increase_rate[y] ;
          printR("cycle change in rates") ;
          printR(delta_rate) ;
          
          // Get accumulated change over time in i state:
          printR("time in state") ;
          printR(m_t_states.col(i)) ;
          arma::vec v_dR_t = 1 + m_t_states.col(i) * delta_rate ;
          printR("cumulative change in rate") ;
          printR(v_dR_t) ;
          
          // Convert the rates back to probabilities after adjusting them:
          arma::colvec v_P_t = 1 - arma::exp(- v_R_t % v_dR_t) ;
          printR("time adjusted probabilities") ;
          printR(v_P_t) ;
    
          c_P_t.slice(i).col(state_to) =  v_P_t ; 
          printR("slice") ;
          printR(c_P_t.slice(i)) ;
        
        }
      }
      
      // Calculate complement and update the target column
      int target = v_states_comp[index] - 1 ;
      printR("target") ;
      printR(target) ;
      c_P_t.slice(i).col(target) = 1 - (arma::sum(c_P_t.slice(i), 1) - c_P_t.slice(i).col(target)) ;
      
      printR("slice corrected") ;
      printR(c_P_t.slice(i)) ;
    }
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
  code12, code13, code14
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
testR <- Probs(
  M_it = v.M_1
)
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
## Advanced ProbsV:
m_t_states <- matrix(
  data = 0,
  nrow = n.i,
  ncol = n.s
)
test9_1 <- ProbsV_Cpp9(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p,
  m_t_states = m_t_states,
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2)
)
test9_2 <- ProbsV_Cpp9(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p,
  m_t_states = m_t_states,
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector()
)
v_M_1_2 <- sample(
  x = 1:n.s, 
  prob = rep(0.25, 4), 
  replace = TRUE, 
  size = n.i
)
m_t_states2 <- matrix(
  data = 0,
  nrow = n.i,
  ncol = n.s
)
time_in_state(
  m_t_states = m_t_states2,
  v_S_t = v_M_1_2,
  v_tracked_states = c(2, 3)
)
test9_3 <- ProbsV_Cpp9(
  v_S_t = v_M_1_2,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p,
  m_t_states = m_t_states2,
  v_states_from = c(2),
  v_states_to = c(4),
  v_states_comp = c(1),
  v_increase_rate = c(rp.S1S2)
)
test9_4 <- ProbsV_Cpp9(
  v_S_t = v_M_1_2,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p,
  m_t_states = m_t_states2,
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2)
)
test9_5 <- ProbsV_Cpp9(
  v_S_t = v_M_1_2,
  n_I = n.i,
  n_S = n.s,
  t_P = m_t_p,
  m_t_states = m_t_states2,
  v_states_from = c(2, 2, 3),
  v_states_to = c(3, 4, 4),
  v_states_comp = c(2, 2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2, rp.S1S2)
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "ProbsV_R" = Probs(
    M_it = v.M_1
  ),
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
  ),
  "ProbsV_Cpp9_1" = ProbsV_Cpp9(
    v_S_t = v_M_1_2,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p,
    m_t_states = m_t_states2,
    v_states_from = c(2),
    v_states_to = c(4),
    v_states_comp = c(1),
    v_increase_rate = c(rp.S1S2)
  ),
  "ProbsV_Cpp9_2" = ProbsV_Cpp9(
    v_S_t = v_M_1_2,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p,
    m_t_states = m_t_states2,
    v_states_from = c(2, 3),
    v_states_to = c(4, 4),
    v_states_comp = c(2, 3),
    v_increase_rate = c(rp.S1S2, rp.S1S2)
  ),
  "ProbsV_Cpp9_3" = ProbsV_Cpp9(
    v_S_t = v_M_1_2,
    n_I = n.i,
    n_S = n.s,
    t_P = m_t_p,
    m_t_states = m_t_states2,
    v_states_from = c(2, 2, 3),
    v_states_to = c(3, 4, 4),
    v_states_comp = c(2, 2, 3),
    v_increase_rate = c(rp.S1S2, rp.S1S2, rp.S1S2)
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)