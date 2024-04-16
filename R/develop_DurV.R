# Parameters:----
n.i   <- 1e5                    # number of simulated individuals
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
ru.S1S2 <- 0.01               # decrease in utility of treated sick individuals with every additional year being sick/sicker

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, 0, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.D", "c.Trt")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, 0, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.D", "u.Trt")

# Create a vector containing dis-utilities parameters:
ru_vec <- c(0, ru.S1S2, ru.S1S2, 0)
names(ru_vec) <- c("ru.H", "ru.S1", "ru.S2", "ru.D")

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
## Code 1 - time_in_states:----
code1 <- 
  'arma::mat time_in_state( arma::mat& m_t_states,
                         arma::colvec& v_S_t,
                         std::vector<int>& v_tracked_states) {
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
    // Increment the time for all individuals in the column.
    // tmp_col += 1 ;
    arma::vec tmp_col = m_t_states.col(i) + 1 ;
    
    // Reset time in other states to 0 using .elem():
    // First, ensure other_states has valid indices within the range of tmp_col
    if (!other_states.empty()) {
      tmp_col.elem(other_states).fill(0) ;
    }
    
    // Record time back to states time matrix:
    m_t_states.col(i) = tmp_col ;
  }
  
  return(m_t_states) ;
}'
## Code 2 - time_in_states:----
code2 <- 
  'arma::mat time_in_state2( arma::mat& m_t_states,
                         arma::colvec& v_S_t,
                         std::vector<int>& v_tracked_states) {
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
    arma::vec tmp_col = m_t_states.col(i) ;
    
    // Increment the time for all individuals in the column.
    tmp_col += 1 ;
    
    // Ensure other_states has valid indices within the range of tmp_col
    if (!other_states.empty()) {
      // Reset time in other states to 0 using .elem():
      tmp_col.elem(other_states).fill(0) ;
    }
    
    // Record time back to states time matrix:
    m_t_states.col(i) = tmp_col ;
  }
  
  return(m_t_states) ;
}'
## Code 3 - time_in_states:----
code3 <- 
  'void time_in_state3( arma::mat& m_t_states,
                         arma::colvec& v_S_t,
                         std::vector<int>& v_tracked_states) {
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
    arma::vec tmp_col = m_t_states.col(i) ;
    
    // Increment the time for all individuals in the column.
    tmp_col += 1 ;
    
    // Ensure other_states has valid indices within the range of tmp_col
    if (!other_states.empty()) {
      // Reset time in other states to 0 using .elem():
      tmp_col.elem(other_states).fill(0) ;
    }
    
    // Record time back to states time matrix:
    m_t_states.col(i) = tmp_col ;
  }
}'

## Code 4 - time_in_states:----
code4 <- 
  'void time_in_state4( arma::mat& m_t_states,
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
## Code 5 - time_in_states:----
### Returns the matrix manipulated by the function
code5 <- 
  'arma::mat time_in_state5( arma::mat& m_t_states,
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
  
  return(m_t_states);
}'

## Code 6 - time_in_states:----
### Returns the matrix passed to the function
code6 <- 
  'arma::mat time_in_state6( arma::mat m_t_states,
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
  
  return(m_t_states);
}'
# Compile C++ code:----
cpp_functions_defs <- list(
  code1, code2, code3, code4, code5, code6
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
v_M_1 <- sample(
  x = 1:n.s, 
  prob = rep(0.25, 4), 
  replace = TRUE, 
  size = n.i
)
m_t_states <- matrix(
  data = 0,
  nrow = n.i,
  ncol = n.s
)
m_t_states2 <- matrix(
    data = 0,
    nrow = n.i,
    ncol = n.s
  )
m_t_states3 <- matrix(
    data = 0,
    nrow = n.i,
    ncol = n.s
  )
m_t_states4 <- matrix(
    data = 0,
    nrow = n.i,
    ncol = n.s
  )
m_t_states5 <- matrix(
    data = 0,
    nrow = n.i,
    ncol = n.s
  )
m_t_states6 <- matrix(
    data = 0,
    nrow = n.i,
    ncol = n.s
  )
test1 <- time_in_state(
  m_t_states = m_t_states,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)
test2 <- time_in_state2(
  m_t_states = m_t_states2,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)
time_in_state3(
  m_t_states = m_t_states3,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)
time_in_state4(
  m_t_states = m_t_states4,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)
time_in_state5(
  m_t_states = m_t_states5,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)
test6 <- time_in_state6(
  m_t_states = m_t_states6,
  v_S_t = v_M_1,
  v_tracked_states = c(2, 3)
)

# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "time_in_state" = time_in_state(
    m_t_states = m_t_states,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  ),
  "time_in_state2" = time_in_state2(
    m_t_states = m_t_states2,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  ),
  "time_in_state3" = time_in_state3(
    m_t_states = m_t_states3,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  ),
  "time_in_state4" = time_in_state4( # fastest
    m_t_states = m_t_states4,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  ),
  "time_in_state5" = time_in_state5(
    m_t_states = m_t_states5,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  ),
  "time_in_state6" = time_in_state6(
    m_t_states = m_t_states6,
    v_S_t = v_M_1,
    v_tracked_states = c(2, 3)
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)