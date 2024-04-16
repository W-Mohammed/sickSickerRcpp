# Parameters:----
n.i   <- 100000                # number of simulated individuals
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

# R function:----
Effs <- function (M_it, Trt = FALSE, cl = 1) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[M_it == "H"]  <- u.H      # update the utility if healthy
  u.it[M_it == "S1"] <- Trt * u.Trt + (1 - Trt) * u.S1  # update the utility if sick conditional on treatment
  u.it[M_it == "S2"] <- u.S2     # update the utility if sicker
  QALYs <-  u.it * cl            # calculate the QALYs during cycle t
  return(QALYs)                  # return the QALYs
}
# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
## Code 1 - Sick-Sicker EffsV:----
code1 <- 
  'arma::colvec EffsV_Cpp( arma::colvec v_S_t,
                        int& n_I,
                        int& n_S,
                        NumericVector& v_Utilities,
                        bool b_Trt = false,
                        int cl = 1 ) {
  // v_S_t: vector of health states occupied by individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // v_Utilities: a vector containing utilities values for each health state
  // b_Trt: bool indicating whether treatment costs (default is false).
  // cl: integer variable indicating cycle length (cl) - default is 1
  
  // declaring effects matrix for the individuals at time t and set it to 0:
  arma::mat m_E_t (n_I, 1, fill::zeros) ;
  
  // define state QALYs:
  arma::mat m_H  = arma::repmat( arma::rowvec {(v_Utilities["u.H"] * cl)}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {((b_Trt * v_Utilities["u.Trt"] + (1 - b_Trt) * v_Utilities["u.S1"]) * cl)}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {(v_Utilities["u.S2"] * cl)}, n_I, 1 ) ;
  
  // update m_E_t with the appropriate health effects:
  m_E_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_E_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_E_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  
  return(m_E_t.col(0)) ;
}'
## Code 2 - EffsV:----
code2 <- 
  'arma::colvec EffsV_Cpp2( arma::colvec& v_S_t,
                         arma::vec& v_Utilities,
                         int cycle = 1 ) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Utilities: a vector containing utilities values for each health state.
  // cycle: integer variable indicating cycle length in years - default is 1.

  // Transforming state occupancy to allign with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  // Calculating QALYs, multiplying utilities by the length of the cycle length:
  v_Utilities = v_Utilities * cycle ;

  // Use states indecies to assign utilities correctly:
  arma::colvec v_col_QALYs = v_Utilities.elem(uv_indices) ;

  return(v_col_QALYs) ;
}'
## Code 3 - EffsV:----
code3 <- 
  'arma::colvec EffsV_Cpp3( arma::colvec& v_S_t,
                         arma::vec& v_Utilities,
                         int cycle = 1 ) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Utilities: a vector containing utilities values for each health state.
  // cycle: integer variable indicating cycle length in years - default is 1.

  // Transforming state occupancy to allign with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  // Calculating QALYs, multiplying utilities by the length of the cycle length:
  arma::vec v_QALYs = v_Utilities * cycle ;

  // Use states indecies to assign utilities correctly:
  arma::colvec v_col_QALYs = v_QALYs.elem(uv_indices) ;

  return(v_col_QALYs) ;
}'
## Code 4 - EffsV:----
code4 <- 
  'arma::colvec EffsV_Cpp4( arma::colvec& v_S_t,
                         arma::vec& v_Utilities,
                         arma::vec& v_d_Utilities,
                         arma::mat& m_t_states,
                         int cycle = 1 ) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Utilities: vector containing utilities values for each health state.
  // v_d_Utilities: vector containing dis-utilities for each health state.
  // m_t_states: matrix storing times spent in each health state.
  // cycle: integer variable indicating cycle length in years - default is 1.

  // Transforming state occupancy to allign with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  // Estimating disutilities based on the time spent in each state:
  v_d_Utilities = m_t_states * v_d_Utilities ;
  
  // Use states indecies to assign utilities correctly:
  arma::vec v_col_Utilities = v_Utilities.elem(uv_indices) - v_d_Utilities ;
  
  // Calculating QALYs, multiplying utilities by the length of the cycle length:
  arma::colvec v_col_QALYs = v_col_Utilities  * cycle ;
  
  return(v_col_QALYs) ;
}'
## Code 5 - EffsV:----
code5 <- 
  'arma::colvec EffsV_Cpp5(arma::colvec& v_S_t,
                        arma::vec& v_Utilities,
                        arma::vec* v_d_Utilities, // Use a pointer instead of a reference
                        arma::mat& m_t_states,
                        int cycle = 1) {
                        
  // Transforming state occupancy to align with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  arma::vec v_col_Utilities = v_Utilities.elem(uv_indices) ;

  // Check if v_d_Utilities is not NULL
  if (v_d_Utilities != nullptr) {
  
    // Estimating disutilities based on the time spent in each state if v_d_Utilities is provided:
    arma::vec v_col_d_Utilities = m_t_states * (*v_d_Utilities) ; // Dereference the pointer to use its value
    v_col_Utilities -= v_col_d_Utilities;
  }
  
  // Calculating effects, multiplying utilities by the length of the cycle:
  arma::colvec v_col_effs = v_col_Utilities * cycle ;
  
  return(v_col_effs) ;
}'
# EffsV_Cpp5(v_S_t, v_Utilities, nullptr, m_t_states, cycle);
# EffsV_Cpp5(v_S_t, v_Utilities, &v_d_Utilities, m_t_states, cycle);
## Code 6 - EffsV:----
code6 <- 
  'arma::colvec EffsV_Cpp6(arma::colvec& v_S_t,
                        arma::vec& v_Utilities,
                        Rcpp::Nullable<Rcpp::NumericVector> v_d_Utilities, // Changed to Nullable
                        arma::mat& m_t_states,
                        int cycle = 1) {
  // Transforming state occupancy to align with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  arma::vec v_col_Utilities = v_Utilities.elem(uv_indices) ;

  // Check if v_d_Utilities is not NULL
  if (v_d_Utilities.isNotNull()) {
  
    // Convert Rcpp::NumericVector to arma::vec
    arma::vec v_col_d_Utilities = arma::vec(Rcpp::as<arma::vec>(v_d_Utilities)) ;
    
    // Estimating disutilities based on the time spent in each state if v_d_Utilities is provided:
    v_col_d_Utilities = m_t_states * v_col_d_Utilities; // No need to dereference
    v_col_Utilities -= v_col_d_Utilities;
  }
  
  // Calculating effects, multiplying utilities by the length of the cycle:
  arma::colvec v_col_effs = v_col_Utilities * cycle;
  
  return v_col_effs;
}'
## code 7 - EffsV:----
code7 <- 
  'arma::colvec EffsV_Cpp7(arma::colvec& v_S_t,
                        arma::vec& v_Utilities,
                        Rcpp::Nullable<Rcpp::NumericVector> v_d_Utilities, // Remains Nullable
                        Rcpp::Nullable<Rcpp::NumericMatrix> m_t_states, // Now Nullable
                        int cycle = 1) {
  
  // Transforming state occupancy to align with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1);
  
  arma::vec v_col_Utilities = v_Utilities.elem(uv_indices);

  // Check if both v_d_Utilities and m_t_states are not NULL in the same IF statement
  if (v_d_Utilities.isNotNull() && m_t_states.isNotNull()) {
  
    // Convert Rcpp::NumericVector to arma::vec for disutilities
    arma::vec v_col_d_Utilities = Rcpp::as<arma::vec>(v_d_Utilities);
    
    // Convert Rcpp::NumericMatrix to arma::mat for m_t_states
    arma::mat m_t_states = Rcpp::as<arma::mat>(m_t_states);
    
    // Estimating disutilities based on the time spent in each state:
    v_col_d_Utilities = m_t_states * v_col_d_Utilities;
    v_col_Utilities -= v_col_d_Utilities;
  }
  
  // Calculating effects, multiplying utilities by the length of the cycle:
  arma::colvec v_col_effs = v_col_Utilities * cycle;
  
  return v_col_effs;
}'
# Compile C++ code:----
cpp_functions_defs <- list(
  code1, code2, code3, code4, code6, code7#, code5
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
v.M_1 <- v.n[v_M_1]
time_in_states <- cbind(v_M_1, 0, 0, 0)
time_in_states2 <- cbind(
  0, ifelse(
    v_M_1 == 2, sample(
      x = 1:50, 
      prob = rep(1/50, 50), 
      replace = TRUE, 
      size = n.i
  ), 0), 0, 0)
testR <- Effs(
  M_it = v.M_1,
  Trt = FALSE,
  cl = 1
)
testR2 <- u_vec[v_M_1]
test1 <- EffsV_Cpp(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  v_Utilities = u_vec,
  b_Trt = FALSE,
  cl = 1
)
test2 <- EffsV_Cpp2(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  cycle = 1
)
test3 <- EffsV_Cpp3(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  cycle = 1
)
test4 <- EffsV_Cpp4(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = ru_vec,
  m_t_states = time_in_states,
  cycle = 1
)
test5 <- EffsV_Cpp4(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = ru_vec,
  m_t_states = time_in_states2,
  cycle = 1
)
test6_1 <- EffsV_Cpp6(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = ru_vec,
  m_t_states = time_in_states2,
  cycle = 1
)
test6_2 <- EffsV_Cpp6(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = NULL,
  m_t_states = time_in_states2,
  cycle = 1
)
test7_1 <- EffsV_Cpp7(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = ru_vec,
  m_t_states = time_in_states2,
  cycle = 1
)
test7_2 <- EffsV_Cpp7(
  v_S_t = v_M_1,
  v_Utilities = u_vec,
  v_d_Utilities = NULL,
  m_t_states = NULL,
  cycle = 1
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "Effs_R" = Effs(
    M_it = v.M_1,
    Trt = FALSE,
    cl = 1
  ),
  "Effs_R2" = u_vec[v_M_1],
  "EffsV_Cpp" = EffsV_Cpp(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    v_Utilities = u_vec,
    b_Trt = FALSE,
    cl = 1
  ),
  "EffsV_Cpp2" = EffsV_Cpp2(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    cycle = 1
  ),
  "EffsV_Cpp3" = EffsV_Cpp3(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    cycle = 1
  ),
  "EffsV_Cpp4" = EffsV_Cpp4(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    v_d_Utilities = ru_vec,
    m_t_states = time_in_states2,
    cycle = 1
  ),
  "EffsV_Cpp6_1" = EffsV_Cpp6(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    v_d_Utilities = ru_vec,
    m_t_states = time_in_states2,
    cycle = 1
  ),
  "EffsV_Cpp6_2" = EffsV_Cpp6(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    v_d_Utilities = NULL,
    m_t_states = time_in_states2,
    cycle = 1
  ),
  "EffsV_Cpp7_1" = EffsV_Cpp7(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    v_d_Utilities = ru_vec,
    m_t_states = time_in_states2,
    cycle = 1
  ),
  "EffsV_Cpp7_2" = EffsV_Cpp7(
    v_S_t = v_M_1,
    v_Utilities = u_vec,
    v_d_Utilities = NULL,
    m_t_states = NULL,
    cycle = 1
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)