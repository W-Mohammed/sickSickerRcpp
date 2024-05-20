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

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, 0, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.D", "c.Trt")

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
Costs <- function (M_it, Trt = FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual being treated? (default is FALSE) 
  
  c.it <- 0                                  # by default the cost for everyone is zero 
  c.it[M_it == "H"]  <- c.H                  # update the cost if healthy
  c.it[M_it == "S1"] <- c.S1 + c.Trt * Trt   # update the cost if sick conditional on treatment
  c.it[M_it == "S2"] <- c.S2 + c.Trt * Trt   # update the cost if sicker conditional on treatment
  c.it[M_it == "D"]  <- 0                    # update the cost if dead
  
  return(c.it)        		                   # return the costs
}
# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
## Code 1 - Sick-Sicker CostsV:----
code1 <- 
  'arma::colvec CostsV_Cpp( arma::colvec v_S_t,
                         int& n_I,
                         int& n_S,
                         NumericVector& v_Costs,
                         bool b_Trt = false ) {
  // v_S_t: vector of health states occupied by individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  // v_Costs: a vector containing cost parameters.
  // b_Trt: a bool indicating whether treatment costs (default is false).
  
  // declaring costs for the individuals at time t and set it to 0:
  arma::mat m_C_t (n_I, 1, fill::zeros) ;
  
  // define state costs:
  arma::mat m_H  = arma::repmat( arma::rowvec {(v_Costs["c.H"])}, n_I, 1 ) ;
  arma::mat m_S1 = arma::repmat( arma::rowvec {(v_Costs["c.S1"] + (v_Costs["c.Trt"] * b_Trt))}, n_I, 1 ) ;
  arma::mat m_S2 = arma::repmat( arma::rowvec {(v_Costs["c.S2"] + (v_Costs["c.Trt"] * b_Trt))}, n_I, 1 ) ;
  
  // update m_C_t with the appropriate costs:
  m_C_t.rows ( arma::find( v_S_t == 1 ) ) = m_H.rows  ( arma::find( v_S_t == 1 ) ) ;
  m_C_t.rows ( arma::find( v_S_t == 2 ) ) = m_S1.rows ( arma::find( v_S_t == 2 ) ) ;
  m_C_t.rows ( arma::find( v_S_t == 3 ) ) = m_S2.rows ( arma::find( v_S_t == 3 ) ) ;
  
  return(m_C_t.col(0)) ;
}'
## Code 2 - CostsV:----
code2 <- 
  'arma::colvec CostsV_Cpp2( arma::colvec& v_S_t,
                         arma::vec& v_Costs) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Costs: a vector containing cost parameters.

  // Transforming state occupancy to allign with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;

  // Use states indecies to assign costs correctly:
  arma::colvec v_col_costs = v_Costs.elem(uv_indices) ;

  return(v_col_costs) ;
}'

# Compile C++ code:----
cpp_functions_defs <- list(
  code1, code2
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
testR <- Costs(
  M_it = v.M_1,
  Trt = FALSE
)
test1 <- CostsV_Cpp(
  v_S_t = v_M_1,
  n_I = n.i,
  n_S = n.s,
  v_Costs = c_vec,
  b_Trt = FALSE
)
test2 <- CostsV_Cpp2(
  v_S_t = v_M_1,
  v_Costs = c_vec
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "Costs_R" = Costs(
    M_it = v.M_1,
    Trt = FALSE
  ),
  "Costs_Cpp" = CostsV_Cpp(
    v_S_t = v_M_1,
    n_I = n.i,
    n_S = n.s,
    v_Costs = c_vec,
    b_Trt = FALSE
  ),
  "Costs_Cpp2" = CostsV_Cpp2(
    v_S_t = v_M_1,
    v_Costs = c_vec
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)