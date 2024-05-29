# clear R working global environment
rm(list = ls())

# source calc_discount_wtsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "calc_discount_wtsC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the calc_discount_wts functions:
calc_discount_wts <- function(
    discount_rate,
    num_cycles,
    cycle_length) {
  
  # calculate discount weights based on the number & length (in years) of cycles
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  
  return(v_discount_wts)
}

#------------------------------------------------------------------------------#

# define function inputs:
## General parameters
seed                <- 1234    # random number generator state
num_cycles          <- 30      # time horizon
num_i               <- 1e6     # number of simulated individuals
discount_rate_costs <- 0.03    # annual discount rate for costs
discount_rate_QALYs <- 0.015   # annual discount rate for health outcomes
cycle_length        <- 1       # length of cycle in years

#------------------------------------------------------------------------------#

# run the calc_discount_wts function:
R_results1 <- calc_discount_wts(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)
R_results2 <- calc_discount_wts(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)
# run the calc_discount_wtsC function:
C_results1 <- calc_discount_wtsC1(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)
C_results2 <- calc_discount_wtsC1(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)
C_results3 <- calc_discount_wtsC2(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)[, 1]
C_results4 <- calc_discount_wtsC2(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)[, 1]
C_results5 <- calc_discount_wtsC3(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)[, 1]
C_results6 <- calc_discount_wtsC3(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)[, 1]
C_results7 <- calc_discount_wtsC4(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)
C_results8 <- calc_discount_wtsC4(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)
C_results9 <- calc_discount_wtsC5(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = cycle_length
)
C_results10 <- calc_discount_wtsC5(
  discount_rate = discount_rate_costs,
  num_cycles    = num_cycles,
  cycle_length  = 0.5
)
# check results
identical(R_results1, C_results1)
identical(R_results2, C_results2)
identical(R_results1, C_results3)
identical(R_results2, C_results4)
identical(R_results1, C_results5)
identical(R_results2, C_results6)
identical(R_results1, C_results7)
identical(R_results2, C_results8)
identical(R_results1, C_results9)
identical(R_results2, C_results10)

#------------------------------------------------------------------------------#

# benchmark the functions
calc_discount_wts_RvC <- bench::mark(
  "R_1" = calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_1" = calc_discount_wtsC1(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_2" = calc_discount_wtsC2(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_3" = calc_discount_wtsC3(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_4" = calc_discount_wtsC4(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_5" = calc_discount_wtsC5(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  check = FALSE
)

calc_discount_wts_RvC[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

# Use the microbenchmark::microbenchmark()
calc_discount_wts_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_1" = calc_discount_wtsC1(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_2" = calc_discount_wtsC2(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_3" = calc_discount_wtsC3(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_4" = calc_discount_wtsC4(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  ),
  "C_5" = calc_discount_wtsC5(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
)

calc_discount_wts_RvC2
plot(calc_discount_wts_RvC2)