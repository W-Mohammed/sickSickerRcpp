# clear R working global environment
rm(list = ls())

# source calc_costsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "calc_costsC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the calc_costsV functions:
calc_costsV1 <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # estimate costs based on occupied state
  v_state_costs                           <- rep(0, length(v_occupied_state))                              # by default the cost for everyone is zero
  v_state_costs[v_occupied_state == "H"]  <- v_states_costs["H"]                                           # update the cost if healthy
  v_state_costs[v_occupied_state == "S1"] <- v_states_costs["S1"] + v_indi_costs[v_occupied_state == "S1"] # update the cost if sick
  v_state_costs[v_occupied_state == "S2"] <- v_states_costs["S2"] + v_indi_costs[v_occupied_state == "S2"] # update the cost if sicker
  v_state_costs[v_occupied_state == "D"]  <- v_states_costs["D"]                                           # update the cost if dead
  
  return(v_state_costs)                                                                                    # return the costs
}

calc_costsV2 <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # estimate costs based on occupied state
  v_state_costs                        <- rep(0, length(v_occupied_state))                        # by default the cost for everyone is zero
  v_state_costs[v_occupied_state == 1] <- v_states_costs[1]                                       # update the cost if healthy
  v_state_costs[v_occupied_state == 2] <- v_states_costs[2] + v_indi_costs[v_occupied_state == 2] # update the cost if sick
  v_state_costs[v_occupied_state == 3] <- v_states_costs[3] + v_indi_costs[v_occupied_state == 3] # update the cost if sicker
  v_state_costs[v_occupied_state == 4] <- v_states_costs[4]                                       # update the cost if dead
  
  return(v_state_costs)                                                                           # return the costs
}

#------------------------------------------------------------------------------#

# define function inputs:
## General parameters
seed    <- 1234                            # random number generator state
num_i   <- 1e6                             # number of simulated individuals
v_states_names <- c("H","S1", "S2", "D")   # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_states_index <- 1:length(v_states_names) # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
set.seed(seed)                             # set the seed to ensure reproducible samples below
v_occupied_state <- sample(                # sample current health state
  x       = v_states_names,                # from the four health states
  size    = num_i,                         # sample a state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)         # sample 75% as healthy, 20% as sick, and 5% as sicker
)
set.seed(seed)                             # set the seed to ensure reproducible samples below
v_occupied_state2 <- sample(               # sample current health state
  x       = v_states_index,                # from the four health states
  size    = num_i,                         # sample a state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)         # sample 75% as healthy, 20% as sick, and 5% as sicker
)
set.seed(seed)                             # set the seed to ensure reproducible samples below
v_time_in_state  <- sample(                # sample time in current health state
  x       = 1:10,                          # for values between 1 and 10
  size    = num_i,                         # sample time in state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = rep(                           # values 1 to 10 will have equal chances of being sampled
    x     = 1/length(1:10),                # the chance of sampling any the values is 1/the number of values
    times = length(1:10)                   # repeat the probabilities as many times as there as values
  )
)
mean_age            <- 50                  # mean age in the simulated population
sd_age              <- 3                   # standard deviation of the age in the simulated population
prop_females        <- 0.6                 # proportion of females in the simulated population
prop_males          <- 1 - prop_females    # proportion of males in the simulated population
set.seed(seed)                             # set a seed to ensure reproducible samples
v_occupied_state <- sample(                # sample current health state
  x       = v_states_names,                # from the four health states
  size    = num_i,                         # sample a state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)         # sample 75% as healthy, 20% as sick, and 5% as sicker
)
v_time_in_state  <- sample(                # sample time in current health state
  x       = 1:10,                          # for values between 1 and 10
  size    = num_i,                         # sample time in state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = rep(                           # values 1 to 10 will have equal chances of being sampled
    x     = 1/length(1:10),                # the chance of sampling any the values is 1/the number of values
    times = length(1:10)                   # repeat the probabilities as many times as there as values
  )
)
m_indi_features  <- cbind(                 # simulate individuals characteristics
  "age" = rnorm(                           # get random samples for 'age' from a normal distribution
    n = num_i,
    mean = mean_age,
    sd = sd_age
  ),
  "sex" = sample(                          # get random samples for 'sex' based on sex distribution
    x = c(0, 1),
    size = num_i,
    replace = TRUE,
    prob = c(prop_females, prop_males)
  )
)

## Cost inputs
c_H       <- 2000                             # cost of remaining one cycle healthy
c_S1      <- 4000                             # cost of remaining one cycle sick
c_S2      <- 15000                            # cost of remaining one cycle sicker
c_D       <- 0                                # cost associated with being dead
c_age_cof <- 11.5                             # cost age coefficient
c_sex_cof <- 300                              # cost sex coefficient, where 0 is female and 1 is male

v_states_costs <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_cost_coeffs <- c(                           # pack the cost regression coefficients in a vector
  "age" = c_age_cof, "sex" = c_sex_cof
)

#------------------------------------------------------------------------------#

# run the calc_costsV function:
R_results1 <- calc_costsV1(
  v_occupied_state = v_occupied_state,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
R_results2 <- calc_costsV2(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
# run the calc_costsC function:
C_results1 <- calc_costsC1(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
C_results2 <- calc_costsC2(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
C_results3 <- calc_costsC3(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
C_results4 <- calc_costsC4(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
C_results5 <- calc_costsC5(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
C_results6 <- calc_costsC6(
  v_occupied_state = v_occupied_state2,
  v_states_costs = v_states_costs,
  m_indi_features = m_indi_features,
  v_cost_coeffs = v_cost_coeffs
)
# check results
identical(R_results1, R_results1)
identical(R_results2, C_results1)
identical(C_results1, C_results2[, 1])
identical(C_results2, C_results3)
identical(C_results3, C_results4)
identical(C_results4, C_results5)
identical(C_results5, C_results6)

#------------------------------------------------------------------------------#

# benchmark the functions

calc_costs_RvC <- bench::mark(
  "R_1" = calc_costsV1(
    v_occupied_state = v_occupied_state,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "R_2" = calc_costsV2(
    v_occupied_state = v_occupied_state,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_1" = calc_costsC1(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_2" = calc_costsC2(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_3" = calc_costsC3(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_4" = calc_costsC4(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_5" = calc_costsC5(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_6" = calc_costsC6(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  check = FALSE
)

calc_costs_RvC[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

calc_costs_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = calc_costsV1(
    v_occupied_state = v_occupied_state,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "R_2" = calc_costsV2(
    v_occupied_state = v_occupied_state,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_1" = calc_costsC1(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_2" = calc_costsC2(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_3" = calc_costsC3(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_4" = calc_costsC4(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_5" = calc_costsC5(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  ),
  "C_6" = calc_costsC6(
    v_occupied_state = v_occupied_state2,
    v_states_costs = v_states_costs,
    m_indi_features = m_indi_features,
    v_cost_coeffs = v_cost_coeffs
  )
)

calc_costs_RvC2
plot(calc_costs_RvC2)