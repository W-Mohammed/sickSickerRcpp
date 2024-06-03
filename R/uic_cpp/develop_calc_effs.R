# clear R working global environment
rm(list = ls())

# source calc_effsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "calc_effsC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the calc_effsV functions:
calc_effsV1 <- function (
    v_occupied_state,
    v_states_utilities,
    m_indi_features,
    v_util_coeffs,
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- (m_indi_features %*% v_util_coeffs)[,1]
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == "S1"] <- v_util_t_decs["S1"] * v_time_in_state[v_occupied_state == "S1"]
  time_decrement[v_occupied_state == "S2"] <- v_util_t_decs["S2"] * v_time_in_state[v_occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- v_ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  v_state_utility                           <- rep(0, length(v_occupied_state))                               # by default the utility for everyone is zero
  v_state_utility[v_occupied_state == "H"]  <- v_states_utilities["H"]  + decrement[v_occupied_state == "H"]  # update the utility if healthy
  v_state_utility[v_occupied_state == "S1"] <- v_states_utilities["S1"] + decrement[v_occupied_state == "S1"] # update the utility if sick
  v_state_utility[v_occupied_state == "S2"] <- v_states_utilities["S2"] + decrement[v_occupied_state == "S2"] # update the utility if sicker
  v_state_utility[v_occupied_state == "D"]  <- v_states_utilities["D"]                                        # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                                    # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                                                               # return the QALYs
}

calc_effsV2 <- function (
    v_occupied_state,
    v_states_utilities,
    m_indi_features,
    v_util_coeffs,
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- (m_indi_features %*% v_util_coeffs)[,1]
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == 2] <- v_util_t_decs[1] * v_time_in_state[v_occupied_state == 2]
  time_decrement[v_occupied_state == 3] <- v_util_t_decs[2] * v_time_in_state[v_occupied_state == 3]
  
  # estimate total decrements
  decrement <- v_ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  v_state_utility                        <- rep(0, length(v_occupied_state))                         # by default the utility for everyone is zero
  v_state_utility[v_occupied_state == 1] <- v_states_utilities[1] + decrement[v_occupied_state == 1] # update the utility if healthy
  v_state_utility[v_occupied_state == 2] <- v_states_utilities[2] + decrement[v_occupied_state == 2] # update the utility if sick
  v_state_utility[v_occupied_state == 3] <- v_states_utilities[3] + decrement[v_occupied_state == 3] # update the utility if sicker
  v_state_utility[v_occupied_state == 4] <- v_states_utilities[4]                                    # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                           # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                                                      # return the QALYs
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

## Utility inputs
u_H       <- 1            # utility when healthy
u_S1      <- 0.75         # utility when sick, untreated
u_S2      <- 0.5          # utility when sicker, untreated
u_D       <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
u_sex_cof <- -0.015       # utility sex coefficient, where 0 is female and 1 is male
ru_S1     <- -0.0015      # change in utility of individuals with every additional year being sick
ru_S2     <- -0.0020      # change in utility of individuals with every additional year being sicker

v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)
v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof, "sex" = u_sex_cof
)
v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "S1" = ru_S1, "S2" = ru_S2
)

#------------------------------------------------------------------------------#

# run the calc_effsV function:
R_results1 <- calc_effsV1(
  v_occupied_state = v_occupied_state,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
R_results2 <- calc_effsV2(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
# run the calc_effsC function:
C_results0 <- calc_effsC0(
  v_occupied_state = v_occupied_state,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results1 <- calc_effsC1(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results2 <- calc_effsC2(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results3 <- calc_effsC3(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results4 <- calc_effsC4(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results5 <- calc_effsC5(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results6 <- calc_effsC6(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
C_results7 <- calc_effsC7(
  v_occupied_state = v_occupied_state2,
  v_states_utilities = v_states_utilities,
  cycle_length = 1,
  m_indi_features = m_indi_features,
  v_util_coeffs = v_util_coeffs,
  v_util_t_decs = v_util_t_decs,
  v_time_in_state = v_time_in_state
)
# check results
identical(R_results1, C_results0)
identical(R_results2, C_results1)
identical(C_results1, C_results2[, 1])
identical(C_results2[, 1], C_results3[, 1])
identical(C_results3[, 1], C_results4[, 1])
identical(C_results4[, 1], C_results5[, 1])
identical(C_results5[, 1], C_results6[, 1])
identical(C_results6[, 1], C_results7[, 1])

#------------------------------------------------------------------------------#

# benchmark the functions

calc_effs_RvC <- bench::mark(
  "R_1" = calc_effsV1(
    v_occupied_state = v_occupied_state,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "R_2" = calc_effsV2(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_0" = calc_effsC0(
    v_occupied_state = v_occupied_state,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_1" = calc_effsC1(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_2" = calc_effsC1(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_3" = calc_effsC3(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_4" = calc_effsC4(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_5" = calc_effsC5(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_6" = calc_effsC6(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_7" = calc_effsC7(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  check = FALSE
)

calc_effs_RvC[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

calc_effs_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = calc_effsV1(
    v_occupied_state = v_occupied_state,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "R_2" = calc_effsV2(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_0" = calc_effsC0(
    v_occupied_state = v_occupied_state,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_1" = calc_effsC1(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_2" = calc_effsC1(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_3" = calc_effsC3(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_4" = calc_effsC4(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_5" = calc_effsC5(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_6" = calc_effsC6(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  ),
  "C_7" = calc_effsC7(
    v_occupied_state = v_occupied_state2,
    v_states_utilities = v_states_utilities,
    cycle_length = 1,
    m_indi_features = m_indi_features,
    v_util_coeffs = v_util_coeffs,
    v_util_t_decs = v_util_t_decs,
    v_time_in_state = v_time_in_state
  )
)

calc_effs_RvC2
plot(calc_effs_RvC2)
saveRDS(object = calc_effs_RvC2, file = "calc_effs_RvC2")