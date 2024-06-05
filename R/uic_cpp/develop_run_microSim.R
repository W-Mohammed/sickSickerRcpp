# clear R working global environment
rm(list = ls())

# source calc_discount_wtsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "run_microSimC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the run_microSim() and its nested functions:
# update_probsV
update_probsV <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs,
    v_time_in_state) {
  
  with(
    data = l_trans_probs,
    expr = {
      
      # update probabilities of death after first converting them to rates and applying the rate ratio
      r_S1D <-  - log(1 - p_S1D)
      r_S2D <-  - log(1 - p_S2D)
      # calculate p_S1D/p_S2D conditional on current state and duration of being sick/sicker
      p_S1D  <- 1 - exp(- r_S1D * (1 + v_time_in_state[v_occupied_state == "S1"] * rp_S1))
      p_S2D  <- 1 - exp(- r_S2D * (1 + v_time_in_state[v_occupied_state == "S2"] * rp_S2))
      
      p_S1S1 <- 1 - p_S1S2 - p_S1H - p_S1D
      p_HD   <- rep(p_HD, length(which(v_occupied_state == "H")))
      p_DD   <- rep(1, length(v_occupied_state[v_occupied_state == "D"]))
      # Create a state transition probabilities matrix
      m_probs <- matrix(
        nrow = length(v_time_in_state), # a row for each individual
        ncol = length(v_states_names),  # a column for each state
        dimnames = list(
          v_occupied_state,             # name each row based on the occupied state
          v_states_names                # give each column one of the states names
        )
      )
      
      # update m_probs with the appropriate probabilities
      m_probs[v_occupied_state == "H", ]  <- cbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD) # transition probabilities when healthy
      m_probs[v_occupied_state == "S1", ] <- cbind(p_S1H, p_S1S1, p_S1S2, p_S1D)     # transition probabilities when sick
      m_probs[v_occupied_state == "S2", ] <- cbind(0, 0, 1 - p_S2D, p_S2D)           # transition probabilities when sicker
      m_probs[v_occupied_state == "D", ]  <- cbind(0, 0, 0, p_DD)                    # transition probabilities when dead
      
      # sanity check
      ifelse(
        test = all(abs(rowSums(m_probs) - 1) < 1e-12), # check if the transition probabilities add up to 1
        yes = return(m_probs),                         # return the transition probabilities
        no = stop("Probabilities do not sum to 1")     # or produce an error
      )
    }
  )
}
# sampleV
sampleV <- function(
    m_trans_probs,
    v_states_names) {
  
  # create an upper triangular matrix of ones
  m_upper_tri <- upper.tri(
    x = diag(ncol(m_trans_probs)),
    diag = TRUE
  )
  
  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  colnames(m_cum_probs) <- v_states_names
  
  # ensure that the maximum cumulative probabilities are equal to 1
  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.000000)) {
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }
  
  # sample random values from Uniform standard distribution for each individual
  v_rand_values <- runif(n = nrow(m_trans_probs))
  
  # repeat each sampled value to have as many copies as the number of states
  m_rand_values <- matrix(
    data  = rep(
      x = v_rand_values,
      each = length(v_states_names)
    ),
    nrow  = nrow(m_trans_probs),
    ncol  = length(v_states_names),
    byrow = TRUE
  )
  
  # identify transitions, compare random samples to cumulative probabilities
  m_transitions <- m_rand_values > m_cum_probs # transitions from first state
  
  # sum transitions to identify health state in next cycle
  v_transitions <- rowSums(m_transitions)
  
  # identify health state to which each individual is transitioning
  v_health_states <- v_states_names[1 + v_transitions]
  
  return(v_health_states)
}
# calc_costsV
calc_costsV <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # estimate costs based on occupied state
  v_state_costs                           <- rep(NA, length(v_occupied_state))
  v_state_costs[v_occupied_state == "H"]  <- v_states_costs["H"]                                           # update the cost if healthy
  v_state_costs[v_occupied_state == "S1"] <- v_states_costs["S1"] + v_indi_costs[v_occupied_state == "S1"] # update the cost if sick
  v_state_costs[v_occupied_state == "S2"] <- v_states_costs["S2"] + v_indi_costs[v_occupied_state == "S2"] # update the cost if sicker
  v_state_costs[v_occupied_state == "D"]  <- v_states_costs["D"]                                           # update the cost if dead
  
  return(v_state_costs)                                                                                    # return the costs
}
# calc_effsV
calc_effsV <- function (
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
  v_state_utility                           <- rep(NA, length(v_occupied_state))
  v_state_utility[v_occupied_state == "H"]  <- v_states_utilities["H"]  + decrement[v_occupied_state == "H"]  # update the utility if healthy
  v_state_utility[v_occupied_state == "S1"] <- v_states_utilities["S1"] + decrement[v_occupied_state == "S1"] # update the utility if sick
  v_state_utility[v_occupied_state == "S2"] <- v_states_utilities["S2"] + decrement[v_occupied_state == "S2"] # update the utility if sicker
  v_state_utility[v_occupied_state == "D"]  <- v_states_utilities["D"]                                        # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                                    # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                                                               # return the QALYs
}
# calc_discount_wts
calc_discount_wts <- function(
    discount_rate,
    num_cycles,
    cycle_length) {
  
  # calculate discount weights based on the number & length (in years) of cycles
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  
  return(v_discount_wts)
}

# run_microSimV
run_microSimV <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs,
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs,
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 1,
    starting_seed = 1) {
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1
  )
  
  # set the seed for every individual for the random number generator
  set.seed(starting_seed)
  
  # initialize parameter tracking time in current state
  v_time_in_state <- rep(1, times = num_i)
  
  # get the initial health state
  m_States[, 1] <- v_starting_states
  
  # calculate the costs incurred in their starting health state
  m_Costs[, 1]  <- calc_costsV(
    v_occupied_state = m_States[, 1],
    v_states_costs   = v_states_costs,
    m_indi_features  = m_indi_features,
    v_cost_coeffs    = v_cost_coeffs
  )
  
  # calculate the QALYs accrued in their starting health state
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )
  
  # for each 't' of the 'num_cycles' cycles:
  for (t in 1:num_cycles) {
    # update the transition probabilities at cycle 't'
    m_trans_probs     <- update_probsV(
      v_states_names   = v_states_names,
      v_occupied_state = m_States[, t],
      l_trans_probs    = l_trans_probs,
      v_time_in_state  = v_time_in_state
    )
    
    # sample the health state at 't + 1'
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names
    )
    
    # keep track of time in state at 't + 1'
    stayed                   <- m_States[, t] == m_States[, t + 1] # check if remains in current state at 't + 1'
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        # increment time spent in state
    v_time_in_state[!stayed] <- 1                                  # reset time once transitioned
    
    # keep track of time in the model
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1
    
    # calculate the costs incurred in their 't + 1' health state
    m_Costs[, t + 1]  <- calc_costsV(
      v_occupied_state = m_States[, t + 1],
      v_states_costs   = v_states_costs,
      m_indi_features  = m_indi_features,
      v_cost_coeffs    = v_cost_coeffs
    )
    
    # calculate the QALYs accrued in their 't + 1' health state
    m_Effs[, t + 1]   <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      m_indi_features    = m_indi_features,
      v_util_coeffs      = v_util_coeffs,
      v_util_t_decs      = v_util_t_decs,
      v_time_in_state    = v_time_in_state,
      cycle_length       = cycle_length
    )
    
  } # close the loop for the cycles 't'
  
  # Calculate discount weights for both outcomes:
  v_c_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  v_e_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_QALYs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  # Compute costs and QALYs:
  v_total_costs <- rowSums(m_Costs)         # calculate total costs per individual
  v_total_qalys <- rowSums(m_Effs)          # calculate total QALYs per individual
  mean_costs    <- mean(v_total_costs)      # calculate average costs
  mean_qalys    <- mean(v_total_qalys)      # calculate average QALYs
  
  # Compute discounted costs and QALYs:
  v_total_Dcosts <- m_Costs %*% v_c_dsc_wts # calculate total discounted costs per individual
  v_total_Dqalys <- m_Effs  %*% v_e_dsc_wts # calculate total discounted QALYs per individual
  mean_Dcosts    <- mean(v_total_Dcosts)    # calculate average discounted costs
  mean_Dqalys    <- mean(v_total_Dqalys)    # calculate average discounted QALYs
  
  # store the results in a list:
  results <- list(
    m_States       = m_States,
    m_Costs        = m_Costs,
    m_Effs         = m_Effs,
    v_total_costs  = v_total_costs,
    v_total_qalys  = v_total_qalys,
    v_total_Dcosts = v_total_Dcosts,
    v_total_Dqalys = v_total_Dqalys,
    mean_costs     = mean_costs,
    mean_qalys     = mean_qalys,
    mean_Dcosts    = mean_Dcosts,
    mean_Dqalys    = mean_Dqalys
  )
  
  # return the results
  return(results)
}

# run_microSimVC
# The 'run_microSimVC()' function is a combined R and C++ implementation of the 
# run_microSim() function.
run_microSimVC <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_index,
    v_states_costs,
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs,
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 1,
    starting_seed = 1) {
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1
  )
  m_Costs <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1
  )
  m_Effs <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1
  )
  
  # set the seed for every individual for the random number generator
  set.seed(starting_seed)
  
  # initialize parameter tracking time in current state
  v_time_in_state <- rep(1, times = num_i)
  
  # get the initial health state
  m_States[, 1] <- v_starting_states
  
  # calculate the costs incurred in their starting health state
  m_Costs[, 1]  <- calc_costsC(
    v_occupied_state = m_States[, 1],
    v_states_costs   = v_states_costs,
    m_indi_features  = m_indi_features,
    v_cost_coeffs    = v_cost_coeffs
  )
  
  # calculate the QALYs accrued in their starting health state
  m_Effs[, 1]   <- calc_effsC(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )
  
  # for each 't' of the 'num_cycles' cycles:
  for (t in 1:num_cycles) {
    # update the transition probabilities at cycle 't'
    m_trans_probs     <- update_probsC(
      v_states_index   = v_states_index,
      v_occupied_state = m_States[, t],
      l_trans_probs    = l_trans_probs,
      v_time_in_state  = v_time_in_state
    )
    
    # sample the health state at 't + 1'
    m_States[, t + 1] <- sampleC(
      m_trans_probs    = m_trans_probs
    )
    
    # keep track of time in state at 't + 1'
    stayed                   <- m_States[, t] == m_States[, t + 1] # check if remains in current state at 't + 1'
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        # increment time spent in state
    v_time_in_state[!stayed] <- 1                                  # reset time once transitioned
    
    # keep track of time in the model
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1
    
    # calculate the costs incurred in their 't + 1' health state
    m_Costs[, t + 1]  <- calc_costsC(
      v_occupied_state = m_States[, t + 1],
      v_states_costs   = v_states_costs,
      m_indi_features  = m_indi_features,
      v_cost_coeffs    = v_cost_coeffs
    )
    
    # calculate the QALYs accrued in their 't + 1' health state
    m_Effs[, t + 1]   <- calc_effsC(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      m_indi_features    = m_indi_features,
      v_util_coeffs      = v_util_coeffs,
      v_util_t_decs      = v_util_t_decs,
      v_time_in_state    = v_time_in_state,
      cycle_length       = cycle_length
    )
    
  } # close the loop for the cycles 't'
  
  # Calculate discount weights for both outcomes:
  v_c_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  v_e_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_QALYs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  # Compute costs and QALYs:
  v_total_costs <- rowSums(m_Costs)         # calculate total costs per individual
  v_total_qalys <- rowSums(m_Effs)          # calculate total QALYs per individual
  mean_costs    <- mean(v_total_costs)      # calculate average costs
  mean_qalys    <- mean(v_total_qalys)      # calculate average QALYs
  
  # Compute discounted costs and QALYs:
  v_total_Dcosts <- m_Costs %*% v_c_dsc_wts # calculate total discounted costs per individual
  v_total_Dqalys <- m_Effs  %*% v_e_dsc_wts # calculate total discounted QALYs per individual
  mean_Dcosts    <- mean(v_total_Dcosts)    # calculate average discounted costs
  mean_Dqalys    <- mean(v_total_Dqalys)    # calculate average discounted QALYs
  
  # store the results in a list:
  results <- list(
    m_States       = m_States,
    m_Costs        = m_Costs,
    m_Effs         = m_Effs,
    v_total_costs  = v_total_costs,
    v_total_qalys  = v_total_qalys,
    v_total_Dcosts = v_total_Dcosts,
    v_total_Dqalys = v_total_Dqalys,
    mean_costs     = mean_costs,
    mean_qalys     = mean_qalys,
    mean_Dcosts    = mean_Dcosts,
    mean_Dqalys    = mean_Dqalys
  )
  
  # return the results
  return(results)
}

#------------------------------------------------------------------------------#

# define function inputs:
## General parameters
num_i               <- 1e6               # number of simulated individuals
num_cycles          <- 31                # time horizon if each cycle is a year long
cycle_length        <- 1                 # length of cycle (in years)
seed                <- 1234              # random number generator state
wtp                 <- 30000             # Willingness to pay for each QALY ($)
discount_rate_costs <- 0.03              # annual discount rate for costs
discount_rate_QALYs <- 0.015             # annual discount rate for health outcomes

## Population characteristics/features
mean_age            <- 50                # mean age in the simulated population
sd_age              <- 3                 # standard deviation of the age in the simulated population
prop_females        <- 0.6               # proportion of females in the simulated population
prop_males          <- 1 - prop_females  # proportion of males in the simulated population

generate_indi_features <- function(num_i) {
  set.seed(seed)                           # set a seed to ensure reproducible samples
  cbind(            # simulate individuals characteristics
    "age" = rnorm(                         # get random samples for 'age' from a normal distribution
      n = num_i,
      mean = mean_age,
      sd = sd_age
    ),
    "sex" = sample(                        # get random samples for 'sex' based on sex distribution
      x = c(0, 1),
      size = num_i,
      replace = TRUE,
      prob = c(prop_females, prop_males)
    )
  )
  
}
m_indi_features     <- generate_indi_features(num_i = num_i)
m_indi_features1    <- generate_indi_features(num_i = num_i)
m_indi_features2    <- generate_indi_features(num_i = num_i)

## Health states
v_states_names <- c("H","S1", "S2", "D")   # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_states_index <- 1:length(v_states_names) # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_starting_states  <- rep("H", num_i)      # everyone begins in the healthy state
v_starting_statesC <- rep(1, num_i)        # everyone begins in the healthy state

## Transition probabilities (per cycle)
p_HD   <- 0.005                          # probability to die when healthy
p_HS1  <- 0.15                           # probability to become sick when healthy
p_S1H  <- 0.5                            # probability to become healthy when sick
p_S1S2 <- 0.105         	               # probability to become sicker when sick
rr_S1  <- 3                              # rate ratio of death in sick vs healthy
rr_S2  <- 10            	               # rate ratio of death in sicker vs healthy
r_HD   <- -log(1 - p_HD)                 # rate of death in healthy
r_S1D  <- rr_S1 * r_HD                   # rate of death in sick
r_S2D  <- rr_S2 * r_HD  	               # rate of death in sicker
p_S1D  <- 1 - exp(- r_S1D)               # probability to die in sick
p_S2D  <- 1 - exp(- r_S2D)               # probability to die in sicker
rp_S1  <- 0.2                            # increase in mortality rate with every additional year being sick
rp_S2  <- 0.29                           # increase in mortality rate with every additional year being sicker

l_trans_probs <- list(                   # pack the transition probabilities and rates in a list
  "p_HD"   = p_HD,
  "p_HS1"  = p_HS1,
  "p_S1H"  = p_S1H,
  "p_S1S2" = p_S1S2,
  "p_S1D"  = p_S1D,
  "p_S2D"  = p_S2D,
  "rp_S1"  = rp_S1,
  "rp_S2"  = rp_S2
)

## Cost and utility inputs
c_H       <- 2000         # cost of remaining one cycle healthy
c_S1      <- 4000         # cost of remaining one cycle sick
c_S2      <- 15000        # cost of remaining one cycle sicker
c_S1_Trt1 <- c_S1 + 12000 # cost of remaining one cycle sick under treatment 1
c_S2_Trt1 <- c_S2 + 12000 # cost of remaining one cycle sicker under treatment 1
c_S1_Trt2 <- c_S1 + 11350 # cost of remaining one cycle sick under treatment 2
c_S2_Trt2 <- c_S2 + 11350 # cost of remaining one cycle sicker under treatment 2
c_D       <- 0            # cost associated with being dead
c_age_cof <- 11.5         # cost age coefficient
c_sex_cof <- 300          # cost sex coefficient, where 0 is female and 1 is male

v_cost_coeffs <- c(       # pack the cost regression coefficients in a vector
  "age" = c_age_cof, "sex" = c_sex_cof
)

u_H       <- 1            # utility when healthy
u_S1      <- 0.75         # utility when sick, untreated
u_S2      <- 0.5          # utility when sicker, untreated
u_S1_Trt1 <- u_S1 + 0.2   # utility when sick, treatment 1
u_S2_Trt1 <- u_S2         # utility when sicker, treatment 1
u_S1_Trt2 <- u_S1 + 0.15  # utility when sick, treatment 2
u_S2_Trt2 <- u_S2 + 0.05  # utility when sicker, treatment 2
u_D       <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
u_sex_cof <- -0.015       # utility sex coefficient, where 0 is female and 1 is male
ru_S1     <- -0.0015      # change in utility of individuals with every additional year being sick
ru_S2     <- -0.0020      # change in utility of individuals with every additional year being sicker

v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof, "sex" = u_sex_cof
)
v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "S1" = ru_S1, "S2" = ru_S2
)
### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)
### Payoffs - treatment 1
v_states_costs1     <- c("H" = c_H, "S1" = c_S1_Trt1, "S2" = c_S2_Trt1, "D" = c_D)
v_states_utilities1 <- c("H" = u_H, "S1" = u_S1_Trt1, "S2" = u_S2_Trt1, "D" = u_D)
### Payoffs - treatment 2
v_states_costs2     <- c("H" = c_H, "S1" = c_S1_Trt2, "S2" = c_S2_Trt2, "D" = c_D)
v_states_utilities2 <- c("H" = u_H, "S1" = u_S1_Trt2, "S2" = u_S2_Trt2, "D" = u_D)

#------------------------------------------------------------------------------#

# run the run_microSimV function:
R_results1 <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
# R_results2 <- run_microSimV(
#   v_starting_states   = v_starting_states,
#   num_i               = num_i,
#   num_cycles          = num_cycles,
#   m_indi_features     = m_indi_features,
#   v_states_names      = v_states_names,
#   v_states_costs      = v_states_costs1,
#   v_cost_coeffs       = v_cost_coeffs,
#   v_states_utilities  = v_states_utilities1,
#   v_util_coeffs       = v_util_coeffs,
#   v_util_t_decs       = v_util_t_decs,
#   l_trans_probs       = l_trans_probs,
#   discount_rate_costs = discount_rate_costs,
#   discount_rate_QALYs = discount_rate_QALYs,
#   cycle_length        = cycle_length,
#   starting_seed       = seed
# )
# R_results3 <- run_microSimV(
#   v_starting_states   = v_starting_states,
#   num_i               = num_i,
#   num_cycles          = num_cycles,
#   m_indi_features     = m_indi_features,
#   v_states_names      = v_states_names,
#   v_states_costs      = v_states_costs2,
#   v_cost_coeffs       = v_cost_coeffs,
#   v_states_utilities  = v_states_utilities2,
#   v_util_coeffs       = v_util_coeffs,
#   v_util_t_decs       = v_util_t_decs,
#   l_trans_probs       = l_trans_probs,
#   discount_rate_costs = discount_rate_costs,
#   discount_rate_QALYs = discount_rate_QALYs,
#   cycle_length        = cycle_length,
#   starting_seed       = seed
# )
# run the run_microSimVC function:
RC_results1 <- run_microSimVC(
  v_starting_states   = v_starting_statesC,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_index      = v_states_index,
  v_states_costs      = v_states_costs,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)

# run the run_microSimC function:
C_results0 <- run_microSimC0(
  v_starting_states   = v_starting_states,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features1,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
C_results1 <- run_microSimC1(
  v_starting_states   = v_starting_statesC,
  num_i               = num_i,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features2,
  v_states_index      = v_states_index,
  v_states_costs      = v_states_costs,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  discount_rate_costs = discount_rate_costs,
  discount_rate_QALYs = discount_rate_QALYs,
  cycle_length        = cycle_length,
  starting_seed       = seed,
  age_column_index    = which(colnames(m_indi_features) == "age") - 1
)
# check results
## Combined (RC):
identical(R_results1$mean_costs, RC_results1$mean_costs)
identical(R_results1$v_total_costs, RC_results1$v_total_costs)
testthat::expect_equal(R_results1$mean_costs, RC_results1$mean_costs)
testthat::expect_equal(R_results1$v_total_costs, RC_results1$v_total_costs)
R_results1$mean_costs - RC_results1$mean_costs
identical(R_results1$mean_qalys, RC_results1$mean_qalys)
identical(R_results1$v_total_qalys, RC_results1$v_total_qalys)
testthat::expect_equal(R_results1$mean_qalys, RC_results1$mean_qalys)
testthat::expect_equal(R_results1$v_total_qalys, RC_results1$v_total_qalys)
R_results1$mean_qalys - RC_results1$mean_qalys

identical(R_results1$mean_Dcosts, RC_results1$mean_Dcosts)
identical(R_results1$v_total_Dcosts, RC_results1$v_total_Dcosts)
testthat::expect_equal(R_results1$mean_Dcosts, RC_results1$mean_Dcosts)
testthat::expect_equal(R_results1$v_total_Dcosts, RC_results1$v_total_Dcosts)
R_results1$mean_Dcosts - RC_results1$mean_Dcosts
identical(R_results1$mean_Dqalys, RC_results1$mean_Dqalys)
identical(R_results1$v_total_Dqalys, RC_results1$v_total_Dqalys)
testthat::expect_equal(R_results1$mean_Dqalys, RC_results1$mean_Dqalys)
testthat::expect_equal(R_results1$v_total_Dqalys, RC_results1$v_total_Dqalys)
R_results1$mean_Dqalys - RC_results1$mean_Dqalys

## Rcpp (C0):
identical(R_results1$mean_costs, C_results0$mean_costs)
identical(R_results1$v_total_costs, C_results0$v_total_costs)
testthat::expect_equal(R_results1$mean_costs, C_results0$mean_costs)
testthat::expect_equal(R_results1$v_total_costs, C_results0$v_total_costs)
R_results1$mean_costs - C_results0$mean_costs
identical(R_results1$mean_qalys, C_results0$mean_qalys)
identical(R_results1$v_total_qalys, C_results0$v_total_qalys)
testthat::expect_equal(R_results1$mean_qalys, C_results0$mean_qalys)
testthat::expect_equal(R_results1$v_total_qalys, C_results0$v_total_qalys)
R_results1$mean_qalys - C_results0$mean_qalys

identical(R_results1$mean_Dcosts, C_results0$mean_Dcosts)
identical(R_results1$v_total_Dcosts[,1], C_results0$v_total_Dcosts)
testthat::expect_equal(R_results1$mean_Dcosts, C_results0$mean_Dcosts)
testthat::expect_equal(R_results1$v_total_Dcosts[,1], C_results0$v_total_Dcosts)
R_results1$mean_Dcosts - C_results0$mean_Dcosts
identical(R_results1$mean_Dqalys, C_results0$mean_Dqalys)
identical(R_results1$v_total_Dqalys[,1], C_results0$v_total_Dqalys)
testthat::expect_equal(R_results1$mean_Dqalys, C_results0$mean_Dqalys)
testthat::expect_equal(R_results1$v_total_Dqalys[,1], C_results0$v_total_Dqalys)
R_results1$mean_Dqalys - C_results0$mean_Dqalys

## Arma (C1):
identical(R_results1$mean_costs, C_results1$mean_costs)
identical(R_results1$v_total_costs, C_results1$v_total_costs[,1])
testthat::expect_equal(R_results1$mean_costs, C_results1$mean_costs)
testthat::expect_equal(R_results1$v_total_costs, C_results1$v_total_costs[,1])
R_results1$mean_costs - C_results1$mean_costs
identical(R_results1$mean_qalys, C_results1$mean_qalys)
identical(R_results1$v_total_qalys, C_results1$v_total_qalys[,1])
testthat::expect_equal(R_results1$mean_qalys, C_results1$mean_qalys)
testthat::expect_equal(R_results1$v_total_qalys, C_results1$v_total_qalys[,1])
R_results1$mean_qalys - C_results1$mean_qalys

identical(R_results1$mean_Dcosts, C_results1$mean_Dcosts)
identical(R_results1$v_total_Dcosts, C_results1$v_total_Dcosts)
testthat::expect_equal(R_results1$mean_Dcosts, C_results1$mean_Dcosts)
testthat::expect_equal(R_results1$v_total_Dcosts, C_results1$v_total_Dcosts)
R_results1$mean_Dcosts - C_results1$mean_Dcosts
identical(R_results1$mean_Dqalys, C_results1$mean_Dqalys)
identical(R_results1$v_total_Dqalys, C_results1$v_total_Dqalys)
testthat::expect_equal(R_results1$mean_Dqalys, C_results1$mean_Dqalys)
testthat::expect_equal(R_results1$v_total_Dqalys, C_results1$v_total_Dqalys)
R_results1$mean_Dqalys - C_results1$mean_Dqalys

#------------------------------------------------------------------------------#

# benchmark the functions
run_microSim_RvC <- bench::press(
  num_i = seq(1e4, 1e6, length.out = 10),
  {
    m_indi_features     <- generate_indi_features(num_i = num_i)
    m_indi_features1    <- generate_indi_features(num_i = num_i)
    m_indi_features2    <- generate_indi_features(num_i = num_i)
    v_starting_states   <- rep("H", num_i)
    v_starting_statesC  <- rep(1, num_i)
    
    bench::mark(
      "R_1" = run_microSimV(
        v_starting_states   = v_starting_states,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features,
        v_states_names      = v_states_names,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      "RC" = run_microSimVC(
        v_starting_states   = v_starting_statesC,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features,
        v_states_index      = v_states_index,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      "C_0" = run_microSimC0(
        v_starting_states   = v_starting_states,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features1,
        v_states_names      = v_states_names,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      "C_1" = run_microSimC1(
        v_starting_states   = v_starting_statesC,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features2,
        v_states_index      = v_states_index,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed,
        age_column_index    = which(colnames(m_indi_features) == "age") - 1
      ),
      check = FALSE
    )
  }
)

run_microSim_RvC |>
  summary() |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = as.numeric(num_i),
      y = as.numeric(median),
      color = as.character(expression),
      group = as.character(expression)
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::labs(
    x = "Number of simulated individuals",
    y = "Median (s)",
    color = "Implementation"
  ) +
  ggplot2::theme(legend.position = "top")

run_microSim_RvC1 <- bench::mark(
  "R_1" = run_microSimV(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "RC" = run_microSimVC(
    v_starting_states   = v_starting_statesC,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_index      = v_states_index,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "C_0" = run_microSimC0(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features1,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "C_1" = run_microSimC1(
    v_starting_states   = v_starting_statesC,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features2,
    v_states_index      = v_states_index,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed,
    age_column_index    = which(colnames(m_indi_features2) == "age") - 1
  ),
  check = FALSE
)

run_microSim_RvC1[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

# Use the microbenchmark::microbenchmark()
run_microSim_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = run_microSimV(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "RC" = run_microSimVC(
    v_starting_states   = v_starting_statesC,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_index      = v_states_index,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "C_0" = run_microSimC0(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features1,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  "C_1" = run_microSimC1(
    v_starting_states   = v_starting_statesC,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features2,
    v_states_index      = v_states_index,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed,
    age_column_index    = which(colnames(m_indi_features2) == "age") - 1
  )
)

run_microSim_RvC2
plot(run_microSim_RvC2)
saveRDS(object = calc_effs_RvC2, file = "calc_effs_RvC2")