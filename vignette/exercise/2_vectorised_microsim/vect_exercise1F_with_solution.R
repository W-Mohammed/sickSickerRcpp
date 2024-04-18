## @knitr solution_microSimV_function

# Healthy-Sick-Sicker-Dead microsimulation model - vectorised version

# clear R's session memory (Global Environment) 
rm(list = ls())

# 0. Execute the lines between the #---# below.
#------------------------------------------------------------------------------#

                        ### Defining model functions ###

# Define model functions
## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by each individual at cycle 't' and the time spent
### in the states
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
        test = rowSums(m_probs) == 1,               # check if the transition probabilities add up to 1
        yes = return(m_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )   
    }
  )
}

## Sample Health State function
### This function identifies the health state each individual will transition
### to in the next model cycle
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

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by the individuals at cycle 't' and the individual characteristics
calc_costsV <- function (
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

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by each individuals at cycle 't',
### time spent in the states and the cycle_length (measured in years)
calc_effsV <- function (
    v_occupied_state, 
    v_states_utilities,
    m_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- m_indi_features %*% v_util_coeffs
  
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

## Run Microsimulation function
### This function runs the microsimulation function of the Healthy-Sick-Dead 
### model
run_microSimV <- function(
    v_starting_states, 
    num_sim, 
    num_cycles, 
    m_indi_features,
    v_states_names,
    v_states_costs,
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs,
    cycle_length = 1,
    starting_seed = 1) {
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(
    nrow = num_sim, 
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_sim,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )
  
  # set the seed for every individual for the random number generator
  set.seed(starting_seed)
  
  # initialize parameter tracking time in current state
  v_time_in_state <- rep(0, times = num_sim)
  
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
    
    # keep track of time in current health state
    stayed                   <- m_States[, t] == m_States[, t + 1] # check if remains in current state at 't + 1'
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        # increment time spent in state
    v_time_in_state[!stayed] <- 0                                  # reset time once transitioned
    
  } # close the loop for the cycles 't' 

  # Compute costs and QALys:
  v_total_costs <- rowSums(m_Costs)        # calculate total cost per individual
  v_total_qalys <- rowSums(m_Effs)         # calculate total QALYs per individual
  total_costs   <- mean(v_total_costs)     # calculate average discounted cost 
  total_qalys   <- mean(v_total_qalys)     # calculate average discounted QALYs
  
  # store the results in a list:
  results <- list(
    m_States      = m_States, 
    m_Costs       = m_Costs, 
    m_Effs        = m_Effs, 
    v_total_costs = v_total_costs, 
    v_total_qalys = v_total_qalys, 
    total_costs   = total_costs, 
    total_qalys   = total_qalys
  )
  
  # return the results
  return(results)
}


                        ### Defining model parameters ###

# Define model inputs
## General parameters
num_sim             <- 100               # number of simulated individuals
num_cycles          <- 30                # time horizon if each cycle is a year long
cycle_length        <- 1                 # length of cycle (in years)
seed                <- 1234              # random number generator state
wtp                 <- 30000             # Willingness to pay for each QALY ($)

## Population characteristics/features
mean_age            <- 50                # mean age in the simulated population
sd_age              <- 3                 # standard deviation of the age in the simulated population
prop_females        <- 0.6               # proportion of females in the simulated population
prop_males          <- 1 - prop_females  # proportion of males in the simulated population
set.seed(seed)                           # set a seed to ensure reproducible samples
m_indi_features     <- cbind(            # simulate individuals characteristics
  "age" = rnorm(                         # get random samples for 'age' from a normal distribution
    n = num_sim, 
    mean = mean_age, 
    sd = sd_age
  ),
  "sex" = sample(                        # get random samples for 'sex' based on sex distribution
    x = c(0, 1), 
    size = num_sim, 
    replace = TRUE, 
    prob = c(prop_females, prop_males)
  )
)

## Health states
v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_starting_states <- rep("H", num_sim)   # everyone begins in the healthy state 

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

#------------------------------------------------------------------------------#

### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)
# 1. Considering how the payoffs of the "no treatment" option defined in the 2
# lines above, define the payoffs for "treatment 1" and "treatment 2".

### Payoffs - treatment 1
v_states_costs1     <- c("H" = c_H, "S1" = c_S1_Trt1, "S2" = c_S2_Trt1, "D" = c_D)
v_states_utilities1 <- c("H" = u_H, "S1" = u_S1_Trt1, "S2" = u_S2_Trt1, "D" = u_D)
### Payoffs - treatment 2
v_states_costs2     <- c("H" = c_H, "S1" = c_S1_Trt2, "S2" = c_S2_Trt2, "D" = c_D)
v_states_utilities2 <- c("H" = u_H, "S1" = u_S1_Trt2, "S2" = u_S2_Trt2, "D" = u_D)




# 2. Print all payoffs vectors to the console. Are the states rewards defined
# correctly? i.e. are the payoffs set with the correct states names?
# HINT: review the "Cost and utility inputs" section above if in doubt

v_states_costs
v_states_costs1
v_states_costs2
v_states_utilities
v_states_utilities1
v_states_utilities2
# Yes, all states' rewards are defined correctly




# 3. Print the probabilities list. Are the probabilities set out correctly?

l_trans_probs
# Yes, the probabilities are set out correctly.




# 4. What are the arguments (parameters) of the 'run_microSimV()' function?
# HINT: Check the definition of the 'run_microSimV()' function above.
# HINT: The function 'formals()' can be used to identify functions' arguments.
# Run '?formals' in the console to access the 'formals()' function help file.
# Run 'formals(run_microSimV)' in the console for R to list the argumnets.

# v_starting_states, num_sim, num_cycles, m_indi_features, v_states_names,
# v_states_costs, v_cost_coeffs, v_states_utilities, v_util_coeffs, 
# v_util_t_decs, l_trans_probs, cycle_length, and starting_seed 




# 5. The code below runs the simulation for the "no treatment" option. 
res_no_Trt <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_sim             = num_sim,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
# 5.1. Run the simulation for the two treatments ("treatment 1" and "treatment 2"),
# and assign the outputs to objects named "res_Trt1" and "res_Trt2", 
# respectively.
# HINT: set the correct payoffs before calling the 'run_microSimV()' function
# for each of the treatments.

## For being treated with treatment 1
res_Trt1 <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_sim             = num_sim,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs1,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities1,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
## For being treated with treatment 2
res_Trt2 <- run_microSimV(
  v_starting_states   = v_starting_states,
  num_sim             = num_sim,
  num_cycles          = num_cycles,
  m_indi_features     = m_indi_features,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs2,
  v_cost_coeffs       = v_cost_coeffs,
  v_states_utilities  = v_states_utilities2,
  v_util_coeffs       = v_util_coeffs,
  v_util_t_decs       = v_util_t_decs,
  l_trans_probs       = l_trans_probs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)




# 5.2. Print to the console the following data:
# - the outcomes (average costs/QALYs) associated with the "no treatment" choice
# - the average costs associated with "treatment 1", and
# - the average QALYs associated with "treatment 2"
# HINT: use the 'names()' function to retreive the names of the objects in any 
# of the results "res_" objects. E.g. 'names(res_Trt2)'
# HINT: use the 'View()' function to view the contents of the microsimulations
# results' objects.

# The outcomes (average costs/QALYs) associated with the "no treatment" choice
res_no_Trt$total_costs
res_no_Trt$total_qalys
# The average costs associated with "treatment 1"
res_Trt1$total_costs
# The average QALYs associated with "treatment 2"
res_Trt2$total_qalys




# 5.3. Compute the net monetary benefits (NMB) at $30,000 per QALY. Which of the
# three choices (no treatment, treatment 1 and treatment 2) is the 
# cost-effective one.
# HINT: the NMB = (QALYs * 30000) - Costs

nb_no_Trt <- res_no_Trt$total_qalys * wtp - res_no_Trt$total_costs
nb_Trt1   <- res_Trt1$total_qalys   * wtp - res_Trt1$total_costs
nb_Trt2   <- res_Trt2$total_qalys   * wtp - res_Trt2$total_costs
which.max(c("no treatment" = nb_no_Trt, "treatment 1" = nb_Trt1, "treatment 2" = nb_Trt2))



