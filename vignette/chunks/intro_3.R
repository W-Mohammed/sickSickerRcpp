## @knitr demo_microsim_surv

# Healthy-Sick-Sicker-Dead microsimulation model

# clear R's session memory (Global Environment) 
rm(list = ls())

#------------------------------------------------------------------------------#

### Defining model functions ###

# Define model functions
## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by individual 'i' at cycle 't'
update_probs <- function(
    v_states_names,
    occupied_state,
    l_surv_params, 
    time_in_state,
    cycle_length = 1) { 
  
  with(
    data = l_surv_params,
    expr = {
      
      # Calculate delta cumulative hazards
      h_HS1 <- h_HD <- h_S1D <- 0
      ## H to S1: change in the cumulative hazard during one cycle being healthy
      h_HS1[occupied_state == "H"] <-
        pweibull(
          q = (time_in_state[occupied_state == "H"] - 1) * cycle_length,
          scale = scale_HS1,
          shape = shape_HS1,
          lower.tail = FALSE,
          log.p = TRUE
        ) -
        pweibull(
          q = time_in_state[occupied_state == "H"] * cycle_length,
          scale = scale_HS1,
          shape = shape_HS1,
          lower.tail = FALSE,
          log.p = TRUE
        )
      ## H to D: change in the cumulative hazard during one cycle being healthy
      h_HD[occupied_state == "H"] <-
        pweibull(
          q = (time_in_state[occupied_state == "H"] - 1) * cycle_length,
          scale = scale_HD,
          shape = shape_HD,
          lower.tail = FALSE,
          log.p = TRUE
        ) -
        pweibull(
          q = time_in_state[occupied_state == "H"] * cycle_length,
          scale = scale_HD,
          shape = shape_HD,
          lower.tail = FALSE,
          log.p = TRUE
        )
      ## S1 to D: change in the cumulative hazard during one cycle being sick
      h_S1D[occupied_state == "S1"] <-
        pweibull(
          q = (time_in_state[occupied_state == "S1"] - 1) * cycle_length,
          scale = scale_S1D,
          shape = shape_S1D,
          lower.tail = FALSE,
          log.p = TRUE
        ) -
        pweibull(
          q = time_in_state[occupied_state == "S1"] * cycle_length,
          scale = scale_S1D,
          shape = shape_S1D,
          lower.tail = FALSE,
          log.p = TRUE
        )
      
      # probabilities leaving H
      h_H   <- h_HS1 + h_HD
      p_HH  <- exp(-h_H)
      p_HS1 <- h_HS1/h_H * (1 - exp(-h_H))
      p_HD  <- h_HD/h_H  * (1 - exp(-h_H))
      
      # probabilities leaving S1
      p_S1S1 <- exp(-h_S1D)
      p_S1D  <- 1 - exp(-h_S1D)
      
      # create vector of state transition probabilities
      v_probs        <- rep(NA, length(v_states_names))
      names(v_probs) <- v_states_names
      
      # update v_probs with the appropriate probabilities   
      v_probs[occupied_state == "H"]  <- c(p_HH, p_HS1, p_HD) # transition probabilities when healthy
      v_probs[occupied_state == "S1"] <- c(0, p_S1S1, p_S1D)  # transition probabilities when sick
      v_probs[occupied_state == "D"]  <- c(0, 0, 1 )          # transition probabilities when dead   
      
      # sanity check
      ifelse(
        test = sum(v_probs) == 1,                   # check if the transition probabilities add up to 1
        yes = return(v_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )   
    }
  )
}

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by individual 'i' at cycle 't'
calc_costs <- function (
    occupied_state,
    v_states_costs,
    v_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  indi_costs <- 0
  for (variable in names(v_cost_coeffs)) {
    indi_costs <- indi_costs + (v_indi_features[[variable]] * v_cost_coeffs[[variable]])
  }                    
  
  # estimate costs based on occupied state
  state_costs <- 0                                                         # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]               # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] + indi_costs # update the cost if sick

  return(state_costs)                                         # return the costs
}

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by individual 'i' at cycle 't' and
### the cycle_length (measured in years)
calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    v_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  ind_decrement <- 0
  for (variable in names(v_util_coeffs)) {
    ind_decrement <- ind_decrement + (v_indi_features[[variable]] * v_util_coeffs[[variable]])
  }
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- 0
  time_decrement[occupied_state == "S1"] <-  v_util_t_decs["S1"] * time_in_state[occupied_state == "S1"]

  # estimate total decrements
  decrement <- ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  state_utility <- 0                                                            # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  + decrement # update the utility if healthy
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] + decrement # update the utility if sick

  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}

## Run Microsimulation function
### This function runs the microsimulation function of the Healthy-Sick-Dead model
run_microSim <- function(
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
    l_surv_params,
    cycle_length = 1,
    starting_seed = 1) {
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(
    nrow = num_sim, 
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_sim,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )
  
  # for each 'i' of the 'num_sim' simulated individual:
  for (i in 1:num_sim) {       
    # set the seed for every individual for the random number generator
    set.seed(starting_seed + i)
    
    # initialize parameter tracking time in current state
    time_in_state <- 1
    
    # get individual features
    v_indi_features <- m_indi_features[i,] 
    
    # get the initial health state
    m_States[i, 1] <- v_starting_states[i]
    
    # calculate the costs incurred in their starting health state
    m_Costs[i, 1]  <- calc_costs(          
      occupied_state  = m_States[i, 1],
      v_states_costs  = v_states_costs,
      v_indi_features = v_indi_features,
      v_cost_coeffs   = v_cost_coeffs
    )
    
    # calculate the QALYs accrued in their starting health state
    m_Effs[i, 1]   <- calc_effs(
      occupied_state     = m_States[i, 1], 
      v_states_utilities = v_states_utilities,
      v_indi_features    = v_indi_features,
      v_util_coeffs      = v_util_coeffs, 
      v_util_t_decs      = v_util_t_decs,
      time_in_state      = time_in_state,
      cycle_length       = cycle_length
    )
    
    # for each 't' of the 'num_cycles' cycles:
    for (t in 1:num_cycles) {       
      # update the transition probabilities at cycle 't'
      v_trans_probs <- update_probs(
        v_states_names = v_states_names,
        occupied_state = m_States[i, t],
        l_surv_params  = l_surv_params,
        time_in_state  = time_in_state,
        cycle_length   = cycle_length
      ) 
      
      # sample the health state at 't + 1' 
      m_States[i, t + 1] <- sample(
        x = v_states_names, 
        prob = v_trans_probs, 
        size = 1
      )  
      
      # calculate the costs incurred in their 't + 1' health state
      m_Costs[i, t + 1]  <- calc_costs(
        occupied_state  = m_States[i, t + 1],
        v_states_costs  = v_states_costs,
        v_indi_features = v_indi_features,
        v_cost_coeffs   = v_cost_coeffs
      )
      
      # calculate the QALYs accrued in their 't + 1' health state
      m_Effs[i, t + 1]   <- calc_effs(
        occupied_state     = m_States[i, t + 1], 
        v_states_utilities = v_states_utilities,
        v_indi_features    = v_indi_features,
        v_util_coeffs      = v_util_coeffs, 
        v_util_t_decs      = v_util_t_decs,
        time_in_state      = time_in_state,
        cycle_length       = cycle_length
      )
      
      # keep track of time in current health state
      stayed <- m_States[i, t] == m_States[i, t + 1]      # check if remains in current state at 't + 1'
      time_in_state[stayed]  <- time_in_state[stayed] + 1 # increment time spent in state
      time_in_state[!stayed] <- 1                         # reset time once transitioned
      
    } # close the loop for the cycles 't' 
    # display the progress of the simulation
    if(i/100 == round(i/100,0)) {
      cat('\r', paste(i/num_sim * 100, "% done", sep = " "))
    }
    
  } # close the loop for the individuals 'i' 
  
  # Step 5:
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

#------------------------------------------------------------------------------#

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
v_states_names <- c("H","S1", "D")       # the model states: Healthy (H), Sick (S1), Dead (D)
v_starting_states <- rep("H", num_sim)   # everyone begins in the healthy state 

## Transition probabilities (per cycle)
scale_HS1   <- 5                         # scale parameters for H -> S1 Weibull distribution
shape_HS1   <- 2.5                       # shape parameters for H -> S1 Weibull distribution
scale_HD    <- 10                        # scale parameters for H -> D Weibull distribution
shape_HD    <- 2.5                       # shape parameters for H -> D Weibull distribution
scale_S1D   <- 7.5                       # scale parameters for S1 -> D Weibull distribution
shape_S1D   <- 2.5                       # shape parameters for S1 -> D Weibull distribution

l_surv_params <- list(                   # pack the transition probabilities and rates in a list
  "scale_HS1" = scale_HS1, 
  "shape_HS1" = shape_HS1, 
  "scale_HD"  = scale_HD, 
  "shape_HD"  = shape_HD, 
  "scale_S1D" = scale_S1D, 
  "shape_S1D" = shape_S1D
)

## Cost and utility inputs 
c_H       <- 2000         # cost of remaining one cycle healthy
c_S1      <- 4000         # cost of remaining one cycle sick
c_S2      <- 15000        # cost of remaining one cycle sicker
c_S1_Trt1 <- c_S1 + 12000 # cost of remaining one cycle sick under treatment 1
c_S1_Trt2 <- c_S1 + 11350 # cost of remaining one cycle sick under treatment 2
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
u_S1_Trt2 <- u_S1 + 0.15  # utility when sick, treatment 2
u_D       <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
u_sex_cof <- -0.015       # utility sex coefficient, where 0 is female and 1 is male
ru_S1     <- -0.0015      # change in utility of individuals with every additional year being sick
ru_S2     <- -0.0020      # change in utility of individuals with every additional year being sicker

v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof, "sex" = u_sex_cof
)
v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "S1" = ru_S1
)

### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "D" = u_D)
### Payoffs - treatment 1
v_states_costs1     <- c("H" = c_H, "S1" = c_S1_Trt1, "D" = c_D)
v_states_utilities1 <- c("H" = u_H, "S1" = u_S1_Trt1, "D" = u_D)
### Payoffs - treatment 2
v_states_costs2     <- c("H" = c_H, "S1" = c_S1_Trt2, "D" = c_D)
v_states_utilities2 <- c("H" = u_H, "S1" = u_S1_Trt2, "D" = u_D)

#------------------------------------------------------------------------------#

### Running the simulation ###

# Run the simulation:
## For the not treated
res_no_Trt <- run_microSim(
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
  l_surv_params       = l_surv_params,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
## For being treated with treatment 1
res_Trt1 <- run_microSim(
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
  l_surv_params       = l_surv_params,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
## For being treated with treatment 2
res_Trt2 <- run_microSim(
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
  l_surv_params       = l_surv_params,
  cycle_length        = cycle_length,
  starting_seed       = seed
)

nb_no_Trt <- res_no_Trt$total_qalys * wtp - res_no_Trt$total_costs
nb_Trt1   <- res_Trt1$total_qalys   * wtp - res_Trt1$total_costs
nb_Trt2   <- res_Trt2$total_qalys   * wtp - res_Trt2$total_costs
which.max(c("no treatment" = nb_no_Trt, "treatment 1" = nb_Trt1, "treatment 2" = nb_Trt2))



