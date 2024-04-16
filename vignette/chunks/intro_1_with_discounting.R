## @knitr demo_microsim_disc

# Healthy-Sick-Dead microsimulation model demo

# clear R's session memory (Global Environment) 
rm(list = ls())

#------------------------------------------------------------------------------#

                        ### Defining model functions ###

# Define model functions
## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by individual 'i' at cycle 't'
update_probs <- function(
    occupied_state,
    m_trans_probs) { 
  
  v_probs <- rep(NA, num_states)   # create vector of state transition probabilities
  names(v_probs) <- v_states_names # name the vector
  
  # update v_probs with the appropriate probabilities   
  v_probs[occupied_state == "H"]  <- m_trans_probs["H",]  # transition probabilities when healthy
  v_probs[occupied_state == "S1"] <- m_trans_probs["S1",] # transition probabilities when sick
  v_probs[occupied_state == "D"]  <- m_trans_probs["D",]  # transition probabilities when dead   
  
  # sanity check
  ifelse(
    test = sum(v_probs) == 1,                   # check if the transition probabilities add up to 1
    yes = return(v_probs),                      # return the transition probabilities
    no = print("Probabilities do not sum to 1") # or produce an error
  )
}

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by individual 'i' at cycle 't'
calc_costs <- function (
    occupied_state,
    v_states_costs) {
  
  state_costs <- 0                                            # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]  # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] # update the cost if sick conditional
  
  return(state_costs)                                         # return the costs
}

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by individual 'i' at cycle 't' and
### the cycle_length (measured in years)
calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    cycle_length = 1) {

  state_utility <- 0                                                # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  # update the utility if healthy
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] # update the utility if sick
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}

## Run Microsimulation function
### This function runs the microsimulation function of the Healthy-Sick-Dead model
run_microSim <- function(
    v_starting_states, 
    num_sim, 
    num_cycles, 
    v_states_names,
    v_states_costs,
    v_states_utilities,
    m_trans_probs,
    discount_rate,
    cycle_length = 1,
    starting_seed = 1) {

  # calculate discount weights accounting for the number and length (in years) of cycles
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(nrow = num_sim, ncol = num_cycles + 1)  
  
  for (i in 1:num_sim) {                   # for each 'i' of the 'num_sim' simulated individual:
    set.seed(starting_seed + i)            # set the seed for every individual for the random number generator
    
    # Step 1:
    m_States[i, 1] <- v_starting_states[i] # indicate the initial health state
    m_Costs[i, 1]  <- calc_costs(          # calculate the costs incurred in their starting health state
      occupied_state = m_States[i, 1],
      v_states_costs = v_states_costs
    )
    m_Effs[i, 1]   <- calc_effs(           # calculate the QALYs accrued in their starting health state
      occupied_state = m_States[i, 1], 
      v_states_utilities = v_states_utilities,
      cycle_length = cycle_length
    )
    
    for (t in 1:num_cycles) {              # for each 't' of the 'num_cycles' cycles:
      # Step 2:
      v_trans_probs <- update_probs(       # update the transition probabilities at cycle 't'
        occupied_state = m_States[i, t],
        m_trans_probs = m_trans_probs
      ) 
      
      # Step 3:
      m_States[i, t + 1] <- sample(        # sample the health state at 't + 1' 
        x = v_states_names, 
        prob = v_trans_probs, 
        size = 1
      )  
      
      # Step 4:
      m_Costs[i, t + 1]  <- calc_costs(    # calculate the costs incurred in their 't + 1' health state
        occupied_state = m_States[i, t + 1],
        v_states_costs = v_states_costs
      )
      m_Effs[i, t + 1]   <- calc_effs(     # calculate the QALYs accrued in their 't + 1' health state
        occupied_state = m_States[i, t + 1], 
        v_states_utilities = v_states_utilities,
        cycle_length = cycle_length
      )
      
    } # close the loop for the cycles 't' 
    
    if(i/100 == round(i/100,0)) {          # display the progress of the simulation
      cat('\r', paste(i/num_sim * 100, "% done", sep = " "))
    }
    
  } # close the loop for the individuals 'i' 
  
  # Step 5:
  v_total_costs <- rowSums(m_Costs)             # calculate total costs per individual
  v_total_qalys <- rowSums(m_Effs)              # calculate total QALYs per individual
  total_costs   <- mean(v_total_costs)          # calculate average costs
  total_qalys   <- mean(v_total_qalys)          # calculate average QALYs
  v_total_D_costs <- m_Costs %*% v_discount_wts # calculate total discounted costs per individual
  v_total_D_qalys <- m_Effs%*% v_discount_wts   # calculate total discounted QALYs per individual
  total_D_costs   <- mean(v_total_D_costs)      # calculate average discounted costs 
  total_D_qalys   <- mean(v_total_D_qalys)      # calculate average discounted QALYs
  
  # store the results in a list:
  results <- list(
    m_States = m_States, 
    m_Costs = m_Costs, 
    m_Effs = m_Effs, 
    v_total_costs = v_total_costs, 
    v_total_qalys = v_total_qalys, 
    total_costs = total_costs, 
    total_qalys = total_qalys, 
    v_total_D_costs = v_total_D_costs, 
    v_total_D_qalys = v_total_D_qalys, 
    total_D_costs = total_D_costs, 
    total_D_qalys = total_D_qalys
  )
  
  # return the results
  return(results)
}

#------------------------------------------------------------------------------#

                        ### Defining model parameters ###

# Define model inputs
## General parameters
num_sim <- 100000                     # number of simulated individuals
num_cycles <- 30                      # time horizon
discount_rate <- 0.03                 # equal discounting of costs and health outcomes

## Health states
v_states_names <- c("H","S1","D")      # the model states: Healthy (H), Sick (S1), Dead (D)
num_states <- length(v_states_names)   # the number of health states
v_starting_states <- rep("H", num_sim) # everyone begins in the healthy state 

## Transition probabilities (per cycle)
p_HD  <- 0.005                        # probability to die when healthy
p_HS1 <- 0.15                         # probability to become sick when healthy
p_S1H <- 0.5                          # probability to become healthy when sick
rr_S1 <- 3                            # rate ratio of death in sick vs healthy
r_HD  <- -log(1 - p_HD)               # rate of death in healthy 
r_S1D <- rr_S1 * r_HD                 # rate of death in sick
p_S1D <- 1 - exp(- r_S1D)             # probability to die in sick

m_trans_probs <- matrix(                        # create a transition probability matrix
  data = c(
    1 - p_HS1 - p_HD,             p_HS1,  p_HD, # transition probabilities when healthy
               p_S1H, 1 - p_S1H - p_S1D, p_S1D, # transition probabilities when sick
                   0,                 0,     1  # transition probabilities when dead
  ),
  nrow = num_states,
  byrow = TRUE,
  dimnames = list(v_states_names, v_states_names)
)

## Cost and utility inputs 
c_H <- 2000                           # cost of remaining one cycle healthy
c_S1 <- 4000                          # cost of remaining one cycle sick 
u_H <- 1                              # utility when healthy 
u_S1 <- 0.75                          # utility when sick 

v_states_costs     <- c("H" = c_H, "S1" = c_S1, "D" = 0) # named costs vector
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "D" = 0) # named utilities vector

#------------------------------------------------------------------------------#

                            ### Running the simulation ###

# Run the simulation:
microsim_results <- run_microSim(
  v_starting_states = v_starting_states,
  num_sim = num_sim,
  num_cycles = num_cycles,
  v_states_names = v_states_names,
  v_states_costs = v_states_costs,
  v_states_utilities = v_states_utilities,
  m_trans_probs = m_trans_probs,
  discount_rate = 0.03,
  cycle_length = 1,
  starting_seed = 1
)

# View the results:
str(microsim_results)

microsim_results$v_total_D_costs[1:10]
microsim_results$total_D_costs

microsim_results$v_total_D_qalys[1:10]
microsim_results$total_D_qalys