## @knitr exercise_run_Trts_microSim

# Healthy-Sick-Sicker-Dead microsimulation model

# clear R's session memory (Global Environment) 
rm(list = ls())

# 1. Execute the lines between the #---# below.
#------------------------------------------------------------------------------#

                        ### Defining model functions ###

# Define model functions
## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by individual 'i' at cycle 't'
update_probs <- function(
    occupied_state,
    m_trans_probs) { 
  
  v_probs <- rep(NA, nrow(m_trans_probs))                 # create state transition probabilities vector
  
  # update v_probs with the appropriate probabilities   
  v_probs[occupied_state == "H"]  <- m_trans_probs["H",]  # transition probabilities when healthy
  v_probs[occupied_state == "S1"] <- m_trans_probs["S1",] # transition probabilities when sick
  v_probs[occupied_state == "S2"] <- m_trans_probs["S2",] # transition probabilities when sicker
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
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] # update the cost if sick
  state_costs[occupied_state == "S2"] <- v_states_costs["S2"] # update the cost if sicker
  
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
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] # update the utility if sicker

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
    cycle_length = 1,
    starting_seed = 1) {
  
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
num_sim             <- 100                # number of simulated individuals
num_cycles          <- 30                 # time horizon if each cycle is a year long
cycle_length        <- 1                  # length of cycle (in years)
seed                <- 1234               # random number generator state
wtp                 <- 30000              # Willingness to pay for each QALY ($)

## Health states
v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
num_states <- length(v_states_names)     # the number of health states
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

m_trans_probs <- matrix(                              # create a transition probability matrix
  data = c(
    1 - p_HS1 - p_HD, p_HS1, 0, p_HD,                 # transition probabilities when healthy
    p_S1H, 1 - p_S1S2 - p_S1H - p_S1D, p_S1S2, p_S1D, # transition probabilities when sick
    0, 0, 1 - p_S2D, p_S2D,                           # transition probabilities when sicker
    0, 0, 0, 1                                        # transition probabilities when dead
  ),
  nrow = num_states,
  byrow = TRUE,
  dimnames = list(v_states_names, v_states_names)
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

u_H       <- 1            # utility when healthy 
u_S1      <- 0.75         # utility when sick 
u_S2      <- 0.5          # utility when sicker
u_S1_Trt1 <- u_S1 + 0.2   # utility when sick under treatment 1
u_S2_Trt1 <- u_S2         # utility when sicker under treatment 1
u_S1_Trt2 <- u_S1 + 0.15  # utility when sick under treatment 2
u_S2_Trt2 <- u_S2 + 0.05  # utility when sicker under treatment 2
u_D       <- 0            # utility when dead

### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)

#------------------------------------------------------------------------------#

### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)
# 2. Considering how the payoffs of the "no treatment" option defined in the 2
# lines above, define the payoffs for "treatment 1" and "treatment 2".




# 3. Print all payoffs vectors to the console. Are the states rewards defined
# correctly? i.e. are the payoffs set with the correct states names?
# HINT: review the "Cost and utility inputs" section above if in doubt




# 4. Print the transition matrix. Are the probabilities set out correctly?
# HINT: check if the row sums (using the function 'rowSums()') are equal to 1.




# 5. What are the arguments (parameters) of the 'run_microSim()' function?
# HINT: Check the definition of the 'run_microSim()' function above.
# HINT: The function 'formals()' can be used to identify functions' arguments.
# Run '?formals' in the console to access the 'formals()' function help file.
# Run 'formals(run_microSim)' in the console for R to list the argumnets.




## For the no treatment option
res_no_Trt <- run_microSim(
  v_starting_states   = v_starting_states,
  num_sim             = num_sim,
  num_cycles          = num_cycles,
  v_states_names      = v_states_names,
  v_states_costs      = v_states_costs,
  v_states_utilities  = v_states_utilities,
  m_trans_probs       = m_trans_probs,
  cycle_length        = cycle_length,
  starting_seed       = seed
)
# 6. The code above runs the simulation for the "no treatment" option. 
# 6.1. Run the simulation for the two treatments ("treatment 1" and "treatment 2"),
# and assign the outputs to objects named "res_Trt1" and "res_Trt2", 
# respectively.
# HINT: set the correct payoffs before calling the 'run_microSim()' function
# for each of the treatments.




# 6.2. Print to the console the following data:
# - the outcomes (average costs/QALYs) associated with the "no treatment" choice
# - the average costs associated with "treatment 1", and
# - the average QALYs associated with "treatment 2"
# HINT: use the 'names()' function to retreive the names of the objects in any 
# of the results "res_" objects. E.g. 'names(res_Trt2)'
# HINT: use the 'View()' function to view the contents of the microsimulations
# results' objects.




# 6.3. Compute the net monetary benefits (NMB) at $30,000 per QALY. Which of the
# three choices (no treatment, treatment 1 and treatment 2) is the 
# cost-effective one.
# HINT: the NMB = (QALYs * 30000) - Costs



