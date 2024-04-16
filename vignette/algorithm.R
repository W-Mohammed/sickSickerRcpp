## @knitr microSim_algorithm

# create matrices to capture the simulation of 1,000 individuals over 30 cycles
m_States <- matrix(nrow = 1000, ncol = 30 + 1) # to store state name
m_Costs  <- matrix(nrow = 1000, ncol = 30 + 1) # to store state-related calc_costs
m_Effs   <- matrix(nrow = 1000, ncol = 30 + 1) # to store state-related health outcomes (QALYs) 

for (i in 1:1000) {                            # for each 'i' of the 1,000 simulated individual:
  
  # Step 1:
  m_States[i, 1] <- v_starting_states[i]       # indicate the initial health state
  m_Costs[i, 1]  <- calc_costs(m_States[i, 1]) # calculate the costs incurred in their starting health state
  m_Effs[i, 1]   <- calc_effs(m_States[i, 1])  # calculate the QALYs accrued in their starting health state
  
  for (t in 1:30) {                               # for each 't' of the 30 cycles:
    # Step 2:
    v_trans_probs <- update_porbs(m_States[i, t]) # update the transition probabilities at cycle 't'
    
    # Step 3:
    m_States[i, t + 1] <- sample(state_names, prob = v_trans_probs, size = 1)  # sample the health state at 't + 1' 
    
    # Step 4:
    m_Costs[i, t + 1]  <- calc_costs(m_States[i, t + 1]) # calculate the costs incurred in their 't + 1' health state
    m_Effs[i, t + 1]   <- calc_effs( m_States[i, t + 1]) # calculate the QALYs accrued in their 't + 1' health state
  } # close the loop for the cycles 't' 
  
} # close the loop for the individuals 'i'

# Step 5:
v_total_costs <- rowSums(m_Costs)           # calculate total cost per individual
v_total_qalys <- rowSums(m_Effs)            # calculate total QALYs per individual
total_costs   <- mean(v_total_costs)        # calculate average discounted cost 
total_qalys   <- mean(v_total_qalys)        # calculate average discounted QALYs