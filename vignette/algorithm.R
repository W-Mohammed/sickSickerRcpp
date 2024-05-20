## @knitr microSim_diagram

DiagrammeR::grViz("
digraph flowchart {
  node [fontname = Helvetica, fontsize=20, shape = box]
  edge [fontname = Helvetica]

  Start [label = 'Start simulation', shape = ellipse, style=filled, fillcolor='palegreen']
  ForIndividuals [label = 'For each simulated individual i:', shape = diamond, style=filled, fillcolor='grey', fontsize=24, fontname='Helvetica-Bold']
  SetInitialState [label = 'Identify the starting health state', style=filled, fillcolor='orange']
  CalcInitialOutcomes [label = 'Calculate the costs (QALYs) incurred (accrued) in the initial health state', style=filled, fillcolor='orange']
  ForCycles [label = 'In each cycle t:', shape = diamond, style=filled, fillcolor='skyblue', fontsize=24, fontname='Helvetica-Bold']
  UpdateProbs [label = 'Update the transition probabilities based on the current health state', style=filled, fillcolor='skyblue']
  SampleState [label = 'Identify the health state to be occupied in the next cycle', style=filled, fillcolor='skyblue']
  CalcOutcomes [label = 'Calculate the costs (QALYs) to be incurred (accrued) in the next cycle', style=filled, fillcolor='skyblue']
  CheckEndCycles [label = 'Was this the last cycle?', shape = diamond, style=filled, fillcolor='skyblue', fontsize=24, fontname='Helvetica-Bold']
  CheckEndIndividuals [label = 'Was this the last individual?', shape = diamond, style=filled, fillcolor='grey', fontsize=24, fontname='Helvetica-Bold']
  CalcTotOutcomes [label = 'Calculate the total costs (QALYs) incurred (accrued) per individual']
  End [label = 'End simulation', shape = ellipse, style=filled, fillcolor='red']

  Start -> ForIndividuals -> SetInitialState -> CalcInitialOutcomes
  CalcInitialOutcomes -> ForCycles -> UpdateProbs -> SampleState -> CalcOutcomes -> CheckEndCycles
  CheckEndCycles -> ForCycles [label = 'No\nprocess next cycle', fontcolor='skyblue', fontsize=18]
  CheckEndCycles -> CheckEndIndividuals [label = 'Yes']
  CheckEndIndividuals -> ForIndividuals [label = 'No\nprocess next individual', fontcolor='orange', fontsize=18]
  CheckEndIndividuals -> CalcTotOutcomes [label = 'Yes']
  CalcTotOutcomes -> End
}
")

## @knitr microSim_algorithm

# create matrices to capture the simulation of num_i individuals over num_cycles cycles
m_States <- matrix(nrow = num_i, ncol = num_cycles + 1) # to store state name
m_Costs  <- matrix(nrow = num_i, ncol = num_cycles + 1) # to store state-related calc_costs
m_Effs   <- matrix(nrow = num_i, ncol = num_cycles + 1) # to store state-related health outcomes (QALYs) 

for (i in 1:num_i) {                              # for each 'i' of the num_i simulated individual:
  
  # Step 1:
  m_States[i, 1] <- v_starting_states[i]          # indicate the initial health state
  m_Costs[i, 1]  <- calc_costs(m_States[i, 1])    # calculate the costs incurred in the starting health state
  m_Effs[i, 1]   <- calc_effs(m_States[i, 1])     # calculate the QALYs accrued in the starting health state
  
  for (t in 1:num_cycles) {                       # for each 't' of the num_cycles cycles:
    # Step 2:
    v_trans_probs <- update_porbs(m_States[i, t]) # update the transition probabilities at 't'
    
    # Step 3:
    m_States[i, t + 1] <- sample(state_names, prob = v_trans_probs, size = 1)  # sample the state at 't + 1' 
    
    # Step 4:
    m_Costs[i, t + 1]  <- calc_costs(m_States[i, t + 1]) # calculate the costs incurred at 't + 1'
    m_Effs[i, t + 1]   <- calc_effs( m_States[i, t + 1]) # calculate the QALYs accrued at 't + 1'
  } # close the loop for the cycles 't' 
  
} # close the loop for the individuals 'i'

# Step 5:
v_total_costs <- rowSums(m_Costs)           # calculate total cost per individual
v_total_qalys <- rowSums(m_Effs)            # calculate total QALYs per individual
mean_costs    <- mean(v_total_costs)        # calculate average discounted cost 
mean_qalys    <- mean(v_total_qalys)        # calculate average discounted QALYs
