## @knitr solution_costs_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

## Cost inputs 
c_H  <- 2000         # cost of remaining one cycle healthy
c_S1 <- 4000         # cost of remaining one cycle sick
c_S2 <- 15000        # cost of remaining one cycle sicker
c_D  <- 0            # cost associated with being dead

v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
occupied_state <- "H"                    # current health state
state_costs <- 0                         # by default the cost for everyone is zero 

# 1. The code line below is intended to create a vector 'v_states_costs'. It is
# a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_states_costs'.
# HINT: the first missing element is 'c_H'.

v_states_costs <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)




# 2. Print the 'v_states_costs' object on the console. Is the data in the printed
# object as expected? i.e. are the values associated with the correct names?

v_states_costs
# Yes!




# 3. Complete the code below to extract the costs associated with being healthy
# from the 'v_states_costs' vector
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

state_costs[occupied_state == "H"]  <- v_states_costs["H"]




# 4. Complete the code below to retrieve the costs if a simulated individual was
# sick (occupying the sick state) from the 'v_states_costs' vector
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

state_costs[occupied_state == "S1"] <- v_states_costs["S1"] 




# 5. Based on tasks (3) and (4), add two lines of code to extract the costs 
# associated with being sicker and dead from the 'v_states_costs'
# HINT: Review and compare code in tasks (3) and (4)

state_costs[occupied_state == "S2"] <- v_states_costs["S2"]
state_costs[occupied_state == "D"]  <- v_states_costs["D"] 




# 6. Review the 'calc_costs' function below. 
calc_costs <- function (
    occupied_state,
    v_states_costs) {
  
  state_costs <- 0                                            # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]  # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] # update the cost if sick
  state_costs[occupied_state == "D"]  <- v_states_costs["D"]  # update the cost if dead
  
  return(state_costs)                                         # return the costs
}
# 6.1. What are arguments (inputs) of this function?

# 'occupied_state' and 'v_states_costs'




# 6.2. Is this function suitable to extract the costs associated with each of 
# the four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "S2"
# health state.




# 6.3. Create a copy the 'calc_costs' function below and add the missing line.
# HINT: The required code was one of the lines from the outputs of task (5).

calc_costs <- function (
    occupied_state,
    v_states_costs) {
  
  state_costs <- 0                                            # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]  # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] # update the cost if sick
  state_costs[occupied_state == "S2"] <- v_states_costs["S2"] # update the cost if sicker
  state_costs[occupied_state == "D"]  <- v_states_costs["D"]  # update the cost if dead
  
  return(state_costs)                                         # return the costs
}




