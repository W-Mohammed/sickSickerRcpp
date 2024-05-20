## @knitr solution_effects_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

## Utility inputs 
u_H  <- 1                 # utility when healthy 
u_S1 <- 0.75              # utility when sick 
u_S2 <- 0.5               # utility when sicker
u_D  <- 0                 # utility when dead

v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
occupied_state <- "H"                    # current health state
state_utility <- 0                       # by default the utility for everyone is zero

# 1. The code line below is intended to create a vector 'v_states_utilities'.
# It is a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_states_utilities'.
# HINT: the first missing element is 'u_H'.

v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)




# 2. Print the 'v_states_utilities' object on the console. Is the data in the 
# printed object as expected? i.e. are the values associated with the correct 
# names?

v_states_utilities
# Yes!




# 3. Complete the code below to extract the utility associated with being 
# healthy from the 'v_states_utilities' vector
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

state_utility[occupied_state == "H"]  <- v_states_utilities["H"]




# 4. Complete the code below to retrieve the utility if a simulated individual
# was sick (occupying the sick state) from the the 'v_states_utilities' vector
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] 




# 5. Based on tasks (3) and (4), add two lines of code to extract the utilities 
# associated with being sicker and dead from the 'v_states_utilities'
# HINT: Review and compare code in tasks (3) and (4)

state_utility[occupied_state == "S2"] <- v_states_utilities["S2"]
state_utility[occupied_state == "D"]  <- v_states_utilities["D"] 




# 6. Review the 'calc_effs' function below. 
calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    cycle_length = 1) {
  
  state_utility <- 0                                                # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  # update the utility if healthy
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] # update the utility if sicker
  state_utility[occupied_state == "D"]  <- v_states_utilities["D"]  # update the utility if dead
  
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}
# 6.1. What are arguments (inputs) of this function?

# 'occupied_state', 'v_states_utilities' and 'cycle_length'




# 6.2. Is this function suitable to compute the health outcomes (QALYs) 
# associated with each of the four health states in the healthy-sick-sicker-dead
# model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "S1"
# health state.




# 6.3. Create a copy the 'calc_effs' function below and add the missing line.
# HINT: The required code was one of the lines from the outputs of task (5).

calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    cycle_length = 1) {
  
  state_utility <- 0                                                # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  # update the utility if healthy
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] # update the utility if sick
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] # update the utility if sicker
  state_utility[occupied_state == "D"]  <- v_states_utilities["D"]  # update the utility if dead
  
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}




