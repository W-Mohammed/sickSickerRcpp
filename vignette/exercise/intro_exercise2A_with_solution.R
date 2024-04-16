## @knitr solution_probs2_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

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
rp_S1  <- 0.02                           # increase in mortality rate with every additional year being sick
rp_S2  <- 0.05                           # increase in mortality rate with every additional year being sicker

v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
occupied_state <- "S2"                   # current health state
time_in_state  <- 3                      # time in current health state

# 1. Replace the "..." below to pack the transition probabilities and rates of 
# increase in mortality in the list 'l_trans_probs'.
# HINT: the replacements are objects and/or their names, hence, use quotation marks
# "" carefully. Only names of objects need to be wrapped with "".
# HINT: for example, the first "..." should be replaced with 'p_HD' not '"p_HD"'
# HINT: 'rr_S1' is used to compute 'r_HD' then 'r_S1D' then 'p_S1D'; hence, only
# 'p_S1D' is packed into 'l_trans_probs'.
# HINT: pack all probabilities and the rates that are yet to be used
 
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




# 2. Print the list to the console. 

l_trans_probs




# 3. Are all values numeric? i.e. none of the numbers is surrounded by quotation
# marks? 

# Yes! the probabilities are packed as numeric.




# 4. Update the probabilities of death for sick and sicker states.
r_S1D <-  - log(1 - p_S1D)
# 4.1. In the line above, the probability of death for those in the sick (S1)
# state is converted into the rate of death for those in the sick state. Update
# the line of code blow to compute the same for those in the sicker (S2) state.

r_S2D <-  - log(1 - p_S2D)




time_in_state[occupied_state == "S1"]
# 4.2. The line of code above checks the time spent in the current state if and 
# only if the current state was "S1". 
# 4.2.1 Write a line of code to check the time spent in the current state if the
# current state was "S2". Run the code and notice how the output is different
# from that of "S1" above. 

time_in_state[occupied_state == "S2"]




# 4.2.2. Since 'occupied_state' = "S2", running the code 'occupied_state == "S1"'
# in the console returns 'FALSE'. What does R return when the value to subset is
# 'FALSE'?
# HINT: execute 'time_in_state[occupied_state == "S1"]'
 
# A numeric object of length 0, i.e. an object that contains no data.




p_S1D[occupied_state == "S1"] <- 1 - exp(- r_S1D * (1 + time_in_state[occupied_state == "S1"] * rp_S1))
# 4.3. In the line above the probability of death for those in the sick (S1) 
# state is updated based on how long an individual spent in that health state.
# 4.3.1. Execute the portion 'p_S1D[occupied_state == "S1"]' in the console.
# What does it mean to assign the outputs of the code to the right of the '<-'
# (assignment operator) to numeric(0)? i.e. will anything be saved in 'p_S1D'?
# HINT: print the object 'p_S1D', then print 'p_S1D[occupied_state == "S1"]'.
# HINT: the '[occupied_state == "S1"]' refers to no where in the object 'p_S1D'.

# Nothing will be saved in 'p_S1D'. This code is intended to only update, or
# override, the default value of 'p_S1D' if the current health state was (S1).

 


# 4.3.2. Printing 'p_S2D[occupied_state == "S2"]' to console shows the default
# value of 'p_S2D'. Why?
# HINT: execute 'occupied_state == "S2"' (without the '')
# HINT: what is the difference between 'p_S2D[TRUE]' and 'p_S2D'?

# The object 'occupied_state', which stores the name of the current state has
# the value '"S2"'. Therefore, 'occupied_state == "S2"' evaluates to 'TRUE' 
# and 'p_S2D[TRUE]' is equivalent to 'p_S2D'.




# 4.3.3. Write a line of code to update the probability of death for those in 
# the sicker (S2) state based on how long an individual spent in that health 
# state.
# HINT: copy then amend the line of code provided for (right above) task 4.3.

p_S2D[occupied_state == "S2"] <- 1 - exp(- r_S2D * (1 + time_in_state[occupied_state == "S2"] * rp_S2))




# 5. Review the 'update_probs' function below. 
update_probs <- function(
    v_states_names,
    occupied_state,
    l_trans_probs, 
    time_in_state) { 
  
  with(
    data = time_in_state,
    expr = {
      
      # update probabilities of death after first converting them to rates and applying the rate ratio
      r_S1D <-  - log(1 - p_S1D)
      r_S2D <-  - log(1 - p_S2D)
      # calculate p_S1D/p_S2D conditional on current state and duration of being sick/sicker
      p_S1D[occupied_state == "S1"] <- 1 - exp(- r_S1D * (1 + time_in_state[occupied_state == "S1"] * rp_S1))
      p_S2D[occupied_state == "S2"] <- 1 - exp(- r_S2D * (1 + time_in_state[occupied_state == "S2"] * rp_S2))
      
      # create vector of state transition probabilities
      v_probs        <- rep(NA, length(v_states_names))
      names(v_probs) <- v_states_names
      
      # update v_probs with the appropriate probabilities   
      v_probs[occupied_state == "H"]  <- c( 1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                # transition probabilities when healthy
      v_probs[occupied_state == "S1"] <- c(p_S1H, 1 - p_S1S2 - p_S1H - p_S1D, p_S1S2, p_S1D) # transition probabilities when sick
      v_probs[occupied_state == "S2"] <- c(0, 0, 1 - p_S2D, p_S2D)                           # transition probabilities when sicker
      v_probs[occupied_state == "D"]  <- c(0, 0, 0, 1 )                                      # transition probabilities when dead   
      
      # sanity check
      ifelse(
        test = sum(v_probs) == 1,                   # check if the transition probabilities add up to 1
        yes = return(v_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )   
    }
  )
}
# 5.1. What are arguments (inputs) of this function?

# 'v_states_names', 'occupied_state', 'l_trans_probs', and 'time_in_state'




# 5.2. Is this function suitable to extract the transition probabilities from 
# all four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "D"
# health state.




# 5.3. In addition to the missing code for the fourth health state, the wrong
# object was used in the 'with()' function. Identify the mistake and suggest a
# fix.
# HINT: you can learn more about the 'with()' function by running 'help("with")'
# HINT: the 'with()' function makes the code (or expression) called within it 
# aware of the data objects passed by the user.
# HINT: what parameters (data objects) do the lines of code within the 'with()'
# function refer to?
# HINT: where are those parameters (data) stored?
# HINT: what function argument (input) is yet to be referred to in the function?
# HINT: run 'rm(p_HD, p_HS1)' and test the function, what is the error message?

# The mistake was setting the argument "data" of the 'with()' function to the
# wrong value, 'data = time_in_state'. Most of the data structures referred to
# in the code within the 'with()' function are stored in the list 
# 'l_trans_probs'.




# 5.4. Create a copy the 'update_probs' function below, fix the mistake from
# task (5.3) and add the missing line from task (5.2).
# HINT: the missing code relates to the transitions from the dead "D" state.
# HINT: "D" state is an absorbing state, hence, all transitions from it are 0.
# test the function, are there any error messages?

update_probs <- function(
    v_states_names,
    occupied_state,
    l_trans_probs, 
    time_in_state) { 
  
  with(
    data = l_trans_probs,
    expr = {
      
      # update probabilities of death after first converting them to rates and applying the rate ratio
      r_S1D <-  - log(1 - p_S1D)
      r_S2D <-  - log(1 - p_S2D)
      # calculate p_S1D/p_S2D conditional on current state and duration of being sick/sicker
      p_S1D[occupied_state == "S1"] <- 1 - exp(- r_S1D * (1 + time_in_state[occupied_state == "S1"] * rp_S1))
      p_S2D[occupied_state == "S2"] <- 1 - exp(- r_S2D * (1 + time_in_state[occupied_state == "S2"] * rp_S2))
      
      # create vector of state transition probabilities
      v_probs        <- rep(NA, length(v_states_names))
      names(v_probs) <- v_states_names
      
      # update v_probs with the appropriate probabilities   
      v_probs[occupied_state == "H"]  <- c( 1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                # transition probabilities when healthy
      v_probs[occupied_state == "S1"] <- c(p_S1H, 1 - p_S1S2 - p_S1H - p_S1D, p_S1S2, p_S1D) # transition probabilities when sick
      v_probs[occupied_state == "S2"] <- c(0, 0, 1 - p_S2D, p_S2D)                           # transition probabilities when sicker
      v_probs[occupied_state == "D"]  <- c(0, 0, 0, 1 )                                      # transition probabilities when dead   
      
      # sanity check
      ifelse(
        test = sum(v_probs) == 1,                   # check if the transition probabilities add up to 1
        yes = return(v_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )   
    }
  )
}



