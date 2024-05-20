## @knitr exercise_probs_function

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

v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
num_states     <- length(v_states_names) # the number of health states
v_probs <- rep(NA, num_states)           # create vector of state transition probabilities
occupied_state <- "H"                    # current health state

# 1. Complete the transition matrix 'm_trans_probs'. Replace "..." with the 
# appropriate transition probabilities from the ones defined above.
# HINT: the four lines represents the probability of moving from healthy, sick,
# sicker and dead states, respectively. For instance, the probabilities of,
# moving from the dead state are [0, 0, 0, 1] indicating that individuals remain
# in the dead state. 
# HINT: the probability missing from the first line is 'p_HD'
 
m_trans_probs <- matrix(                              # create a transition probability matrix
  data = c(
    1 - p_HS1 - p_HD, p_HS1, 0, "...",                # transition probabilities when healthy
    "...", 1 - p_S1S2 - p_S1H - p_S1D, p_S1S2, p_S1D, # transition probabilities when sick
    0, 0, 1 - p_S2D, "...",                           # transition probabilities when sicker
    0, 0, 0, 1                                        # transition probabilities when dead
  ),
  nrow = num_states,
  byrow = TRUE,
  dimnames = list(v_states_names, v_states_names)
)




# 2. Print the updated matrix to the console. 




# 3. Are the probabilities in the correct order? 




# 4. Check if the probabilities of 
# moving from each health state add up to 1?
# HINT: Each row in 'm_trans_probs' gives the probabilities of moving from a
# health state.
# HINT: R provides sum functions for columns 'colSums()' and rows 'rowSums()'.
# HINT: To calculate the column sums we would run 'colSums(m_trans_probs)'.
# However, we want row sums.
# HINT: We need to use the equality operator '==' to check if two sets of values
# (the row sums and 1) are equal




# 5. Complete the code below to extract the transition probabilities when 
# healthy from the 'm_trans_probs' matrix
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

v_probs[occupied_state == "H"] <- m_trans_probs["...",]




# 6. Complete the code below to retrieve the transition probabilities if a 
# simulated individual was sick (occupying the sick state) from the 
# 'm_trans_probs' matrix
# HINT: Replace "..." with the name of the health state
# HINT: Print the 'v_states_names' object to see the states names

v_probs[occupied_state == "..."] <- m_trans_probs["S1",]




# 7. What would probabilities would are retrieved from the 'm_trans_probs' 
# matrix using the command 'm_trans_probs[,"S1"]'
# HINT, Notice the change in the position of the comma inside the square-brackets
# HINT, What do the columns represent in the transition matrix 'm_trans_probs'




# 8. Based on tasks (5) and (6), add two lines of code to extract the transition 
# probabilities when sicker, and dead from the 'm_trans_probs'
# HINT: Compare (5) and (6)




# 9. Review the 'update_probs' function below. 
update_probs <- function(
    occupied_state,
    m_trans_probs) { 
  
  v_probs <- rep(NA, nrow(m_trans_probs))                 # create state transition probabilities vector
  
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
# 9.1. What are arguments (inputs) of this function?



# 9.2. Is this function suitable to extract the transition probabilities from 
# all four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.



# 9.3. Create a copy the 'update_probs' function below and add the missing line.
# HINT: The required code was one of the lines from the outputs of task (8).



