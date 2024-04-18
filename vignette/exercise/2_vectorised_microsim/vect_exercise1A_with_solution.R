## @knitr solution_probsV_function

# 0. Execute the lines below.
## clear R's session memory (Global Environment) 
rm(list = ls())

## General parameters
seed    <- 1234                          # random number generator state
num_sim <- 100                           # number of simulated individuals
v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
set.seed(seed)                           # set the seed to ensure reproducible samples below
v_occupied_state <- sample(              # sample current health state
  x       = v_states_names,              # from the four health states
  size    = num_sim,                     # sample a state for each of the simulated individuals
  replace = TRUE,                        # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)       # sample 75% as healthy, 20% as sick, and 5% as sicker
)
v_time_in_state  <- sample(              # sample time in current health state
  x       = 1:10,                        # for values between 1 and 10
  size    = num_sim,                     # sample time in state for each of the simulated individuals
  replace = TRUE,                        # allow sampled states to be re-sampled
  prob    = rep(                         # values 1 to 10 will have equal chances of being sampled
    x     = 1/length(1:10),              # the chance of sampling any the values is 1/the number of values
    times = length(1:10)                 # repeat the probabilities as many times as there as values
  )
)

## Transition probabilities (per cycle)
p_HD    <- 0.005                         # probability to die when healthy
p_HS1   <- 0.15                          # probability to become sick when healthy
p_S1H   <- 0.5                           # probability to become healthy when sick
p_S1S2  <- 0.105         	               # probability to become sicker when sick
rr_S1   <- 3                             # rate ratio of death in sick vs healthy
rr_S2   <- 10            	               # rate ratio of death in sicker vs healthy 
r_HD    <- -log(1 - p_HD)                # rate of death in healthy 
r_S1D   <- rr_S1 * r_HD                  # rate of death in sick
r_S2D   <- rr_S2 * r_HD  	               # rate of death in sicker
p_S1D   <- 1 - exp(- r_S1D)              # probability to die in sick
p_S2D   <- 1 - exp(- r_S2D)              # probability to die in sicker
rp_S1   <- 0.02                          # increase in mortality rate with every additional year being sick
rp_S2   <- 0.05                          # increase in mortality rate with every additional year being sicker

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
# 1. Subset healthy individual
v_time_in_state[v_occupied_state == "H"]
# 1.1. The code above identifies the time each simulated individual spent in the
# "H" state. Copy and the amend the code above to subset and print the time 
# spent by individuals in S1 and S2 health states.
 
v_time_in_state[v_occupied_state == "S1"]
v_time_in_state[v_occupied_state == "S2"]




# 1.2. Write a line of code to identify long did the third individual occupying
# the "S1" state remained in it?
# HINT: subset the time in state "S1" for (1.1) using the square brackets[]

v_time_in_state[v_occupied_state == "S1"][3]




# 1.3. Change the time spent in "S1" for the third simulated individual to 4.
# HINT: use the code from task (1.2) in addition to '<-' to change the value

v_time_in_state[v_occupied_state == "S1"][3] <- 4




# 1.4. Write a single line of code to double the time spent in "H", printing the 
# results to the console. 
# - DO NOT USE A LOOP. 
# - compare the scaled values to the original ones; print time spent in "H".
# HINT: run the code '4 * c(1, 3)'.
# HINT: R multiplies each element in c(1, 3) by 4 resulting in c(4, 12).
# HINT: multiply the time spent for everyone in "H" by 2.

v_time_in_state[v_occupied_state == "H"] * 2
v_time_in_state[v_occupied_state == "H"]




# 2. Compute probabilities based on time spent in each state.
r_S1D <-  - log(1 - p_S1D)
r_S2D <-  - log(1 - p_S2D)
1 - exp(- r_S1D * (1 + v_time_in_state[v_occupied_state == "S1"] * rp_S1))
# 2.1. Print the baseline 'p_S1D'

p_S1D
# '0.01492512'




# 2.2. The line below computes the probability of dying for those in "S1" based
# on the time each individual remained in the "S1" state.
# '1 - exp(- r_S1D * (1 + v_time_in_state[v_occupied_state == "S1"] * rp_S1))'
# - Write a line of code to assign the estimated values to 'p_S1D' for those in
# "S1" health state.
# HINT: subset 'p_S1D' values to be replaced for only those in "S1", see (1.1) 

p_S1D[v_occupied_state == "S1"] <- 1 - exp(- r_S1D * (1 + v_time_in_state[v_occupied_state == "S1"] * rp_S1))




# 2.3. Print 'p_S1D' and 'v_occupied_state == "S1"'.
# - See how the number of element is different from that in task (2.1)
# - See how the values in 'p_S1D' relate to 'v_occupied_state == "S1"'
# HINT: print both 'p_S1D' and 'v_occupied_state == "S1"'

p_S1D
print(v_occupied_state == "S1")
# All values where 'v_occupied_state == "S1"' is 'TRUE', except the first one is 
# set to some probability value. 
# The first value remains in place because 'v_occupied_state == "S1"' was 'FALSE'
# for the first element and therefore R left it unchanged. The first value
# will still be ignored by the remaining code in the 'update_prob()' function.




# 2.4. Print the 'p_S2D'

p_S2D
# '0.04888987'




# 2.5. Write a line of code to compute the probability of dying for those in
# "S2" based on the time each individual remained in "S2".
# HINT: the required line of code is similar to the output from task (2.2)

p_S2D[v_occupied_state == "S2"] <- 1 - exp(- r_S2D * (1 + v_time_in_state[v_occupied_state == "S2"] * rp_S2))




# 2.6. Print 'p_S2D' and 'v_occupied_state == "S1"'.
# - See how the number of element is different from that in task (2.1)
# - See how the values in 'p_S2D' relate to 'v_occupied_state == "S1"'
# HINT: print both 'p_S2D' and 'v_occupied_state == "S1"'

p_S2D
print(v_occupied_state == "S2")
# The only value where 'v_occupied_state == "S2"' is 'TRUE', except the first one
# is set to some probability value. 
# The first value remains in place because 'v_occupied_state == "S2"' was 'FALSE'
# for the first element and therefore R left it unchanged. The first value
# will still be ignored by the remaining code in the 'update_prob()' function.




# 3. Execute the code below to create a state transition probabilities matrix
m_probs <- matrix(
  nrow = length(v_time_in_state), # a row for each individual
  ncol = length(v_states_names),  # a column for each state
  dimnames = list(
    v_occupied_state,             # name each row based on the occupied state
    v_states_names                # give each column one of the states names
  )
)
# 3.1. Use the 'head()' function to print part of the 'm_probs' matrix

head(m_probs)




# 3.2. Run the code below to populate the matrix at the indices of the healthy
# simulated individuals.
# HINT: The indices of the simulated individuals in any of the four state can be
# identified by using 'v_occupied_state == "H"'  

m_probs[v_occupied_state == "H",]  <- c( 1 - p_HS1 - p_HD, p_HS1, 0, p_HD)




# 3.3 Use the 'head()' function twice, once to print part of the 'm_probs' 
# matrix, and once to print part of the 'm_probs' matrix for those in "H" state
# HINT: 'm_probs' is a matrix; use the square brackets with a comma to indicate
# the rows and columns. For example, matrix_name[row, column]
# HINT: since we want all columns, use matrix_name[row, ]
# HINT: use the logical expression (the one with '==') to subset targeted rows

head(m_probs)
head(m_probs[v_occupied_state == "H",])




# 3.4. The transitions from "S1" depends on 'p_S1H', 'p_S1S2', and 'p_S1D'.
# However, 'p_S1D' is not of equal length to the others. Execute the line below
# to construct a matrix where all probabilities are set correctly. Investigate
# the results.

cbind(p_S1H, 1 - p_S1S2 - p_S1H - p_S1D[v_occupied_state == "S1"], p_S1S2, p_S1D[v_occupied_state == "S1"])




# 3.5. Similar to the code in task (3.2), write three lines of code to populate
# the 'm_probs' for individuals in "S1", "S2" and "D" states.
# - "S1" probabilities: cbind(p_S1H, 1 - p_S1S2 - p_S1H - p_S1D[v_occupied_state == "S1"], p_S1S2, p_S1D[v_occupied_state == "S1"])
# - "S2" probabilities: cbind(0, 0, 1 - p_S2D[v_occupied_state == "S2"], p_S2D[v_occupied_state == "S2"])
# - "D"  probabilities: c(0, 0, 0, 1) 

m_probs[v_occupied_state == "S1",] <- cbind(p_S1H, 1 - p_S1S2 - p_S1H - p_S1D[v_occupied_state == "S1"], p_S1S2, p_S1D[v_occupied_state == "S1"])
m_probs[v_occupied_state == "S2",] <- cbind(0, 0, 1 - p_S2D[v_occupied_state == "S2"], p_S2D[v_occupied_state == "S2"])                          
m_probs[v_occupied_state == "D",]  <- c(0, 0, 0, 1 )                                        




# 3.5. Print the 'm_probs' function and investigate the results.

print(m_probs)




# 4. Review the 'update_probsV' function below. 
update_probsV <- function(
    v_states_names,
    v_occupied_state,
    l_trans_probs, 
    v_time_in_state) { 
  
  with(
    data = v_occupied_state,
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

      # sanity check
      ifelse(
        test = rowSums(m_probs) == 1,               # check if the transition probabilities add up to 1
        yes = return(m_probs),                      # return the transition probabilities
        no = print("Probabilities do not sum to 1") # or produce an error
      )   
    }
  )
}
# 4.1. What are arguments (inputs) of this function?

# 'v_states_names', 'v_occupied_state', 'l_trans_probs', and 'v_time_in_state'




# 4.2. Is this function suitable to extract the transition probabilities from 
# all four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "S2"
# and "D" health state.




# 4.3. In addition to the missing code for the fourth health state, the wrong
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
# wrong value, 'data = v_occupied_state'. Most of the data structures referred to
# in the code within the 'with()' function are stored in the list 
# 'l_trans_probs'.




# 4.4. Create a copy the 'update_probsV' function below, fix the mistake from
# task (4.3) and add the missing line from task (4.2).
# HINT: the missing code relates to the transitions from the dead "D" state.
# HINT: "D" state is an absorbing state, hence, all transitions from it are 0.
# test the function, are there any error messages?

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



