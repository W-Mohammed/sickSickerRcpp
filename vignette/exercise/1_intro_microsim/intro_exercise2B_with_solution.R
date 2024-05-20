## @knitr solution_costs2_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

## Population characteristics/features
seed                <- 1234                # random number generator state
num_sim             <- 100                 # number of simulated individuals
mean_age            <- 50                  # mean age in the simulated population
sd_age              <- 3                   # standard deviation of the age in the simulated population
prop_females        <- 0.6                 # proportion of females in the simulated population
prop_males          <- 1 - prop_females    # proportion of males in the simulated population
set.seed(seed)                             # set a seed to ensure reproducible samples
m_indi_features     <- cbind(              # simulate individuals characteristics
  "age" = rnorm(                           # get random samples for 'age' from a normal distribution
    n = num_sim, 
    mean = mean_age, 
    sd = sd_age
  ),
  "sex" = sample(                          # get random samples for 'sex' based on sex distribution
    x = c(0, 1), 
    size = num_sim, 
    replace = TRUE, 
    prob = c(prop_females, prop_males)
  )
)
v_indi_features     <- m_indi_features[1,] # get the characteristics of the first individual 

## Cost inputs 
c_H       <- 2000                          # cost of remaining one cycle healthy
c_S1      <- 4000                          # cost of remaining one cycle sick
c_S2      <- 15000                         # cost of remaining one cycle sicker
c_D       <- 0                             # cost associated with being dead
c_age_cof <- 11.5                          # cost age coefficient
c_sex_cof <- 300                           # cost sex coefficient, where 0 is female and 1 is male

v_states_costs <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)

v_states_names <- c("H","S1", "S2", "D")   # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
occupied_state <- "H"                      # current health state
state_costs    <- 0                        # by default the cost for everyone is zero 
indi_costs     <- 0                        # individual-specific sick/sicker costs
 
# 1. The code line below is intended to create a vector 'v_cost_coeffs'. It is
# a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_cost_coeffs'.
# HINT: the vector 'v_cost_coeffs' should contain the cost regression 
# coefficients for the age and sex variables.

v_cost_coeffs <- c(       # pack the cost regression coefficients in a vector
  "age" = c_age_cof, "sex" = c_sex_cof
)




# 2. Print the 'v_cost_coeffs' object on the console. Is the data in the printed
# object as expected? i.e. are the values associated with the correct names?

v_cost_coeffs
# Yes!




# 3. Write a line of code to extract the coefficient associated with being a 
# male from the vector 'v_cost_coeffs'.
# HINT: use the angular brackets '[[]]'
# HINT: print 'v_cost_coeffs' if you need to remember the names of the elements
# in the vector 'v_cost_coeffs'.

v_cost_coeffs[["sex"]]




# 4. Print the 'v_indi_features' vector and investigate the printed elements.

v_indi_features 




# 5. Write a line of code to extract the sex of the individual from the vector 
# 'v_indi_features'.

v_indi_features[["sex"]]




# 6. Calculate the total individual-specific costs based on their age and sex, 
# assigning the outputs to 'indi_costs'.
# HINT: sex-associated costs can be calculated by multiplying 
# 'v_indi_features[["sex"]]' by the 'v_cost_coeffs[["sex"]]'
# HINT: sum sex-related and age-related sick/sicker costs to get the total 
# individual-specific sick/sicker costs

indi_costs <- v_indi_features[["sex"]] * v_cost_coeffs[["sex"]] + 
  v_indi_features[["age"]] * v_cost_coeffs[["age"]]




# 7. The code you wrote in task (6) is not flexible; it would need to be updated
# if the names of the characteristics change (e.g. sex to gender) or if the 
# model was updated to use more or other individual-level characteristics. Use a
# loop to make this code flexible.
# HINT: use a for loop.
# HINT: loop over all the elements in 'v_cost_coeffs' using their names; for
# example 'names(v_cost_coeffs)'
# HINT: add up the results of each multiplication to 'indi_costs'

indi_costs <- 0
for (variable in names(v_cost_coeffs)) {
  indi_costs <- indi_costs + (v_indi_features[[variable]] * v_cost_coeffs[[variable]])
}




# 8. Review the 'calc_costs' function below. 
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
  state_costs                         <- 0                                 # by default the cost for everyone is zero 
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] + indi_costs # update the cost if sick
  state_costs[occupied_state == "S2"] <- v_states_costs["S2"] + indi_costs # update the cost if sicker
  state_costs[occupied_state == "D"]  <- v_states_costs["D"]               # update the cost if dead
  
  return(state_costs)                                                      # return the costs
}
# 8.1. What are arguments (inputs) of this function?

# 'occupied_state', 'v_states_costs', v_indi_features, and v_cost_coeffs




# 8.2. Is this function suitable to extract the costs associated with each of 
# the four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "H"
# health state.




# 8.3. Create a copy the 'calc_costs' function below and add the missing line.

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
  state_costs                         <- 0                                 # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]               # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] + indi_costs # update the cost if sick
  state_costs[occupied_state == "S2"] <- v_states_costs["S2"] + indi_costs # update the cost if sicker
  state_costs[occupied_state == "D"]  <- v_states_costs["D"]               # update the cost if dead
  
  return(state_costs)                                                      # return the costs
}




# Extra: 1. There are many ways to improve the efficiency of the 'calc_costs()'.
# Suggest and implement one way to improve the efficiency of the code 
# estimating the 'indi_costs'.
# HINT: the code can do without a loop.
# HINT: if you multiply two vectors of equal sizes, R will perform element
# multiplication.
# HINT: what is the result of multiplying the two vectors, what else remains to
# get the total individual-specific costs?
# HINT: use the 'sum()' function.

calc_costs <- function (
    occupied_state,
    v_states_costs,
    v_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  indi_costs <- sum(v_indi_features * v_cost_coeffs)
  
  # estimate costs based on occupied state
  state_costs                         <- 0                                 # by default the cost for everyone is zero 
  state_costs[occupied_state == "H"]  <- v_states_costs["H"]               # update the cost if healthy
  state_costs[occupied_state == "S1"] <- v_states_costs["S1"] + indi_costs # update the cost if sick
  state_costs[occupied_state == "S2"] <- v_states_costs["S2"] + indi_costs # update the cost if sicker
  state_costs[occupied_state == "D"]  <- v_states_costs["D"]               # update the cost if dead
  
  return(state_costs)                                                      # return the costs
}



