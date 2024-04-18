## @knitr solution_effects2_function

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

## Utility inputs 
u_H       <- 1            # utility when healthy 
u_S1      <- 0.75         # utility when sick, untreated
u_S2      <- 0.5          # utility when sicker, untreated
u_D       <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
u_sex_cof <- -0.015       # utility sex coefficient, where 0 is female and 1 is male
ru_S1     <- -0.0015      # change in utility of individuals with every additional year being sick
ru_S2     <- -0.0020      # change in utility of individuals with every additional year being sicker

v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)

v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
occupied_state <- "S2"                    # current health state
time_in_state  <- 3                      # time in current health state
state_utility  <- 0                      # by default the utility for everyone is zero
ind_decrement  <- 0                      # individual-specific sick/sicker utility
 
# 1. The code line below is intended to create a vector 'v_util_coeffs'. It is
# a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_util_coeffs'.
# HINT: the vector 'v_util_coeffs' should contain the utility regression 
# coefficients for the age and sex variables.

v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof, "sex" = u_sex_cof
)




# 2. Print the 'v_util_coeffs' object on the console. Is the data in the printed
# object as expected? i.e. are the values associated with the correct names?

v_util_coeffs
# Yes!




# 3. The code line below is intended to create a vector 'v_util_t_decs'. It is
# a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_util_t_decs'.
# HINT: the vector 'v_util_t_decs' should contain the utility regression 
# coefficients for the age and sex variables.

v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "S1" = ru_S1, "S2" = ru_S2
)



# 4. Print the 'v_util_t_decs' vector and investigate the printed elements.

v_util_t_decs 




# 5. Write a line of code to extract the decrements associated with being sick 
# for one cycle from the vector 'v_util_t_decs'.

v_util_t_decs[["S1"]]




# 6. Calculate the total individual-specific utility decrements based on their
# age and sex, assigning the outputs to 'ind_decrement'.
# HINT: sex-associated decrements can be calculated by multiplying 
# 'v_indi_features[["sex"]]' by the 'v_util_coeffs[["sex"]]'
# HINT: sum sex-related and age-related sick/sicker decrements to get the total 
# individual-specific sick/sicker utility decrements

ind_decrement <- v_indi_features[["sex"]] * v_util_coeffs[["sex"]] + 
  v_indi_features[["age"]] * v_util_coeffs[["age"]]




# 7. The code you wrote in task (6) is not flexible; it would need to be updated
# if the names of the characteristics change (e.g. sex to gender) or if the 
# model was updated to use more or other individual-level characteristics. Use a
# loop to make this code flexible.
# HINT: use a for loop.
# HINT: loop over all the elements in 'v_util_coeffs' using their names; for
# example 'names(v_util_coeffs)'
# HINT: add up the results of each multiplication to 'ind_decrement'

ind_decrement <- 0
for (variable in names(v_util_coeffs)) {
  ind_decrement <- ind_decrement + (v_indi_features[[variable]] * v_util_coeffs[[variable]])
}




# 8. Review the 'calc_effs' function below. 
calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    v_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  ind_decrement <- 0
  for (variable in names(v_util_coeffs)) {
    ind_decrement <- ind_decrement + (v_indi_features[[variable]] * v_util_coeffs[[variable]])
  }
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- 0
  time_decrement[occupied_state == "S1"] <-  v_util_t_decs["S1"] * time_in_state[occupied_state == "S1"]
  time_decrement[occupied_state == "S2"] <-  v_util_t_decs["S2"] * time_in_state[occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  state_utility <- 0                                                            # by default the utility for everyone is zero
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] + decrement # update the utility if sick
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] + decrement # update the utility if sicker
  state_utility[occupied_state == "D"]  <- v_states_utilities["D"]              # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}
# 8.1. What are arguments (inputs) of this function?

# 'occupied_state', 'v_states_utilities', v_indi_features, v_util_coeffs, 
# v_util_t_decs, time_in_state, and cycle_length




# 8.2. Is this function suitable to extract the decrements associated with each
# of the four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "H"
# health state.




# 8.3. Create a copy the 'calc_effs' function below and add the missing line.

calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    v_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  ind_decrement <- 0
  for (variable in names(v_util_coeffs)) {
    ind_decrement <- ind_decrement + (v_indi_features[[variable]] * v_util_coeffs[[variable]])
  }
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- 0
  time_decrement[occupied_state == "S1"] <-  v_util_t_decs["S1"] * time_in_state[occupied_state == "S1"]
  time_decrement[occupied_state == "S2"] <-  v_util_t_decs["S2"] * time_in_state[occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  state_utility <- 0                                                            # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  + decrement # update the utility if healthy
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] + decrement # update the utility if sick
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] + decrement # update the utility if sicker
  state_utility[occupied_state == "D"]  <- v_states_utilities["D"]              # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}




# Extra: 1. There are many ways to improve the efficiency of the 'calc_effs()'.
# Suggest and implement one way to improve the efficiency of the code 
# estimating the 'ind_decrement'.
# HINT: the code can do without a loop.
# HINT: if you multiply two vectors of equal sizes, R will perform element
# multiplication.
# HINT: what is the result of multiplying the two vectors, what else remains to
# get the total individual-specific decrements?
# HINT: use the 'sum()' function.

calc_effs <- function (
    occupied_state, 
    v_states_utilities,
    v_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  ind_decrement <- 0
  ind_decrement  <- sum(v_indi_features * v_util_coeffs)   
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- 0
  time_decrement[occupied_state == "S1"] <- v_util_t_decs["S1"] * time_in_state[occupied_state == "S1"]
  time_decrement[occupied_state == "S2"] <- v_util_t_decs["S2"] * time_in_state[occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  state_utility <- 0                                                            # by default the utility for everyone is zero
  state_utility[occupied_state == "H"]  <- v_states_utilities["H"]  + decrement # update the utility if healthy
  state_utility[occupied_state == "S1"] <- v_states_utilities["S1"] + decrement # update the utility if sick
  state_utility[occupied_state == "S2"] <- v_states_utilities["S2"] + decrement # update the utility if sicker
  state_utility[occupied_state == "D"]  <- v_states_utilities["D"]              # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  state_utility * cycle_length                            # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                     # return the QALYs
}



