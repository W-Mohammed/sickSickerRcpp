## @knitr solution_effectsV_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

## General parameters
## Population characteristics/features
seed                <- 1234                   # random number generator state
num_sim             <- 100                    # number of simulated individuals
v_states_names      <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
mean_age            <- 50                     # mean age in the simulated population
sd_age              <- 3                      # standard deviation of the age in the simulated population
prop_females        <- 0.6                    # proportion of females in the simulated population
prop_males          <- 1 - prop_females       # proportion of males in the simulated population
set.seed(seed)                                # set a seed to ensure reproducible samples
v_occupied_state <- sample(                   # sample current health state
  x       = v_states_names,                   # from the four health states
  size    = num_sim,                          # sample a state for each of the simulated individuals
  replace = TRUE,                             # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)            # sample 75% as healthy, 20% as sick, and 5% as sicker
)
v_time_in_state  <- sample(                   # sample time in current health state
  x       = 1:10,                             # for values between 1 and 10
  size    = num_sim,                          # sample time in state for each of the simulated individuals
  replace = TRUE,                             # allow sampled states to be re-sampled
  prob    = rep(                              # values 1 to 10 will have equal chances of being sampled
    x     = 1/length(1:10),                   # the chance of sampling any the values is 1/the number of values
    times = length(1:10)                      # repeat the probabilities as many times as there as values
  )
)
m_indi_features  <- cbind(                    # simulate individuals characteristics
  "age" = rnorm(                              # get random samples for 'age' from a normal distribution
    n = num_sim, 
    mean = mean_age, 
    sd = sd_age
  ),
  "sex" = sample(                             # get random samples for 'sex' based on sex distribution
    x = c(0, 1), 
    size = num_sim, 
    replace = TRUE, 
    prob = c(prop_females, prop_males)
  )
)

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

v_state_utility  <- rep(x = 0, times = num_sim) # by default the utility for everyone is zero
v_ind_decrement  <- rep(x = 0, times = num_sim) # individual-specific sick/sicker utility
 
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
# age and sex, assigning the outputs to 'v_ind_decrement'.
# HINT: sex-associated decrements can be calculated by multiplying 
# 'm_indi_features[, "sex"]' by the 'v_util_coeffs[["sex"]]'
# HINT: sum sex-related and age-related sick/sicker decrements to get the total 
# individual-specific sick/sicker utility decrements

v_ind_decrement <- m_indi_features[, "sex"] * v_util_coeffs[["sex"]] + 
  m_indi_features[, "age"] * v_util_coeffs[["age"]]




# 7. The code you wrote in task (6) is not flexible; it would need to be updated
# if the names of the characteristics change (e.g. sex to gender) or if the 
# model was updated to use more or other individual-level characteristics.Use 
# matrix multiplication to compute the total individual-specific utility 
# decrements.
# HINT: number of columns in the first matrix must be equal to the number of 
# rows in the second matrix in matrix multiplication.
# HINT: by default, R will treat the vector as either row-matrix or 
# column-matrix as needed by the matrix multiplication.
# HINT: the 'm_indi_features' matrix has two columns which is equal to the 
# number of rows in the column-vector 'v_ind_decrement'

v_ind_decrement <- m_indi_features %*% v_util_coeffs




# 8. Review the 'calc_effsV' function below. 
calc_effsV <- function (
    v_occupied_state, 
    v_states_utilities,
    m_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- m_indi_features %*% v_util_coeffs
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == "S1"] <- v_util_t_decs["S1"] * v_time_in_state[v_occupied_state == "S1"]
  time_decrement[v_occupied_state == "S2"] <- v_util_t_decs["S2"] * v_time_in_state[v_occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- v_ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  v_state_utility                           <- rep(0, length(v_occupied_state))                               # by default the utility for everyone is zero
  v_state_utility[v_occupied_state == "S1"] <- v_states_utilities["S1"] + decrement[v_occupied_state == "S1"] # update the utility if sick
  v_state_utility[v_occupied_state == "S2"] <- v_states_utilities["S2"] + decrement[v_occupied_state == "S2"] # update the utility if sicker
  v_state_utility[v_occupied_state == "D"]  <- v_states_utilities["D"]                                        # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                                    # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                                                               # return the QALYs
}
# 8.1. What are arguments (inputs) of this function?

# v_occupied_state, v_states_utilities, m_indi_features, v_util_coeffs, 
# v_util_t_decs, v_time_in_state, and cycle_length




# 8.2. Is this function suitable to extract the decrements associated with each
# of the four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "H"
# health state.




# 8.3. Create a copy the 'calc_effsV' function below and add the missing line.

calc_effsV <- function (
    v_occupied_state, 
    v_states_utilities,
    m_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- m_indi_features %*% v_util_coeffs
  
  # calculate time-dependent state-specific utility decrements
  time_decrement <- rep(0, length(v_occupied_state))
  time_decrement[v_occupied_state == "S1"] <- v_util_t_decs["S1"] * v_time_in_state[v_occupied_state == "S1"]
  time_decrement[v_occupied_state == "S2"] <- v_util_t_decs["S2"] * v_time_in_state[v_occupied_state == "S2"]
  
  # estimate total decrements
  decrement <- v_ind_decrement + time_decrement
  
  # estimate utilities based on occupied state
  v_state_utility                           <- rep(0, length(v_occupied_state))                               # by default the utility for everyone is zero
  v_state_utility[v_occupied_state == "H"]  <- v_states_utilities["H"]  + decrement[v_occupied_state == "H"]  # update the utility if healthy
  v_state_utility[v_occupied_state == "S1"] <- v_states_utilities["S1"] + decrement[v_occupied_state == "S1"] # update the utility if sick
  v_state_utility[v_occupied_state == "S2"] <- v_states_utilities["S2"] + decrement[v_occupied_state == "S2"] # update the utility if sicker
  v_state_utility[v_occupied_state == "D"]  <- v_states_utilities["D"]                                        # update the utility if dead
  
  # calculate Quality Adjusted Life Years (QALYs)
  QALYs <-  v_state_utility * cycle_length                                                                    # calculate the QALYs during cycle `t`
  
  return(QALYs)                                                                                               # return the QALYs
}



