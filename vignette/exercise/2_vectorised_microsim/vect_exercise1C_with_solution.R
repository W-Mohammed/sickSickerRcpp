## @knitr solution_costsV_function

# 0. Execute the lines below.
# clear R's session memory (Global Environment) 
rm(list = ls())

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
m_indi_features     <- cbind(                 # simulate individuals characteristics
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

## Cost inputs 
c_H       <- 2000                             # cost of remaining one cycle healthy
c_S1      <- 4000                             # cost of remaining one cycle sick
c_S2      <- 15000                            # cost of remaining one cycle sicker
c_D       <- 0                                # cost associated with being dead
c_age_cof <- 11.5                             # cost age coefficient
c_sex_cof <- 300                              # cost sex coefficient, where 0 is female and 1 is male

v_states_costs <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)

v_state_costs  <- rep(x = 0, times = num_sim) # by default the cost for everyone is zero 
v_indi_costs   <- rep(x = 0, times = num_sim) # individual-specific sick/sicker costs

# 1. The code line below is intended to create a vector 'v_cost_coeffs'. It is
# a named vector, i.e. the elements in the vector has names. Complete the 
# code replacing "..." with the correct values. Once done, execute the line to 
# create the object 'v_cost_coeffs'.
# HINT: the vector 'v_cost_coeffs' should contain the cost regression 
# coefficients for the age and sex variables.

v_cost_coeffs <- c(                           # pack the cost regression coefficients in a vector
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




# 4. Print the 'm_indi_features' matrix and investigate the printed elements.

m_indi_features 




# 5. Write a line of code to extract the sex of the individuals from the matrix 
# 'm_indi_features'.
# HINT: use the square brackets '[]' and a comma ',' to subset the "sex" column

m_indi_features[, "sex"]




# 6. Calculate the total individual-specific costs based on their age and sex, 
# assigning the outputs to 'v_indi_costs'.
# HINT: sex-associated costs can be calculated by multiplying 
# 'm_indi_features[, "sex"]' by the 'v_cost_coeffs[["sex"]]'
# HINT: sum sex-related and age-related sick/sicker costs to get the total 
# individual-specific sick/sicker costs

v_indi_costs <- m_indi_features[, "sex"] * v_cost_coeffs[["sex"]] + 
  m_indi_features[, "age"] * v_cost_coeffs[["age"]]




# 7. The code you wrote in task (6) is not flexible; it would need to be updated
# if the names of the characteristics change (e.g. sex to gender) or if the 
# model was updated to use more or other individual-level characteristics. Use 
# matrix multiplication to compute the total individual-specific costs.
# HINT: number of columns in the first matrix must be equal to the number of 
# rows in the second matrix in matrix multiplication.
# HINT: by default, R will treat the vector as either row-matrix or 
# column-matrix as needed by the matrix multiplication.
# HINT: the 'm_indi_features' matrix has two columns which is equal to the 
# number of rows in the column-vector 'v_cost_coeffs'

v_indi_costs <- m_indi_features %*% v_cost_coeffs




# 8. Review the 'calc_costsV' function below. 
calc_costsV <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # estimate costs based on occupied state
  v_state_costs                           <- rep(0, length(v_occupied_state))                              # by default the cost for everyone is zero 
  v_state_costs[v_occupied_state == "S1"] <- v_states_costs["S1"] + v_indi_costs[v_occupied_state == "S1"] # update the cost if sick
  v_state_costs[v_occupied_state == "S2"] <- v_states_costs["S2"] + v_indi_costs[v_occupied_state == "S2"] # update the cost if sicker
  v_state_costs[v_occupied_state == "D"]  <- v_states_costs["D"]                                           # update the cost if dead
  
  return(v_state_costs)                                                                                    # return the costs
}
# 8.1. What are arguments (inputs) of this function?

# 'v_occupied_state', 'v_states_costs', v_indi_features, and v_cost_coeffs




# 8.2. Is this function suitable to extract the costs associated with each of 
# the four health states in the healthy-sick-sicker-dead model?
# HINT: Print the 'v_states_names' object to see the states names.

# No. 
# In its current state, the function does not take into consideration the "H"
# health state.




# 8.3. Create a copy the 'calc_costsV' function below and add the missing line.

calc_costsV <- function (
    v_occupied_state,
    v_states_costs,
    m_indi_features,
    v_cost_coeffs) {
  
  # calculate individual-specific costs based on costs regression coefficients
  v_indi_costs <- m_indi_features %*% v_cost_coeffs
  
  # estimate costs based on occupied state
  v_state_costs                           <- rep(0, length(v_occupied_state))                              # by default the cost for everyone is zero 
  v_state_costs[v_occupied_state == "H"]  <- v_states_costs["H"]                                           # update the cost if healthy
  v_state_costs[v_occupied_state == "S1"] <- v_states_costs["S1"] + v_indi_costs[v_occupied_state == "S1"] # update the cost if sick
  v_state_costs[v_occupied_state == "S2"] <- v_states_costs["S2"] + v_indi_costs[v_occupied_state == "S2"] # update the cost if sicker
  v_state_costs[v_occupied_state == "D"]  <- v_states_costs["D"]                                           # update the cost if dead
  
  return(v_state_costs)                                                                                    # return the costs
}



