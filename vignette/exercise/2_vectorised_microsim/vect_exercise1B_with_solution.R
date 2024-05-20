## @knitr solution_sampleV_function

# 0. Execute the lines below.
## clear R's session memory (Global Environment) 
rm(list = ls())

## General parameters
seed    <- 1234                          # random number generator state
num_sim <- 1e6                           # number of simulated individuals
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

## Generate transition probabilities matrix
update_probsV <- function(              # define the update_probsV function to generate m_trans_probs
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
m_trans_probs <- update_probsV(        # generate transition probabilities matrix
  v_states_names   = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs    = l_trans_probs,
  v_time_in_state  = v_time_in_state
) 

# 1. The 'sampleV' function is a vectorised version of R's 'sample' function.
# - Run the 'sample' function in the code below. Can 'sample()' sample a state
# for more than one individual? Remember each individual has their own 
# transition probabilities values.

set.seed(seed)
sample(x = v_states_names, size = 1, replace = TRUE, prob = m_trans_probs[1,])

# H
# No, 'sample()' only accepts one vector of probabilities in the argument 'prob'




# 2. The 'sampleV()' function utilizes multinomial sampling. It requires the
# cumulative probabilities of transitioning from current state for each of the
# simulated individuals.
# 2.1. Print the head of (using 'head()') the 'm_trans_probs' matrix. 
# Investigate the printed values.

head(m_trans_probs)



# 2.2. Each row in 'm_trans_probs' represents each individual's probabilities of
# transitioning to states in the next cycle.
# - Execute the lines of code below (the for loop) to compute the cumulative
# transition probabilities for the first individual
# - Investigate the cumulative values. Are they correctly estimated? 
# HINT: compare the printed values of 'v_indi_cum_probs' to 'm_trans_probs[1, ]'

v_indi_cum_probs <- m_trans_probs[1, ]
for (i in 2:length(v_indi_cum_probs)) {
  v_indi_cum_probs[i] <- v_indi_cum_probs[i - 1] + v_indi_cum_probs[i] 
}
v_indi_cum_probs




# 2.3. Based on the code (the 'for' loop) in task (2.2):
# - compute the cumulative probabilities for each of the simulated individuals. 
# - assign the outputs to 'm_cum_probs', instead of 'v_indi_cum_probs'.
# HINT: In the above code 'm_trans_probs[1, ]' was used to subset the transition
# probabilities for the first individual. Amend this code to take all 
# individuals transition probabilities.
# HINT: In the 'for' loop, amend the code to subset the i-th columns, for 
# example '[, i]', instead of the probability of the i-th state '[i]'.
# HINT: In the 'for' loop amend the code loop through the second to the last
# columns in the data structure.

m_cum_probs2 <- m_trans_probs
for (i in 2:ncol(m_cum_probs2)) {
  m_cum_probs2[, i] <- m_cum_probs2[, i - 1] + m_cum_probs2[, i] 
}
head(m_cum_probs2)




# 3. The 'for' loop in (2.2) and (2.3) works fine but is not the most efficient.
# In R, matrix multiplication is a highly vectorized approach. There are a few
# steps involved in computing cumulative probabilities using matrix 
# multiplication. Follow the sub-tasks below to cover the necessary steps.
# 3.1. First, run the code below to construct a diagonal matrix. 
# - Investigate the the printed matrix.

diag(nrow = 3)




# 3.2. Execute the line of code below to perform element-wise multiplication on
# the matrix created by the 'diag()' function.
# - Study the resulting matrix.

diag(nrow = 3) * 2




# 3.3. Run the code lines below to create an upper triangular matrix of ones. 
# - Investigate the printed matrix.
# - Do you see any ones?
# - Multiply the upper matrix by 2

upper.tri(
  x = diag(nrow = 3),
  diag = TRUE
)




# 3.4. Copy the code from (3.2) to recreate an upper triangular matrix and:
# - multiply (element-wise) the upper triangular matrix by 2.
# - do not assign the results so that R prints the results of the operation.
# - notice how the logical (Boolean) values 'TRUE' are equivalent to 1 in R when
# they are used in numerical contexts.
# Note: In R, the Boolean TRUE and FALSE are effectively treated as 1 and 0 
# respectively when used in numerical calculations. 

upper.tri(
  x = diag(nrow = 3),
  diag = TRUE
) * 2




# 3.5. Execute the lines below to:
# - create an upper triangular matrix,
# - print the created triangular matrix,
# - subset the 'm_trans_probs' to the top five individuals,
# - compute the cumulative probabilities and assign it to 'm_cum_probs',
# - print the 'm_cum_probs' object,
# - compare 'm_cum_probs' to 'm_cum_probs2[1:5, ]', and
# - compute the cumulative probabilities for all individuals and assign it to
# 'm_cum_probs',

m_upper_tri  <- upper.tri(
  x = diag(nrow = ncol(m_trans_probs)),
  diag = TRUE
)
m_upper_tri
m_cum_probs <- m_trans_probs[1:5, ] %*% m_upper_tri
m_cum_probs
identical(m_cum_probs         |> `dimnames<-`(NULL), 
          m_cum_probs2[1:5, ] |> `dimnames<-`(NULL))
m_cum_probs <- m_trans_probs %*% m_upper_tri




# 4. The implementation of multinomial sampling in the 'sampleV()' function
# depends on comparing the cumulative distributions of transitioning 
# probabilities to randomly sampled values (sampled from a standard Uniform 
# distribution U ~ (0, 1)).
# 4.1. Run the lines of code below to draw 5 values from a standard Uniform
# distribution.

set.seed(seed)
v_rand_values <- runif(n = 5)
v_rand_values
# 0.1137034 0.6222994 0.6092747 0.6233794 0.8609154




# 4.2. Execute the following lines of code to transform the vector of randomly
# drawn values in 'v_rand_values' to a matrix.
# - Notice how the "each" argument is used within the 'rep()' function.
# - Confirm that each row contains one sampled value repeated over the columns.

m_rand_values <- matrix(
  data = rep(
    x = v_rand_values, 
    each = ncol(m_cum_probs)),
  ncol = ncol(m_cum_probs),
  byrow = TRUE
)
m_rand_values




# 4.3. Run the lines of code below to compare the random values in 
# 'm_rand_values' to the cumulative distribution probabilities of the top five
# individuals in 'm_cum_probs[1:5, ]'.

m_transitions <- m_rand_values > m_cum_probs[1:5, ]
m_transitions



# 4.4. Execute the following lines of code to get the row sums of the matrix 
# 'm_transitions'. These sums represent the health states each individual would
# transition to or remain in in the next cycle, relative to the first health 
# state. The first health state is the health state corresponding to the 
# probabilities in the first column of the 'm_trans_probs' matrix.

v_transitions <- rowSums(m_transitions)
v_transitions
# H  H  H  H S1 
# 0  0  0  0  1




# 4.5. Run the lines of code below to identify the names or labels of the health
# states to be occupied by the five individuals in the health state.
# - Notice how we added 1 to the 'v_transitions'. This is because the values in
# 'v_transitions' represent transitions relative to the first health state.
# - Do these differ from the states in the current cycled?

v_states_names[1 + v_transitions] # identify the states for the next cycle
v_occupied_state[1:5]

# 5. Review the 'sampleV()' function defined below.
sampleV <- function(
    m_trans_probs,
    v_states_names) {
  
  # create an upper triangular matrix of ones
  m_upper_tri <- upper.tri(
    x = diag(ncol(m_trans_probs)),
    diag = TRUE
  )
  
  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  colnames(m_cum_probs) <- v_states_names
  
  # ensure that the maximum cumulative probabilities are equal to 1
  if (any(m_cum_probs[, ncol(m_cum_probs)] > 1.000000)) {
    stop("Error in multinomial sampling: probabilities do not sum to 1")
  }
  
  # sample random values from Uniform standard distribution for each individual
  v_rand_values <- runif(n = nrow(m_trans_probs))
  
  # repeat each sampled value to have as many copies as the number of states 
  m_rand_values <- matrix(
    data  = rep(
      x = v_rand_values,
      each = length(v_states_names)
    ), 
    nrow  = nrow(m_trans_probs), 
    ncol  = length(v_states_names), 
    byrow = TRUE
  )
  
  # identify transitions, compare random samples to cumulative probabilities
  m_transitions <- m_rand_values > m_cum_probs # transitions from first state
  
  # sum transitions to identify health state in next cycle
  v_transitions <- rowSums(m_transitions)
  
  # identify health state to which each individual is transitioning
  v_health_states <- v_states_names[1 + v_transitions]
  
  return(v_health_states)
}
