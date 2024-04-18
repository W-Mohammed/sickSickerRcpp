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

# Functions:
mat_mult_func <- function(m_trans_probs = m_trans_probs) {
  set.seed(1)
  # create an upper triangular matrix of ones
  m_upper_tri <- upper.tri(
    x = diag(ncol(m_trans_probs)),
    diag = TRUE
  )
  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_upper_tri
  
  return(m_cum_probs)
}

mat_mult_func2 <- function(m_trans_probs = m_trans_probs) {
  set.seed(1)
  # create a lower triangular matrix of ones
  m_lower_tri <- matrix(
    lower.tri(
      x = diag(ncol(m_trans_probs)),
      diag = TRUE
    ), 
    ncol = ncol(m_trans_probs), 
    byrow = TRUE
  )
  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_lower_tri
  
  return(m_cum_probs)
}

loop_func <- function(m_trans_probs = m_trans_probs) {
  set.seed(1)
  # loop through columns to estimate cumulative transition probabilities
  m_cum_probs <- m_trans_probs
  for(i in 2:ncol(m_cum_probs)){
    m_cum_probs[, i] <- m_cum_probs[, i] + m_cum_probs[, i - 1]
  }
  
  return(m_cum_probs)
}

apply_func <- function(m_trans_probs = m_trans_probs) {
  set.seed(1)
  # use apply to estimate cumulative transition probabilities
  m_cum_probs <- t(apply(m_trans_probs, 1, cumsum))
  
  return(m_cum_probs)
}

eline_func <- function (probs = m_trans_probs, m = 1) {
  set.seed(1)
  d <- dim(probs)
  k <- d[2]
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  t(U)
}

identical(loop_func(m_trans_probs = m_trans_probs),
          eline_func(probs = m_trans_probs))
identical(apply_func(m_trans_probs = m_trans_probs),
          eline_func(probs = m_trans_probs))
identical(mat_mult_func(m_trans_probs = m_trans_probs),
          loop_func(m_trans_probs = m_trans_probs) |> `colnames<-`(NULL))
identical(mat_mult_func(m_trans_probs = m_trans_probs),
          mat_mult_func2(m_trans_probs = m_trans_probs))

results <- microbenchmark::microbenchmark(
  "mat" = mat_mult_func(m_trans_probs = m_trans_probs),
  "mat2" = mat_mult_func2(m_trans_probs = m_trans_probs),
  "loop" = loop_func(m_trans_probs = m_trans_probs),
  "eline" = eline_func(probs = m_trans_probs),
  "apply" = apply_func(m_trans_probs = m_trans_probs)
)
results
# Unit: milliseconds
# expr        min        lq       mean     median         uq       max neval
# mat     21.9105   26.1567   45.46324   29.27955   35.32365  603.3019   100
# mat2    22.6374   25.9278   35.14190   28.26885   37.12145  215.3022   100
# loop    67.5821   92.5681  167.10640  111.04605  155.73250  637.9308   100
# eline   96.9682  136.7572  209.37509  150.63835  195.68390  695.3226   100
# apply 2430.5633 3007.2595 3483.54759 3634.91000 3895.88205 4924.6142   100
plot(results)