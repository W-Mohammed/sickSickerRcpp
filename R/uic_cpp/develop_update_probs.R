# clear R working global environment
rm(list = ls())

# source update_probsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "update_probsC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the update_probsV functions:
update_probsV1 <- function(
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
        test = all(abs(rowSums(m_probs) - 1) < 1e-12), # check if the transition probabilities add up to 1
        yes = return(m_probs),                         # return the transition probabilities
        no = stop("Probabilities do not sum to 1")     # or produce an error
      )
    }
  )
}

update_probsV2 <- function(
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
        ncol = length(v_states_names)  # a column for each state
      )
      
      # update m_probs with the appropriate probabilities
      m_probs[v_occupied_state == "H", ]  <- cbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD) # transition probabilities when healthy
      m_probs[v_occupied_state == "S1", ] <- cbind(p_S1H, p_S1S1, p_S1S2, p_S1D)     # transition probabilities when sick
      m_probs[v_occupied_state == "S2", ] <- cbind(0, 0, 1 - p_S2D, p_S2D)           # transition probabilities when sicker
      m_probs[v_occupied_state == "D", ]  <- cbind(0, 0, 0, p_DD)                    # transition probabilities when dead
      
      # sanity check
      ifelse(
        test = all(abs(rowSums(m_probs) - 1) < 1e-12), # check if the transition probabilities add up to 1
        yes = return(m_probs),                         # return the transition probabilities
        no = stop("Probabilities do not sum to 1")     # or produce an error
      )
    }
  )
}

update_probsV3 <- function(
    v_states_index,
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
      p_S1D  <- 1 - exp(- r_S1D * (1 + v_time_in_state[v_occupied_state == 2] * rp_S1))
      p_S2D  <- 1 - exp(- r_S2D * (1 + v_time_in_state[v_occupied_state == 3] * rp_S2))
      
      p_S1S1 <- 1 - p_S1S2 - p_S1H - p_S1D
      p_HD   <- rep(p_HD, length(which(v_occupied_state == 1)))
      p_DD   <- rep(1, length(v_occupied_state[v_occupied_state == 4]))
      # Create a state transition probabilities matrix
      m_probs <- matrix(
        nrow = length(v_time_in_state), # a row for each individual
        ncol = length(v_states_names)  # a column for each state
      )
      
      # update m_probs with the appropriate probabilities
      m_probs[v_occupied_state == 1, ]  <- cbind(1 - p_HS1 - p_HD, p_HS1, 0, p_HD) # transition probabilities when healthy
      m_probs[v_occupied_state == 2, ] <- cbind(p_S1H, p_S1S1, p_S1S2, p_S1D)     # transition probabilities when sick
      m_probs[v_occupied_state == 3, ] <- cbind(0, 0, 1 - p_S2D, p_S2D)           # transition probabilities when sicker
      m_probs[v_occupied_state == 4, ]  <- cbind(0, 0, 0, p_DD)                    # transition probabilities when dead
      
      # sanity check
      ifelse(
        test = all(abs(rowSums(m_probs) - 1) < 1e-12), # check if the transition probabilities add up to 1
        yes = return(m_probs),                         # return the transition probabilities
        no = stop("Probabilities do not sum to 1")     # or produce an error
      )
    }
  )
}

#------------------------------------------------------------------------------#

# define function inputs:
## General parameters
seed    <- 1234                            # random number generator state
num_i   <- 1e6                             # number of simulated individuals
v_states_names <- c("H","S1", "S2", "D")   # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_states_index <- 1:length(v_states_names) # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
set.seed(seed)                             # set the seed to ensure reproducible samples below
v_occupied_state <- sample(                # sample current health state
  x       = v_states_names,                # from the four health states
  size    = num_i,                         # sample a state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)         # sample 75% as healthy, 20% as sick, and 5% as sicker
)
set.seed(seed)                             # set the seed to ensure reproducible samples below
v_occupied_state2 <- sample(               # sample current health state
  x       = v_states_index,                # from the four health states
  size    = num_i,                         # sample a state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = c(0.75, 0.20, 0.05, 0)         # sample 75% as healthy, 20% as sick, and 5% as sicker
)
v_time_in_state  <- sample(                # sample time in current health state
  x       = 1:10,                          # for values between 1 and 10
  size    = num_i,                         # sample time in state for each of the simulated individuals
  replace = TRUE,                          # allow sampled states to be re-sampled
  prob    = rep(                           # values 1 to 10 will have equal chances of being sampled
    x     = 1/length(1:10),                # the chance of sampling any the values is 1/the number of values
    times = length(1:10)                   # repeat the probabilities as many times as there as values
  )
)

## Transition probabilities (per cycle)
p_HD    <- 0.005            # probability to die when healthy
p_HS1   <- 0.15             # probability to become sick when healthy
p_S1H   <- 0.5              # probability to become healthy when sick
p_S1S2  <- 0.105            # probability to become sicker when sick
rr_S1   <- 3                # rate ratio of death in sick vs healthy
rr_S2   <- 10               # rate ratio of death in sicker vs healthy
r_HD    <- -log(1 - p_HD)   # rate of death in healthy
r_S1D   <- rr_S1 * r_HD     # rate of death in sick
r_S2D   <- rr_S2 * r_HD     # rate of death in sicker
p_S1D   <- 1 - exp(- r_S1D) # probability to die in sick
p_S2D   <- 1 - exp(- r_S2D) # probability to die in sicker
rp_S1   <- 0.02             # increase in mortality rate with every additional year being sick
rp_S2   <- 0.05             # increase in mortality rate with every additional year being sicker

l_trans_probs <- list(      # pack the transition probabilities and rates in a list
  "p_HD"   = p_HD,
  "p_HS1"  = p_HS1,
  "p_S1H"  = p_S1H,
  "p_S1S2" = p_S1S2,
  "p_S1D"  = p_S1D,
  "p_S2D"  = p_S2D,
  "rp_S1"  = rp_S1,
  "rp_S2"  = rp_S2
)

#------------------------------------------------------------------------------#

# run the update_probsV function:
R_results1 <- update_probsV1(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
R_results2 <- update_probsV2(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
R_results3 <- update_probsV3(
  v_states_index = v_states_index,
  v_occupied_state = v_occupied_state2,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
# run the update_probsC function:
C_results1 <- update_probsC1(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
C_results2 <- update_probsC2(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
C_results3 <- update_probsC3(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
C_results4 <- update_probsC4(
  v_states_index = v_states_index,
  v_occupied_state = v_occupied_state2,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
C_results5 <- update_probsC5(
  v_states_index = v_states_index,
  v_occupied_state = v_occupied_state2,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
C_results6 <- update_probsC6(
  v_states_index = v_states_index,
  v_occupied_state = v_occupied_state2,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
# check results
identical(R_results1 |> `dimnames<-`(NULL), R_results2)
identical(R_results2, R_results3)
identical(R_results1, C_results1)
identical(C_results2, C_results2)
identical(C_results2, C_results3)
identical(C_results3, C_results4)
identical(C_results4, C_results5)
identical(C_results5, C_results6)

#------------------------------------------------------------------------------#

# benchmark the functions

update_probs_RvC <- bench::mark(
  "R_1" = update_probsV1(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "R_2" = update_probsV2(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "R_3" = update_probsV3(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_1" = update_probsC1(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_2" = update_probsC2(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_3" = update_probsC3(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_4" = update_probsC4(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_5" = update_probsC5(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_6" = update_probsC6(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  check = FALSE
)

update_probs_RvC[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

update_probs_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = update_probsV1(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "R_2" = update_probsV2(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "R_3" = update_probsV3(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_1" = update_probsC1(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_2" = update_probsC2(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_3" = update_probsC3(
    v_states_names = v_states_names,
    v_occupied_state = v_occupied_state,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_4" = update_probsC4(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_5" = update_probsC5(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  ),
  "C_6" = update_probsC6(
    v_states_index = v_states_index,
    v_occupied_state = v_occupied_state2,
    l_trans_probs = l_trans_probs,
    v_time_in_state = v_time_in_state
  )
)

update_probs_RvC2
plot(update_probs_RvC2)
saveRDS(object = update_probs_RvC2, file = "update_probs_RvC2")