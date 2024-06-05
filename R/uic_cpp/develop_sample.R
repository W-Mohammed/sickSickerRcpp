# clear R working global environment
rm(list = ls())

# source update_probsC functions
Rcpp::sourceCpp(
  file = file.path(
    here::here(),
    "src",
    "sampleC.cpp"
  )
)

#------------------------------------------------------------------------------#

# define the sampleV functions:
sampleV1 <- function(
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

sampleV2 <- function(
    m_trans_probs,
    v_states_index) {
  
  # create an upper triangular matrix of ones
  m_upper_tri <- upper.tri(
    x = diag(ncol(m_trans_probs)),
    diag = TRUE
  )
  
  # create matrix with row-wise cumulative transition probabilities
  m_cum_probs <- m_trans_probs %*% m_upper_tri

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
      each = length(v_states_index)
    ),
    nrow  = nrow(m_trans_probs),
    ncol  = length(v_states_index),
    byrow = TRUE
  )
  
  # identify transitions, compare random samples to cumulative probabilities
  m_transitions <- m_rand_values > m_cum_probs # transitions from first state
  
  # sum transitions to identify health state in next cycle
  v_transitions <- rowSums(m_transitions)
  
  # identify health state to which each individual is transitioning
  v_health_states <- 1 + v_transitions
  
  return(v_health_states)
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
set.seed(seed)
m_trans_probs1 <- update_probsV2(
  v_states_names = v_states_names,
  v_occupied_state = v_occupied_state,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)
set.seed(seed)
m_trans_probs2 <- update_probsV3(
  v_states_index = v_states_index,
  v_occupied_state = v_occupied_state2,
  l_trans_probs = l_trans_probs,
  v_time_in_state = v_time_in_state
)

#------------------------------------------------------------------------------#

# run the sampleV function:
set.seed(seed)
R_results1 <- sampleV1(
  m_trans_probs = m_trans_probs1,
  v_states_names = v_states_names
)
set.seed(seed)
R_results2 <- sampleV2(
  m_trans_probs = m_trans_probs2,
  v_states_index = v_states_index
)
# run the sampleC function:
set.seed(seed)
C_results0 <- sampleC0(
  m_trans_probs = m_trans_probs1,
  v_states_names = v_states_names
)
set.seed(seed)
C_results1 <- sampleC1(
  m_trans_probs = m_trans_probs2
) |> 
  as.numeric()
set.seed(seed)
C_results2 <- sampleC2(
  m_trans_probs = m_trans_probs2
)
set.seed(seed)
C_results3 <- sampleC3(
  m_trans_probs = m_trans_probs2
)
set.seed(seed)
C_results4 <- sampleC4(
  m_trans_probs = m_trans_probs2
)
set.seed(seed)
C_results5 <- sampleC5(
  m_trans_probs = m_trans_probs2
)
# check results
identical(R_results1, C_results0)
R_results1.1 <- R_results1
R_results1.1[R_results1 == "H"]  <- 1
R_results1.1[R_results1 == "S1"] <- 2
R_results1.1[R_results1 == "S2"] <- 3
R_results1.1[R_results1 == "D"]  <- 4
identical(R_results1.1 |> as.numeric(), R_results2)
identical(R_results2, C_results1 |> as.numeric())
identical(R_results2, C_results2[, 1])
identical(C_results2[, 1], C_results3[, 1])
identical(C_results3[, 1], C_results4[, 1])
identical(C_results4[, 1], C_results5[, 1]) # integer vs numeric
testthat::expect_equal(C_results4, C_results5)

#------------------------------------------------------------------------------#

# benchmark the functions
sample_RvC <- bench::mark(
  "R_1" = sampleV1(
    m_trans_probs = m_trans_probs1,
    v_states_names = v_states_names
  ),
  "R_2" = sampleV2(
    m_trans_probs = m_trans_probs2,
    v_states_index = v_states_index
  ),
  "C0" = sampleC0(
    m_trans_probs = m_trans_probs1,
    v_states_names = v_states_names
  ),
  "C_1" = sampleC1(
    m_trans_probs = m_trans_probs2
  ),
  "C_2" = sampleC2(
    m_trans_probs = m_trans_probs2
  ),
  "C_3" = sampleC3(
    m_trans_probs = m_trans_probs2
  ),
  "C_4" = sampleC4(
    m_trans_probs = m_trans_probs2
  ),
  "C_5" = sampleC5(
    m_trans_probs = m_trans_probs2
  ),
  check = FALSE
)

sample_RvC[c("expression", "min", "median", "itr/sec", "n_gc", "mem_alloc")]

sample_RvC2 <- microbenchmark::microbenchmark(
  "R_1" = sampleV1(
    m_trans_probs = m_trans_probs1,
    v_states_names = v_states_names
  ),
  "R_2" = sampleV2(
    m_trans_probs = m_trans_probs2,
    v_states_index = v_states_index
  ),
  "C0" = sampleC0(
    m_trans_probs = m_trans_probs1,
    v_states_names = v_states_names
  ),
  "C_1" = sampleC1(
    m_trans_probs = m_trans_probs2
  ),
  "C_2" = sampleC2(
    m_trans_probs = m_trans_probs2
  ),
  "C_3" = sampleC3(
    m_trans_probs = m_trans_probs2
  ),
  "C_4" = sampleC4(
    m_trans_probs = m_trans_probs2
  ),
  "C_5" = sampleC5(
    m_trans_probs = m_trans_probs2
  )
)

sample_RvC2
plot(sample_RvC2)
saveRDS(object = sample_RvC2, file = "sample_RvC2")