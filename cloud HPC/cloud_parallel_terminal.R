## @knitr micro_parallel_terminal_on_Linode_demo

# This script is adapted to be called from terminal with minimal user inputs. It
# aims to perform PSA via 'multisession' and 'multicore' (if supported by the 
# OS). 
# This script is intended to be demoed on cloud computing, Linode.

# From the terminal, cd to folder containing file, and run:
# Rscript parallel_terminal.R

# clear R's session memory (Global Environment) 
rm(list = ls())

# setting session-wide RNG type for reproducible parallel/serial results
RNGkind("L'Ecuyer-CMRG")

#---

                        ### Defining model functions ###

# Define model functions
## Update Transition Probability function
### This function updates the transition probabilities at every cycle based on
### the health state occupied by each individual at cycle 't' and the time spent
### in the states
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
        test = all(abs(rowSums(m_probs) - 1) < 1e-12), # check if the transition probabilities add up to 1
        yes = return(m_probs),                         # return the transition probabilities
        no = stop("Probabilities do not sum to 1")     # or produce an error
      )
    }
  )
}

## Sample Health States function
### This function identifies the health state each individual will transition
### to in the next model cycle
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

## Calculate Costs function
### This function estimates the costs at every cycle based on the health state
### occupied by the individuals at cycle 't' and the individual characteristics
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

## Calculate Health Outcomes function
### This function estimates the Quality Adjusted Life Years (QALYs) at every
### cycle based on the health state occupied by each individuals at cycle 't',
### time spent in the states and the cycle_length (measured in years)
calc_effsV <- function (
    v_occupied_state, 
    v_states_utilities,
    m_indi_features,
    v_util_coeffs, 
    v_util_t_decs,
    v_time_in_state,
    cycle_length = 1) {
  
  # calculate individual-specific utility decrements based on utilities regression coefficients
  v_ind_decrement <- (m_indi_features %*% v_util_coeffs)[,1]
  
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

## Calculate Discount Weights function
### This function estimates the discount weights to be applied to the outputs of
### each cycle in order to scale their values back to the present. 
calc_discount_wts <- function(
    discount_rate,
    num_cycles,
    cycle_length) {
  
  # calculate discount weights based on the number & length (in years) of cycles
  v_discount_wts <- 1 / (1 + discount_rate) ^ ((0:num_cycles) * cycle_length)
  
  return(v_discount_wts)
}

## Run Microsimulation function
### This function runs the microsimulation function of the Healthy-Sick-Sicker-Dead
### model
run_microSimV <- function(
    v_starting_states,
    num_i,
    num_cycles,
    m_indi_features,
    v_states_names,
    v_states_costs,
    v_cost_coeffs,
    v_states_utilities,
    v_util_coeffs,
    v_util_t_decs,
    l_trans_probs,
    p_HS1 = 0.15,
    p_S1S2 = 0.105,
    discount_rate_costs,
    discount_rate_QALYs,
    cycle_length = 1,
    starting_seed = 1) {
  
  # combine the two probability parameters with the others in 'l_trans_probs'
  l_trans_probs <- c(l_trans_probs, "p_HS1" = p_HS1, "p_S1S2" = p_S1S2)
  
  # create matrices to capture states' names, associated costs and QALYs
  m_States <- m_Costs <- m_Effs <-  matrix(
    nrow = num_i,
    ncol = num_cycles + 1,
    dimnames = list(paste("ind",   1:num_i,    sep ="_"),
                    paste("cycle", 0:num_cycles, sep ="_"))
  )
  
  # set the seed for reproducibility - R parallel processing uses "L'Ecuyer-CMRG"
  set.seed(starting_seed, kind = "L'Ecuyer-CMRG")
  
  # initialize parameter tracking time in current state
  v_time_in_state <- rep(0, times = num_i)
  
  # get the initial health state
  m_States[, 1] <- v_starting_states
  
  # calculate the costs incurred in their starting health state
  m_Costs[, 1]  <- calc_costsV(
    v_occupied_state = m_States[, 1],
    v_states_costs   = v_states_costs,
    m_indi_features  = m_indi_features,
    v_cost_coeffs    = v_cost_coeffs
  )
  
  # calculate the QALYs accrued in their starting health state
  m_Effs[, 1]   <- calc_effsV(
    v_occupied_state   = m_States[, 1],
    v_states_utilities = v_states_utilities,
    m_indi_features    = m_indi_features,
    v_util_coeffs      = v_util_coeffs,
    v_util_t_decs      = v_util_t_decs,
    v_time_in_state    = v_time_in_state,
    cycle_length       = cycle_length
  )
  
  # for each 't' of the 'num_cycles' cycles:
  for (t in 1:num_cycles) {
    # update the transition probabilities at cycle 't'
    m_trans_probs     <- update_probsV(
      v_states_names   = v_states_names,
      v_occupied_state = m_States[, t],
      l_trans_probs    = l_trans_probs,
      v_time_in_state  = v_time_in_state
    )
    
    # sample the health state at 't + 1'
    m_States[, t + 1] <- sampleV(
      m_trans_probs  = m_trans_probs,
      v_states_names = v_states_names
    )
    
    # keep track of time in state at 't + 1'
    stayed                   <- m_States[, t] == m_States[, t + 1] # check if remains in current state at 't + 1'
    v_time_in_state[stayed]  <- v_time_in_state[stayed] + 1        # increment time spent in state
    v_time_in_state[!stayed] <- 1                                  # reset time once transitioned
    
    # keep track of time in the model
    m_indi_features[, "age"] <- m_indi_features[, "age"] + 1
    
    # calculate the costs incurred in their 't + 1' health state
    m_Costs[, t + 1]  <- calc_costsV(
      v_occupied_state = m_States[, t + 1],
      v_states_costs   = v_states_costs,
      m_indi_features  = m_indi_features,
      v_cost_coeffs    = v_cost_coeffs
    )
    
    # calculate the QALYs accrued in their 't + 1' health state
    m_Effs[, t + 1]   <- calc_effsV(
      v_occupied_state   = m_States[, t + 1],
      v_states_utilities = v_states_utilities,
      m_indi_features    = m_indi_features,
      v_util_coeffs      = v_util_coeffs,
      v_util_t_decs      = v_util_t_decs,
      v_time_in_state    = v_time_in_state,
      cycle_length       = cycle_length
    )
    
  } # close the loop for the cycles 't'
  
  # Calculate discount weights for both outcomes:
  v_c_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_costs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  v_e_dsc_wts <- calc_discount_wts(
    discount_rate = discount_rate_QALYs,
    num_cycles    = num_cycles,
    cycle_length  = cycle_length
  )
  # Compute costs and QALYs:
  v_total_costs <- rowSums(m_Costs)         # calculate total costs per individual
  v_total_qalys <- rowSums(m_Effs)          # calculate total QALYs per individual
  mean_costs    <- mean(v_total_costs)      # calculate average costs
  mean_qalys    <- mean(v_total_qalys)      # calculate average QALYs
  
  # Compute discounted costs and QALYs:
  v_total_Dcosts <- m_Costs %*% v_c_dsc_wts # calculate total discounted costs per individual
  v_total_Dqalys <- m_Effs  %*% v_e_dsc_wts # calculate total discounted QALYs per individual
  mean_Dcosts    <- mean(v_total_Dcosts)    # calculate average discounted costs
  mean_Dqalys    <- mean(v_total_Dqalys)    # calculate average discounted QALYs
  
  # store the results in a list:
  results <- list(
    m_States       = m_States,
    m_Costs        = m_Costs,
    m_Effs         = m_Effs,
    v_total_costs  = v_total_costs,
    v_total_qalys  = v_total_qalys,
    v_total_Dcosts = v_total_Dcosts,
    v_total_Dqalys = v_total_Dqalys,
    mean_costs     = mean_costs,
    mean_qalys     = mean_qalys,
    mean_Dcosts    = mean_Dcosts,
    mean_Dqalys    = mean_Dqalys
  )
  
  # return the results
  return(results)
}

## Sample PSA parameters
### This function samples PSA values from parameters prior distributions to be
### used in PSA
sample_psa_data <- function(
    psa_params_names,
    psa_params_dists,
    psa_params_dists_args,
    n_sims,
    starting_seed = 1) {
  
  # Prepare inputs:
  loop_list <- list(
    as.list(psa_params_names),
    as.list(psa_params_dists),
    psa_params_dists_args
  )
  
  # set the seed for reproducibility
  set.seed(starting_seed)
  
  # Generate PSA configurations:
  df_psa <- purrr::pmap(
    .l = loop_list,
    .f = function(name_, dist_, args_) {
      func <- paste0("r", dist_)
      
      # sample from user-specified distribution
      purrr::exec(
        .fn = func,
        n_sims,
        !!!args_
      ) |>
        as.data.frame() |>
        `colnames<-`(name_)
    }
  ) |>
    purrr::list_cbind()
  
  return(df_psa)
}

## Run Microsimulation PSA function
### This function performs PSA for the microsimulation function of the 
### Healthy-Sick-Dead model
run_psa <- function(
    f_model,
    l_model_inputs,
    l_psa_parameters,
    v_extracted_results,
    n_sims,
    starting_seed = 1) {
  # sample PSA parameters' configurations:
  psa_params <- sample_psa_data(
    psa_params_names = l_psa_parameters$psa_params_names,
    psa_params_dists = l_psa_parameters$psa_params_dists,
    psa_params_dists_args = l_psa_parameters$psa_params_dists_args,
    n_sims = n_sims,
    starting_seed = starting_seed
  )
  
  # run the model  
  psa_results <- purrr::pmap(
    .l = as.list(psa_params),
    .f = function(...) {
      # grab a set of sampled values
      v_params <- c(...)
      
      # combine with other (non-PSA) parameters passed to the function
      l_params <- c(
        l_model_inputs,
        v_params
      )
      # run the model given the sampled configuration
      l_model_results <- purrr::exec(
        .fn = f_model,
        !!!l_params
      )
      # extract required results
      df_psa_run_results <- if(is.null(v_extracted_results)) {
        c(v_params, l_model_results)
      } else {
        c(v_params, l_model_results[v_extracted_results])
      } |>
        as.data.frame()
    },
    .progress = list(name = "Running PSA")
  ) |>
    purrr::list_rbind()
  
  return(psa_results)
}

## Set-up Parallel Processing function
### This function allows the user to set-up the parallel computations unit in 
### which the PSA function will be processed.
set_parallel <- function(
    num_workers = NULL,
    parallel_method = "multisession") {
  
  # set-up parallel processing unit
  if(is.null(num_workers) | !is.numeric(num_workers)) {
    future::plan(future::sequential)
  } else {
    if(Sys.info()[['sysname']] == "Windows") {
      future::plan(
        future::multisession, 
        workers = num_workers
      )
    } else if (any(Sys.info()[['sysname']] == c("Linux", "Darwin"))) {
      switch(
        EXPR = parallel_method,
        multisession = future::plan(
          future::multisession, 
          workers = num_workers
        ),
        multicore = future::plan(
          future::multicore, 
          workers = num_workers
        ),
        cluster = future::plan(
          future::cluster, 
          workers = num_workers
        )
      )
    } else {
      stop("Operating system unknow!")
    }
  }
}

## Specify Parallel Processing chunks function
### This function helps split the simulation of the sampled PSA configurations
### between available parallel processing workers.
make_psa_chunks <- function(
    n_sims,
    num_workers   = NULL,
    chunk_size    = NULL,
    chunk_type    = "data_indices",
    df_psa_params = NULL) {
  
  # define parallel processing chunks
  chunks <- if(any(chunk_type != "data_sets" | !is.data.frame(df_psa_params))) {
    if(is.null(num_workers)) {
      # return a vector of all indices
      1:n_sims
    } else {
      if(is.null(chunk_size)) {
        # create num_workers vectors all of equal sizes
        furrr:::make_chunks(
          n_x = n_sims,
          n_workers = num_workers,
          chunk_size = n_sims/num_workers
        )
      } else {
        # create num_workers each of size chunk_size
        furrr:::make_chunks(
          n_x = n_sims,
          n_workers = num_workers,
          chunk_size = chunk_size
        ) 
      }
    }
  } else {
    if(is.null(num_workers)) {
      # return the df_psa_params
      df_psa_params
    } else {
      if(is.null(chunk_size)) {
        # split the df_psa_params to num_workers datasets all of equal sizes 
        split(
          x = df_psa_params, 
          f = rep(
            x = 1:num_workers, 
            each = nrow(df_psa_params) / num_workers, 
            length.out = n_sims
          )
        )
      } else {
        # split the df_psa_params to num_workers datasets each of size chunk_size
        split(
          x = df_psa_params, 
          f = rep(
            x = 1:(num_workers * n_sims), 
            each = chunk_size, 
            length.out = n_sims
          )
        )
      }
    }
  }

  return(chunks)
}

## Run Microsimulation PSA in Parallel or Sequential function
### This function performs PSA for the microsimulation function of the 
### Healthy-Sick-Dead model sequentially or in parallel.
run_psa_parallel <- function(
    f_model,
    l_model_inputs,
    l_psa_parameters,
    v_extracted_results,
    n_sims,
    num_workers = NULL,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1) {
  
  # sample PSA parameters' configurations:
  df_psa_params <- sample_psa_data(
    psa_params_names = l_psa_parameters$psa_params_names,
    psa_params_dists = l_psa_parameters$psa_params_dists,
    psa_params_dists_args = l_psa_parameters$psa_params_dists_args,
    n_sims = n_sims,
    starting_seed = starting_seed
  )
  
  # set-up parallel processing
  set_parallel(
    num_workers = num_workers,
    parallel_method = parallel_method
  )
  
  cat("Number of workers: ", future::nbrOfWorkers(), "\n")
  
  # identify parallel processing chunks
  psa_chunks <- make_psa_chunks(
    n_sims = n_sims,
    num_workers = num_workers,
    chunk_size = chunk_size
  )
  
  # run the PSA in parallel  
  psa_results <- furrr::future_map(
    .x = psa_chunks,
    .f = function(chunk_i) {
      # perform the PSA at the worker level
      purrr::pmap(
        .l = as.list(df_psa_params[chunk_i, ]),
        .f = function(...) {
          # grab a set of sampled values
          v_params <- c(...)
          # combine with other (non-PSA) parameters passed to the function
          l_params <- c(
            l_model_inputs,
            v_params
          )
          # run the model given the sampled configuration
          l_model_results <- purrr::exec(
            .fn = f_model,
            !!!l_params
          )
          # extract required results
          df_psa_run_results <- if(is.null(v_extracted_results)) {
            c(v_params, l_model_results)
          } else {
            c(v_params, l_model_results[v_extracted_results])
          } |>
            as.data.frame()
        }
      ) |>
        purrr::list_rbind()
    },
    # set seed = NULL produces identical results interactively and TRUE otherwise:
    .options = furrr::furrr_options(
      seed = if(base::interactive()) {
        NULL
      } else {
        TRUE
      }
    )
  ) |>
    purrr::list_rbind()
  
  return(psa_results)
}

## Run Microsimulation PSA in Parallel using the parallel package
run_psa_parallel2 <- function(
    f_model,
    l_model_inputs,
    l_psa_parameters,
    v_extracted_results,
    n_sims,
    num_workers = NULL,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1) {
  
  # Sample PSA parameters' configurations
  df_psa_params <- sample_psa_data(
    psa_params_names = l_psa_parameters$psa_params_names,
    psa_params_dists = l_psa_parameters$psa_params_dists,
    psa_params_dists_args = l_psa_parameters$psa_params_dists_args,
    n_sims = n_sims,
    starting_seed = starting_seed
  )
  
  # identify parallel processing chunks
  psa_chunks <- make_psa_chunks(
    n_sims = n_sims,
    num_workers = num_workers,
    chunk_size = chunk_size,
    chunk_type = "data_sets",
    df_psa_params = df_psa_params
  )
  
  # Setup parallel environment
  cl <- parallel::makeCluster(spec = num_workers)
  parallel::clusterSetRNGStream(
    cl = cl, 
    iseed = starting_seed
  )
  if (parallel_method == "multisession") {
    # dynamic identification of global environment functions relevant to f_model
    global_env_functions <- codetools::findGlobals(f_model)[
      codetools::findGlobals(f_model) %in% ls(envir = globalenv())
    ]
    parallel::clusterExport(
      cl = cl, 
      varlist = list("f_model", "l_model_inputs", "v_extracted_results"), 
      envir = environment()
    )
    parallel::clusterExport(
      cl = cl, 
      varlist = as.list(global_env_functions), 
      envir = globalenv()
    )
  }
  
  # Define the worker function to run each simulation
  worker_function <- function(psa_params) {
    # loop over psa_params rows
    purrr::pmap(
      .l = as.list(psa_params),
      .f = function(...) {
        # grab a set of sampled values
        v_params <- c(...)
        
        # combine with other (non-PSA) parameters passed to the function
        l_params <- c(
          l_model_inputs,
          v_params
        )
        # run the model given the sampled configuration
        l_model_results <- purrr::exec(
          .fn = f_model,
          !!!l_params
        )
        # extract required results
        df_psa_run_results <- if(is.null(v_extracted_results)) {
          c(v_params, l_model_results)
        } else {
          c(v_params, l_model_results[v_extracted_results])
        } |>
          as.data.frame()
      }
    ) |>
      purrr::list_rbind()
  }
  
  # Run the PSA simulations in parallel
  psa_results <- if (parallel_method == "multisession") {
    parallel::parLapply(
      cl = cl, 
      X = psa_chunks, 
      worker_function
    )
  } else if (parallel_method == "multicore") {
    parallel::mclapply(
      X = psa_chunks, 
      FUN = worker_function, 
      mc.cores = num_workers
    )
  } else {
    stop("User-defined Parallel method is not supported by the")
  }
  
  # Combine results into a single data frame
  psa_results <- psa_results |>
    purrr::list_rbind()
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  return(psa_results)
}
#------------------------------------------------------------------------------#

                        ### Defining model parameters ###

# Define model inputs
## General parameters
num_i               <- 2e4               # number of simulated individuals
num_cycles          <- 30                # time horizon if each cycle is a year long
cycle_length        <- 1                 # length of cycle (in years)
seed                <- 1234              # random number generator state
wtp                 <- 30000             # Willingness to pay for each QALY ($)
discount_rate_costs <- 0.03              # annual discount rate for costs
discount_rate_QALYs <- 0.015             # annual discount rate for health outcomes

## Population characteristics/features
mean_age            <- 50                # mean age in the simulated population
sd_age              <- 3                 # standard deviation of the age in the simulated population
prop_females        <- 0.6               # proportion of females in the simulated population
prop_males          <- 1 - prop_females  # proportion of males in the simulated population
set.seed(seed)                           # set a seed to ensure reproducible samples
m_indi_features     <- cbind(            # simulate individuals characteristics
  "age" = rnorm(                         # get random samples for 'age' from a normal distribution
    n = num_i, 
    mean = mean_age, 
    sd = sd_age
  ),
  "sex" = sample(                        # get random samples for 'sex' based on sex distribution
    x = c(0, 1), 
    size = num_i, 
    replace = TRUE, 
    prob = c(prop_females, prop_males)
  )
)

## Health states
v_states_names <- c("H","S1", "S2", "D") # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
v_starting_states <- rep("H", num_i)     # everyone begins in the healthy state 

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
rp_S1  <- 0.2                            # increase in mortality rate with every additional year being sick
rp_S2  <- 0.29                           # increase in mortality rate with every additional year being sicker

l_trans_probs <- list(                   # pack the transition probabilities and rates in a list
  "p_HD"   = p_HD, 
  "p_S1H"  = p_S1H, 
  "p_S1D"  = p_S1D, 
  "p_S2D"  = p_S2D, 
  "rp_S1"  = rp_S1, 
  "rp_S2"  = rp_S2
)

## Cost and utility inputs 
c_H       <- 2000         # cost of remaining one cycle healthy
c_S1      <- 4000         # cost of remaining one cycle sick
c_S2      <- 15000        # cost of remaining one cycle sicker
c_S1_Trt1 <- c_S1 + 12000 # cost of remaining one cycle sick under treatment 1
c_S2_Trt1 <- c_S2 + 12000 # cost of remaining one cycle sicker under treatment 1
c_S1_Trt2 <- c_S1 + 11350 # cost of remaining one cycle sick under treatment 2
c_S2_Trt2 <- c_S2 + 11350 # cost of remaining one cycle sicker under treatment 2
c_D       <- 0            # cost associated with being dead
c_age_cof <- 11.5         # cost age coefficient
c_sex_cof <- 300          # cost sex coefficient, where 0 is female and 1 is male

v_cost_coeffs <- c(       # pack the cost regression coefficients in a vector
  "age" = c_age_cof, "sex" = c_sex_cof
)

u_H       <- 1            # utility when healthy 
u_S1      <- 0.75         # utility when sick, untreated
u_S2      <- 0.5          # utility when sicker, untreated
u_S1_Trt1 <- u_S1 + 0.2   # utility when sick, treatment 1
u_S2_Trt1 <- u_S2         # utility when sicker, treatment 1
u_S1_Trt2 <- u_S1 + 0.15  # utility when sick, treatment 2
u_S2_Trt2 <- u_S2 + 0.05  # utility when sicker, treatment 2
u_D       <- 0            # utility when dead
u_age_cof <- -0.0018      # utility age coefficient
u_sex_cof <- -0.015       # utility sex coefficient, where 0 is female and 1 is male
ru_S1     <- -0.0015      # change in utility of individuals with every additional year being sick
ru_S2     <- -0.0020      # change in utility of individuals with every additional year being sicker

v_util_coeffs <- c(       # pack the utility regression coefficients in a vector
  "age" = u_age_cof, "sex" = u_sex_cof
)
v_util_t_decs <- c(       # pack the state-specific utility decrements in a vector
  "S1" = ru_S1, "S2" = ru_S2
)

### Payoffs - no treatment
v_states_costs     <- c("H" = c_H, "S1" = c_S1, "S2" = c_S2, "D" = c_D)
v_states_utilities <- c("H" = u_H, "S1" = u_S1, "S2" = u_S2, "D" = u_D)

# PSA parameters
l_psa_parameters <- list(
  "psa_params_names" = c("p_HS1", "p_S1S2"),
  "psa_params_dists" = c("p_HS1" = "unif", "p_S1S2" = "beta"),
  "psa_params_dists_args" = list(
    "p_HS1" = list(
      "min" = 0.11,
      "max" = 0.19
    ),
    "p_S1S2" = list(
      "shape1" = 1.91,
      "shape2" = 16.27
    )
  )
)

# ---------------------------------------------------------------------------- #

# Compare results:
print("Comparing results of run_PSA() and all 3 run_psa_parallel() varients:")
## run_psa():
r0 <- run_psa(
  f_model             = run_microSimV,
  l_model_inputs      = list(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  l_psa_parameters    = l_psa_parameters,
  v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
  n_sims              = 20,
  starting_seed       = seed
)

## run_psa_parallel:
### sequential 
r1 <- run_psa_parallel(
  f_model             = run_microSimV,
  l_model_inputs      = list(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  l_psa_parameters    = l_psa_parameters,
  v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
  n_sims              = 20,
  starting_seed       = seed
)
### multisession:
r2 <- run_psa_parallel(
  f_model             = run_microSimV,
  l_model_inputs      = list(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  l_psa_parameters    = l_psa_parameters,
  v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
  n_sims              = 20,
  num_workers         = 2,
  parallel_method     = "multisession",
  starting_seed       = seed
)
### multicore:
if(all(!base::interactive(), Sys.info()[['sysname']] != "Windows")) {
  r3 <- run_psa_parallel(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    num_workers         = 2,
    parallel_method     = "multicore",
    starting_seed       = seed
  )
}

print("Check if all four results, run_microSim() and run_psa_parallel() veriants:")
identical(r0, r1)
identical(r1, r2)
if(all(!base::interactive(), Sys.info()[['sysname']] != "Windows")) {
  identical(r2, r3)
}

# Testing parallel computing methods:
print("Without future - 20 sims:")
system.time({
  run_psa(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20
  )
})
#   user  system elapsed
# 15.715   0.838  16.633
print("Sequential   - ? workers - 20 sims:")
s0 <- system.time({
  run_psa_parallel(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    starting_seed       = seed
  )
})
s0
#   user  system elapsed
# 16.104   0.817  17.000
print("Multisession - 2 workers - 20 sims:")
s1 <- system.time({
  run_psa_parallel(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    num_workers         = 2,
    parallel_method     = "multisession",
    starting_seed       = seed
  )
})
s1
#  user  system elapsed
# 0.222   0.042   9.959
print("Multisession - 4 workers - 20 sims:")
s3 <- system.time({
  run_psa_parallel(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    num_workers         = 4,
    parallel_method     = "multisession",
    starting_seed       = seed
  )
})
s3
#  user  system elapsed
# 0.239   0.029   6.070
print("Multisession - 4 workers - 100 sims:")
s5 <- system.time({
  run_psa_parallel(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 100,
    num_workers         = 4,
    parallel_method     = "multisession",
    starting_seed       = seed
  )
})
s5
#  user  system elapsed
# 0.683   0.083  27.104

# If on linux or Mac OS, and running the code from the terminal
if(all(!base::interactive(), Sys.info()[['sysname']] != "Windows")) {
  print("Multicore    - 2 workers - 20 sims:")
  s2 <- system.time({
    run_psa_parallel(
      f_model             = run_microSimV,
      l_model_inputs      = list(
        v_starting_states   = v_starting_states,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features,
        v_states_names      = v_states_names,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      l_psa_parameters    = l_psa_parameters,
      v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
      n_sims              = 20,
      num_workers         = 2,
      parallel_method     = "multicore",
      starting_seed       = seed
    )
  })
  print(s2)
  #  user  system elapsed
  # 0.112   0.104   9.331
  print("Multicore    - 4 workers - 20 sims:")
  s4 <- system.time({
    run_psa_parallel(
      f_model             = run_microSimV,
      l_model_inputs      = list(
        v_starting_states   = v_starting_states,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features,
        v_states_names      = v_states_names,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      l_psa_parameters    = l_psa_parameters,
      v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
      n_sims              = 20,
      num_workers         = 4,
      parallel_method     = "multicore",
      starting_seed       = seed
    )
  })
  print(s4)
  #  user  system elapsed
  # 0.085   0.071   5.064
  print("Multicore    - 4 workers - 100 sims:")
  s6 <- system.time({
    run_psa_parallel(
      f_model             = run_microSimV,
      l_model_inputs      = list(
        v_starting_states   = v_starting_states,
        num_i               = num_i,
        num_cycles          = num_cycles,
        m_indi_features     = m_indi_features,
        v_states_names      = v_states_names,
        v_states_costs      = v_states_costs,
        v_cost_coeffs       = v_cost_coeffs,
        v_states_utilities  = v_states_utilities,
        v_util_coeffs       = v_util_coeffs,
        v_util_t_decs       = v_util_t_decs,
        l_trans_probs       = l_trans_probs,
        discount_rate_costs = discount_rate_costs,
        discount_rate_QALYs = discount_rate_QALYs,
        cycle_length        = cycle_length,
        starting_seed       = seed
      ),
      l_psa_parameters    = l_psa_parameters,
      v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
      n_sims              = 100,
      num_workers         = 4,
      parallel_method     = "multicore",
      starting_seed       = seed
    )
  })
  print(s6)
  #  user  system elapsed
  # 0.283   0.256  25.386
}

## Parallel function version
r4 <- run_psa_parallel2(
  f_model             = run_microSimV,
  l_model_inputs      = list(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  l_psa_parameters    = l_psa_parameters,
  v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
  n_sims              = 20,
  num_workers         = 2,
  parallel_method     = "multisession",
  starting_seed       = seed
)
r5 <- run_psa_parallel2(
  f_model             = run_microSimV,
  l_model_inputs      = list(
    v_starting_states   = v_starting_states,
    num_i               = num_i,
    num_cycles          = num_cycles,
    m_indi_features     = m_indi_features,
    v_states_names      = v_states_names,
    v_states_costs      = v_states_costs,
    v_cost_coeffs       = v_cost_coeffs,
    v_states_utilities  = v_states_utilities,
    v_util_coeffs       = v_util_coeffs,
    v_util_t_decs       = v_util_t_decs,
    l_trans_probs       = l_trans_probs,
    discount_rate_costs = discount_rate_costs,
    discount_rate_QALYs = discount_rate_QALYs,
    cycle_length        = cycle_length,
    starting_seed       = seed
  ),
  l_psa_parameters    = l_psa_parameters,
  v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
  n_sims              = 20,
  num_workers         = 2,
  parallel_method     = "multicore",
  starting_seed       = seed
)
print("Compare Parallel package execution results with run_psa():")
identical(r0, r4)
identical(r4, r5)

# Testing parallel computing methods:
print("multisession - 2 workers - 20 sims - Parallel package:")
s7 <- system.time({
  run_psa_parallel2(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    num_workers         = 2,
    parallel_method     = "multisession",
    starting_seed       = seed
  )
})
print(s7)
#  user  system elapsed
# 0.026   0.006   9.618
print("Multicore    - 2 workers - 20 sims - Parallel package:")
s8 <- system.time({
  run_psa_parallel2(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 20,
    num_workers         = 2,
    parallel_method     = "multicore",
    starting_seed       = seed
  )
})
print(s8)
#  user  system elapsed
# 8.752   0.087   9.024
print("multisession - 4 workers - 100 sims - Parallel package:")
s9 <- system.time({
  run_psa_parallel2(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 100,
    num_workers         = 4,
    parallel_method     = "multisession",
    starting_seed       = seed
  )
})
print(s9)
print("Multicore    - 4 workers - 100 sims - Parallel package:")
s10 <- system.time({
  run_psa_parallel2(
    f_model             = run_microSimV,
    l_model_inputs      = list(
      v_starting_states   = v_starting_states,
      num_i               = num_i,
      num_cycles          = num_cycles,
      m_indi_features     = m_indi_features,
      v_states_names      = v_states_names,
      v_states_costs      = v_states_costs,
      v_cost_coeffs       = v_cost_coeffs,
      v_states_utilities  = v_states_utilities,
      v_util_coeffs       = v_util_coeffs,
      v_util_t_decs       = v_util_t_decs,
      l_trans_probs       = l_trans_probs,
      discount_rate_costs = discount_rate_costs,
      discount_rate_QALYs = discount_rate_QALYs,
      cycle_length        = cycle_length,
      starting_seed       = seed
    ),
    l_psa_parameters    = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims              = 100,
    num_workers         = 4,
    parallel_method     = "multicore",
    starting_seed       = seed
  )
})
print(s10)