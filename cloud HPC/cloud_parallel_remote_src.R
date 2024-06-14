## @knitr micro_cloud_parallel_remote_on_Linode_demo_src

# This script sources the 
# Source the Healthy-Sick-Sicker-Dead microsimulation model dependencies functions
source(
  file = file.path(here::here(), "code_chunks", "11_microsim", "parallel_terminal_src.R")
)

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
### Healthy-Sick-Sicker-Dead model sequentially or in parallel.
run_psa_parallel <- function(
    f_model = run_microSimV,
    l_model_inputs = list(
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
    l_psa_parameters = l_psa_parameters,
    v_extracted_results = c("mean_Dcosts", "mean_Dqalys"),
    n_sims = 1e4,
    num_workers = NULL,
    parallel_method = "multisession",
    chunk_size = NULL,
    starting_seed = 1,
    cluster = FALSE,
    num_clusters = 4) {
  
  # sample PSA parameters' configurations:
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
    num_workers = ifelse(isTRUE(cluster), num_clusters, num_workers),
    chunk_size = chunk_size,
    chunk_type = ifelse(isTRUE(cluster), "data_sets", "data_indices"),
    df_psa_params = df_psa_params
  )
  
  # Conditional statement based on whether clusters are in play
  psa_results <- if(isTRUE(cluster)) {
    
    # loop over cluster
    furrr::future_map(
      .x = psa_chunks,
      .f = function(df_chunk_i) {
        
        # identify parallel processing chunks
        inner_psa_chunks <- make_psa_chunks(
          n_sims = nrow(df_chunk_i),
          num_workers = num_workers,
          chunk_size = chunk_size,
          chunk_type = "data_indices",
          df_psa_params = df_chunk_i
        )
        
        # run the PSA in parallel, multisession or multicore 
        furrr::future_map(
          .x = inner_psa_chunks,
          .f = function(chunk_i) {
            
            # perform the PSA at the worker level
            purrr::pmap(
              .l = as.list(df_chunk_i[chunk_i, ]),
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
    
  } else {
   
    # run the PSA in parallel, multisession or multicore 
    furrr::future_map(
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
    
  }
  
  return(psa_results)
}

## Connect to Linode instance(s)
### Creates a connection to a Linode cluster to allow code from a local R
### session to be processed remotely
connect_to_linode <- function(public_ip, ssh_private_key_file) {
  future::makeClusterPSOCK(
    
    # Public IP number of Linode instance
    workers = public_ip,
    
    # User name (could be 'ubuntu')
    user = "root",
    
    # Use private SSH key registered with Linode
    rshopts = c(
      "-o", "StrictHostKeyChecking=no",
      "-o", "IdentitiesOnly=yes",
      "-i", ssh_private_key_file
    ),
    
    # Switch this to TRUE to see the code that is run on the workers without
    # making the connection
    dryrun = FALSE
  )
}