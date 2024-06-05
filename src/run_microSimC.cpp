#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

///////////////// run_microSimC1 & run_microSimC2 functions ///////////////////

// update_probsC6 (function name in update_probsC.cpp)
// 
// Compared to the other C++ candidate functions, 'C6' uses arma::cube (array),
// arma::ivec to pass vectors of type integer (to saves memory) and uses state
// indexes instead of names.
//
// [[Rcpp::export]]
arma::mat update_probsC(const arma::ivec& v_states_index,
                        const arma::ivec& v_occupied_state,
                        const Rcpp::List& l_trans_probs,
                        const arma::vec& v_time_in_state) {
  
  // initialize two integers to hold number of individuals and number of states
  int num_i = v_time_in_state.size();
  int num_states = v_states_index.size();
  
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_probs(num_i, num_states, arma::fill::zeros);
  
  // update probabilities of death 
  // first converting probabilities to rates
  double r_S1D = -std::log(1 - Rcpp::as<double>(l_trans_probs["p_S1D"]));
  double r_S2D = -std::log(1 - Rcpp::as<double>(l_trans_probs["p_S2D"]));
  
  // applying the rate ratio and converting rates back to probabilities
  arma::vec p_S1D = 1 - exp(-r_S1D * (1 + v_time_in_state * Rcpp::as<double>(l_trans_probs["rp_S1"])));
  arma::vec p_S2D = 1 - exp(-r_S2D * (1 + v_time_in_state * Rcpp::as<double>(l_trans_probs["rp_S2"])));
  
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube a_trans_probs(num_i, num_states, num_states);
  
  // assign probabilities to the slices
  // probabilities moving from the first health state "H"
  double p_HS1 = Rcpp::as<double>(l_trans_probs["p_HS1"]);
  double p_HD = Rcpp::as<double>(l_trans_probs["p_HD"]);
  double p_HH = 1 - (Rcpp::as<double>(l_trans_probs["p_HS1"]) + 
                     Rcpp::as<double>(l_trans_probs["p_HD"]));
  a_trans_probs.slice(0).col(0).fill(p_HH);
  a_trans_probs.slice(0).col(1).fill(p_HS1);
  a_trans_probs.slice(0).col(3).fill(p_HD);
  
  // probabilities moving from the second health state "S1"
  double p_S1H = Rcpp::as<double>(l_trans_probs["p_S1H"]);
  double p_S1S2 = Rcpp::as<double>(l_trans_probs["p_S1S2"]);
  a_trans_probs.slice(1).col(0).fill(p_S1H);
  a_trans_probs.slice(1).col(2).fill(p_S1S2);
  a_trans_probs.slice(1).col(3) = p_S1D;
  a_trans_probs.slice(1).col(1) = 1 - arma::sum(a_trans_probs.slice(1), 1);
  
  // probabilities moving from the third health state "S2"
  a_trans_probs.slice(2).col(3) = p_S2D;
  a_trans_probs.slice(2).col(2) = 1 - arma::sum(a_trans_probs.slice(2), 1);  
  
  // probabilities moving from the dead health state "D"
  a_trans_probs.slice(3).col(3).fill(1);
  
  // assign probabilities from the slices based on occupied health state 
  for(int i = 0; i < num_states; i++) {
    arma::uvec state_index = arma::find(v_occupied_state == (i + 1));
    m_probs.rows(state_index) = a_trans_probs.slice(i).rows(state_index);
  }
  
  // Sanity check
  arma::vec row_sums = sum(m_probs, 1);
  if (arma::any(arma::abs(row_sums - 1) > 1e-12)) {
    Rcpp::stop("Probabilities do not sum to 1");
  }
  
  return m_probs;
}

// sampleC5 (function name in sampleC.cpp)
// 
// Compared to the other C++ candidate functions, 'C5' returns arma::ivec,
// uses the arma::cumsum() function instead of Rcpp and
// 'arma::trimatu * m_trans_probs' to estimate the cumulative probabilities.
//
// [[Rcpp::export]]
arma::ivec sampleC(const arma::mat& m_trans_probs) {
  
  // Number of individuals and states 
  int num_i = m_trans_probs.n_rows;
  int num_states = m_trans_probs.n_cols;
  
  // Create matrix with row-wise cumulative transition probabilities
  arma::mat m_cum_probs = arma::cumsum(m_trans_probs, 1);
  
  // Ensure that the maximum cumulative probabilities are equal to 1
  if (arma::any(m_cum_probs.col(num_states - 1) > 1.000000)) {
    Rcpp::stop("Error in multinomial sampling: probabilities do not sum to 1");
  }
  
  // Sample random values from Uniform standard distribution for each individual
  arma::vec v_rand_values = arma::randu<arma::vec>(num_i);
  
  // Repeat each sampled value to have as many copies as the number of states
  arma::mat m_rand_values = arma::repmat(v_rand_values, 1, num_states);
  
  // Identify transitions, compare random samples to cumulative probabilities
  arma::umat m_transitions = m_rand_values > m_cum_probs;
  
  // Sum transitions to identify health state in the next cycle
  arma::uvec v_transitions = arma::sum(m_transitions, 1);
  
  // Identify health state to which each individual is transitioning
  arma::ivec v_health_states = 1 + arma::conv_to<arma::ivec>::from(v_transitions);
  
  return v_health_states;
}

// calc_costsC7 (function name in calc_costsC.cpp)
//
// Compared to the other C++ candidate functions, 'C5' removes loop and IF
// statements in assigning costs, does not use the arma::vec '.elem()' method to
// subset vectors of type arma::vec, nor arma::vec '.fill' method to populate a
// vector with the value of a scalar.
//
// [[Rcpp::export]]
arma::vec calc_costsC(const arma::ivec& v_occupied_state,
                      const arma::vec& v_states_costs,
                      const arma::mat& m_indi_features,
                      const arma::vec& v_cost_coeffs) {
  
  // Calculate individual-specific costs based on costs regression coefficients
  arma::vec v_indi_costs = m_indi_features * v_cost_coeffs;
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Estimate costs based on occupied state
  arma::vec v_state_costs(num_i);
  for (int i = 0; i < num_i; ++i) {
    if (v_occupied_state[i] == 1) {
      v_state_costs[i] = v_states_costs[0];
    } else if (v_occupied_state[i] == 2) {
      v_state_costs[i] = v_states_costs[1] + v_indi_costs[i];
    } else if (v_occupied_state[i] == 3) {
      v_state_costs[i] = v_states_costs[2] + v_indi_costs[i];
    } else if (v_occupied_state[i] == 4) {
      v_state_costs[i] = v_states_costs[3];
    }
  }
  
  return v_state_costs;
}

// calc_effsC6 (function name in calc_effsC.cpp)
//
// uses for loop and IF statement to assign costs and a for loop to calculate
// the cost regression estimate for each individual and passing objects by
// reference (&)
//
// [[Rcpp::export]]
arma::vec calc_effsC(const arma::ivec& v_occupied_state,
                     const arma::vec& v_states_utilities,
                     const arma::mat& m_indi_features,
                     const arma::vec& v_util_coeffs,
                     const arma::vec& v_util_t_decs,
                     const arma::vec& v_time_in_state,
                     const double cycle_length = 1.0) {
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Calculate individual-specific utility decrements based on utilities regression coefficients
  arma::vec v_ind_decrement = m_indi_features * v_util_coeffs;
  
  // Calculate time-dependent state-specific utility decrements
  arma::vec time_decrement(num_i);
  for (int i = 0; i < num_i; ++i) {
    if (v_occupied_state[i] == 2) {
      time_decrement[i] = v_util_t_decs[0] * v_time_in_state[i];
    } else if (v_occupied_state[i] == 3) {
      time_decrement[i] = v_util_t_decs[1] * v_time_in_state[i];
    }
  }
  
  // Estimate total decrements
  arma::vec decrement = v_ind_decrement + time_decrement;
  
  // Estimate utilities based on occupied state
  arma::vec v_state_utility(num_i);
  for (int i = 0; i < num_i; ++i) {
    if (v_occupied_state[i] == 1) {
      v_state_utility[i] = v_states_utilities[0] + decrement[i];
    } else if (v_occupied_state[i] == 2) {
      v_state_utility[i] = v_states_utilities[1] + decrement[i];
    } else if (v_occupied_state[i] == 3) {
      v_state_utility[i] = v_states_utilities[2] + decrement[i];
    } else if (v_occupied_state[i] == 4) {
      v_state_utility[i] = v_states_utilities[3];
    }
  }
  
  // Calculate Quality Adjusted Life Years (QALYs)
  arma::vec QALYs = v_state_utility * cycle_length;
  
  return QALYs;
}

// calc_discount_wtsC3 (function name in calc_discount_wtsC.cpp)
// 
// Expands the used data structures and related methods to those supported by
// arma and uses int instead of double for num_cycles
// 
// [[Rcpp::export]]
arma::vec calc_discount_wtsC(const double discount_rate, 
                             const int num_cycles,
                             const double cycle_length) {
  
  // Prepare object of size time_horizon + 1:
  arma::vec weights = arma::linspace(0, (num_cycles * cycle_length), 
                                     static_cast<size_t>(num_cycles) + 1);
  
  // Calculate discounting weights using element-wise operation
  weights.transform([discount_rate](double time) -> double {
    return 1.0 / std::pow(1.0 + discount_rate, time);
  });
  
  return weights;
}

// run_microSimC1
// 
// This implementation of the Sick-Sicker-Dead microsimulation model uses state
// indexes rather than state names. 
//
// [[Rcpp::export]]
Rcpp::List run_microSimC1(const arma::ivec& v_starting_states,
                          const int num_i,
                          const int num_cycles,
                          arma::mat m_indi_features,
                          const arma::ivec& v_states_index,
                          const arma::colvec& v_states_costs,
                          const arma::colvec& v_cost_coeffs,
                          const arma::colvec& v_states_utilities,
                          const arma::colvec& v_util_coeffs,
                          const arma::colvec& v_util_t_decs,
                          const Rcpp::List& l_trans_probs,
                          const double discount_rate_costs,
                          const double discount_rate_QALYs,
                          const double cycle_length = 1,
                          const int starting_seed = 1,
                          const int age_column_index = 0) {
  // declare R's set.seed function to control setting the seed number:
  Rcpp::Function set_seed("set.seed");
  
  // declare R's print and cat functions for updates and debugging:
  Rcpp::Function print("print");
  Rcpp::Function cat("cat");
  
  // declare transitions, costs, and effects matrices:
  arma::imat m_States(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Costs(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Effs(num_i, (num_cycles + 1), arma::fill::zeros);
  
  // set random number generator stream (RNG seed) for reproducibility
  // this is done in R as WM still does not know how to do the same in C++
  set_seed(Rcpp::Named("seed", (starting_seed)));
  
  // initialize parameter tracking time in current state
  arma::colvec v_time_in_state(num_i, arma::fill::ones);
  
  // indicate the initial health state:
  m_States.col(0) = v_starting_states;
  
  // estimate costs for all individuals at the initial health state:
  m_Costs.col(0) = calc_costsC(
    m_States.col(0), 
    v_states_costs, 
    m_indi_features, 
    v_cost_coeffs
  );
  
  // estimate QALYs for all individuals at the initial health state:
  m_Effs.col(0) = calc_effsC(
    m_States.col(0),
    v_states_utilities, 
    m_indi_features, 
    v_util_coeffs, 
    v_util_t_decs, 
    v_time_in_state, 
    cycle_length
  );
  
  for(int t = 0; t < num_cycles; t++) {
    
    // calculate the transition probabilities for cycle (t + 1):
    arma::mat m_trans_probs  = update_probsC(
      v_states_index,
      m_States.col(t),
      l_trans_probs,
      v_time_in_state
    );
    
    // sample the health state at (t+1):
    m_States.col(t + 1) = sampleC(
      m_trans_probs
    );
    
    // Check if remains in current state at 't + 1'
    arma::uvec stayed = arma::find(m_States.col(t) == m_States.col(t + 1));
    arma::uvec transitioned = arma::find(m_States.col(t) != m_States.col(t + 1));
    
    // Increment time spent in state for those who remained in current state
    v_time_in_state(stayed) += 1;
    
    // Reset time to 1 once transitioned
    v_time_in_state(transitioned).ones();
    
    // Assuming "age" is the first column (index 0)
    m_indi_features.col(age_column_index) += 1;    
    
    // estimate costs for all individuals occupying health states at (t+1):
    m_Costs.col(t + 1) = calc_costsC(
      m_States.col(t + 1), 
      v_states_costs, 
      m_indi_features, 
      v_cost_coeffs
    );
    
    // estimate QALYs for all individuals occupying health states at (t+1):
    m_Effs.col(t + 1) = calc_effsC(
      m_States.col(t + 1),
      v_states_utilities, 
      m_indi_features, 
      v_util_coeffs, 
      v_util_t_decs, 
      v_time_in_state, 
      cycle_length
    );
    
  } // end of loops through cycles:
  
  // calculate discount weights for both outcomes:
  arma::colvec v_c_dsc_wts = calc_discount_wtsC(
    discount_rate_costs,
    num_cycles,
    cycle_length
  );
  arma::colvec  v_e_dsc_wts = calc_discount_wtsC(
    discount_rate_QALYs,
    num_cycles,
    cycle_length
  );
  
  // compute summary stats:
  // total and average costs and QALYs (per individual):
  arma::vec v_total_costs = arma::sum(m_Costs, 1);
  arma::vec v_total_qalys = arma::sum(m_Effs, 1);
  double mean_costs = arma::mean(v_total_costs);
  double mean_qalys = arma::mean(v_total_qalys);
  
  // total and average discounted costs and QALYs (per individual):
  arma::vec v_total_Dcosts = m_Costs * v_c_dsc_wts;
  arma::vec v_total_Dqalys = m_Effs * v_e_dsc_wts;
  double mean_Dcosts = arma::mean(v_total_Dcosts);
  double mean_Dqalys = arma::mean(v_total_Dqalys);
  
  // save the results in a list:
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_("m_States") = m_States,
    Rcpp::_["m_Costs"] = m_Costs, 
    Rcpp::_["m_Effs"] = m_Effs,
    Rcpp::_["v_total_costs"] = v_total_costs, 
    Rcpp::_["v_total_qalys"] = v_total_qalys, 
    Rcpp::_["v_total_Dcosts"] = v_total_Dcosts, 
    Rcpp::_["v_total_Dqalys"] = v_total_Dqalys, 
    Rcpp::_["mean_costs"] = mean_costs,
    Rcpp::_["mean_qalys"] = mean_qalys,
    Rcpp::_["mean_Dcosts"] = mean_Dcosts,
    Rcpp::_["mean_Dqalys"] = mean_Dqalys
  );
  
  return(results) ;
}

// run_microSimC2
// 
// Similar to 'run_microSimC1' but only returns discounted and un-discounted
// outcomes
//
// [[Rcpp::export]]
Rcpp::List run_microSimC2(const arma::ivec& v_starting_states,
                          const int num_i,
                          const int num_cycles,
                          arma::mat m_indi_features,
                          const arma::ivec& v_states_index,
                          const arma::colvec& v_states_costs,
                          const arma::colvec& v_cost_coeffs,
                          const arma::colvec& v_states_utilities,
                          const arma::colvec& v_util_coeffs,
                          const arma::colvec& v_util_t_decs,
                          const Rcpp::List& l_trans_probs,
                          const double discount_rate_costs,
                          const double discount_rate_QALYs,
                          const double cycle_length = 1,
                          const int starting_seed = 1,
                          const int age_column_index = 0) {
  // declare R's set.seed function to control setting the seed number:
  Rcpp::Function set_seed("set.seed");
  
  // declare R's print and cat functions for updates and debugging:
  Rcpp::Function print("print");
  Rcpp::Function cat("cat");
  
  // declare transitions, costs, and effects matrices:
  arma::imat m_States(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Costs(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Effs(num_i, (num_cycles + 1), arma::fill::zeros);
  
  // set random number generator stream (RNG seed) for reproducibility
  // this is done in R as WM still does not know how to do the same in C++
  set_seed(Rcpp::Named("seed", (starting_seed)));
  
  // initialize parameter tracking time in current state
  arma::colvec v_time_in_state(num_i, arma::fill::ones);
  
  // indicate the initial health state:
  m_States.col(0) = v_starting_states;
  
  // estimate costs for all individuals at the initial health state:
  m_Costs.col(0) = calc_costsC(
    m_States.col(0), 
    v_states_costs, 
    m_indi_features, 
    v_cost_coeffs
  );
  
  // estimate QALYs for all individuals at the initial health state:
  m_Effs.col(0) = calc_effsC(
    m_States.col(0),
    v_states_utilities, 
    m_indi_features, 
    v_util_coeffs, 
    v_util_t_decs, 
    v_time_in_state, 
    cycle_length
  );
  
  for(int t = 0; t < num_cycles; t++) {
    
    // calculate the transition probabilities for cycle (t + 1):
    arma::mat m_trans_probs  = update_probsC(
      v_states_index,
      m_States.col(t),
      l_trans_probs,
      v_time_in_state
    );
    
    // sample the health state at (t+1):
    m_States.col(t + 1) = sampleC(
      m_trans_probs
    );
    
    // Check if remains in current state at 't + 1'
    arma::uvec stayed = arma::find(m_States.col(t) == m_States.col(t + 1));
    arma::uvec transitioned = arma::find(m_States.col(t) != m_States.col(t + 1));
    
    // Increment time spent in state for those who remained in current state
    v_time_in_state(stayed) += 1;
    
    // Reset time to 1 once transitioned
    v_time_in_state(transitioned).ones();
    
    // Assuming "age" is the first column (index 0)
    m_indi_features.col(age_column_index) += 1;    
    
    // estimate costs for all individuals occupying health states at (t+1):
    m_Costs.col(t + 1) = calc_costsC(
      m_States.col(t + 1), 
      v_states_costs, 
      m_indi_features, 
      v_cost_coeffs
    );
    
    // estimate QALYs for all individuals occupying health states at (t+1):
    m_Effs.col(t + 1) = calc_effsC(
      m_States.col(t + 1),
      v_states_utilities, 
      m_indi_features, 
      v_util_coeffs, 
      v_util_t_decs, 
      v_time_in_state, 
      cycle_length
    );
    
  } // end of loops through cycles:
  
  // calculate discount weights for both outcomes:
  arma::colvec v_c_dsc_wts = calc_discount_wtsC(
    discount_rate_costs,
    num_cycles,
    cycle_length
  );
  arma::colvec  v_e_dsc_wts = calc_discount_wtsC(
    discount_rate_QALYs,
    num_cycles,
    cycle_length
  );
  
  // compute summary stats:
  // total and average costs and QALYs (per individual):
  arma::vec v_total_costs = arma::sum(m_Costs, 1);
  arma::vec v_total_qalys = arma::sum(m_Effs, 1);
  double mean_costs = arma::mean(v_total_costs);
  double mean_qalys = arma::mean(v_total_qalys);
  
  // total and average discounted costs and QALYs (per individual):
  arma::vec v_total_Dcosts = m_Costs * v_c_dsc_wts;
  arma::vec v_total_Dqalys = m_Effs * v_e_dsc_wts;
  double mean_Dcosts = arma::mean(v_total_Dcosts);
  double mean_Dqalys = arma::mean(v_total_Dqalys);
  
  // save the results in a list:
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_["mean_costs"] = mean_costs,
    Rcpp::_["mean_qalys"] = mean_qalys,
    Rcpp::_["mean_Dcosts"] = mean_Dcosts,
    Rcpp::_["mean_Dqalys"] = mean_Dqalys
  );
  
  return(results) ;
}

////////////////////// run_microSimC0 functions ///////////////////////////

// update_probsC1 (function name in update_probsC.cpp)
//
// Like R's 'V' version of the function, 'C1' utilizes named vector
//
Rcpp::NumericMatrix update_probsC0(const Rcpp::CharacterVector& v_states_names,
                                   const Rcpp::CharacterVector& v_occupied_state,
                                   const Rcpp::List& l_trans_probs,
                                   const Rcpp::NumericVector& v_time_in_state) {
  
  int n = v_time_in_state.size();
  Rcpp::NumericMatrix m_probs(n, v_states_names.size());
  Rcpp::rownames(m_probs) = v_occupied_state;
  Rcpp::colnames(m_probs) = v_states_names;
  
  double r_S1D = -log(1 - Rcpp::as<double>(l_trans_probs["p_S1D"]));
  double r_S2D = -log(1 - Rcpp::as<double>(l_trans_probs["p_S2D"]));
  Rcpp::NumericVector p_S1D(n), p_S2D(n), p_HD(n);
  
  for (int i = 0; i < n; ++i) {
    if (v_occupied_state[i] == "S1") {
      p_S1D[i] = 1 - exp(-r_S1D * (1 + v_time_in_state[i] * Rcpp::as<double>(l_trans_probs["rp_S1"])));
    } else if (v_occupied_state[i] == "S2") {
      p_S2D[i] = 1 - exp(-r_S2D * (1 + v_time_in_state[i] * Rcpp::as<double>(l_trans_probs["rp_S2"])));
    } else if (v_occupied_state[i] == "H") {
      p_HD[i] = Rcpp::as<double>(l_trans_probs["p_HD"]);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    if (v_occupied_state[i] == "H") {
      m_probs(i, Rcpp::_) = Rcpp::NumericVector::create(1 - Rcpp::as<double>(l_trans_probs["p_HS1"]) - p_HD[i], Rcpp::as<double>(l_trans_probs["p_HS1"]), 0, p_HD[i]);
    } else if (v_occupied_state[i] == "S1") {
      m_probs(i, Rcpp::_) = Rcpp::NumericVector::create(Rcpp::as<double>(l_trans_probs["p_S1H"]), 1 - Rcpp::as<double>(l_trans_probs["p_S1S2"]) - Rcpp::as<double>(l_trans_probs["p_S1H"]) - p_S1D[i], Rcpp::as<double>(l_trans_probs["p_S1S2"]), p_S1D[i]);
    } else if (v_occupied_state[i] == "S2") {
      m_probs(i, Rcpp::_) = Rcpp::NumericVector::create(0, 0, 1 - p_S2D[i], p_S2D[i]);
    } else if (v_occupied_state[i] == "D") {
      m_probs(i, Rcpp::_) = Rcpp::NumericVector::create(0, 0, 0, 1);
    }
  }
  
  // Sanity check
  for (int i = 0; i < n; ++i) {
    if (fabs(sum(m_probs(i, Rcpp::_)) - 1) > 1e-12) {
      Rcpp::stop("Probabilities do not sum to 1");
    }
  }
  
  return m_probs;
}

// sampleC0 (function name in sampleC.cpp)
//
// Like R's 'V' version of the function, 'C0' utilizes named vector
//
Rcpp::CharacterVector sampleC0(const Rcpp::NumericMatrix& m_trans_probs,
                               const Rcpp::CharacterVector& v_states_names) {
  
  // Number of individuals and states 
  int num_i = m_trans_probs.nrow();
  int num_states = m_trans_probs.ncol();
  
  // Create an upper triangular matrix of ones
  Rcpp::NumericMatrix m_upper_tri(num_states, num_states);
  for (int i = 0; i < num_states; i++) {
    for (int j = i; j < num_states; j++) {
      m_upper_tri(i, j) = 1;
    }
  }
  
  // Create matrix with row-wise cumulative transition probabilities
  Rcpp::NumericMatrix m_cum_probs(num_i, num_states);
  for (int i = 0; i < num_i; i++) {
    for (int j = 0; j < num_states; j++) {
      for (int k = 0; k < num_states; k++) {
        m_cum_probs(i, j) += m_trans_probs(i, k) * m_upper_tri(k, j);
      }
    }
  }
  
  // Ensure that the maximum cumulative probabilities are equal to 1
  for (int i = 0; i < num_i; i++) {
    if (m_cum_probs(i, num_states - 1) > 1.000000) {
      Rcpp::stop("Error in multinomial sampling: probabilities do not sum to 1");
    }
  }
  
  // Sample random values from Uniform standard distribution for each individual
  Rcpp::NumericVector v_rand_values = Rcpp::runif(num_i);
  
  // Repeat each sampled value to have as many copies as the number of states
  Rcpp::NumericMatrix m_rand_values(num_i, num_states);
  for (int i = 0; i < num_i; i++) {
    for (int j = 0; j < num_states; j++) {
      m_rand_values(i, j) = v_rand_values[i];
    }
  }
  
  // Identify transitions, compare random samples to cumulative probabilities
  Rcpp::LogicalMatrix m_transitions(num_i, num_states);
  for (int i = 0; i < num_i; i++) {
    for (int j = 0; j < num_states; j++) {
      m_transitions(i, j) = m_rand_values(i, j) > m_cum_probs(i, j);
    }
  }
  
  // Sum transitions to identify health state in next cycle
  Rcpp::IntegerVector v_transitions(num_i);
  for (int i = 0; i < num_i; i++) {
    int sum_transitions = 0;
    for (int j = 0; j < num_states; j++) {
      if (m_transitions(i, j)) {
        sum_transitions++;
      }
    }
    v_transitions[i] = sum_transitions;
  }
  
  // Identify health state to which each individual is transitioning
  Rcpp::CharacterVector v_health_states(num_i);
  for (int i = 0; i < num_i; i++) {
    v_health_states[i] = v_states_names[v_transitions[i]];
  }
  
  return v_health_states;
}

// calc_costsC0 (function name in calc_costsC.cpp)
//
// Like R's 'V' version of the function, 'C0' utilizes named vector
//
Rcpp::NumericVector calc_costsC0(const Rcpp::CharacterVector v_occupied_state,
                                 const Rcpp::NumericVector v_states_costs,
                                 const Rcpp::NumericMatrix m_indi_features,
                                 const Rcpp::NumericVector v_cost_coeffs) {
  
  int n_rows = m_indi_features.nrow();
  int n_cols = m_indi_features.ncol();
  Rcpp::NumericVector v_indi_costs(n_rows);
  
  // Manually perform matrix-vector multiplication
  for (int i = 0; i < n_rows; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n_cols; ++j) {
      sum += m_indi_features(i, j) * v_cost_coeffs[j];
    }
    v_indi_costs[i] = sum;
  }
  
  // Initialize the cost vector
  Rcpp::NumericVector v_state_costs(v_occupied_state.size(), NA_REAL);
  
  // Estimate costs based on occupied state
  for (int i = 0; i < v_occupied_state.size(); ++i) {
    std::string state = Rcpp::as<std::string>(v_occupied_state[i]);
    if (state == "H") {
      v_state_costs[i] = v_states_costs["H"];
    } else if (state == "S1") {
      v_state_costs[i] = v_states_costs["S1"] + v_indi_costs[i];
    } else if (state == "S2") {
      v_state_costs[i] = v_states_costs["S2"] + v_indi_costs[i];
    } else if (state == "D") {
      v_state_costs[i] = v_states_costs["D"];
    }
  }
  
  return v_state_costs;
}

// calc_effsC0 (function name in calc_effsC.cpp)
//
// Like R's 'V' version of the function, 'C0' utilizes named vector
//
Rcpp::NumericVector calc_effsC0(const Rcpp::CharacterVector& v_occupied_state,
                                const Rcpp::NumericVector& v_states_utilities,
                                const Rcpp::NumericMatrix& m_indi_features,
                                const Rcpp::NumericVector& v_util_coeffs,
                                const Rcpp::NumericVector& v_util_t_decs,
                                const Rcpp::NumericVector& v_time_in_state,
                                const double cycle_length = 1) {
  
  int n_rows = m_indi_features.nrow();
  int n_cols = m_indi_features.ncol();
  Rcpp::NumericVector v_ind_decrement(n_rows);
  
  // Manually perform matrix-vector multiplication to calculate individual-specific utility decrements
  for (int i = 0; i < n_rows; ++i) {
    double sum = 0.0;
    for (int j = 0; j < n_cols; ++j) {
      sum += m_indi_features(i, j) * v_util_coeffs[j];
    }
    v_ind_decrement[i] = sum;
  }
  
  // Calculate time-dependent state-specific utility decrements
  Rcpp::NumericVector time_decrement(v_occupied_state.size(), 0.0);
  for (int i = 0; i < v_occupied_state.size(); ++i) {
    std::string state = Rcpp::as<std::string>(v_occupied_state[i]);
    if (state == "S1") {
      time_decrement[i] = v_util_t_decs["S1"] * v_time_in_state[i];
    } else if (state == "S2") {
      time_decrement[i] = v_util_t_decs["S2"] * v_time_in_state[i];
    }
  }
  
  // Estimate total decrements
  Rcpp::NumericVector decrement = v_ind_decrement + time_decrement;
  
  // Estimate utilities based on occupied state
  Rcpp::NumericVector v_state_utility(v_occupied_state.size(), 0.0);
  for (int i = 0; i < v_occupied_state.size(); ++i) {
    std::string state = Rcpp::as<std::string>(v_occupied_state[i]);
    if (state == "H") {
      v_state_utility[i] = v_states_utilities["H"] + decrement[i];
    } else if (state == "S1") {
      v_state_utility[i] = v_states_utilities["S1"] + decrement[i];
    } else if (state == "S2") {
      v_state_utility[i] = v_states_utilities["S2"] + decrement[i];
    } else if (state == "D") {
      v_state_utility[i] = v_states_utilities["D"];
    }
  }
  
  // Calculate Quality Adjusted Life Years (QALYs)
  Rcpp::NumericVector QALYs = v_state_utility * cycle_length;
  
  return QALYs;
}

// calc_discount_wtsC6 (function name in calc_discount_wtsC.cpp)
// 
// Expands the used data structures and related methods to those supported by
// Rcpp and uses int instead of double for num_cycles
//
Rcpp::NumericVector calc_discount_wtsC0(const double discount_rate,
                                        const int num_cycles,
                                        const double cycle_length) {
  // Prepare object of size num_cycles + 1
  int time = num_cycles + 1;
  Rcpp::NumericVector weights(time); // Properly initialize the vector
  
  // Calculate discounting weights using element-wise operation
  for (int i = 0; i < time; ++i) {
    weights[i] = 1.0 / pow(1.0 + discount_rate, i * cycle_length);
  }
  
  return weights;
}

// run_microSimC0
// 
// This implementation of the Sick-Sicker-Dead microsimulation model uses state
// indexes rather than state names. 
//
// [[Rcpp::export]]
Rcpp::List run_microSimC0(const Rcpp::CharacterVector& v_starting_states,
                          const int num_i,
                          const int num_cycles,
                          Rcpp::NumericMatrix m_indi_features,
                          const Rcpp::CharacterVector& v_states_names,
                          const Rcpp::NumericVector& v_states_costs,
                          const Rcpp::NumericVector& v_cost_coeffs,
                          const Rcpp::NumericVector& v_states_utilities,
                          const Rcpp::NumericVector& v_util_coeffs,
                          const Rcpp::NumericVector& v_util_t_decs,
                          const Rcpp::List& l_trans_probs,
                          const double discount_rate_costs,
                          const double discount_rate_QALYs,
                          const double cycle_length = 1,
                          const int starting_seed = 1) {
  // declare R's set.seed function to control setting the seed number:
  Rcpp::Function set_seed("set.seed");
  
  // declare R's print and cat functions for updates and debugging:
  Rcpp::Function cat("cat");
  Rcpp::Function print("print");
  
  // declare transitions, costs, and effects matrices:
  // Rcpp::CharacterMatrix matrixName(row, column)
  // Rcpp::NumericMatrix matrixName(row, column)
  Rcpp::CharacterMatrix m_States(num_i, (num_cycles + 1));
  Rcpp::NumericMatrix m_Costs(num_i, (num_cycles + 1));
  Rcpp::NumericMatrix m_Effs(num_i, (num_cycles + 1));
  
  // set random number generator stream (RNG seed) for reproducibility
  // this is done in R as WM still does not know how to do the same in C++
  set_seed(Rcpp::Named("seed", (starting_seed)));
  
  // initialize parameter tracking time in current state
  // Rcpp::NumericVector vectorName(length, default value)
  Rcpp::NumericVector v_time_in_state(num_i, 1.0);
  
  // indicate the initial health state:
  m_States(Rcpp::_, 0) = v_starting_states;
  
  // estimate costs for all individuals at the initial health state:
  m_Costs(Rcpp::_, 0) = calc_costsC0(
    m_States(Rcpp::_, 0), 
    v_states_costs, 
    m_indi_features, 
    v_cost_coeffs
  );
  
  // estimate QALYs for all individuals at the initial health state:
  m_Effs(Rcpp::_, 0) = calc_effsC0(
    m_States(Rcpp::_, 0),
    v_states_utilities, 
    m_indi_features, 
    v_util_coeffs, 
    v_util_t_decs, 
    v_time_in_state, 
    cycle_length
  );
  
  for(int t = 0; t < num_cycles; t++) {
    
    // calculate the transition probabilities for cycle (t + 1):
    Rcpp::NumericMatrix m_trans_probs  = update_probsC0(
      v_states_names,
      m_States(Rcpp::_, t),
      l_trans_probs,
      v_time_in_state
    );
    
    // sample the health state at (t+1):
    m_States(Rcpp::_, t + 1) = sampleC0(
      m_trans_probs,
      v_states_names
    );
    
    // Check if remains in current state at 't + 1'
    // Create logical vectors for those who stayed or transitioned
    Rcpp::LogicalVector stayed = (m_States(Rcpp::_, t) == m_States(Rcpp::_, t + 1));
    Rcpp::LogicalVector transitioned = !stayed;

    // Increment time spent in state for those who remained in current state
    v_time_in_state = Rcpp::ifelse(stayed, v_time_in_state + 1, v_time_in_state);
    
    // Reset time to 1 for those who transitioned
    v_time_in_state = Rcpp::ifelse(transitioned, 1, v_time_in_state);
    
    // Assuming "age" is the first column (index 0)
    m_indi_features(Rcpp::_, 0) = m_indi_features(Rcpp::_, 0) + 1;    
    
    // estimate costs for all individuals occupying health states at (t+1):
    m_Costs(Rcpp::_, t + 1) = calc_costsC0(
      m_States(Rcpp::_, t + 1), 
      v_states_costs, 
      m_indi_features, 
      v_cost_coeffs
    );
    
    // estimate QALYs for all individuals occupying health states at (t+1):
    m_Effs(Rcpp::_, t + 1) = calc_effsC0(
      m_States(Rcpp::_, t + 1),
      v_states_utilities, 
      m_indi_features, 
      v_util_coeffs, 
      v_util_t_decs, 
      v_time_in_state, 
      cycle_length
    );
  } // end of loops through cycles:
  
  // calculate discount weights for both outcomes:
  Rcpp::NumericVector v_c_dsc_wts = calc_discount_wtsC0(
    discount_rate_costs,
    num_cycles,
    cycle_length
  );
  Rcpp::NumericVector v_e_dsc_wts = calc_discount_wtsC0(
    discount_rate_QALYs,
    num_cycles,
    cycle_length
  );
  
  // compute summary stats:
  // total and average costs and QALYs (per individual):
  Rcpp::NumericVector v_total_costs = Rcpp::rowSums(m_Costs);
  Rcpp::NumericVector v_total_qalys = Rcpp::rowSums(m_Effs);
  double mean_costs = Rcpp::mean(v_total_costs);
  double mean_qalys = Rcpp::mean(v_total_qalys);
  
  // total and average discounted costs and QALYs (per individual):
  Rcpp::NumericMatrix m_total_Dcosts(num_i, num_cycles + 1);
  Rcpp::NumericMatrix m_total_Dqalys(num_i, num_cycles + 1);
  for(int i = 0; i <= num_cycles ; i++) {
    m_total_Dcosts(Rcpp::_, i) = m_Costs(Rcpp::_, i) * v_c_dsc_wts[i];
    m_total_Dqalys(Rcpp::_, i) = m_Effs(Rcpp::_, i) * v_e_dsc_wts[i];
  }
  Rcpp::NumericVector v_total_Dcosts = Rcpp::rowSums(m_total_Dcosts);
  Rcpp::NumericVector v_total_Dqalys = Rcpp::rowSums(m_total_Dqalys);
  double mean_Dcosts = Rcpp::mean(v_total_Dcosts);
  double mean_Dqalys = Rcpp::mean(v_total_Dqalys);

  // save the results in a list:
  Rcpp::List results = Rcpp::List::create(
    Rcpp::_("m_States") = m_States,
    Rcpp::_["m_Costs"] = m_Costs, 
    Rcpp::_["m_Effs"] = m_Effs,
    Rcpp::_["v_total_costs"] = v_total_costs, 
    Rcpp::_["v_total_qalys"] = v_total_qalys, 
    Rcpp::_["v_total_Dcosts"] = v_total_Dcosts, 
    Rcpp::_["v_total_Dqalys"] = v_total_Dqalys, 
    Rcpp::_["mean_costs"] = mean_costs,
    Rcpp::_["mean_qalys"] = mean_qalys,
    Rcpp::_["mean_Dcosts"] = mean_Dcosts,
    Rcpp::_["mean_Dqalys"] = mean_Dqalys
  );
  
  return(results) ;
}