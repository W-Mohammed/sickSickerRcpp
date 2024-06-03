#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// update_probsC6:
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

// sampleC5
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

// calc_costsC7
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

// calc_effsC6
//
// uses for loop and IF statement to assign costs and a for loop to calculate
// the cost regression estimate for each individual and passing objects by
// reference (&)
//
// [[Rcpp::export]]
arma::vec calc_effsC(const arma::ivec v_occupied_state,
                     const arma::vec v_states_utilities,
                     const arma::mat m_indi_features,
                     const arma::vec v_util_coeffs,
                     const arma::vec v_util_t_decs,
                     const arma::vec v_time_in_state,
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

// calc_discount_wtsC3
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

// run_microSimC
//
// [[Rcpp::export]]
Rcpp::List run_microSimC1(arma::ivec& v_starting_states,
                          int num_i,
                          int num_cycles,
                          arma::mat& m_indi_features,
                          arma::ivec& v_states_index,
                          arma::colvec& v_states_costs,
                          arma::colvec& v_cost_coeffs,
                          arma::colvec& v_states_utilities,
                          arma::colvec& v_util_coeffs,
                          arma::colvec& v_util_t_decs,
                          Rcpp::List l_trans_probs,
                          double discount_rate_costs,
                          double discount_rate_QALYs,
                          double cycle_length = 1,
                          int starting_seed = 1,
                          int age_column_index = 0) {
  // declare R's set.seed function to control setting the seed number:
  Rcpp::Function set_seed("set.seed");
  
  // declare R's print and cat functions for updates and debugging:
  Rcpp::Function cat("cat");
  Rcpp::Function paste("paste");
  
  // declare transitions, costs, and effects matrices:
  arma::imat m_States(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Costs(num_i, (num_cycles + 1), arma::fill::zeros);
  arma::mat m_Effs(num_i, (num_cycles + 1), arma::fill::zeros);
  
  // set random number generator stream (RNG seed) for reproducibility
  // this is done in R as WM still does not know how to do the same in C++
  set_seed(Rcpp::Named("seed", (starting_seed)));
  
  // initialize parameter tracking time in current state
  arma::colvec v_time_in_state(num_i, arma::fill::zeros);
  
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