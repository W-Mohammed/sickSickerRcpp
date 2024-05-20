#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::vec estimate_discounting_weights( const double discount_rate, 
                                        const double time_horizon) {
  // discount_rate: the discount rate. 
  // time_horizon: the model time horizon.
  
  // Prepare object of size time_horizon + 1:
  arma::vec weights = arma::linspace(0, time_horizon, static_cast<size_t>(time_horizon) + 1);
  
  // Calculate discounting weights using element-wise operation
  weights.transform([discount_rate](double time) -> double {
    return 1.0 / std::pow(1.0 + discount_rate, time);
  });
  
  return(weights) ;
}

// [[Rcpp::export]]
void time_in_state( arma::mat& m_t_states,
                    const arma::colvec& v_S_t,
                    const std::vector<int>& v_tracked_states) {
  // m_t_states: matrix storing times spent in each health state.
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_tracked_states: vector of health states in which to track time.
  // This function tracks the time (in number of cycles) each individual spends
  // in tracked states. Individuals exiting states gets their times reset to 0.
  // This version of the function returns nothing but makes changes to the 
  // function in place.
  
  // Track time if in tracked states:
  for(int s : v_tracked_states) {
    // Re-scale states numerical identifier (s) to C++ index (i):
    int i = s - 1 ;
    
    // Create a logical vector where state matches s:
    arma::uvec other_states = arma::find(v_S_t != s) ;
    
    // Copy column data to a temporary vector:
    arma::vec state_col = m_t_states.col(i) ;
    
    // Increment the time for all individuals in the column.
    state_col += 1 ;
    
    // Ensure other_states has valid indices within the range of tmp_col
    if (!other_states.empty()) {
      // Reset time in other states to 0 using .elem():
      state_col.elem(other_states).fill(0) ;
    }
    
    // Record time back to states time matrix:
    m_t_states.col(i) = state_col ;
  }
}
// 
// // [[Rcpp::export]]
// arma::mat time_in_state( arma::mat m_t_states,
//                          const arma::colvec& v_S_t,
//                          const std::vector<int>& v_tracked_states) {
//   // m_t_states: matrix storing times spent in each health state.
//   // v_S_t: vector of health states occupied by individuals at cycle t.
//   // v_tracked_states: vector of health states for which to track time.
//   // This function tracks the time (in number of cycles) each individual spends
//   // in tracked states. Individuals exiting states gets their times reset to 0.
//   // This version of the function returns a matrix.
//   
//   // Track time if in tracked states:
//   for(int s : v_tracked_states) {
//     // Re-scale states numerical identifier (s) to C++ index (i):
//     int i = s - 1 ;
//     
//     // Create a logical vector where state matches s:
//     arma::uvec other_states = arma::find(v_S_t != s) ;
//     
//     // Copy column data to a temporary vector:
//     arma::vec state_col = m_t_states.col(i) ;
//     
//     // Increment the time for all individuals in the column.
//     state_col += 1 ;
//     
//     // Ensure other_states has valid indices within the range of tmp_col
//     if (!other_states.empty()) {
//       // Reset time in other states to 0 using .elem():
//       state_col.elem(other_states).fill(0) ;
//     }
//     
//     // Record time back to states time matrix:
//     m_t_states.col(i) = state_col ;
//   }
//   
//   return(m_t_states);
// }
// 
// // [[Rcpp::export]]
// arma::mat ProbsV_Cpp( arma::rowvec v_S_t,
//                       int& n_I,
//                       int& n_S,
//                       arma::mat& t_P ) {
//   // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
//   // n_I: number of simulated individuals
//   // n_S: number of health states
//   // t_P: transition probabilities matrix.
//   
//   // create a matrix for the state transition probabilities at time t (m_P_t):
//   // mat(n_rows, n_cols) initiated with zeros
//   arma::mat m_P_t(n_I, n_S) ;
//   
//   // create a cube/array for the probabilities of transitioning from current states:
//   // cube(n_rows, n_cols, n_slices) initiated with zeros
//   arma::cube c_P_t(n_I, n_S, n_S) ;
//   
//   // assign probabilities from t_P to the slices
//   for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
//     
//     // Get the transition probabilities from health state i as a row vector:
//     // Replicate t_P.row(i) n_I times vertically, and 1 time horizontally:
//     c_P_t.slice(i) = arma::repmat(t_P.row(i), n_I, 1) ;
//   }
//   
//   // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
//   for(int i = 0; i < n_S; i++) {
//     // Identify individuals occupying state i at time t:
//     // Adjust states occupancy to allign with C++ indexing, which starts at 0:
//     arma::uvec v_O_t = arma::find(v_S_t == (i + 1)) ;
//     
//     // Grab transition probabilities from state i
//     // Assign probabilities for state i occupancy according to v_O_t:
//     m_P_t.rows(v_O_t) = c_P_t.slice(i).rows(v_O_t) ;
//   }
//   
//   // check if vector of probabilities sum to 1
//   // need to round up to the 1e-6, otherwise it breaks
//   // Rounding to the 6th decimal place
//   int n = 6 ;
//   arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / 
//     std::pow(10, n) ;
//   bool notSumToOne = arma::any(sum_row != 1.000000) ;
//   
//   if(notSumToOne) {
//     stop("Probabilities do not sum to 1!") ;
//   }
//   else {
//     return(m_P_t) ;
//   }
// }

// [[Rcpp::export]]
arma::mat ProbsV_Cpp( arma::rowvec v_S_t,
                      int& n_I,
                      int& n_S,
                      arma::mat& t_P,
                      arma::mat& m_t_states, 
                      std::vector<int>& v_states_from, 
                      std::vector<int>& v_states_to,
                      std::vector<int>& v_states_comp,
                      std::vector<double>& v_increase_rate) {
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // n_I: number of simulated individuals
  // n_S: number of health states
  // t_P: transition probabilities matrix.
  // m_t_states: matrix containing time spent in tracked states.
  // v_states_from: vector containing the states from which transitions occur
  // v_states_to: vector containing the states to which transitions occur
  // v_states_comp: vector containing the states complementing the transition probabilities
  // v_increase_rate: change in transition rates with every additional cycle in tracked states
  
  // create a matrix for the state transition probabilities at time t (m_P_t):
  // mat(n_rows, n_cols) initiated with zeros
  arma::mat m_P_t(n_I, n_S) ;
  
  // create a cube/array for the probabilities of transitioning from current states:
  // cube(n_rows, n_cols, n_slices) initiated with zeros
  arma::cube c_P_t(n_I, n_S, n_S) ;
  
  // assign probabilities from t_P to the slices
  for(arma::uword i = 0; i < c_P_t.n_slices; i++) {
    
    // Get the transition probabilities from health state i as a row vector:
    // Replicate t_P.row(i) n_I times vertically, and 1 time horizontally:
    c_P_t.slice(i) = arma::repmat(t_P.row(i), n_I, 1) ;
    
    // Get the health state to which the transition occur, re-scale to C++:
    int state_from = i + 1 ;
    
    // Use std::find to check if i is in v_states_from
    auto it = std::find(v_states_from.begin(), v_states_from.end(), state_from) ;
    if (it != v_states_from.end()) {
      // Found state_from in v_states_from, get the index:
      auto index = std::distance(v_states_from.begin(), it) ;
      
      // Check if state_from appears more than once in the v_states_from:
      for(size_t y = index; y <  v_states_from.size(); y++) {
        if(state_from == v_states_from[y]) {
          
          // Get the health state to which the transition occur, re-scale to C++:
          int state_to = v_states_to[y] - 1 ;
          
          // Convert the probabilities in the target column to rates:
          arma::vec v_R_t = - arma::log(1 - c_P_t.slice(i).col(state_to)) ;
          
          // Get cycle change in rates:
          double delta_rate = v_increase_rate[y] ;
          
          // Get accumulated change over time in i state:
          arma::vec v_dR_t = 1 + m_t_states.col(i) * delta_rate ;
          
          // Convert the rates back to probabilities after adjusting them:
          arma::colvec v_P_t = 1 - arma::exp(- v_R_t % v_dR_t) ;
          
          c_P_t.slice(i).col(state_to) =  v_P_t ; 
        }
      }
      
      // Calculate complement and update the target column
      int target = v_states_comp[index] - 1 ;
      c_P_t.slice(i).col(target) = 1 - (arma::sum(c_P_t.slice(i), 1) - c_P_t.slice(i).col(target)) ;
    }
  }
  
  // Add probabilities to m_P_t based on time t and state i or S as in v_S_t:
  for(int i = 0; i < n_S; i++) {
    // Identify individuals occupying state i at time t:
    // Adjust states occupancy to allign with C++ indexing, which starts at 0:
    arma::uvec v_O_t = arma::find(v_S_t == (i + 1)) ;
    
    // Grab transition probabilities from state i
    // Assign probabilities for state i occupancy according to v_O_t:
    m_P_t.rows(v_O_t) = c_P_t.slice(i).rows(v_O_t) ;
  }
  
  // check if vector of probabilities sum to 1
  // need to round up to the 1e-6, otherwise it breaks
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec sum_row = arma::round(arma::sum(m_P_t, 1) * std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = arma::any(sum_row != 1.000000) ;
  
  if(notSumToOne) {
    stop("Probabilities do not sum to 1!") ;
  }
  else {
    return(m_P_t) ;
  }
}

// [[Rcpp::export]]
arma::colvec SampleV_Cpp( arma::mat& m_P_t,
                          int& n_I,
                          int& n_S) {
  // m_P_t: numeric matrix containing the transition probabilities for individuals at cycle t
  // n_I: number of simulated individuals.
  // n_S: number of health states.
  
  // create m_CumProb_t matrix with row-wise cumulative transition probabilities:
  arma::mat m_CumProb_t = arma::cumsum(m_P_t, 1) ;
  
  // create a column vector for sampled health states (v_s_t):
  arma::colvec v_s_t = arma::ones<arma::colvec>(n_I) ;
  
  // recheck the probabilities:
  // need to round up to the 1e-6, otherwise it breaks
  // Rounding to the 6th decimal place
  int n = 6 ;
  arma::colvec v_CumProb_t = arma::round(m_CumProb_t.col(m_CumProb_t.n_cols - 1) * 
    std::pow(10, n)) / std::pow(10, n) ;
  bool notSumToOne = arma::any(v_CumProb_t != 1.000000) ;
  
  if(notSumToOne) {
    stop("error in multinom: probabilities do not sum to 1") ;
  }
  
  // Initialise a matrix to save values sampled from the U~(0, 1)  
  arma::mat m_U(n_I, n_S, arma::fill::ones) ;
  // in each row, sample one random value for n_I individuals and repeat that value n_S times
  m_U = arma::repmat( arma::randu<arma::colvec>(n_I), 1, n_S ) ;
  
  // identify the first individual-specific health-state with a cumulative probability higher than the their corresponding sampled value
  // using a logical (true/false or 1/0), matrix get the value to be added to 1 (the starting health-state)
  // one plus the sum of each row of the logical matrix gives the health-state for the corresponding individuals at time t + 1
  v_s_t = v_s_t + arma::sum( (m_U > m_CumProb_t), 1 ) ;
  
  return(v_s_t) ;
}

// [[Rcpp::export]]
arma::colvec CostsV_Cpp( arma::colvec v_S_t,
                         arma::vec& v_Costs) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Costs: a vector containing cost parameters.
  
  // Transforming state occupancy to allign with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  // Use states indecies to assign costs correctly:
  arma::colvec v_col_costs = v_Costs.elem(uv_indices) ;
  
  return(v_col_costs) ;
}

// [[Rcpp::export]]
arma::colvec EffsV_Cpp( arma::colvec v_S_t,
                        arma::vec& v_Utilities,
                        double cycle = 1 ) {
  // v_S_t: vector of health states occupied by individuals at cycle t.
  // v_Utilities: a vector containing utilities values for each health state.
  // cycle: integer variable indicating cycle length in years - default is 1.
  
  // Transforming state occupancy to align with C++ indexing:
  arma::uvec uv_indices = arma::conv_to<arma::uvec>::from(v_S_t - 1) ;
  
  // Calculating QALYs, multiplying utilities by the length of the cycle length:
  arma::vec v_QALYs = v_Utilities * cycle ;
  
  // Use states indecies to assign utilities correctly:
  arma::colvec v_col_QALYs = v_QALYs.elem(uv_indices) ;
  
  return(v_col_QALYs) ;
}

/************* Micro-simulation ***************/

// [[Rcpp::export]]
// Micro-simulation generic function. This version of the function can track the
// time spent in user-defined states and employ this information to affect 
// probabilities of transitioning to other health states, costs and/or
// consequences associated with the tracked health states.
// v_S_t numeric vector containing the health states occupied by the individuals at cycle t.
// m_t_P matrix containing transition probabilities.
// v_Costs numeric vector containing costs.
// v_Utilities numeric vector containing utilities.
// n_I integer scalar specifying the number of simulated individuals.
// v_tracked_states integer vector of health states for which to track time.
// v_states_from integer vector containing the states from which transitions occur.
// v_states_to integer vector containing the states to which transitions occur.
// v_states_comp integer vector containing the states complementing the transition probabilities.
// v_increase_rate double vector specifying the change in transition rates with every additional cycle in tracked states.
// n_S integer scalar number of health-states. Default is 4.
// n_T integer scalar defining the total number of cycles, model's time-horizon. Default is 31.
// n_Cycle_length double scalar defining, in years, the length of model cycle. Default is 1 year.
// d_dC double scalar defining the discount rate for the costs. Default is 3.5%.
// d_dE double scalar defining the discount rate for health outcome, e.g. QALYs. Default is 3.5%.
// n_Seed integer scalar setting the starting seed number for random number generator. Default is 1.
Rcpp::List MicroSimV_Cpp_3( arma::colvec& v_S_t,
                            arma::mat& m_t_P,
                            arma::vec& v_Costs,
                            arma::vec& v_Utilities,
                            int n_I,
                            std::vector<int>& v_tracked_states, 
                            std::vector<int>& v_states_from, 
                            std::vector<int>& v_states_to,
                            std::vector<int>& v_states_comp,
                            std::vector<double>& v_increase_rate,
                            int n_S = 4,
                            int n_T = 30,
                            double n_Cycle_length = 1,
                            double d_dC = 0.03,
                            double d_dE = 0.03,
                            int n_Seed = 1 ) {
  // Arguments:
  // v_S_t: numeric vector containing the health states occupied by the individuals at cycle t
  // m_t_P: matrix containing transition probabilities
  // v_Costs: vector containing costs
  // v_Utilities: vector containing utilities
  // n_I: number of individuals
  // v_tracked_states: vector of health states for which to track time
  // v_states_from: vector containing the states from which transitions occur
  // v_states_to: vector containing the states to which transitions occur
  // v_states_comp: vector containing the states complementing the transition probabilities
  // v_increase_rate: change in transition rates with every additional cycle in tracked states
  // n_S: number of health-states
  // n_T: total number of cycles to run the model
  // n_Cycle_length: length of model cycle (default is 1 year)
  // d_dC: discount rate for costs
  // d_dE: discount rate for health outcome (QALYs)
  // n_Seed: starting seed number for random number generator (default is 1)
  
  // Makes use of:
  // estimate_discounting_weights: calculates discounting weights
  // ProbsV_Cpp:  estimates the transition probabilities at times t
  // SampleV_Cpp: samples the health states that will be occupied by individuals at times t+1
  // CostsV_Cpp:  estimates the costs associated with state occupied at times t
  // EffsV_Cpp:   estimates the health effects associated with state occupied at times t
  
  // declare R's set.seed function to control setting the seed number:
  Rcpp::Function set_seed("set.seed") ;
  
  // declare other R functions:
  Rcpp::Function cat("cat") ;
  Rcpp::Function paste("paste") ;
  Rcpp::Function print("print") ;
  
  // declare discounting vectors and estimate them:
  arma::vec v_dwc = estimate_discounting_weights(d_dC, n_T) ;
  arma::vec v_dwe = estimate_discounting_weights(d_dE, n_T) ;
  
  // declare transitions, costs, and effects matrices:
  arma::mat m_M(n_I, (n_T + 1), arma::fill::zeros) ;
  arma::mat m_C(n_I, (n_T + 1), arma::fill::zeros) ;
  arma::mat m_E(n_I, (n_T + 1), arma::fill::zeros) ;
  arma::mat m_t_states(n_I, n_S, arma::fill::zeros) ;
  
  // indicate the initial health state:
  m_M.col(0) = v_S_t ;
  
  // set the seed for every individual for the random number generator
  // this is done in R as WM still does not know how to do the same in R
  set_seed(Named("seed", (n_Seed))) ;
  
  // estimate costs for all individuals at the initial health state conditional on treatment:
  m_C.col(0) = CostsV_Cpp(m_M.col(0), v_Costs) ;
  
  // estimate QALYs for all individuals at the initial health state conditional on treatment:
  m_E.col(0) = EffsV_Cpp (m_M.col(0), v_Utilities, n_Cycle_length) ;
  
  for(int t = 0; t < n_T; t++) {
    
    // track time in states (calling the void version of function):
    time_in_state(m_t_states, m_M.col(t), v_tracked_states) ;
    // calculate time-dependent transition probabilities for cycle (t + 1):
    arma::mat m_P  = ProbsV_Cpp (arma::trans(m_M.col(t)), n_I, n_S, m_t_P, m_t_states, v_states_from, v_states_to, v_states_comp, v_increase_rate) ;
    // sample the health state at (t+1):
    m_M.col(t + 1) = SampleV_Cpp(m_P, n_I, n_S) ;
    // estimate costs for all individuals occupying health states at (t+1) conditional on treatment:
    m_C.col(t + 1) = CostsV_Cpp(m_M.col(t + 1), v_Costs) ;
    // estimate QALYs for all individuals occupying health states at (t+1) conditional on treatment:
    m_E.col(t + 1) = EffsV_Cpp (m_M.col(t + 1), v_Utilities, n_Cycle_length) ;
    
    // print the progress of the simulation:
    cat(Named("...", '\r')) ;
    cat(Named("...", 
              paste(Named("...", (round(((double)(t + 1)/n_T) * 100 * 10) / 10)),
                    Named("...", "% done "), Named("sep", " ")))) ;
  } // end of loops through cycles:
  
  // compute summary stats:
  // total discounted costs and QALYs (per individual):
  arma::vec tc = m_C * v_dwc ;
  arma::vec te = m_E * v_dwe ;
  
  // average discounted costs and average discounted QALYs:
  double tc_hat = arma::mean(tc) ;
  double te_hat = arma::mean(te) ;
  
  // save the results in a list:
  Rcpp::List results = List::create(
    Named("m.M") = m_M, _["m.C"] = m_C, _["m.E"] = m_E, _["tc"] = tc, 
                       _["te"] = te, _["tc_hat"] = tc_hat, _["te_hat"] = te_hat) ;
  
  return(results) ;
}
