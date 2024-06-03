#include <RcppArmadillo.h>

// calc_costsC0
//
// Like the R's 'V' version of the function, 'C0' utilizes named vector
//
// [[Rcpp::export]]
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

// calc_costsC1
//
// Limits the data structures and methods to those supported by Rcpp
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector calc_costsC1(Rcpp::IntegerVector v_occupied_state,
                                 Rcpp::NumericVector v_states_costs,
                                 Rcpp::NumericMatrix m_indi_features,
                                 Rcpp::NumericVector v_cost_coeffs) {
  
  // Number of individuals
  int num_i = v_occupied_state.size();
  
  // Number of features
  int num_features = m_indi_features.ncol();
  
  // Calculate individual-specific costs based on costs regression coefficients
  // Rcpp::NumericVector v_indi_costs(num_i, 0.0); if to be initialized 0s.
  Rcpp::NumericVector v_indi_costs(num_i);

  for (int i = 0; i < num_i; ++i) {
    for (int j = 0; j < num_features; ++j) {
      v_indi_costs[i] += m_indi_features(i, j) * v_cost_coeffs[j];
    }
  }
  
  // Estimate costs based on occupied state
  Rcpp::NumericVector v_state_costs(num_i);
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

// calc_costsC2
//
// Moves to using arma instead of Rcpp 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC2(arma::ivec& v_occupied_state,
                       arma::vec& v_states_costs,
                       arma::mat& m_indi_features,
                       arma::vec& v_cost_coeffs) {
  
  // Calculate individual-specific costs based on costs regression coefficients
  arma::vec v_indi_costs = m_indi_features * v_cost_coeffs;
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Estimate costs based on occupied state
  // arma::vec v_state_costs(num_i, arma::fill::zeros); if to be initialized 0s.
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

// calc_costsC3
//
// Suspends using ivec (integer vector) for the argument 'v_occupied_state'
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC3(arma::vec& v_occupied_state,
                       arma::vec& v_states_costs,
                       arma::mat& m_indi_features,
                       arma::vec& v_cost_coeffs) {
  
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

// calc_costsC4
//
// Passing constant values (objects) compared to the function 'C3'.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC4(const arma::ivec& v_occupied_state,
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

// calc_costsC5
//
// Removes loop and IF statements in assigning costs
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC5(const arma::ivec& v_occupied_state,
                       const arma::vec& v_states_costs,
                       const arma::mat& m_indi_features,
                       const arma::vec& v_cost_coeffs) {
  
  // Calculate individual-specific costs based on costs regression coefficients
  arma::vec v_indi_costs = m_indi_features * v_cost_coeffs;
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Estimate costs based on occupied state
  arma::vec v_state_costs(num_i);
  
  // Efficiently handle costs based on state
  v_state_costs.elem(find(v_occupied_state == 1)).fill(v_states_costs[0]);
  v_state_costs.elem(find(v_occupied_state == 2)) = v_states_costs[1] + v_indi_costs.elem(find(v_occupied_state == 2));
  v_state_costs.elem(find(v_occupied_state == 3)) = v_states_costs[2] + v_indi_costs.elem(find(v_occupied_state == 3));
  v_state_costs.elem(find(v_occupied_state == 4)).fill(v_states_costs[3]);
  
  return v_state_costs;
}

// calc_costsC6
//
// Replaces ivec with vec v_occupied_state
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC6(const arma::vec& v_occupied_state,
                       const arma::vec& v_states_costs,
                       const arma::mat& m_indi_features,
                       const arma::vec& v_cost_coeffs) {
  
  // Calculate individual-specific costs based on costs regression coefficients
  arma::vec v_indi_costs = m_indi_features * v_cost_coeffs;
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Estimate costs based on occupied state
  arma::vec v_state_costs(num_i);
  
  // Efficiently handle costs based on state
  v_state_costs.elem(find(v_occupied_state == 1)).fill(v_states_costs[0]);
  v_state_costs.elem(find(v_occupied_state == 2)) = v_states_costs[1] + v_indi_costs.elem(find(v_occupied_state == 2));
  v_state_costs.elem(find(v_occupied_state == 3)) = v_states_costs[2] + v_indi_costs.elem(find(v_occupied_state == 3));
  v_state_costs.elem(find(v_occupied_state == 4)).fill(v_states_costs[3]);
  
  return v_state_costs;
}

// calc_costsC7
//
// Similar to 'C6' but use a loop instead of arma::find to assign individual-
// specific costs 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_costsC7(const arma::vec& v_occupied_state,
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