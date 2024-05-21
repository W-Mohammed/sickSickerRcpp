#include <RcppArmadillo.h>

// calc_costsC1
// [[Rcpp::depends(RcppArmadillo)]]
//
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
  Rcpp::NumericVector v_indi_costs(num_i, 0.0);
  for (int i = 0; i < num_i; ++i) {
    for (int j = 0; j < num_features; ++j) {
      v_indi_costs[i] += m_indi_features(i, j) * v_cost_coeffs[j];
    }
  }
  
  // Estimate costs based on occupied state
  Rcpp::NumericVector v_state_costs(num_i, 0.0);
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
// arma::ivec is a vector of type integer, it saves memory 
// [[Rcpp::depends(RcppArmadillo)]]
//
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
  arma::vec v_state_costs(num_i, arma::fill::zeros);
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
// [[Rcpp::depends(RcppArmadillo)]]
//
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
  arma::vec v_state_costs(num_i, arma::fill::zeros);
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
// [[Rcpp::depends(RcppArmadillo)]]
//
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
  arma::vec v_state_costs(num_i, arma::fill::zeros);  // By default, the cost for everyone is zero
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
// [[Rcpp::depends(RcppArmadillo)]]
//
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
  arma::vec v_state_costs(num_i, arma::fill::zeros);
  
  // Efficiently handle costs based on state
  v_state_costs.elem(find(v_occupied_state == 1)).fill(v_states_costs[0]);
  v_state_costs.elem(find(v_occupied_state == 2)) = v_states_costs[1] + v_indi_costs.elem(find(v_occupied_state == 2));
  v_state_costs.elem(find(v_occupied_state == 3)) = v_states_costs[2] + v_indi_costs.elem(find(v_occupied_state == 3));
  v_state_costs.elem(find(v_occupied_state == 4)).fill(v_states_costs[3]);
  
  return v_state_costs;
}