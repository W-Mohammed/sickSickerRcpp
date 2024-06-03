#include <RcppArmadillo.h>

// calc_effsC0
//
// Like the R's 'V' version of the function, 'C0' utilizes named vector
//
// [[Rcpp::export]]
Rcpp::NumericVector calc_effsC0(
    Rcpp::CharacterVector v_occupied_state,
    Rcpp::NumericVector v_states_utilities,
    Rcpp::NumericMatrix m_indi_features,
    Rcpp::NumericVector v_util_coeffs,
    Rcpp::NumericVector v_util_t_decs,
    Rcpp::NumericVector v_time_in_state,
    double cycle_length = 1) {
  
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

// calc_effsC1
//
// Limits data structures and functions used to those defined by Rcpp
//
// [[Rcpp::export]]
Rcpp::NumericVector calc_effsC1(Rcpp::IntegerVector v_occupied_state,
                                Rcpp::NumericVector v_states_utilities,
                                Rcpp::NumericMatrix m_indi_features,
                                Rcpp::NumericVector v_util_coeffs,
                                Rcpp::NumericVector v_util_t_decs,
                                Rcpp::NumericVector v_time_in_state,
                                double cycle_length = 1.0) {
  
  // Number of individuals
  int num_i = v_occupied_state.size();
  
  // Number of features
  int num_features = m_indi_features.ncol();
  
  // Calculate individual-specific utility decrements based on utilities regression coefficients
  // Rcpp::NumericVector v_ind_decrement(num_i, 0.0);
  Rcpp::NumericVector v_ind_decrement(num_i);
  for (int i = 0; i < num_i; ++i) {
    for (int j = 0; j < num_features; ++j) {
      v_ind_decrement[i] += m_indi_features(i, j) * v_util_coeffs[j];
    }
  }
  
  // Calculate time-dependent state-specific utility decrements
  Rcpp::NumericVector time_decrement(num_i, 0.0);
  for (int i = 0; i < num_i; ++i) {
    if (v_occupied_state[i] == 2) {
      time_decrement[i] = v_util_t_decs[0] * v_time_in_state[i];
    } else if (v_occupied_state[i] == 3) {
      time_decrement[i] = v_util_t_decs[1] * v_time_in_state[i];
    }
  }
  
  // Estimate total decrements
  Rcpp::NumericVector decrement = v_ind_decrement + time_decrement;
  
  // Estimate utilities based on occupied state
  Rcpp::NumericVector v_state_utility(num_i, 0.0);
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
  Rcpp::NumericVector QALYs = v_state_utility * cycle_length;
  
  return QALYs;
}

// calc_effsC2
//
// uses arma instead of Rcpp
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC2(arma::vec v_occupied_state,
                      arma::vec v_states_utilities,
                      arma::mat m_indi_features,
                      arma::vec v_util_coeffs,
                      arma::vec v_util_t_decs,
                      arma::vec v_time_in_state,
                      double cycle_length = 1.0) {
  
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
  arma::vec v_state_utility(num_i, arma::fill::zeros);
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

// calc_effsC3
//
//
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC3(const arma::vec& v_occupied_state,
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
  arma::vec v_state_utility(num_i, arma::fill::zeros);
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


// calc_effsC4
//
// arma::ivec is a vector of type integer, it saves memory 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC4(const arma::ivec& v_occupied_state,
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
  // use find() to get the indices where v_occupied_state matches specific values (2 and 3)
  // apply the corresponding decrements directly using elem()
  arma::vec time_decrement(num_i);
  time_decrement.elem(find(v_occupied_state == 2)) = v_util_t_decs[0] * v_time_in_state.elem(find(v_occupied_state == 2));
  time_decrement.elem(find(v_occupied_state == 3)) = v_util_t_decs[1] * v_time_in_state.elem(find(v_occupied_state == 3));
  
  // Estimate total decrements
  arma::vec decrement = v_ind_decrement + time_decrement;
  
  // Estimate utilities based on occupied state
  arma::vec v_state_utility = arma::zeros<arma::vec>(num_i);
  v_state_utility.elem(find(v_occupied_state == 1)) = v_states_utilities[0] + decrement.elem(find(v_occupied_state == 1));
  v_state_utility.elem(find(v_occupied_state == 2)) = v_states_utilities[1] + decrement.elem(find(v_occupied_state == 2));
  v_state_utility.elem(find(v_occupied_state == 3)) = v_states_utilities[2] + decrement.elem(find(v_occupied_state == 3));
  // v_state_utility.elem(find(v_occupied_state == 4)) is a subview and cannot be assigned a scalar directly, use fill()
  v_state_utility.elem(find(v_occupied_state == 4)).fill(v_states_utilities[3]);
  
  // Calculate Quality Adjusted Life Years (QALYs)
  arma::vec QALYs = v_state_utility * cycle_length;
  
  return QALYs;
}

// calc_effsC5
//
// Similar to C2, uses for loop and IF statement to assign costs but uses
// arma::find to calculate the cost regression estimate for each individual
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC5(arma::ivec v_occupied_state,
                      arma::vec v_states_utilities,
                      arma::mat m_indi_features,
                      arma::vec v_util_coeffs,
                      arma::vec v_util_t_decs,
                      arma::vec v_time_in_state,
                      double cycle_length = 1.0) {
  
  // Number of individuals
  int num_i = v_occupied_state.n_elem;
  
  // Calculate individual-specific utility decrements based on utilities regression coefficients
  arma::vec v_ind_decrement = m_indi_features * v_util_coeffs;
  
  // Calculate time-dependent state-specific utility decrements
  // use find() to get the indices where v_occupied_state matches specific values (2 and 3)
  // apply the corresponding decrements directly using elem()
  // arma::vec time_decrement = arma::zeros<arma::vec>(num_i);
  arma::vec time_decrement(num_i);
  time_decrement.elem(arma::find(v_occupied_state == 2)) = v_util_t_decs[0] * v_time_in_state.elem(arma::find(v_occupied_state == 2));
  time_decrement.elem(arma::find(v_occupied_state == 3)) = v_util_t_decs[1] * v_time_in_state.elem(arma::find(v_occupied_state == 3));
  
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

// calc_effsC6
//
// Similar to C2 but passes an ivec instead of vec for 'v_occupied_state' and
// passing objects by reference (&)
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC6(const arma::ivec& v_occupied_state,
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

// calc_effsC7
//
// Similar to C6 but uses arma::find() to assign individual-specific utilities. 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_effsC7(const arma::ivec& v_occupied_state,
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
  arma::vec v_state_utility(num_i, arma::fill::zeros);
  // Estimate utilities based on occupied state
  v_state_utility.elem(find(v_occupied_state == 1)) = v_states_utilities[0] + decrement.elem(find(v_occupied_state == 1));
  v_state_utility.elem(find(v_occupied_state == 2)) = v_states_utilities[1] + decrement.elem(find(v_occupied_state == 2));
  v_state_utility.elem(find(v_occupied_state == 3)) = v_states_utilities[2] + decrement.elem(find(v_occupied_state == 3));
  v_state_utility.elem(find(v_occupied_state == 4)).fill(v_states_utilities[3]);
  
  // Calculate Quality Adjusted Life Years (QALYs)
  arma::vec QALYs = v_state_utility * cycle_length;
  
  return QALYs;
}
