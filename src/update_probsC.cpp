#include <RcppArmadillo.h>

// update_probsC1
//
// Like the R's 'V' version of the function, 'C1' utilizes named vector
//
// [[Rcpp::export]]
Rcpp::NumericMatrix update_probsC1(const Rcpp::CharacterVector& v_states_names,
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

// update_probsC2
//
// Similar to 'C1' but stops using column and row names
//
// [[Rcpp::export]]
Rcpp::NumericMatrix update_probsC2(Rcpp::CharacterVector v_states_names,
                                   Rcpp::CharacterVector v_occupied_state,
                                   Rcpp::List l_trans_probs,
                                   Rcpp::NumericVector v_time_in_state) {
  
  int n = v_time_in_state.size();
  Rcpp::NumericMatrix m_probs(n, v_states_names.size());
  
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

// update_probsC3
// 
// Uses arma in addition to Rcpp.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_probsC3(Rcpp::CharacterVector v_states_names,
                         Rcpp::CharacterVector v_occupied_state,
                         Rcpp::List l_trans_probs,
                         Rcpp::NumericVector v_time_in_state) {
  
  int n = v_time_in_state.size();
  arma::mat m_probs(n, v_states_names.size());
  m_probs.zeros();
  
  double r_S1D = -std::log(1 - Rcpp::as<double>(l_trans_probs["p_S1D"]));
  double r_S2D = -std::log(1 - Rcpp::as<double>(l_trans_probs["p_S2D"]));
  arma::vec p_S1D = arma::zeros<arma::vec>(n);
  arma::vec p_S2D = arma::zeros<arma::vec>(n);
  
  for (int i = 0; i < n; ++i) {
    if (v_occupied_state[i] == "S1") {
      p_S1D(i) = 1 - std::exp(-r_S1D * (1 + v_time_in_state[i] * Rcpp::as<double>(l_trans_probs["rp_S1"])));
    } else if (v_occupied_state[i] == "S2") {
      p_S2D(i) = 1 - std::exp(-r_S2D * (1 + v_time_in_state[i] * Rcpp::as<double>(l_trans_probs["rp_S2"])));
    }
  }
  
  for (int i = 0; i < n; ++i) {
    if (v_occupied_state[i] == "H") {
      m_probs(i, 0) = 1 - Rcpp::as<double>(l_trans_probs["p_HS1"]) - Rcpp::as<double>(l_trans_probs["p_HD"]);
      m_probs(i, 1) = Rcpp::as<double>(l_trans_probs["p_HS1"]);
      m_probs(i, 3) = Rcpp::as<double>(l_trans_probs["p_HD"]);
    } else if (v_occupied_state[i] == "S1") {
      m_probs(i, 0) = Rcpp::as<double>(l_trans_probs["p_S1H"]);
      m_probs(i, 1) = 1 - Rcpp::as<double>(l_trans_probs["p_S1S2"]) - Rcpp::as<double>(l_trans_probs["p_S1H"]) - p_S1D(i);
      m_probs(i, 2) = Rcpp::as<double>(l_trans_probs["p_S1S2"]);
      m_probs(i, 3) = p_S1D(i);
    } else if (v_occupied_state[i] == "S2") {
      m_probs(i, 2) = 1 - p_S2D(i);
      m_probs(i, 3) = p_S2D(i);
    } else if (v_occupied_state[i] == "D") {
      m_probs(i, 3) = 1;
    }
  }
  
  // Sanity check
  for (int i = 0; i < n; ++i) {
    if (std::abs(arma::accu(m_probs.row(i)) - 1) > 1e-12) {
      Rcpp::stop("Probabilities do not sum to 1");
    }
  }
  
  return m_probs;
}

// update_probsC4
// 
// Switches from using characters and strings, v_states_names, to using numbers,
// v_states_index
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_probsC4(arma::vec v_states_index,
                         arma::vec v_occupied_state,
                         Rcpp::List l_trans_probs,
                         arma::vec v_time_in_state) {
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

// update_probsC5
//
// Passes const values or arguments compared to 'C4'.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_probsC5(const arma::vec& v_states_index,
                         const arma::vec& v_occupied_state,
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

// update_probsC6
//
// uses arma::ivec, a vector of type integer, which saves memory 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_probsC6(const arma::ivec& v_states_index,
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
