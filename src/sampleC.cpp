#include <RcppArmadillo.h>

// sampleC0
//
// Like the R's 'V' version of the function, 'C0' utilizes named vector
//
// [[Rcpp::export]]
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

// sampleC1
//
// Limits data structures and supported methods to those provided by 'Rcpp'.
// Stops using named vectors
//
// [[Rcpp::export]]
Rcpp::IntegerVector sampleC1(Rcpp::NumericMatrix m_trans_probs) {
  
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
  Rcpp::IntegerVector v_health_states(num_i);
  for (int i = 0; i < num_i; i++) {
    v_health_states[i] = 1 + v_transitions[i];
  }
  
  return v_health_states;
}

// sampleC2
//
// Uses data structures and functions defined by Armadillo (arma)
// Uses arma::trimatu to create an upper triangle and cumulative probabilities 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec sampleC2(const arma::mat& m_trans_probs) {
  
  // Number of individuals and states 
  int num_i = m_trans_probs.n_rows;
  int num_states = m_trans_probs.n_cols;
  
  // Create an upper triangular matrix of ones, 'trimatu' copies the upper triangle
  arma::mat m_upper_tri = arma::trimatu(arma::ones<arma::mat>(num_states, num_states));
  
  // Create matrix with row-wise cumulative transition probabilities
  arma::mat m_cum_probs = m_trans_probs * m_upper_tri;
  
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
  arma::uvec v_health_states = 1 + v_transitions;
  
  return v_health_states;
}

// sampleC3
//
// Returns arma::colv instead of arma's unsigned vector arma::uvec
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec sampleC3(const arma::mat& m_trans_probs) {
  
  // Number of individuals and states 
  int num_i = m_trans_probs.n_rows;
  int num_states = m_trans_probs.n_cols;
  
  // Create an upper triangular matrix of ones, 'trimatu' copies the upper triangle
  arma::mat m_upper_tri = arma::trimatu(arma::ones<arma::mat>(num_states, num_states));
  
  // Create matrix with row-wise cumulative transition probabilities
  arma::mat m_cum_probs = m_trans_probs * m_upper_tri;
  
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
  arma::colvec v_health_states = 1 + arma::conv_to<arma::vec>::from(v_transitions);
  
  return v_health_states;
}

// sampleC4
//
// Uses arma::cumsum(m_trans_probs) instead of 'arma::trimatu * m_trans_probs' 
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec sampleC4(const arma::mat& m_trans_probs) {
  
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
  arma::colvec v_health_states = 1 + arma::conv_to<arma::vec>::from(v_transitions);
  
  return v_health_states;
}

// sampleC5
//
// Returns arma::ivec instead of arma::colvec
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::ivec sampleC5(const arma::mat& m_trans_probs) {
  
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