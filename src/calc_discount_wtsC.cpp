#include <RcppArmadillo.h>

// calc_discount_wtsC1
// 
// Limited the data structures and methods to those supported by Rcpp, std to be
// specific here
// 
// [[Rcpp::export]]
std::vector<double> calc_discount_wtsC1(const double discount_rate, 
                                        const double num_cycles,
                                        const double cycle_length) {
  
  // Calculate the number of time points
  std::size_t size = static_cast<std::size_t>(num_cycles) + 1;
  
  // Prepare the vector of size num_cycles + 1
  std::vector<double> weights(size);
  
  // Fill the vector with time values and calculate discounting weights
  for (std::size_t i = 0; i < size; ++i) {
    double time = i * cycle_length;
    weights[i] = 1.0 / std::pow(1.0 + discount_rate, time);
  }
  
  return weights;
}

// calc_discount_wtsC2
// 
// Expands the used data structures and related methods to those supported by
// arma
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_discount_wtsC2(const double discount_rate, 
                              const double num_cycles,
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

// calc_discount_wtsC3
// 
// Uses int instead of double for num_cycles
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calc_discount_wtsC3(const double discount_rate, 
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

// calc_discount_wtsC4
// 
// Testing the effects of switching cycle_length to int from double. This change
// might be illogical since cycle_length could be a year, shorter or even longer
// 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> calc_discount_wtsC4(const double discount_rate, 
                                        const int num_cycles,
                                        const int cycle_length) {
  // Calculate the number of time points
  std::size_t size = static_cast<std::size_t>(num_cycles) + 1;
  
  // Prepare the vector of size time_horizon + 1
  std::vector<double> weights(size);
  
  // Pr+ecompute the discount factor
  double discount_factor = 1.0 / (1.0 + discount_rate);
  
  // Fill the vector with time values
  std::vector<double> times(size);
  for (std::size_t i = 0; i < size; ++i) {
    times[i] = i * cycle_length;
  }
  
  // Calculate discounting weights using a lambda function and std::transform
  std::transform(times.begin(), times.end(), weights.begin(),
                 [discount_factor](double time) -> double {
                   return std::pow(discount_factor, time);
                 });
  
  return weights;
}

// calc_discount_wtsC5
// 
// Similar to 'C1' but allows the limited use of arma
// 
// [[Rcpp::export]]
arma::vec calc_discount_wtsC5(const double discount_rate, 
                                        const double num_cycles,
                                        const double cycle_length) {
  
  // Calculate the number of time points
  std::size_t size = static_cast<std::size_t>(num_cycles) + 1;
  
  // Prepare the vector of size num_cycles + 1
  arma::vec weights(size);
  
  // Fill the vector with time values and calculate discounting weights
  for (std::size_t i = 0; i < size; ++i) {
    double time = i * cycle_length;
    weights[i] = 1.0 / std::pow(1.0 + discount_rate, time);
  }
  
  return weights;
}
