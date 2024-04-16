# Parameters:----
n.i   <- 100000                # number of simulated individuals
# Utilities equation coefficients: intercept, age/10, (age/10)2, female
Utilities_coefficients <- c(0.8502, -0.0080, -0.0007, 0.0206)
# Individuals Characteristics:
age_mean <- 50
set.seed(1234)
individuals_ages <- rnorm(n.i, mean = age_mean, sd = age_mean)
female_proportion <- 0.4
set.seed(1234)
individuals_female <- sample(
  x = 0:1, 
  size = n.i, 
  replace = TRUE, 
  prob = c(1 - female_proportion, female_proportion)
)
Individuals_characteristics <- cbind(
  "intercept" = 1,
  "age_b10" = individuals_ages / 10,
  "age_b10_sq" = (individuals_ages / 10)^2,
  "female" = individuals_female
)
# R function:----
UtilitiesV <- function (Indi_characteristics, Utils_coefficients) {
  # Indi_characteristics: Individual characteristics
  # Utils_coefficients: Utilities' equation coefficients
  
  utilities <- Indi_characteristics %*% Utils_coefficients
  
  return(utilities)
}
# Code dependencies:----
depends <- c("RcppArmadillo")
plugins <- c("cpp11", "cpp17")
includes <- c(
  'using namespace Rcpp;
  using namespace arma;'
)

# C++ functions:----
## Code 1 - UtilitiesV:----
code1 <- 
  'arma::colvec UtilitiesV_Cpp( arma::mat m_I_characteristics,
                        arma::colvec v_U_coefficients) {
  // Indi_characteristics: Individual characteristics
  // Utils_coefficients: Utilities equation coefficients
  
  // Calculate individual level utilities by matrix-multiplication:
  arma::colvec m_I_Utilities = m_I_characteristics * v_U_coefficients ;
  
  return(m_I_Utilities) ;
}'

# Compile C++ code:----
cpp_functions_defs <- list(
  code1
)

for (code in cpp_functions_defs) {
  Rcpp::cppFunction(
    code = code,
    depends = depends,
    includes = includes#,
    #plugins = plugins
  )
}

# Test Rcpp functions:----
testR <- UtilitiesV(
  Indi_characteristics = Individuals_characteristics,
  Utils_coefficients = Utilities_coefficients
)
test1 <- UtilitiesV_Cpp(
  m_I_characteristics = Individuals_characteristics,
  v_U_coefficients = Utilities_coefficients
)
# Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "UtilitiesV_R" = UtilitiesV(
    Indi_characteristics = Individuals_characteristics,
    Utils_coefficients = Utilities_coefficients
  ),
  "UtilitiesV_Cpp" = UtilitiesV_Cpp(
    m_I_characteristics = Individuals_characteristics,
    v_U_coefficients = Utilities_coefficients
  )
)

# Print the results
print(results)

# For a visual comparison
boxplot(results)