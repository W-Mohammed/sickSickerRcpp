# MicroSim:----
## Parameters:----
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker
rp.S1S2 <- 0.2                 # increase of the mortality rate with every additional year being sick

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated
ru.S1S2 <- 0.01                # decrease in utility of treated sick individuals with every additional year being sick/sicker

# Define starting health state, using numbers instead of characters to identify the health states:
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p <- c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec <- c(c.H, c.S1, c.S2, 0, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.D", "c.Trt")
c_T_vec <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, 0)
names(c_T_vec) <- c("c.H", "c.S1", "c.S2", "c.D")

# Create a vector containing utilities parameters:
u_vec <- c(u.H, u.S1, u.S2, 0, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.D", "u.Trt")
u_T_vec <- c(u.H, u.Trt, u.S2, 0)
names(u_T_vec) <- c("u.H", "u.S1", "u.S2", "u.D")

# Create a vector containing dis-utilities parameters:
ru_vec <- c(0, ru.S1S2, ru.S1S2, 0)
names(ru_vec) <- c("ru.H", "ru.S1", "ru.S2", "ru.D")

# Transition matrix:
m_t_p <- matrix(
  data = c(1 - t_p["p.HS1"] - t_p["p.HD"], t_p["p.HS1"], 0, t_p["p.HD"],
           t_p["p.S1H"], 1 - t_p["p.S1H"] - t_p["p.S1S2"] - t_p["p.S1D"], 
           t_p["p.S1S2"],  t_p["p.S1D"], 0, 0, 1 - t_p["p.S2D"], t_p["p.S2D"], 0,
           0, 0, 1),
  nrow = n.s,
  ncol = n.s,
  byrow = TRUE
)

## Compile C++ function:----
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_1.0.cpp"
  )
)
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_2.0.cpp"
  )
)
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_3.0.cpp"
  )
)
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src/MicroSim_4.0.cpp"
  )
)
## Test function:----
results_1 <- MicroSimV_Cpp_1(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  b_Trt = FALSE,
  n_Seed = 1
)
results_1_T <- MicroSimV_Cpp_1(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  b_Trt = TRUE,
  n_Seed = 1
)
results_2 <- MicroSimV_Cpp_2(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_2_T <- MicroSimV_Cpp_2(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_T <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_2 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_T_2 <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_3_3_LYs <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0,
  d_dE = 0,
  n_Seed = 1
)
results_3_4_dLYs <- MicroSimV_Cpp_3(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_4 <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  v_d_Utilities = ru_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_Seed = 1
)
results_4_T <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  v_d_Utilities = ru_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  v_tracked_states = vector(),
  v_states_from = vector(),
  v_states_to = vector(),
  v_states_comp = vector(),
  v_increase_rate = vector(),
  n_Seed = 1
)
results_4_2 <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_4_T_2 <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_T_vec,
  v_Utilities = u_T_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
results_4_3_LYs <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0,
  d_dE = 0,
  n_Seed = 1
)
results_4_4_dLYs <- MicroSimV_Cpp_4(
  v_S_t = v_M_1,
  m_t_P = m_t_p,
  v_Costs = c_vec,
  v_Utilities = c(rep(1, n.s - 1), 0),
  n_I = n.i,
  v_tracked_states = c(2, 3),
  v_states_from = c(2, 3),
  v_states_to = c(4, 4),
  v_states_comp = c(2, 3),
  v_increase_rate = c(rp.S1S2, rp.S1S2),
  n_S = n.s,
  n_T = n.t,
  n_Cycle_length = 1,
  d_dC = 0.03,
  d_dE = 0.03,
  n_Seed = 1
)
## Compare functions:----
results <- microbenchmark::microbenchmark(
  times = 1000,
  "MicroSim_Cpp_1" = MicroSimV_Cpp_1( # 908ms median 955ms mean
    v_S_t = v_M_1,
    t_P = t_p,
    v_C = c_vec,
    v_U = u_vec,
    n_I = n.i,
    n_S = n.s,
    n_T = n.t,
    n_Cl = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    b_Trt = FALSE,
    n_Seed = 1
  ),
  "MicroSim_Cpp_2" = MicroSimV_Cpp_2( # 542ms median 577ms mean
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  ),
  "MicroSimV_Cpp_3_1" = MicroSimV_Cpp_3(
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    v_tracked_states = vector(),
    v_states_from = vector(),
    v_states_to = vector(),
    v_states_comp = vector(),
    v_increase_rate = vector(),
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  ),
  "MicroSimV_Cpp_3_2" = MicroSimV_Cpp_3(
    v_S_t = v_M_1,
    m_t_P = m_t_p,
    v_Costs = c_vec,
    v_Utilities = u_vec,
    n_I = n.i,
    v_tracked_states = c(2, 3),
    v_states_from = c(2, 3),
    v_states_to = c(4, 4),
    v_states_comp = c(2, 3),
    v_increase_rate = c(rp.S1S2, rp.S1S2),
    n_S = n.s,
    n_T = n.t,
    n_Cycle_length = 1,
    d_dC = 0.03,
    d_dE = 0.03,
    n_Seed = 1
  )
)
# Print the results
print(results)
# For a visual comparison
boxplot(results)