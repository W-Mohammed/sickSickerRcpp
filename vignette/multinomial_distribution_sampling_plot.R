library(ggplot2)
library(reshape2)

# Example transition probabilities
v_states_names <- c("H", "S1", "S2", "D")
m_trans_probs <- matrix(c(0.8, 0.15, 0.05, 0, 
                          0.1, 0.7, 0.1, 0.1,
                          0.05, 0.05, 0.8, 0.1,
                          0, 0, 0, 1), nrow = 4, byrow = TRUE)
colnames(m_trans_probs) <- v_states_names
rownames(m_trans_probs) <- v_states_names

# Calculate cumulative probabilities
m_cum_prob <- t(apply(m_trans_probs, 1, cumsum))

# Convert the matrix to data frame for plotting
df_cum_prob <- as.data.frame(m_cum_prob)
df_cum_prob$State <- rownames(m_cum_prob)

# Melt the data frame to long format for ggplot2
df_long <- melt(df_cum_prob, id.vars = "State", variable.name = "NextState", value.name = "CumProb")

# Random draws for demonstration
set.seed(123)  # For reproducibility
random_draws <- runif(length(unique(df_long$State)))  # One random draw per state
df_random_draws <- data.frame(State = unique(df_long$State), RandomDraw = random_draws)

# Create the plot
p <- ggplot() +
  geom_bar(data = df_long, aes(x = NextState, y = CumProb, fill = NextState), stat = "identity", position = "dodge") +
  geom_text(data = df_long, aes(x = NextState, y = CumProb, label = sprintf("%.2f", CumProb)), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(title = "Cumulative Transition Probabilities",
       x = "Next State",
       y = "Cumulative Probability",
       fill = "States") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add random draw lines
p <- p + geom_segment(data = df_random_draws, aes(x = 0, xend = State, y = RandomDraw, yend = RandomDraw),
                      linetype = "dashed", color = "black") +
  geom_segment(data = df_random_draws, aes(x = State, xend = State, y = RandomDraw, yend = 0),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "red")  +
  geom_text(data = df_random_draws, aes(x = 0, y = RandomDraw, label = sprintf("%.2f", RandomDraw)), hjust = 0.01, vjust = -0.5, color = "blue") +
  coord_cartesian(ylim = c(0, 1.2)) +
  facet_wrap(~State, scales = "free_y")

# Print the plot
print(p)
