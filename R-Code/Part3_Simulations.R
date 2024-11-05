# Parameters for the AR(2)-ARMA(2,2) model
phi1 <- 0.1          # First AR coefficient for Phi_t
phi2 <- 0.2         # Second AR coefficient for Phi_t
theta1 <- 0.3        # First MA coefficient
theta2 <- 0.4        # Second MA coefficient
mu <- 0.5            # Mean of Phi_t
sigma_zeta <- 0.1    # Standard deviation of zeta_t
sigma_epsilon <- 0.2 # Standard deviation of epsilon_t

# Number of observations
N <- 200             # Length of the time series

# Pre-allocate vectors
Phi <- numeric(N)
X <- numeric(N)
epsilon <- numeric(N)

# Initial values
Phi[1] <- mu
Phi[2] <- mu
X[1] <- 0
X[2] <- 0
epsilon[1] <- 0
epsilon[2] <- 0

# Compute delta
delta <- mu * (1 - phi1 - phi2)
delta
set.seed(123)  # For reproducibility

# Generate noise terms
zeta <- rnorm(N, mean = 0, sd = sigma_zeta)         # State noise for Phi_t
epsilon_t <- rnorm(N, mean = 0, sd = sigma_epsilon)  # Observation noise

# Simulate the time series
for (t in 3:N) {
  # Update Phi_t with AR(2)-ARMA(2,2) structure
  Phi[t] <- delta + phi1 * Phi[t - 1] + phi2 * Phi[t - 2] +
    zeta[t] + theta1 * zeta[t - 1] + theta2 * zeta[t - 2]
  
  # Update X_t with AR(2)-MA(2) structure
  X[t] <- Phi[t] * X[t - 1] + Phi[t - 1] * X[t - 2] +
    theta1 * epsilon[t - 1] + theta2 * epsilon[t - 2] + epsilon_t[t]
  
  # Store epsilon_t for use in the next iteration
  epsilon[t] <- epsilon_t[t]
}

# Plot X_t
par(mfrow = c(1, 1))
plot(X, type = 'l', col = 'blue', main = 'Simulated Time Series X_t', ylab = expression(X[t]), xlab = 'Time', ylim=c(-5, 5))

# Plot Phi_t
par(mfrow = c(1, 1))
plot(Phi, type = 'l', col = 'red', main = 'Time-varying AR Coefficient Phi_t', ylab = expression(Phi[t]), xlab = 'Time', ylim=c(0, 0.8))
abline(h = mu, col = 'blue', lty = 2, lwd = 1)

# Combined plot of X_t and Phi_t
par(mfrow = c(2, 1))
plot(X, type = 'l', col = 'blue', main = 'Simulated Time Series X_t', ylab = expression(X[t]), xlab = 'Time', ylim=c(-5, 5))
plot(Phi, type = 'l', col = 'red', main = 'Time-varying AR Coefficient Phi_t', ylab = expression(Phi[t]), xlab = 'Time', ylim=c(0, 0.8))
abline(h = mu, col = 'blue', lty = 2, lwd = 1)




# Combined plot of X_t and Phi_t with dual y-axes

par(mfrow = c(1, 1))

# Plot X_t with its own y-axis (left axis)
plot(X, type = 'l', col = 'blue', main = 'Simulated Time Series X_t (observed) and Phi_t (hidden)',
     ylab = expression(X[t]), xlab = 'Time', ylim = c(min(X), max(X)))

# Add Phi_t on a secondary y-axis (right axis)
par(new = TRUE)
plot(Phi, type = 'l', col = 'red', axes = FALSE, xlab = '', ylab = '',
     ylim = c(min(Phi), max(Phi)))

# Add right y-axis for Phi_t
axis(4, col = 'red', col.axis = 'red')
mtext(expression(Phi[t]), side = 4, line = 3, col = 'red')

# Add legend to differentiate the lines
legend("topright",
       legend = c(bquote(X[t] ~ "(observed)"), bquote(Phi[t] ~ "(hidden)")),
       col = c("blue", "red"),
       lty = 1)

