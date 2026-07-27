# Define target mean and a chosen standard deviation for the lognormal curve
target_mean <- 40
target_std <- 1 # You can adjust this spread

# Calculate the mu and sigma parameters
variance <- target_std^2
sigma2 <- log(variance / (target_mean^2) + 1)
sdlog <- sqrt(sigma2)
meanlog <- log(target_mean) - (sigma2 / 2)

# check
c(
  "mean" = meanlog,
  "sd" = sdlog
)

# Generate 10,000 random samples
samples <- rlnorm(n=10000, meanlog=meanlog, sdlog=sdlog)

# Verify the actual mean
cat("Empirical mean of samples:", mean(samples), "\n")

# Plot the distribution
hist(samples, breaks=50, col="lightblue", main="Log-Normal Distribution (Mean = 110)", xlab="Value")
