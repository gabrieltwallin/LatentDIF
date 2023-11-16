

rm(list=ls())
graphics.off()

library(LaplacesDemon)
library(varhandle)


#--------------------------------------
# Multigroup latent DIF implementation
#--------------------------------------


N = 1000
J = 50
p = 20 # p=10 when J=25, and p=20 when J=50

pi = c(0.5, 0.3, 0.2) 

pi2 = 0.3
pi3 = 0.2

step <- 1
lambda <- 0.1

# Respondent class threshold:
class_thre <- 0.5

# Number of classes: 
K <- 3


# Number of replications
sim <- 100

bias_mean1 <- c()
bias_mean2 <- c()
bias_sd1 <- c()
bias_sd2 <- c()
bias_pi2 <- c()
bias_pi3 <- c()
bias_dif1 <- matrix(0, nrow = J, ncol = sim)
bias_dif2 <- matrix(0, nrow = J, ncol = sim)
bias_d.vec <- matrix(0, nrow = J, ncol = sim)
bias_a.vec <- matrix(0, nrow = J, ncol = sim)

mse_mean1 <- c()
mse_mean2 <- c()
mse_sd1 <- c()
mse_sd2 <- c()
mse_pi2 <- c()
mse_pi3 <- c()
mse_dif1 <- matrix(0, nrow = J, ncol = sim)
mse_dif2 <- matrix(0, nrow = J, ncol = sim)
mse_d.vec <- matrix(0, nrow = J, ncol = sim)
mse_a.vec <- matrix(0, nrow = J, ncol = sim)

est_mean1 <- c()
est_mean2 <- c()
est_sd1 <- c()
est_sd2 <- c()
est_pi2 <- c()
est_pi3 <- c()
est_d.vec <- matrix(0, nrow = J, ncol = sim)
est_a.vec <- matrix(0, nrow = J, ncol = sim)


true_positives_delta1 <- rep(0, sim)
true_positives_delta2 <- rep(0, sim)

false_positives_delta1 <- rep(0, sim)
false_positives_delta2 <- rep(0, sim)


true_positives_xi1.true <- rep(0, sim)
true_positives_xi2.true <-rep(0, sim)

false_positives_xi1.true <- rep(0, sim)
false_positives_xi2.true <- rep(0, sim)

classification_error <- rep(0, sim)
classification_error_true <- rep(0, sim)

# Initialize list to store ROC curves for each iteration
roc_curves_est_class2 <- list()
roc_curves_est_class3 <- list()

roc_curves_true_class2 <- list()
roc_curves_true_class3 <- list()


calculate_roc_curve <- function(class_assign, xi_vec, thresholds) {
  fpr <- numeric(length(thresholds))
  tpr <- numeric(length(thresholds))
  
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    predicted <- ifelse(class_assign > threshold, 1, 0)
    tp <- sum(predicted == 1 & xi_vec == 1)
    fp <- sum(predicted == 1 & xi_vec == 0)
    fn <- sum(predicted == 0 & xi_vec == 1)
    tn <- sum(predicted == 0 & xi_vec == 0)
    fpr[i] <- fp / (fp + tn)
    tpr[i] <- tp / (tp + fn)
  }
  
  return(data.frame(fpr, tpr))
}


set.seed(776155)

# Discrimination parameter 
a.vec <- runif(J, 0.5, 1.5)

# Item easiness (intercept parameter)
d.vec <- runif(J, -2, 2)

# DIF effect parameter (only the first p=5 items have DIF):
K = 3 # Number of latent classes, incl. baseline 
delta.mat = matrix(cbind(rep(0, J), # Baseline class
                         c(runif(p, 0.5, 1), rep(0, J-p)), # Latent class 1
                         c(runif(p, 1, 1.5), rep(0, J-p))), # Latent class 2
                   nrow = J, ncol = K)

colnames(delta.mat) = c("Baseline", "LC1", "LC2")

delta1.vec = delta.mat[,2]
delta2.vec = delta.mat[,3]

delta1.vec_binary <- rep(0, J)
delta2.vec_binary <- rep(0, J)

delta1.vec_binary[delta1.vec > 0] <- 1
delta2.vec_binary[delta2.vec > 0] <- 1

# Latent classes (0 is baseline group)
xi.vec = rcat(n=N, p=pi)
xi.vec = xi.vec - 1 

# Ability parameter in the baseline group:
theta.vec <- rnorm(N, 0, 1)

# Mean and standard deviation of ability in latent class 1
mu1 = 1
sigma1 = 1.2 

# Mean and standard deviation of ability in latent class 2
mu2 = 0.5
sigma2 = 0.8

# Ability parameter in the latent classes
theta.vec[xi.vec==1] = rnorm(sum(xi.vec==1), mu1, sigma1)
theta.vec[xi.vec==2] = rnorm(sum(xi.vec==2), mu2, sigma2)

# Transform LCs to 0/1 columns, with ncol = K
xi.mat <- to.dummy(xi.vec, "LC")

# Row constant
row_cons = rep(1, N)

# Gather row constant, theta and the LCs
theta.mat = cbind(row_cons, theta.vec, xi.mat)

# Column constant
col_cons = rep(1, J)

# Gather item easiness, column constant and DIF effects for each LC
A.mat = cbind(d.vec, col_cons, delta.mat)

# Probability of correct answer according to the Rasch model
temp = theta.mat %*% t(A.mat)
prob = exp(temp) / (1 + exp(temp))


for(j in 1:sim){
  # Generate 0/1 data according to the specified Rasch model:
  data <- matrix(0, N, J)
  data[] <- rbinom(N*J, 1, prob)
  
  # Specify grid for the Gauss-Hermite quadrature of theta integration: 
  grid <- seq(-4, 4, length.out = 31)
  
  # Weights for the grid:
  weight = dnorm(grid)
  
  # Normalize:
  weight = weight/sum(weight)
  
  
  # Expected complete data log-likelihood
  llik <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2){
    # log-likelihood for baseline group
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31))
    prob = exp(temp) / (1 + exp(temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = exp(temp1) %*% weight
    
    # log-likelihood for latent class 1
    grid.1 = grid * sigma1 + mu1
    temp = a.vec %*% t(grid.1) + (d.vec+delta1.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = exp(temp2) %*% weight
    
    # log-likelihood for latent class 3
    grid.2 = grid * sigma2 + mu2
    temp = a.vec %*% t(grid.2) + (d.vec+delta2.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp3 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp3 = exp(temp3) %*% weight
    
    sum(log(temp1 * (1-(pi2+pi3)) + temp2 * pi2 + temp3 * pi3))
  }
  
  
  # Matrix of posterior probabilities for class allocation:  
  post.matr <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2){
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31))
    prob = exp(temp) / (1 + exp(temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = t(t(exp(temp1)) * weight) * (1-(pi2 + pi3))
    
    grid.1 = grid * sigma1 + mu1
    temp = a.vec %*% t(grid.1) + (d.vec+delta1.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = t(t(exp(temp2)) * weight) * pi2
    
    grid.2 = grid * sigma2 + mu2
    temp = a.vec %*% t(grid.2) + (d.vec+delta2.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp3 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp3 = t(t(exp(temp3)) * weight) * pi3
    
    post <- cbind(temp1, temp2, temp3)
    
    post = post/rowSums(post)
    post
  }
  
  
  
  classify_respondents <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) {
    # Calculate probability of belonging to first class
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = t(t(exp(temp1)) * weight) * (1-(pi2 + pi3))
    
    # Calculate probability of belonging to second class
    grid.1 = grid * sigma1 + mu1
    temp = a.vec %*% t(grid.1) + (d.vec+delta1.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = t(t(exp(temp2)) * weight) * pi2
    
    
    # Calculate probability of belonging to third class
    grid.2 = grid * sigma2 + mu2
    temp = a.vec %*% t(grid.2) + (d.vec+delta2.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp3 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp3 = t(t(exp(temp3)) * weight) * pi3
    
    # Combine probabilities
    post <- cbind(temp1, temp2, temp3)
    
    # Normalize probabilities
    post = post/rowSums(post)
    
    # Assign class based on the highest probability
    #class_assignments <- ifelse(rowSums(post[,1:31]) > rowSums(post[,32:62]), 0, 1)
    #return(class_assignments)
    return(list(class1 = rowSums(post[,32:62]), class2 = rowSums(post[,63:93])))
  }
  
  class.assign_tmp_true <- classify_respondents(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2)
  
  class1.assign_TruParVal <- ifelse(class.assign_tmp_true[[1]] > class_thre, 1, 0)
  
  class2.assign_TruParVal <- ifelse(class.assign_tmp_true[[2]] > class_thre, 1, 0)
  
  
  predicted_labels_true <- ifelse(class1.assign_TruParVal == 1, 1, ifelse(class2.assign_TruParVal == 1, 2, 0))
  
  # Calculate the number of incorrect predictions
  incorrect_predictions_true <- sum(xi.vec != predicted_labels_true)
  
  # Calculate the total number of predictions
  total_predictions <- length(xi.vec)
  
  # Calculate the classification error
  classification_error_true[j] <- incorrect_predictions_true / total_predictions
  
  
  
  # Calculate true positive rate for xi for the true model
  true_positives_xi1.true[j] <- sum((xi.vec == 1) & (class1.assign_TruParVal == 1)) / sum(xi.vec == 1)
  
  true_positives_xi2.true[j] <- sum((xi.vec == 2) & (class2.assign_TruParVal == 1)) / sum(xi.vec == 2)
  
  # Calculate false positive rate for xi for the true model
  false_positives_xi1.true[j] <- sum((xi.vec == 0) & (class1.assign_TruParVal == 1)) / sum(xi.vec == 0)
  
  # Calculate false positive rate for xi for the true model
  false_positives_xi2.true[j] <- sum((xi.vec == 0) & (class2.assign_TruParVal == 1)) / sum(xi.vec == 0)
  
  
  
  # Calculating the class proportion:
  Update.pi2 <- function(post){
    
    temp1 = colMeans(post[,1:31])
    temp2 = colMeans(post[,32:62])
    temp3 = colMeans(post[,63:93])
    
    return(sum(temp2))
  }
  
  # Calculating the class proportion:
  Update.pi3 <- function(post){
    
    temp1 = colMeans(post[,1:31])
    temp2 = colMeans(post[,32:62])
    temp3 = colMeans(post[,63:93])
    
    return(sum(temp3))
  }
  
  
  # Soft thresholding function:
  soft_thre <- function(x, lambda){
    temp = x
    temp[abs(x) <= lambda] = 0
    temp[x > lambda] = temp[x > lambda] - lambda
    temp[x < -lambda] = temp[x < -lambda] + lambda
    temp
  }
  
  # Proximal gradient descent
  prox.grad <- function(x, grad, lambda, step){
    soft_thre(x - step * grad, lambda * step)
  }
  
  
  line.search <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, 
                          mu1, mu2, sigma1, sigma2, grad.a, grad.d, 
                          grad.delta1, grad.delta2, grad.sigma1, grad.sigma2, 
                          grad.mu1, grad.mu2, lambda, step=1){
    
    
    # Set the initial step size to the input step size
    step.size <- step
    
    # Set the initial value of the objective function to the current value
    obj.val <- -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) + lambda * sum(abs(delta1.vec)) + lambda*sum(abs(delta2.vec))
    
    # Set the tolerance for the line search algorithm
    tol <- 1e-6
    
    # Set the maximum number of iterations for the line search algorithm
    max.iter <- 100
    
    # Set the initial number of iterations to 0
    iter <- 0
    
    # Iteratively adjust the step size and evaluate the performance of the algorithm
    while (iter < max.iter) {
      # Update the parameters using the current step size
      a.new <- a.vec - step.size * grad.a  
      d.new <- d.vec - step.size * grad.d 
      delta1.new <- prox.grad(delta1.vec, grad.delta1, lambda, step.size)
      delta2.new <- prox.grad(delta2.vec, grad.delta2, lambda, step.size)
      sigma1.new <- sigma1 - step.size * grad.sigma1
      sigma2.new <- sigma2 - step.size * grad.sigma2
      mu1.new <- mu1 - step.size * grad.mu1 
      mu2.new <- mu2 - step.size * grad.mu2
      
      
      # Calculate the new value of the objective function
      new.obj.val <- -llik(data, a.new, d.new, delta1.new, delta2.new, pi2, pi3, mu1.new, mu2.new, sigma1.new, sigma2.new) + 
        lambda * sum(abs(delta1.new)) + lambda*sum(abs(delta2.new))
      
      # Check if the change in the objective function is within the tolerance
      if (abs(new.obj.val - obj.val) < tol) {
        # If the change is within the tolerance, return the current step size
        return(step.size)
      } else {
        # If the change is not within the tolerance, adjust the step size
        # and update the value of the objective function
        step.size <- step.size / 2
        obj.val <- new.obj.val
        iter <- iter + 1
      }
    }
    # If the maximum number of iterations is reached, return the current step size
    return(step.size)
  }
  
  
  
  line.search.conf <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, 
                               mu1, mu2, sigma1, sigma2, grad.a, grad.d, grad.delta1, 
                               grad.delta2, grad.sigma1, grad.sigma2, grad.mu1, grad.mu2, step=1) {
    
    # Set the initial step size to the input step size
    step.size <- step
    
    # Set the initial value of the objective function to the current value
    obj.val <- -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
    
    # Set the tolerance for the line search algorithm
    tol <- 1e-6
    
    # Set the maximum number of iterations for the line search algorithm
    max.iter <- 100
    
    # Set the initial number of iterations to 0
    iter <- 0
    
    # Iteratively adjust the step size and evaluate the performance of the algorithm
    while (iter < max.iter) {
      # Update the parameters using the current step size
      a.new <- a.vec - step.size * grad.a  
      d.new <- d.vec - step.size * grad.d 
      delta1.new <- delta1.vec - step.size * grad.delta1
      delta2.new <- delta2.vec - step.size * grad.delta2
      sigma1.new <- sigma1 - step.size * grad.sigma1 
      sigma2.new <- sigma2 - step.size * grad.sigma2
      mu1.new <- mu1 - step.size * grad.mu1 
      mu2.new <- mu2 - step.size * grad.mu2 
      
      # Calculate the new value of the objective function
      new.obj.val <- -llik(data, a.new, d.new, delta1.new, delta2.new, pi2, pi3, mu1.new, mu2.new, sigma1.new, sigma2.new)
      
      # Check if the change in the objective function is within the tolerance
      if (abs(new.obj.val - obj.val) < tol) {
        # If the change is within the tolerance, return the current step size
        return(step.size)
      } else {
        # If the change is not within the tolerance, adjust the step size
        # and update the value of the objective function
        step.size <- step.size / 2
        obj.val <- new.obj.val
        iter <- iter + 1
      }
    }
    # If the maximum number of iterations is reached, return the current step size
    return(step.size)
  }
  
  
  
  # Update the parameter estimates through gradient descent and proximal gradient descent:
  Update.par <- function(post, data, a.vec, d.vec, delta1.vec, delta2.vec, 
                         mu1, mu2, sigma1, sigma2, lambda, step){
    
    # Posterior probabilities, baseline group  
    temp1 = post[,1:31]
    
    # Posterior probabilities, LC1
    temp2 = post[,32:62]
    
    # Posterior probabilities, LC2
    temp3 = post[,63:93]
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31))
    prob = exp(temp) / (1 + exp(temp))
    
    # log-likelihood for LC1
    grid.1 = grid * sigma1 + mu1
    temp = a.vec %*% t(grid.1) + (d.vec+delta1.vec) %*% t(rep(1, 31)) 
    prob2 = 1/(1+exp(-temp))
    
    # log-likelihood for LC3
    grid.2 = grid * sigma2 + mu2
    temp = a.vec %*% t(grid.2) + (d.vec+delta2.vec) %*% t(rep(1, 31)) 
    prob3 = 1/(1+exp(-temp))
    
    
    
    # Gradients of the parameters
    grad.d = rowSums(t(data) %*% temp1) - prob %*% (colSums(temp1) ) + 
      rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2)) + 
      rowSums(t(data) %*% temp3) - prob3 %*% (colSums(temp3))
    grad.d = -(grad.d/N)
    
    
    grad.a = t(data) %*% temp1 %*% grid - prob %*% (colSums(temp1) * grid) +  
      t(data) %*% temp2 %*% grid.1 - prob2 %*% (colSums(temp2) * grid.1) + 
      t(data) %*% temp3 %*% grid.2 - prob3 %*% (colSums(temp3) * grid.2)
    grad.a = -(grad.a/N)
    
    
    grad.delta1 = rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.delta1 = -(grad.delta1/N)
    
    grad.delta2 = rowSums(t(data) %*% temp3) - prob3 %*% (colSums(temp3))
    grad.delta2 = -(grad.delta2/N)
    
    
    grad.sigma1 = t(a.vec) %*% (t(data) %*% temp2 %*% grid) - t(a.vec)%*% prob2 %*% (colSums(temp2)*grid) 
    grad.sigma1 = -grad.sigma1/N
    
    grad.sigma2 = t(a.vec) %*% (t(data) %*% temp3 %*% grid) - t(a.vec)%*% prob3 %*% (colSums(temp3)*grid) 
    grad.sigma2 = -grad.sigma2/N
    
    
    grad.mu1 = sum(t(a.vec) %*% (t(data) %*% temp2)) - t(a.vec)%*% prob2 %*% (colSums(temp2)) 
    grad.mu1 = -grad.mu1/N
    
    grad.mu2 = sum(t(a.vec) %*% (t(data) %*% temp3)) - t(a.vec)%*% prob3 %*% (colSums(temp3)) 
    grad.mu2 = -grad.mu2/N
    
    
    # Select step size 
    step <- line.search(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, 
                        mu1, mu2, sigma1, sigma2, grad.a, grad.d, grad.delta1, 
                        grad.delta2, grad.sigma1, grad.sigma2, grad.mu1, grad.mu2, 
                        lambda, step)
    
    
    # Gradient descent to update parameters:
    d.new = d.vec - step*grad.d
    
    a.new = a.vec - step*grad.a
    
    delta1.new = prox.grad(delta1.vec, grad.delta1, lambda, step)
    
    delta2.new = prox.grad(delta2.vec, grad.delta2, lambda, step)
    
    mu1.new = mu1 - step*grad.mu1
    
    mu2.new = mu2 - step*grad.mu2
    
    sigma1.new = sigma1 - step*grad.sigma1
    
    sigma2.new = sigma2 - step*grad.sigma2
    
    
    
    list(a.vec = a.new, d.vec = d.new, delta1.vec = delta1.new, delta2.vec = delta2.new,
         mu1 = mu1.new[1], mu2 = mu2.new[1], sigma1 = sigma1.new[1], sigma2 = sigma2.new[1])  
  }
  
  
  
  Update.par_conf <- function(post, data, a.vec, d.vec, delta1.vec, delta2.vec, mu1, mu2, sigma1, sigma2, step){
    
    # Set delta for the first number of items to 0
    delta1.vec <- ifelse(delta1.vec != 0, delta1.vec, 0)
    delta2.vec <- ifelse(delta2.vec != 0, delta2.vec, 0)
    
    
    # Posterior probabilities, baseline group  
    temp1 = post[,1:31]
    
    # Posterior probabilities, LC1
    temp2 = post[,32:62]
    
    # Posterior probabilities, LC2
    temp3 = post[,63:93]
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    
    # log-likelihood for LC1
    grid.1 = grid * sigma1 + mu1
    temp = a.vec %*% t(grid.1) + (d.vec+delta1.vec) %*% t(rep(1, 31)) 
    prob2 = 1/(1+exp(-temp))
    
    # log-likelihood for LC3
    grid.2 = grid * sigma2 + mu2
    temp = a.vec %*% t(grid.2) + (d.vec+delta2.vec) %*% t(rep(1, 31)) 
    prob3 = 1/(1+exp(-temp))
    
    
    
    # Gradients of the parameters
    grad.d = rowSums(t(data) %*% temp1) - prob %*% (colSums(temp1) ) + 
      rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2)) + 
      rowSums(t(data) %*% temp3) - prob3 %*% (colSums(temp3))
    grad.d = -(grad.d/N)
    
    
    grad.a = t(data) %*% temp1 %*% grid - prob %*% (colSums(temp1) * grid) +  
      t(data) %*% temp2 %*% grid.1 - prob2 %*% (colSums(temp2) * grid.1) + 
      t(data) %*% temp3 %*% grid.2 - prob3 %*% (colSums(temp3) * grid.2)
    grad.a = -(grad.a/N)
    
    grad.delta1 = rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.delta1 = -(grad.delta1/N)
    
    grad.delta2 = rowSums(t(data) %*% temp3) - prob3 %*% (colSums(temp3))
    grad.delta2 = -(grad.delta2/N)
    
    
    grad.sigma1 = t(a.vec) %*% (t(data) %*% temp2 %*% grid) - t(a.vec)%*% prob2 %*% (colSums(temp2)*grid) 
    grad.sigma1 = -grad.sigma1/N
    
    grad.sigma2 = t(a.vec) %*% (t(data) %*% temp3 %*% grid) - t(a.vec)%*% prob3 %*% (colSums(temp3)*grid) 
    grad.sigma2 = -grad.sigma2/N
    
    
    grad.mu1 = sum(t(a.vec) %*% (t(data) %*% temp2)) - t(a.vec)%*% prob2 %*% (colSums(temp2)) 
    grad.mu1 = -grad.mu1/N
    
    grad.mu2 = sum(t(a.vec) %*% (t(data) %*% temp3)) - t(a.vec)%*% prob3 %*% (colSums(temp3)) 
    grad.mu2 = -grad.mu2/N
    
    
    # Select step size 
    step <- line.search.conf(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, 
                             mu1, mu2, sigma1, sigma2, grad.a, grad.d, grad.delta1, 
                             grad.delta2, grad.sigma1, grad.sigma2, grad.mu1, grad.mu2, step)
    
    
    # Gradient descent to update parameters:
    d.new = d.vec - step*grad.d
    
    a.new = a.vec - step*grad.a
    
    #delta1.new = delta1.vec - step*grad.delta1
    delta1.new <- ifelse(delta1.vec != 0, delta1.vec - step * grad.delta1, 0)
    
    delta2.new <- ifelse(delta2.vec != 0, delta2.vec - step * grad.delta2, 0)
    
    # delta2.new = delta2.vec - step*grad.delta2
    
    mu1.new = mu1 - step*grad.mu1 
    
    mu2.new = mu2 - step*grad.mu2
    
    sigma1.new = sigma1 - step*grad.sigma1
    
    sigma2.new = sigma2 - step*grad.sigma2  
    
    
    list(a.vec = a.new, d.vec = d.new, delta1.vec = delta1.new, delta2.vec = delta2.new,
         mu1 = mu1.new[1], mu2 = mu2.new[1], sigma1 = sigma1.new[1], sigma2 = sigma2.new[1])  
  }
  
  
  # EM algorithm:
  EM <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2, lambda, step){
    
    # Objective function: 
    obj0 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) + 
      lambda * sum(abs(delta1.vec)) + lambda * sum(abs(delta2.vec)) 
    
    # Posterior probabilities class membership: 
    post = post.matr(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
    
    pi2 = Update.pi2(post)
    pi3 = Update.pi3(post)
    
    # Run gradient descent 10 steps
    for(i in 1:10){
      res = Update.par(post, data, a.vec, d.vec, delta1.vec, delta2.vec, mu1, mu2, sigma1, sigma2, lambda, step)  
      a.vec = res$a.vec
      d.vec = res$d.vec
      delta1.vec = res$delta1.vec
      delta2.vec = res$delta2.vec
      mu1 = res$mu1
      mu2 = res$mu2
      sigma1 = res$sigma1
      sigma2 = res$sigma2
    }
    
    # Update value of objective function:
    obj1 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) + 
      lambda * sum(abs(delta1.vec)) + lambda * sum(abs(delta2.vec))
    
    while(obj0-obj1 > 1e-4){
      
      obj0=obj1
      
      post = post.matr(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
      
      pi2 = Update.pi2(post)
      pi3 = Update.pi3(post)
      
      # Run gradient descent 10 steps
      for(i in 1:10){
        res = Update.par(post, data, a.vec, d.vec, delta1.vec, delta2.vec, mu1, mu2, sigma1, sigma2, lambda, step)  
        a.vec = res$a.vec
        d.vec = res$d.vec
        delta1.vec = res$delta1.vec
        delta2.vec = res$delta2.vec
        mu1 = res$mu1
        mu2 = res$mu2
        sigma1 = res$sigma1
        sigma2 = res$sigma2
      }
      
      obj1 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) + 
        lambda * sum(abs(delta1.vec)) + lambda * sum(abs(delta2.vec))
      print(obj1)
    }
    list(a.vec = a.vec, d.vec = d.vec, delta1.vec = delta1.vec, delta2.vec = delta2.vec, 
         mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, pi2 = pi2, pi3 = pi3, post = post)
  }
  
  
  # EM algorithm:
  EM_conf <- function(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2, step){
    
    # Objective function: 
    obj0 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
    
    # Posterior probabilities class membership: 
    post = post.matr(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
    
    pi2 = Update.pi2(post)
    pi3 = Update.pi3(post)
    
    # Run gradient descent 10 steps
    for(i in 1:10){
      res = Update.par_conf(post, data, a.vec, d.vec, delta1.vec, delta2.vec, mu1, mu2, sigma1, sigma2, step=1)  
      a.vec = res$a.vec
      d.vec = res$d.vec
      delta1.vec = res$delta1.vec
      delta2.vec = res$delta2.vec
      mu1 = res$mu1
      mu2 = res$mu2
      sigma1 = res$sigma1
      sigma2 = res$sigma2
    }
    
    # Update value of objective function:
    obj1 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
    
    while(obj0-obj1 > 1e-4){
      
      obj0=obj1
      
      post = post.matr(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
      
      pi2 = Update.pi2(post)
      pi3 = Update.pi3(post)
      
      # Run gradient descent 10 steps
      for(i in 1:10){
        res = Update.par_conf(post, data, a.vec, d.vec, delta1.vec, delta2.vec, mu1, mu2, sigma1, sigma2, step=1)  
        a.vec = res$a.vec
        d.vec = res$d.vec
        delta1.vec = res$delta1.vec
        delta2.vec = res$delta2.vec
        mu1 = res$mu1
        mu2 = res$mu2
        sigma1 = res$sigma1
        sigma2 = res$sigma2
      }
      
      obj1 = -llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) 
      print(obj1)
    }
    
    bic=-2*llik(data, a.vec, d.vec, delta1.vec, delta2.vec, pi2, pi3, mu1, mu2, sigma1, sigma2) + log(N)*(J+J+6 + sum(delta1.vec!=0) + sum(delta2.vec!=0))
    
    list(a.vec = a.vec, d.vec = d.vec, delta1.vec = delta1.vec, delta2.vec = delta2.vec, 
         mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, pi2 = pi2, pi3 = pi3, post = post, bic = bic)
  }
  
  
  a.vec0 = a.vec+rnorm(J, 0, 0.1)
  d.vec0 <- d.vec+rnorm(J, 0, 0.1)
  mu10 <- mu1 + rnorm(1, 0, 0.1)
  mu20 <- mu2 + rnorm(1, 0, 0.1)
  sigma10 <- sigma1 + rnorm(1, 0, 0.1)
  sigma20 <- sigma2 + rnorm(1, 0, 0.1)
  delta1.vec0 = delta1.vec+rnorm(J, 0, 0.1)
  delta2.vec0 = delta2.vec+rnorm(J, 0, 0.1)
  
  lambda.vec <- seq(0.01, 0.0005, by = -0.001)
  
  res = EM(data, a.vec0, d.vec0, delta1.vec0, delta2.vec0, pi2, pi3, mu10, mu20, sigma10, sigma20, lambda = lambda.vec[1], step)
  
  
  
  bic = Inf
  print('bic')
  
  c = 0.5
  
  for(i in 2:length(lambda.vec)){
    res <- EM(data, res$a.vec, res$d.vec, res$delta1.vec, res$delta2.vec, res$pi2, res$pi3, res$mu1, res$mu2, res$sigma1, res$sigma2, lambda = lambda.vec[i], step)
    
    delta1_lambda0 <- res$delta1.vec
    delta1_lambda0[delta1_lambda0 < c] = 0
    delta2_lambda0 <- res$delta2.vec
    delta2_lambda0[delta2_lambda0 < c] = 0
    
    res.conf <- EM_conf(data, res$a.vec, res$d.vec, delta1_lambda0, delta2_lambda0, res$pi2, res$pi3, res$mu1, res$mu2, res$sigma1, res$sigma2, step)
    
    bic_new = res.conf$bic
    
    if(bic_new < bic){
      bic = bic_new
      delta1_bic = res.conf$delta1.vec
      delta2_bic = res.conf$delta2.vec
      #print(delta_bic)
      a.vec_bic = res.conf$a.vec
      d.vec_bic = res.conf$d.vec
      mu1_bic = res.conf$mu1
      mu2_bic = res.conf$mu2
      sigma1_bic = res.conf$sigma1
      sigma2_bic = res.conf$sigma2
      pi2_bic = res.conf$pi2
      pi3_bic = res.conf$pi3
      post_bic = res.conf$post
      lambda = lambda.vec[i]
      print(lambda)
    }
  }
  
  
  #--------------------------------------#
  # Classification error estimated model #
  #--------------------------------------#
  class.assign_tmp <- classify_respondents(data, a.vec_bic, d.vec_bic, delta1_bic, delta2_bic, pi2_bic, pi3_bic, mu1_bic, mu2_bic, sigma1_bic, sigma2_bic)
  
  predicted_labels_class1 <- ifelse(class.assign_tmp[[1]] > class_thre, 1, 0)
  
  predicted_labels_class2 <- ifelse(class.assign_tmp[[2]] > class_thre, 1, 0)
  
  predicted_labels <- ifelse(predicted_labels_class1 == 1, 1, ifelse(predicted_labels_class2 == 1, 2, 0))
  
  # Calculate the number of incorrect predictions
  incorrect_predictions <- sum(xi.vec != predicted_labels)
  
  # Calculate the total number of predictions
  total_predictions <- length(xi.vec)
  
  # Calculate the classification error
  classification_error[j] <- incorrect_predictions / total_predictions
  
  
  
  delta1.bic_binary <- c()
  delta2.bic_binary <- c()
  
  delta1.bic_binary <- ifelse(delta1_bic != 0, 1, 0)
  delta2.bic_binary <- ifelse(delta2_bic != 0, 1, 0)
  
  # Calculate true positive rate for delta
  true_positives_delta1[j] <- sum((delta1.vec_binary == 1) & (delta1.bic_binary == 1)) / sum(delta1.vec_binary == 1)
  true_positives_delta2[j] <- sum((delta2.vec_binary == 1) & (delta2.bic_binary == 1)) / sum(delta2.vec_binary == 1)
  
  # Calculate true positive rate for delta
  false_positives_delta1[j] <- sum((delta1.vec_binary == 0) & (delta1.bic_binary == 1)) / sum(delta1.vec_binary == 0)
  false_positives_delta2[j] <- sum((delta2.vec_binary == 0) & (delta2.bic_binary == 1)) / sum(delta2.vec_binary == 0)
  
  
  # Bias
  bias_mean1[j] <- mu1_bic - mu1
  bias_mean2[j] <- mu2_bic - mu2
  bias_sd1[j] <- sigma1_bic - sigma1
  bias_sd2[j] <- sigma2_bic - sigma2
  bias_pi2[j] <- pi2_bic - pi2
  bias_pi3[j] <- pi3_bic - pi3
  bias_dif1[,j] <- delta1_bic - delta1.vec
  bias_dif2[,j] <- delta2_bic - delta2.vec
  bias_d.vec[,j] <- d.vec_bic - d.vec
  bias_a.vec[,j] <- a.vec_bic - a.vec
  
  # MSE
  mse_mean1[j] <- (mu1_bic - mu1)^2
  mse_mean2[j] <- (mu2_bic - mu2)^2
  mse_sd1[j] <- (sigma1_bic - sigma1)^2
  mse_sd2[j] <- (sigma2_bic - sigma2)^2
  mse_pi2[j] <- (pi2_bic - pi2)^2
  mse_pi3[j] <- (pi3_bic - pi3)^2
  mse_dif1[,j] <- (delta1_bic - delta1.vec)^2
  mse_dif2[,j] <- (delta2_bic - delta2.vec)^2
  mse_d.vec[,j] <- (d.vec_bic - d.vec)^2
  mse_a.vec[,j] <- (a.vec_bic - a.vec)^2
  
  
  
  
  # Calculate the ROC curve based on estimated parameters
  #---------------------------------------------------------------
  thresholds <- seq(0, 1, length.out = 100)
  fpr_class2_est <- numeric(100)
  tpr_class2_est <- numeric(100)
  
  xi.vec_class2 <- ifelse(xi.vec == 1, 1, 0)
  
  for (i in 1:100) {
    threshold <- thresholds[i]
    predicted <- ifelse(class.assign_tmp[[1]] > threshold, 1, 0)
    tp <- sum(predicted == 1 & xi.vec_class2 == 1)
    fp <- sum(predicted == 1 & xi.vec_class2 == 0)
    fn <- sum(predicted == 0 & xi.vec_class2 == 1)
    tn <- sum(predicted == 0 & xi.vec_class2 == 0)
    fpr_class2_est[i] <- fp / (fp + tn)
    tpr_class2_est[i] <- tp / (tp + fn)
  }
  
  roc_curves_est_class2[[j]] <- data.frame(fpr_class2_est, tpr_class2_est)
  
  
  thresholds <- seq(0, 1, length.out = 100)
  fpr_class3_est <- numeric(100)
  tpr_class3_est <- numeric(100)
  
  xi.vec_class3 <- ifelse(xi.vec == 2, 1, 0)
  
  for (i in 1:100) {
    threshold <- thresholds[i]
    predicted <- ifelse(class.assign_tmp[[2]] > threshold, 1, 0)
    tp <- sum(predicted == 1 & xi.vec_class3 == 1)
    fp <- sum(predicted == 1 & xi.vec_class3 == 0)
    fn <- sum(predicted == 0 & xi.vec_class3 == 1)
    tn <- sum(predicted == 0 & xi.vec_class3 == 0)
    fpr_class3_est[i] <- fp / (fp + tn)
    tpr_class3_est[i] <- tp / (tp + fn)
  }
  
  roc_curves_est_class3[[j]] <- data.frame(fpr_class3_est, tpr_class3_est)
  #---------------------------------------------------------------
  
  
  
  # Calculate the ROC curve based on true parameters
  #---------------------------------------------------------------
  thresholds <- seq(0, 1, length.out = 100)
  fpr_class2_true <- numeric(100)
  tpr_class2_true <- numeric(100)
  
  for (i in 1:100) {
    threshold <- thresholds[i]
    predicted <- ifelse(class.assign_tmp_true[[1]] > threshold, 1, 0)
    tp <- sum(predicted == 1 & xi.vec_class2 == 1)
    fp <- sum(predicted == 1 & xi.vec_class2 == 0)
    fn <- sum(predicted == 0 & xi.vec_class2 == 1)
    tn <- sum(predicted == 0 & xi.vec_class2 == 0)
    fpr_class2_true[i] <- fp / (fp + tn)
    tpr_class2_true[i] <- tp / (tp + fn)
  }
  
  roc_curves_true_class2[[j]] <- data.frame(fpr_class2_true, tpr_class2_true)
  
  
  thresholds <- seq(0, 1, length.out = 100)
  fpr_class3_true <- numeric(100)
  tpr_class3_true <- numeric(100)
  
  for (i in 1:100) {
    threshold <- thresholds[i]
    predicted <- ifelse(class.assign_tmp_true[[2]] > threshold, 1, 0)
    tp <- sum(predicted == 1 & xi.vec_class3 == 1)
    fp <- sum(predicted == 1 & xi.vec_class3 == 0)
    fn <- sum(predicted == 0 & xi.vec_class3 == 1)
    tn <- sum(predicted == 0 & xi.vec_class3 == 0)
    fpr_class3_true[i] <- fp / (fp + tn)
    tpr_class3_true[i] <- tp / (tp + fn)
  }
  
  roc_curves_true_class3[[j]] <- data.frame(fpr_class3_true, tpr_class3_true)
  #---------------------------------------------------------------
  
  
} # END OF LOOP

mean(classification_error)
mean(classification_error_true)

# Calculate average true positive rate for delta
mean(true_positives_delta1)
mean(true_positives_delta2)

# Calculate average false positive rate for delta
mean(false_positives_delta1)
mean(false_positives_delta2)



# AUC calculation

n_thresholds <- 100
average_roc <- sapply(1:n_thresholds, function(i) {
  mean_fpr <- mean(sapply(roc_curves_est_class2, function(x) x$fpr[i]))
  mean_tpr <- mean(sapply(roc_curves_est_class2, function(x) x$tpr[i]))
  c(fpr = mean_fpr, tpr = mean_tpr)
})

fpr <- average_roc[1,]
tpr <- average_roc[2,]

ave.roc <- data.frame(fpr = average_roc[1,], tpr = average_roc[2,])



#library(dplyr)
sorted_roc <- ave.roc %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc_class2 <- sum((sorted_roc$fpr[-1] - sorted_roc$fpr[-nrow(sorted_roc)]) * (sorted_roc$tpr[-1] + sorted_roc$tpr[-nrow(sorted_roc)]) / 2)
auc_class2

auc_class2 <- round(auc_class2, 3)

#library(ggplot2)
ROC_class2_est <- ggplot(ave.roc, aes(x = fpr, y = tpr)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", auc_class2)) + 
  ggtitle("ROC for class 2 based on estimated params.") + 
  theme_classic()




n_thresholds <- 100
average_roc <- sapply(1:n_thresholds, function(i) {
  mean_fpr <- mean(sapply(roc_curves_est_class3, function(x) x$fpr[i]))
  mean_tpr <- mean(sapply(roc_curves_est_class3, function(x) x$tpr[i]))
  c(fpr = mean_fpr, tpr = mean_tpr)
})

fpr <- average_roc[1,]
tpr <- average_roc[2,]

ave.roc <- data.frame(fpr = average_roc[1,], tpr = average_roc[2,])



#library(dplyr)
sorted_roc <- ave.roc %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc_class3 <- sum((sorted_roc$fpr[-1] - sorted_roc$fpr[-nrow(sorted_roc)]) * (sorted_roc$tpr[-1] + sorted_roc$tpr[-nrow(sorted_roc)]) / 2)
auc_class3

auc_class3 <- round(auc_class3, 3)

#library(ggplot2)
ROC_class3_est <- ggplot(ave.roc, aes(x = fpr, y = tpr)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", auc_class3)) +
  ggtitle("ROC for class 3 based on estimated params.") + 
  theme_classic()




n_thresholds <- 100
average_roc <- sapply(1:n_thresholds, function(i) {
  mean_fpr <- mean(sapply(roc_curves_true_class2, function(x) x$fpr[i]))
  mean_tpr <- mean(sapply(roc_curves_true_class2, function(x) x$tpr[i]))
  c(fpr = mean_fpr, tpr = mean_tpr)
})

fpr <- average_roc[1,]
tpr <- average_roc[2,]

ave.roc <- data.frame(fpr = average_roc[1,], tpr = average_roc[2,])



#library(dplyr)
sorted_roc <- ave.roc %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc_class2_true <- sum((sorted_roc$fpr[-1] - sorted_roc$fpr[-nrow(sorted_roc)]) * (sorted_roc$tpr[-1] + sorted_roc$tpr[-nrow(sorted_roc)]) / 2)
auc_class2_true

auc_class2_true <- round(auc_class2_true, 3)

#library(ggplot2)
ROC_class2_true <- ggplot(ave.roc, aes(x = fpr, y = tpr)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", auc_class2_true)) +
  ggtitle("ROC for class 2 based on true params.") + 
  theme_classic()



n_thresholds <- 100
average_roc <- sapply(1:n_thresholds, function(i) {
  mean_fpr <- mean(sapply(roc_curves_true_class3, function(x) x$fpr[i]))
  mean_tpr <- mean(sapply(roc_curves_true_class3, function(x) x$tpr[i]))
  c(fpr = mean_fpr, tpr = mean_tpr)
})

fpr <- average_roc[1,]
tpr <- average_roc[2,]

ave.roc <- data.frame(fpr = average_roc[1,], tpr = average_roc[2,])




#library(dplyr)
sorted_roc <- ave.roc %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc_class3_true <- sum((sorted_roc$fpr[-1] - sorted_roc$fpr[-nrow(sorted_roc)]) * (sorted_roc$tpr[-1] + sorted_roc$tpr[-nrow(sorted_roc)]) / 2)
auc_class3_true

auc_class3_true <- round(auc_class3_true, 3)


#library(ggplot2)
ROC_class3_true <- ggplot(ave.roc, aes(x = fpr, y = tpr)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", auc_class3_true)) +
  # geom_text(x = 0.8, y = 0.2, label = paste0("AUC=", auc_class3_true)) +
  ggtitle("ROC for class 3 based on true params.") + 
  theme_classic()



auc_class2
auc_class3

auc_class2_true
auc_class3_true



# Save bias results
save(bias_mean1, bias_mean2, 
     bias_sd1, bias_sd2, 
     bias_pi2, bias_pi3, 
     bias_dif1, bias_dif2, 
     bias_d.vec, bias_a.vec, 
     mse_mean1, mse_mean2, 
     mse_sd1, mse_sd2, 
     mse_pi2, mse_pi3, 
     mse_dif1, mse_dif2, 
     mse_d.vec, mse_a.vec,
     file = "bias_ROC_N5000_J50_pi0.5_0.3_0.2.rda")


# Bias
round(mean(rowMeans(abs(bias_d.vec))),3)
round(mean(rowMeans(abs(bias_a.vec))),3)
round(mean(rowMeans(abs(bias_dif1))),3)
round(mean(rowMeans(abs(bias_dif2))),3)

round(mean(abs(bias_mean1)),3)
round(mean(abs(bias_mean2)),3)
round(mean(abs(bias_sd1)),3)
round(mean(abs(bias_sd2)),3)
round(mean(abs(bias_pi2)),3)
round(mean(abs(bias_pi3)),3)

# RMSE
round((1/J)*sum(sqrt((1/sim)*rowSums(mse_d.vec))),3)
round((1/J)*sum(sqrt((1/sim)*rowSums(mse_a.vec))),3)
round((1/J)*sum(sqrt((1/sim)*rowSums(mse_dif1))),3)
round((1/J)*sum(sqrt((1/sim)*rowSums(mse_dif2))),3)

round(sqrt((1/sim)*sum(mse_mean1)),3)
round(sqrt((1/sim)*sum(mse_mean2)),3)
round(sqrt((1/sim)*sum(mse_sd1)),3)
round(sqrt((1/sim)*sum(mse_sd2)),3)
round(sqrt((1/sim)*sum(mse_pi2)),3)
round(sqrt((1/sim)*sum(mse_pi3)),3)

