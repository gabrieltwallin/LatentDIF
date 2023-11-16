
rm(list=ls())
graphics.off()


#----------------------------------------------------#
# Latent DIF detection for the 2-PL model           #
# with two latent classes (baseline + outlier group) #
#----------------------------------------------------#

N <- 1000
J <- 50
p <- 20 

pi <- 0.5

# Threshhold for delta values to be classified as non-zero.
#delta_thresh <- 0.5 # used to use 0.5

step <- 1
lambda <- 0.1

# Respondent class threshold:
class_thre <- 0.5 # used to be 0.2

# Number of classes: 
K <- 2

# Number of replications
sim <- 100

bias_mean <- c()
bias_sd <- c()
bias_pi <- c()
bias_dif <- matrix(0, nrow = J, ncol = sim)
bias_d.vec <- matrix(0, nrow = J, ncol = sim)
bias_a.vec <- matrix(0, nrow = J, ncol = sim)

rmse_mean <- c()
rmse_sd <- c()
rmse_pi <- c()
rmse_dif <- matrix(0, nrow = J, ncol = sim)
rmse_d.vec <- matrix(0, nrow = J, ncol = sim)
rmse_a.vec <- matrix(0, nrow = J, ncol = sim)

est_mean <- c()
est_sd <- c()
est_pi<- c()
est_d.vec <- matrix(0, nrow = J, ncol = sim)
est_a.vec <- matrix(0, nrow = J, ncol = sim)

true_positives_xi <- rep(0, sim)
true_negatives_xi <- rep(0, sim)

false_positives_xi <- rep(0, sim)
false_positives_xi.true <- rep(0, sim)
false_negatives_xi <- rep(0, sim)

true_positives_xi.true <- rep(0, sim)
true_negatives_xi.true <- rep(0, sim)

true_positives_delta <- rep(0, sim)
false_positives_delta <- rep(0, sim)

true_negatives_delta <- rep(0, sim)

classification_error <- rep(0, sim)
classification_error_true <- rep(0, sim)  



# Initialize list to store ROC curves for each iteration
roc_curves_est <- list()

roc_curves_true <- list()

roc_curves_alt <- list()

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

a.vec <- runif(J, 0.5, 1.5)
d.vec <- runif(J, -2, 2)

delta.vec <- c(rep(0, J-p), runif(p, 0.5, 1.5)) # usually (0.5, 1.5)
delta.vec_binary <- rep(0, J)
delta.vec_binary[delta.vec > 0] <- 1

xi.vec <- rbinom(N, 1, prob = pi)

theta.vec <- rnorm(N)
mu = 0.5
sigma = 1.5
theta.vec[xi.vec==1] = rnorm(sum(xi.vec==1), mu, sigma)



temp = theta.vec %*% t(a.vec) + rep(1, N) %*% t(d.vec) + xi.vec %*% t(delta.vec)

prob = 1/(1+exp(-temp))

for(j in 1:sim){
  
  
  data <- matrix(0, N, J)
  data[] <- rbinom(N*J, 1, prob)
  
  grid <- seq(-4, 4, length.out = 31)
  weight = dnorm(grid)
  weight = weight/sum(weight)
  
  
  
  llik <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma){
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = exp(temp1) %*% weight
    
    grid.1 = grid * sigma + mu
    temp = a.vec %*% t(grid.1) + (d.vec+delta.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = exp(temp2) %*% weight
    
    sum(log(temp1 * (1-pi) + temp2 * pi))
  }
  
  
  post.matr <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma){
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = t(t(exp(temp1)) * weight) * (1-pi)
    
    grid.1 = grid * sigma + mu
    temp = a.vec %*% t(grid.1) + (d.vec+delta.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = t(t(exp(temp2)) * weight) * pi
    
    post <- cbind(temp1, temp2)
    
    post = post/rowSums(post)
    post
  }
  
  classify_respondents <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma) {
    # Calculate probability of belonging to first class
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp1 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp1 = t(t(exp(temp1)) * weight) * (1-pi)
    
    # Calculate probability of belonging to second class
    grid.1 = grid * sigma + mu
    temp = a.vec %*% t(grid.1) + (d.vec+delta.vec) %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    temp2 = data %*% log(prob) + (1-data) %*% log(1-prob)
    temp2 = t(t(exp(temp2)) * weight) * pi
    
    # Combine probabilities
    post <- cbind(temp1, temp2)
    
    # Normalize probabilities
    post = post/rowSums(post)
    
    # Assign class based on the highest probability
    #class_assignments <- ifelse(rowSums(post[,1:31]) > rowSums(post[,32:62]), 0, 1)
    #return(class_assignments)
    return(rowSums(post[,32:62]))
  }
  
  class.assign_tmp_true <- classify_respondents(data, a.vec, d.vec, delta.vec, pi, mu, sigma)
  
  class.assign_TruParVal <- ifelse(class.assign_tmp_true > class_thre, 1, 0)
  
  
  # Calculate true positive rate for xi for the true model
  true_positives_xi.true[j] <- sum((xi.vec == 1) & (class.assign_TruParVal == 1)) / sum(xi.vec == 1)
  
  # Calculate false positive rate for xi for the true model
  false_positives_xi.true[j] <- sum((xi.vec == 0) & (class.assign_TruParVal == 1)) / sum(xi.vec == 0)
  
  # Calculate true negative rate for xi for the true model
  true_negatives_xi.true[j] <- sum((xi.vec == 0) & (class.assign_TruParVal == 0)) / sum(xi.vec == 0)
  
  
  
  # Calculate the number of incorrect predictions
  incorrect_predictions_true <- sum(xi.vec != class.assign_TruParVal)
  
  # Calculate the total number of predictions
  total_predictions <- length(xi.vec)
  
  # Calculate the classification error
  classification_error_true[j] <- incorrect_predictions_true / total_predictions
  
  
  
  Update.pi <- function(post){
    
    temp1 = colMeans(post[,1:31])
    temp2 = colMeans(post[,32:62])
    
    sum(temp2)
  }
  
  soft_thre <- function(x, lambda){
    temp = x
    temp[abs(x) <= lambda] = 0
    temp[x > lambda] = temp[x > lambda] - lambda
    temp[x < -lambda] = temp[x < -lambda] + lambda
    temp
  }
  
  prox.grad <- function(x, grad, lambda, step){
    soft_thre(x - step * grad, lambda * step)
  }
  
  
  
  line.search <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma, grad.a, 
                          grad.d, grad.delta, grad.sigma, grad.mu, lambda, step=1) {
    # Set the initial step size to the input step size
    step.size <- step
    
    
    # Set the initial value of the objective function to the current value
    obj.val <- -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) + lambda * sum(abs(delta.vec))
    
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
      delta.new <- prox.grad(delta.vec, grad.delta, lambda, step.size)
      sigma.new <- sigma - step.size * grad.sigma 
      mu.new <- mu - step.size * grad.mu 
      
      # Calculate the new value of the objective function
      new.obj.val <- -llik(data, a.new, d.new, delta.new, pi, mu.new, sigma.new) + lambda * sum(abs(delta.new))
      
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
  
  
  line.search.conf <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma, grad.a, grad.d, grad.delta, grad.sigma, grad.mu, step=1) {
    
    # Set the initial step size to the input step size
    step.size <- step
    
    # Set the initial value of the objective function to the current value
    obj.val <- -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
    
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
      delta.new <- delta.vec - step.size * grad.delta
      sigma.new <- sigma - step.size * grad.sigma 
      mu.new <- mu - step.size * grad.mu 
      
      # Calculate the new value of the objective function
      new.obj.val <- -llik(data, a.new, d.new, delta.new, pi, mu.new, sigma.new)
      
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
  
  
  
  Update.par <- function(post, data, a.vec, d.vec, delta.vec,
                         mu, sigma, N, grid, pi, lambda, step){
    
    temp1 = post[,1:31]
    temp2 = post[,32:62]
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    grid.1 = grid * sigma + mu
    temp = a.vec %*% t(grid.1) + (d.vec+delta.vec) %*% t(rep(1, 31)) 
    prob2 = 1/(1+exp(-temp))
    
    grad.a = t(data) %*% temp1 %*% grid - prob %*% (colSums(temp1) * grid) +  t(data) %*% temp2 %*% grid.1 - prob2 %*% (colSums(temp2) * grid.1)
    grad.a = -(grad.a/N)
    
    grad.d = rowSums(t(data) %*% temp1) - prob %*% (colSums(temp1) ) +  rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.d = -(grad.d/N)
    
    grad.delta = rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.delta = -(grad.delta/N)
    
    grad.sigma = t(a.vec)%*% (t(data) %*% temp2 %*% grid) - t(a.vec)%*% prob2 %*% (colSums(temp2)*grid) 
    grad.sigma = -grad.sigma/N
    
    grad.mu = sum(t(a.vec)%*% (t(data) %*% temp2)) - t(a.vec)%*% prob2 %*% (colSums(temp2)) 
    grad.mu = -grad.mu/N
    
    
    # Use the line search function to determine the step size
    step <- line.search(data, a.vec, d.vec, delta.vec, pi, mu, sigma, grad.a, grad.d, 
                        grad.delta, grad.sigma, grad.mu, lambda, step)
    
    a.new = a.vec - step * grad.a
    d.new = d.vec - step*grad.d
    delta.new = prox.grad(delta.vec, grad.delta, lambda, step)
    #if (!is.null(delta_lambda)){
    #  delta.new[delta_lambda==0]=0}
    
    sigma.new = sigma - step*grad.sigma
    mu.new = mu - step*grad.mu
    
    list(a.vec = a.new, d.vec = d.new, delta.vec = delta.new, mu = mu.new[1], sigma = sigma.new[1])  
  }
  
  
  
  Update.par_conf <- function(post, data, a.vec, d.vec, delta.vec, 
                              mu, sigma, step){
    
    # Set delta for the first number of items to 0
    delta.vec <- ifelse(delta.vec != 0, delta.vec, 0)
    
    #delta.vec <- c(delta.vec[1:p], rep(0, J-p)) 
    
    
    # Posterior probabilities, baseline group  
    temp1 = post[,1:31]
    
    # Posterior probabilities, outlier group
    temp2 = post[,32:62]
    
    temp = a.vec %*% t(grid) + d.vec %*% t(rep(1, 31)) 
    prob = 1/(1+exp(-temp))
    grid.1 = grid * sigma + mu
    temp = a.vec %*% t(grid.1) + (d.vec+delta.vec) %*% t(rep(1, 31)) 
    prob2 = 1/(1+exp(-temp))
    
    # First derivative of Q w.r.t. slope:
    grad.a = t(data) %*% temp1 %*% grid - prob %*% (colSums(temp1) * grid) +  t(data) %*% temp2 %*% grid.1 - prob2 %*% (colSums(temp2) * grid.1)
    grad.a = -(grad.a/N)
    
    # First derivative of Q w.r.t. intercept:
    grad.d = rowSums(t(data) %*% temp1) - prob %*% (colSums(temp1) ) +  rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.d = -(grad.d/N)
    
    # First derivative of Q w.r.t. to DIF effect:
    grad.delta = rowSums(t(data) %*% temp2) - prob2 %*% (colSums(temp2))
    grad.delta = -(grad.delta/N)
    
    
    # First derivative of Q w.r.t. standard deviation:
    grad.sigma = t(a.vec) %*% (t(data) %*% temp2 %*% grid) - t(a.vec)%*% prob2 %*% (colSums(temp2)*grid) 
    grad.sigma = -grad.sigma/N
    
    
    # Gradient of Q w.r.t. mean:
    grad.mu = sum(t(a.vec) %*% (t(data) %*% temp2)) - t(a.vec)%*% prob2 %*% (colSums(temp2)) 
    grad.mu = -grad.mu/N
    
    
    # Use the line search function to determine the step size
    step <- line.search.conf(data, a.vec, d.vec, delta.vec, pi, mu, sigma, grad.a, grad.d, 
                             grad.delta, grad.sigma, grad.mu, step)
    
    # Gradient descent to update slope parameter:
    a.new = a.vec - step * grad.a
    
    
    # Gradient descent to update intercept:
    d.new = d.vec - step * grad.d
    
    
    # Proximal gradient descent to update DIF effect:
    # Gradient descent to update DIF effect:
    delta.new <- ifelse(delta.vec != 0, delta.vec - step * grad.delta, 0)
    #delta.new[J-p:J] <- 0
    
    #for(q in 1:p){
    #  if(delta.new[q] < 0){  # Sign constraint on the DIF items 
    #    delta.new[q] = 0
    #  }
    #}
    
    
    # Gradient descent to update standard deviation:
    sigma.new = sigma - step*grad.sigma
    
    
    # Gradient descent to update mean:
    mu.new = mu - step*grad.mu 
    
    
    list(a.vec = a.new, d.vec = d.new, delta.vec = delta.new, mu = mu.new[1], sigma = sigma.new[1])  
  }
  
  
  EM <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma, lambda, step){
    
    obj0 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) + lambda * sum(abs(delta.vec))
    
    post = post.matr(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
    
    pi = Update.pi(post)
    
    for(i in 1:10){
      res = Update.par(post, data, a.vec, d.vec, delta.vec, mu, sigma, N, grid, pi,
                       lambda, step)
      
      a.vec = res$a.vec
      d.vec = res$d.vec
      delta.vec = res$delta.vec
      mu = res$mu
      sigma = res$sigma
    }
    
    obj1 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) + lambda * sum(abs(delta.vec))
    
    while(obj0-obj1 > 1e-4){
      obj0=obj1
      post = post.matr(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
      
      pi = Update.pi(post)
      
      for(i in 1:10){
        res = Update.par(post, data, a.vec, d.vec, delta.vec, mu, sigma, N, grid, pi,
                         lambda, step)
        
        a.vec = res$a.vec
        d.vec = res$d.vec
        delta.vec = res$delta.vec
        mu = res$mu
        sigma = res$sigma
      }
      
      obj1 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) + lambda * sum(abs(delta.vec))
      print(obj1)
    }
    
    list(a.vec = a.vec, d.vec = d.vec, delta.vec = delta.vec, 
         mu = mu, sigma = sigma, pi = pi, post = post)
  }
  
  EM_conf <- function(data, a.vec, d.vec, delta.vec, pi, mu, sigma, step){ # I had 'step' as an argument before
    
    # Objective function:
    obj0 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
    
    # Posterior probabilities class membership:
    post = post.matr(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
    
    # Update class proportion
    pi = Update.pi(post)
    
    
    # Run gradient descent 
    for(i in 1:10){ 
      res = Update.par_conf(post, data, a.vec, d.vec, delta.vec, 
                            mu, sigma, step=1) # , N, niter = 5, stop_ep = 10^(-4), t0 = 1  
      a.vec = res$a.vec
      d.vec = res$d.vec
      delta.vec = res$delta.vec
      mu = res$mu
      sigma = res$sigma
    }
    
    # Update value of objective function:
    obj1 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
    
    while(obj0-obj1 > 1e-4){
      
      obj0=obj1
      
      post = post.matr(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
      
      pi = Update.pi(post)
      
      # Run gradient descent 
      for(i in 1:10){ 
        res = Update.par_conf(post, data, a.vec, d.vec, delta.vec, 
                              mu, sigma, step=1) # , N, niter = 5, stop_ep = 10^(-4), t0 = 1 
        a.vec = res$a.vec
        d.vec = res$d.vec
        delta.vec = res$delta.vec
        mu = res$mu
        sigma = res$sigma
      }
      
      obj1 = -llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) 
      print(obj1)
    }
    
    bic=-2*llik(data, a.vec, d.vec, delta.vec, pi, mu, sigma) + log(N)*(J+J+3 + sum(delta.vec!=0))
    
    list(a.vec = a.vec, d.vec = d.vec, delta.vec = delta.vec, mu = mu, sigma = sigma, pi = pi, post = post, bic = bic)
  }
  
  
  delta.vec0 = delta.vec+rnorm(J, 0, 0.1) 
  
  a.vec0 = a.vec+rnorm(J, 0, 0.1)
  
  d.vec0 <- d.vec+rnorm(J, 0, 0.1)
  
  mu0 <- mu + rnorm(1, 0, 0.1)
  
  sigma0 <- sigma + rnorm(1, 0, 0.1)
  
  lambda.vec <- seq(0.01, 0.0005, by = -0.001)
  
  res = EM(data, a.vec0, d.vec0, delta.vec0, pi, mu0, sigma0, lambda = lambda.vec[1], step)
  
  bic = Inf
  print('bic')
  
  c = 0.5
  
  for(i in 2:length(lambda.vec)){
    res <- EM(data, res$a.vec, res$d.vec, res$delta.vec, res$pi, res$mu, res$sigma, lambda = lambda.vec[i], step)
    delta_lambda0 <- res$delta.vec
    delta_lambda0[delta_lambda0 < c] = 0
    res.conf <- EM_conf(data, res$a.vec, res$d.vec, delta_lambda0, res$pi, res$mu, res$sigma, step)
    bic_new = res.conf$bic
    
    if(bic_new < bic){
      bic = bic_new
      delta_bic = res.conf$delta.vec
      #print(delta_bic)
      a.vec_bic = res.conf$a.vec
      d.vec_bic = res.conf$d.vec
      mu_bic = res.conf$mu
      sigma_bic = res.conf$sigma
      pi_bic = res.conf$pi
      post_bic = res.conf$post
      lambda = lambda.vec[i]
      print(lambda)
    }
  }
  
  
  
  # Binary vector indicating the estimated latent class membership of each individual
  
  class.assign_tmp <- classify_respondents(data, a.vec_bic, d.vec_bic, delta_bic, pi_bic, mu_bic, sigma_bic)
  
  predicted_labels <- ifelse(class.assign_tmp > class_thre, 1, 0)
  
  # Calculate the number of incorrect predictions
  incorrect_predictions <- sum(xi.vec != predicted_labels)
  
  # Calculate the total number of predictions
  total_predictions <- length(xi.vec)
  
  # Calculate the classification error
  classification_error[j] <- incorrect_predictions / total_predictions
  
  
  
  delta.bic_binary <- c()
  delta.bic_binary <- ifelse(delta_bic != 0, 1, 0)
  
  
  # Calculate true positive rate for delta
  true_positives_delta[j] <- sum((delta.vec_binary == 1) & (delta.bic_binary == 1)) / sum(delta.vec_binary == 1)
  
  # Calculate true negative rate for delta
  true_negatives_delta[j] <- sum((delta.vec_binary == 0) & (delta.bic_binary == 0)) / sum(delta.vec_binary == 0)
  
  false_positives_delta[j] <- sum((delta.vec_binary == 0) & (delta.bic_binary == 1)) / sum(delta.vec_binary == 0)
  
  
  
  #----------------------------#  
  # MSE and bias of estimators #
  #----------------------------#  
  
  # Bias
  bias_mean[j] <- mu_bic - mu
  bias_sd[j] <- sigma_bic - sigma
  bias_pi[j] <- pi_bic - pi
  bias_dif[,j] <- delta_bic - delta.vec
  bias_d.vec[,j] <- d.vec_bic - d.vec
  bias_a.vec[,j] <- a.vec_bic - a.vec
  
  
  # MSE
  rmse_mean[j] <- (mu_bic - mu)^2
  rmse_sd[j] <- (sigma_bic - sigma)^2
  rmse_pi[j] <- (pi_bic - pi)^2
  rmse_dif[,j] <- (delta_bic - delta.vec)^2
  rmse_d.vec[,j] <- (d.vec_bic - d.vec)^2
  rmse_a.vec[,j] <- (a.vec_bic - a.vec)^2
  
  
  
  # Calculate the ROC curve based on estimated parameters
  
  ########################################
  
  thresholds <- seq(0, 1, length.out = 100)
  
  roc_curves_est[[j]] <- calculate_roc_curve(class.assign_tmp, xi.vec, thresholds)
  roc_curves_true[[j]] <- calculate_roc_curve(class.assign_tmp_true, xi.vec, thresholds)
  
  
} # END OF LOOP

mean(classification_error)
mean(classification_error_true)


# Calculate average true positive rate for delta
mean(true_positives_delta)

# Calculate average false positive rate for delta
mean(false_positives_delta)

# Calculate average true negative rate for delta
mean(true_negatives_delta)



#library(ggplot2)
#library(dplyr)
# Calculate the average ROC curve based on the estimated model:
n_thresholds <- 100
average_roc <- sapply(1:n_thresholds, function(i) {
  mean_fpr <- mean(sapply(roc_curves_est, function(x) x$fpr[i]))
  mean_tpr <- mean(sapply(roc_curves_est, function(x) x$tpr[i]))
  c(fpr = mean_fpr, tpr = mean_tpr)
})

fpr <- average_roc[1,]
tpr <- average_roc[2,]

ave.roc <- data.frame(fpr = average_roc[1,], tpr = average_roc[2,])

#library(dplyr)
sorted_roc <- ave.roc %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc <- sum((sorted_roc$fpr[-1] - sorted_roc$fpr[-nrow(sorted_roc)]) * (sorted_roc$tpr[-1] + sorted_roc$tpr[-nrow(sorted_roc)]) / 2)
auc


#library(ggplot2)
ROC_est <- ggplot(ave.roc, aes(x = fpr, y = tpr)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", round(auc, 3))) + 
  ggtitle("ROC based on estimated params.") + 
  theme_classic()







# Calculate the average ROC curve based on the true model:
n_thresholds_true <- 100
average_roc_true <- sapply(1:n_thresholds_true, function(i) {
  mean_fpr_true <- mean(sapply(roc_curves_true, function(x) x$fpr[i]))
  mean_tpr_true <- mean(sapply(roc_curves_true, function(x) x$tpr[i]))
  c(fpr = mean_fpr_true, tpr = mean_tpr_true)
})

fpr_true <- average_roc_true[1,]
tpr_true <- average_roc_true[2,]

ave.roc_true <- data.frame(fpr = average_roc_true[1,], tpr = average_roc_true[2,])


#library(dplyr)
sorted_roc_true <- ave.roc_true %>% arrange(fpr)

# calculate the area under the curve by adding up the areas of the trapezoids
auc_true <- sum((sorted_roc_true$fpr[-1] - 
                   sorted_roc_true$fpr[-nrow(sorted_roc_true)]) * (sorted_roc_true$tpr[-1] + 
                                                                     sorted_roc_true$tpr[-nrow(sorted_roc_true)]) / 2)
auc_true

#library(ggplot2)
ROC_true <- ggplot(ave.roc_true, aes(x = fpr_true, y = tpr_true)) + 
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") + 
  annotate("text", x=0.8, y=0.2, label=paste0("AUC=", round(auc_true, 3))) +
  ggtitle("ROC based on true params.") + 
  theme_classic()





auc
auc_true




pdf(file = "C:/Users/gawa0005/Dropbox/Psychometrics/VR project/Cheating detection/Code/LatentDIF/Simulation Results/Plots/My Plot.pdf",
    width = 9, # The width of the plot in inches
    height = 5) # The height of the plot in inches

multiplot(ROC_est,
          ROC_true,
          cols = 2)

dev.off()



# Save bias results
save(bias_mean, bias_sd, bias_pi, bias_dif, bias_d.vec, bias_a.vec, 
     ROC_est, ROC_true, file = "results_N5000_J50_pi05.rda")



library(tidyr)
library(gridExtra)


# MSE
MSE <- function(results, parameter){
  rel.mse <- matrix(0, nrow=nrow(results), ncol=ncol(results))
  for (i in 1:nrow(results)){
    rel.mse[i,] <- (results[i,]-parameter[i])^2
  }
  res <- apply(rel.mse, 1, mean)
  return(res) 
}

estimated.d.vec <- bias_d.vec + d.vec
estimated.a.vec <- bias_a.vec + a.vec
estimated.delta.vec <- bias_dif + delta.vec
estimated.mu <- bias_mean + mu
estimated.sd <- bias_sd + sigma

mse_d.vec <- MSE(estimated.d.vec, d.vec)
mse_a.vec <- MSE(estimated.a.vec, a.vec)
mse_delta.vec <- MSE(estimated.delta.vec, delta.vec)

mse_mu <- mean((estimated.mu-mu)^2)
mse_sigma <- mean((estimated.sd-sigma)^2)



#library(ggplot2)

# Combine the vectors into a data frame
df <- data.frame(
  mse_d = mse_d.vec,
  mse_a = mse_a.vec,
  mse_delta = mse_delta.vec,
  index = 1:J
)

# Convert the data frame to long format for plotting
df_long <- tidyr::gather(df, key = "mse_type", value = "mse_value", -index)

# Set the labels for the mse types
mse_labels <- c("Intercept" = "mse_a", "Slope" = "mse_d", "DIF effect" = "mse_delta")

# Reverse the mse_labels vector to map the mse_types to their corresponding labels
label_mapping <- names(mse_labels)

# Define colors and shapes for each mse type
color_palette <- c("mse_a" = "#F8766D", "mse_d" = "#00BFC4", "mse_delta" = "#C77CFF")
shape_palette <- c("mse_a" = 21, "mse_d" = 22, "mse_delta" = 23)

# Plot the data using ggplot2
ggplot(df_long, aes(x = index, y = mse_value, fill = mse_type, shape = mse_type)) +
  geom_point(size = 6, alpha = 0.8, stroke = 0.2) + 
  xlab("Item") +
  ylab("Mean Squared Error") +
  scale_fill_manual(name = "Estimate", values = color_palette, labels = label_mapping) +
  scale_shape_manual(name = "Estimate", values = shape_palette, labels = label_mapping) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    axis.text.x = element_text(size = 12),    # Adjust x-axis text size
    axis.text.y = element_text(size = 12),    # Adjust y-axis text size
    legend.title = element_text(size = 14)    # Adjust legend title size
  )

####################


get_common_ylims <- function(plot1_data, plot2_data) {
  max_y <- max(c(max(plot1_data$mse_value, na.rm = TRUE), max(plot2_data$mse_value, na.rm = TRUE)))
  min_y <- min(c(min(plot1_data$mse_value, na.rm = TRUE), min(plot2_data$mse_value, na.rm = TRUE)))
  
  return(c(min_y, max_y))
}

# Assuming you have different data frames for each new plot
# ylims1 <- get_common_ylims(df_long, new_df1) 
# ylims2 <- get_common_ylims(new_df2, new_df3) 
# ylims3 <- get_common_ylims(new_df4, new_df5) 


main_plot <- ggplot(df_long, aes(x = index, y = mse_value, fill = mse_type, shape = mse_type)) +
  geom_point(size = 6, alpha = 0.8, stroke = 0.2) + 
  xlab("Item") +
  ylab("Mean Squared Error") +
  scale_fill_manual(name = "Estimate", values = color_palette, labels = label_mapping) +
  scale_shape_manual(name = "Estimate", values = shape_palette, labels = label_mapping) +
  theme_classic() +
  coord_cartesian(ylim = ylims1) +
  theme(
    axis.title.x = element_text(size = 14),   
    axis.title.y = element_text(size = 14),   
    axis.text.x = element_text(size = 12),    
    axis.text.y = element_text(size = 12),    
    legend.title = element_text(size = 14),
    legend.position = "none"  # Hide the legend
  )



grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  legends <- lapply(plots, cowplot::get_legend)
  main_plot <- plots[[1]]
  
  combined_plots <- do.call(gridExtra::arrangeGrob, c(plots, ncol=2))
  legend_grob = cowplot::get_legend(main_plot + theme(legend.position="bottom"))
  
  gridExtra::grid.arrange(arrangeGrob(combined_plots, bottom=legend_grob))
}

grid_arrange_shared_legend(main_plot, main_plot, main_plot, 
                           main_plot, main_plot, main_plot)

install.packages("gridExtra")
install.packages("cowplot")

####################


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



###################################################

# Plotting function with conditional axes
plot_custom <- function(df_long, show_x = TRUE, show_y = TRUE) {
  p <- ggplot(df_long, aes(x = index, y = mse_value, fill = mse_type, shape = mse_type)) +
    geom_point(size = 3, alpha = 0.8, stroke = 0.2) +
    scale_fill_manual(name = "Estimate", values = color_palette, labels = label_mapping) +
    scale_shape_manual(name = "Estimate", values = shape_palette, labels = label_mapping) +
    theme_classic() +
    theme(legend.position = "none")
  
  if (!show_x) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (!show_y) {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  return(p)
}

# Create your plots
p1 <- plot_custom(df_long, show_x = FALSE)
p2 <- plot_custom(df_long, show_x = FALSE, show_y = FALSE)
p3 <- plot_custom(df_long, show_x = FALSE)
p4 <- plot_custom(df_long, show_x = FALSE, show_y = FALSE)
p5 <- plot_custom(df_long)
p6 <- plot_custom(df_long, show_y = FALSE)

# Arrange your plots
layout <- (p1 | p2) + (p3 | p4) + (p5 | p6)

# Display with a shared legend
(layout + plot_layout(guides = "collect")) & theme(legend.position = "bottom")

###################################################
