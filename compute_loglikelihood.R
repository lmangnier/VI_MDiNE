library(CholWishart)
library(MCMCpack)

.sigma_epsilon = function(epsilon){
  #Epsilon is a n*J matrix
  return(epsilon/(rowSums(epsilon)+1))
}
.lambda_epsilon = function(epsilon) 0.5*epsilon*(.sigma_epsilon(epsilon)-0.5)


.lb_ph = function(Y, mu_W, sigma_W, epsilon){
  #Function to compute the log-likelihood of the bounded likelihood of the data
  return(sum(colSums(log(.sigma_epsilon(epsilon)) + mu_W*Y - 0.5*mu_W - epsilon*0.5 -
                       .lambda_epsilon(epsilon)*(mu_W^2+sigma_W-epsilon^2))))
  
  
}

.lb_plambda = function(alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0){
  #Alpha, omega, mu and sigma are K*J matrices corresponding to the number of features and species
  #r0 and delta0 are the hyperparameters
  E_q_lambda = function(fn, a,b) {
    integrate(function(lambda) {
      dgamma(lambda, a, b) * fn(lambda)
    }, 0, Inf)$value
  }
  
  return(sum(sapply(1:ncol(omega_lambda), function(x) {
    sapply(1:nrow(omega_lambda), function(y) {
      0.5*(E_q_lambda(function(lambda) log(lambda), alpha_lambda[y,x],omega_lambda[y,x]) -
             (alpha_lambda[y,x]/omega_lambda[y,x])*(mu_beta[y,x]^2 + sigma_beta[y,x])) + 
        (r0-1)*E_q_lambda(function(lambda) log(lambda),alpha_lambda[y,x],omega_lambda[y,x])*(log(delta0)- lgamma(r0)) - delta0*(alpha_lambda[y,x]/omega_lambda[y,x])
    })})))
  
}

.lb_pa = function(alpha_a, beta_a, A){
  
  #alpha and beta are vector of size J corresponding to the number of species
  #A is an scalar hyper parameter
  E_q_a = function(fn, r,t) {
    integrate(function(a) {
      dinvgamma(a, r, t) * fn(a)
    }, 0, Inf)$value
  }
  
  return(sum(sapply(1:length(alpha_a), function(x){
    -0.5*log(A) - lgamma(0.5)-(1.5)*E_q_a(function(a) log(a), alpha_a[x], beta_a[x]) - (1/A)*(alpha_a[x]/beta_a[x])
  })))
  
}

.lb_pw = function(X,z,mu_W, sigma_W, mu_beta, sigma_beta, b_Sigma0, b_Sigma_1, v0, v1,N,M,J, K){


  pw1 = sum(sapply(1:nrow(X), function(i){
    sum(sapply(1:nrow(mu_beta), function(j){
      0.5*((K+1)*(mu_W[i,j]^2 + sigma_W[i,j]))}))}))


  pw2 = sum(sapply(1:nrow(X), function(i){
    sum(sapply(1:nrow(mu_beta), function(j){
      sum(sapply(1:ncol(X), function(k){
        X[i,k]*(mu_beta[j,k]^2 + sigma_beta[j,k]^2) -2*X[i,k]*mu_beta[j,k]*mu_W[i,j]*
          ((1-z[i])*((M+v0+J-1)/b_Sigma0[j,j]) + z[i]* ((N-M+v1+J-1)/b_Sigma_1[j,j]))
      }))
    }))
  }))
  return(pw1+pw2)
}

#The expectation of eta has to be taken here:
.lb_pSigma = function(alpha_a, beta_a, v, B_Sigma, N, J){
  #alpha_a is a vector of size J
  #beta_a is a vector of size J
  #v is an hyper parameter 
  #B_Sigma is the scale matrix 
  #N the number of individuals
  #J the number of species 
  
  eta = diag(alpha_a/beta_a, ncol=J, nrow=J)
  return(-v*sum(diag(eta)*(N+v+J-1)*solve(B_Sigma))- (v+2*J)*(2*v* sum(log(diag(eta))) - J*log(2) -
           sum(sapply(1:J, function(x){
             digamma((v+x-1)/2)
           }))) + v*(v+J-1)*log(det(eta))-((v+J-1)*J/2)*log(2) - lgamma((v+J-1)/2))
}


# .lb_py = function(Y, y_ref, mu_w){
#   #Y is a matrix of count for each species in each individual
#   #y ref is the reference species 
#   #mu_w is a matrix of size N*J
#   
#   return(sum(sapply(1:nrow(Y), function(i){
#     sum(sapply(1:ncol(Y), function(j){
#       (mu_w[i,j] - log(sum(exp(mu_w[i,]) + 1)))*Y[i,j] - y_ref[i]*log(sum(exp(mu_w[i,]) + 1))
#     }))
#   })))
# }

compute_lb = function(Y, y_ref, X, z, epsilon, alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0,
                      alpha_a0, beta_a0, alpha_a1, beta_a1,A, mu_W, sigma_W, b_Sigma0,  eta0, b_Sigma_1,  eta1, v0, v1,N,M,J, K){
  
  
  return(.lb_plambda(alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0) +
           .lb_ph(Y, mu_W, sigma_W, epsilon)+
           .lb_pSigma(alpha_a0, beta_a0,v0, B_Sigma0, M, J)  + .lb_pSigma(alpha_a1, beta_a1, v1, B_Sigma1, N-M, J)  + 
           .lb_pw(X, z, mu_W, sigma_W, mu_beta, sigma_beta, b_Sigma0, b_Sigma_1, v0, v1, N,M,J,K) + 
           .lb_pa(alpha_a0, beta_a0, A) + .lb_pa(alpha_a1, beta_a1,A))
  # .lb_py(Y, y_ref, mu_W)
}