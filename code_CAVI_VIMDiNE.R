library(CholWishart)
library(MCMCpack)

.compute_entropy_lambda = function(alpha_lambda, omega_lambda){
  #Alpha and omega are matrices of size J * K+1
  #corresponding to the number of species and covariates
  return(sum(sapply(1:ncol(alpha_lambda), function(x) sum(alpha_lambda[,x] - log(omega_lambda[,x]) + lgamma(alpha_lambda[,x]) + (1-alpha_lambda[,x])*digamma(alpha_lambda[,x])))))
}


.compute_entropy_W = function(sigma_W, N,J){
  #sigma is a matrix of size N * J
  #corresponding to the number of individuals and species
  return(0.5*(sum(sapply(1:ncol(sigma_W), function(x) sum(log(2*pi*sigma_W[,x])))) + N*J))
}


.compute_entropy_beta = function(sigma_beta, K,J){
  #sigma is a matrix of size J * K+1
  #corresponding to the number of species and covariates
  return(0.5*(sum(sapply(1:ncol(sigma_beta), function(x) sum(log(2*pi*sigma_beta[,x])))) + (K+1)*J))
}

.compute_entropy_a = function(alpha_a, B_a){
  #alpha and B are vectors of size J corresponding to the number of species
  
  return(sum(alpha_a + log(B_a)*lgamma(alpha_a) + (1-alpha_a)*digamma(alpha_a)))
}

.compute_entropy_Sigma = function(v, B_Sigma, N,J){
  #v is the degree of freedom of the inverse wishart
  #B the scale matrix
  return(log(det(B_Sigma))*(N+v) - 0.5*((N+v)*J)*log(2) -  lmvgamma((N+v)*0.5, J) -
           0.5*(v+N+J-1)*(log(det(B_Sigma))-J*log(2) - sum(sapply(1:J, function(x){
             digamma((v+N-J+x)*0.5)
           }))) - (J/(2*(v+N))))
}

compute_entropy = function(alpha_lambda, omega_lambda,sigma_W,sigma_beta, alpha_a0, B_a0, B_Sigma0,alpha_a1, B_a1,B_Sigma1,v0,v1, N,M,K,J){
  #Compute the entropy for the calculation of the ELBO
  .compute_entropy_lambda(alpha_lambda,omega_lambda) + .compute_entropy_W(sigma_W, N,J)+
    .compute_entropy_beta(sigma_beta, K,J)+.compute_entropy_a(alpha_a0,B_a0)+
    .compute_entropy_a(alpha_a1,B_a1)+.compute_entropy_Sigma(v0, B_Sigma0,M,J)+.compute_entropy_Sigma(v1, B_Sigma1,N-M,J)
}

.lb_plambda = function(alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0){
  
  E_q_lambda = function(fn) {
    integrate(function(lambda) {
      dgamma(lambda, alpha_lambda, omega_lambda) * fn(lambda)
    }, 0, Inf)$value
  }
  
  return(0.5*E_q_lambda(function(lambda) log(lambda))-0.5*(alpha_lambda/omega_lambda)*(mu_beta^2 + sigma_beta) + 
    (r0-1)*E_q_lambda(function(lambda) log(lambda))*(log(delta0)- lgamma(r0)) - delta0*(alpha_lambda/omega_lambda))
}

.lb_pa = function(alpha_a, beta_a, A){
  E_q_a = function(fn) {
    integrate(function(a) {
      dinvgamma(a, alpha_a, beta_a) * fn(a)
    }, 0, Inf)$value
  }
  return(-0.5*log(A) - lgamma(0.5)+(0.5+1)*E_q_a(function(a) log(a)) - (1/A)*(alpha_a/beta_a))
}

.lb_pw = function(x,z,mu_W, sigma_W, mu_beta, sigma_beta, b_Sigma0, b_Sigma_1, v0, v1,N,M,J, K){
  return(0.5*((K+1)*(mu_W^2 + sigma_W) + x*(mu_beta^2 + sigma_beta^2) -2*x*mu_beta*mu_W)*
   ((1-z)*((M+v0+J-1)/b_Sigma0) + z* ((N-M+v1+J-1)/b_Sigma_1)) )
}

.lb_pSigma = function(eta, v, B_Sigma, N, J){
  -0.5*sum(diag(eta))*(N+v+J-1)- (v+J-1)*log(det(B_Sigma)) - J*log(2) -
   sum(sapply(1:J, function(x){
     digamma((N+v-J+x)/2)
   })) + 0.5*v*log(det(eta))-(v*J/2)*log(2) - lgamma(v/2)
}

.lb_py = function(y, y_ref, mu_w){
  
  return((mu_w - log(sum(exp(mu_w) + 1)))*y - y_ref*log(sum(exp(mu_w) + 1)))
}

compute_lb = function(y, y_ref, X, Z, alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0,
                      alpha_a0, beta_a0, A0, alpha_a1, beta_a1, A1, mu_W, sigma_W, mu_beta, sigma_beta, b_Sigma0, b_Sigma_1, v0, v1,N,M,J, K,
                      eta0, B_Sigma0, eta1, B_Sigma1){
  
  .lb_plambda() + .lb_py() + .lb_pSigma() + .lb_pw() + .lb_pa()
  
}

compute_ELBO = function(){
  
  return(compute_lb() - compute_entropy())
  
}

.has_converged = function(ELBO_0, ELBO_1, threshold) abs(ELBO1-ELBO0) < threshold

CAVI_MDINE = function(epsilon){
  
  ELBO = c()
  alpha_lambda = r0 + 0.5 #Is not updated through variational inference ! 
  i = 1
  
  ELBO[i] = compute_ELBO()
  
  while(!.has_converged(ELBO[i], ELBO)){
    
  } 
}


