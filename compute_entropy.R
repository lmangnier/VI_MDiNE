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
  return(log(det(B_Sigma))*(N+v+2*J) - 0.5*((N+v+J-1)*J)*log(2) -  lmvgamma((N+v+J-1)*0.5, J) -
           0.5*(v+N+2*J)*(log(det(B_Sigma))-J*log(2) - sum(sapply(1:J, function(x){
             digamma((v+N-1+x)*0.5)
           }))) - 0.5*(J*(v+N+J-1)))
}

compute_entropy = function(alpha_lambda, omega_lambda,sigma_W,sigma_beta, alpha_a0, B_a0, B_Sigma0,alpha_a1, B_a1,B_Sigma1,v0,v1, N,M,K,J){
  #Compute the entropy for the calculation of the ELBO
  .compute_entropy_lambda(alpha_lambda,omega_lambda) + .compute_entropy_W(sigma_W, N,J)+
    .compute_entropy_beta(sigma_beta, K,J)+.compute_entropy_a(alpha_a0,B_a0)+
    .compute_entropy_a(alpha_a1,B_a1)+.compute_entropy_Sigma(v0, B_Sigma0,M,J)+.compute_entropy_Sigma(v1, B_Sigma1,N-M,J)
}