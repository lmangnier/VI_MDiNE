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

.sigma_epsilon = function(epsilon){
  #Epsilon is a n*J matrix
  return(epsilon/(rowSums(epsilon)+1))
}
.lambda_epsilon = function(epsilon) 0.5*epsilon*(.sigma_epsilon(epsilon)-0.5)


.lb_ph = function(Y, mu_W, sigma_W, epsilon){
  
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
    -0.5*log(A) - lgamma(0.5)+(0.5+1)*E_q_a(function(a) log(a), alpha_a[x], beta_a[x]) - (1/A)*(alpha_a[x]/beta_a[x])
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


.lb_pSigma = function(eta, v, B_Sigma, N, J){
  #eta is a diagonal matrix
  #v is an hyper parameter 
  #B_Sigma is the scale matrix 
  #N the number of individuals
  #J the number of species 
  return(-0.5*sum(diag(eta)*(N+v+J-1)*solve(B_Sigma))- (v+J-1)*log(det(B_Sigma)) - J*log(2) -
   sum(sapply(1:J, function(x){
     digamma((N+v-J+x)/2)
   })) + 0.5*v*log(det(eta))-(v*J/2)*log(2) - lgamma(v/2))
}


.lb_py = function(Y, y_ref, mu_w){
  #Y is a matrix of count for each species in each individual
  #y ref is the reference species 
  #mu_w is a matrix of size N*J
  
  return(sum(sapply(1:nrow(Y), function(i){
    sum(sapply(1:ncol(Y), function(j){
      (mu_w[i,j] - log(sum(exp(mu_w[i,]) + 1)))*Y[i,j] - y_ref[i]*log(sum(exp(mu_w[i,]) + 1))
    }))
  })))
}

compute_lb = function(Y, y_ref, X, z, epsilon, alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0,
                      alpha_a0, beta_a0, alpha_a1, beta_a1,A, mu_W, sigma_W, b_Sigma0,  eta0, b_Sigma_1,  eta1, v0, v1,N,M,J, K){
  
  
  return(.lb_plambda(alpha_lambda, omega_lambda, mu_beta, sigma_beta, r0, delta0) +
    .lb_py(Y, y_ref, mu_W) + .lb_ph(Y, mu_W, sigma_W, epsilon)+
    .lb_pSigma(v0, eta0, B_Sigma0, M, J)  + .lb_pSigma(v1, eta1, B_Sigma1, N-M, J)  + 
    .lb_pw(X, z, mu_W, sigma_W, mu_beta, sigma_beta, b_Sigma0, b_Sigma_1, v0, v1, N,M,J,K) + 
    .lb_pa(alpha_a0, beta_a0, A) + .lb_pa(alpha_a1, beta_a1,A))
  
}

compute_ELBO = function(Y, y_ref, X, z, epsilon, v0, v1, A, delta0, r0,
                        alpha_lambda, omega_lambda,
                        mu_beta, sigma_beta,
                        mu_W, sigma_W,
                        eta0, b_Sigma0,
                        eta1, b_Sigma1,
                        alpha_a0, beta_a0,
                        alpha_a1, beta_a1, N, M, J, K){
  
  return(compute_lb(Y, y_ref, X, z, epsilon, alpha_lambda, omega_lambda,
                    mu_beta, sigma_beta, r0, delta0,
                    alpha_a0, beta_a0,
                    alpha_a1, beta_a1, A,mu_W, sigma_W,
                    b_Sigma0, eta0, b_Sigma1,  v0,v1,N,M,J,K) - compute_entropy(alpha_lambda, omega_lambda,
                                                                                sigma_W, sigma_beta, 
                                                                                alpha_a0, beta_a0,b_Sigma0,
                                                                                alpha_a1, beta_a1,b_Sigma1,
                                                                                v0, v1,N,M,K,J))
  
}

.has_converged = function(ELBO_0, ELBO_1, threshold) abs(ELBO1-ELBO0) < threshold

CAVI_MDINE = function(Y, y_ref, X, z, v0, v1, A, delta0, r0,
                      init_omega_lambda, 
                      init_mu_beta, init_sigma_beta, 
                      init_mu_W, init_sigma_W,
                      init_mu_eta0, init_b_Sigma0,
                      init_mu_eta1, init_b_Sigma1,
                      init_beta_a0,
                      init_beta_a1, init_epsilon, threshold=1e-06){
  
  N = length(z)
  J = ncol(Y)
  K = ncol(X)
  M = sum(z==0)
  
  
  if(length(z) != nrow(X) | length(z) != nrow(Y) | nrow(Y) != nrow(X)){
    stop("z, Y and X do not have the same number of individuals, please check")
  }
  
  alpha_lambda = matrix(rep(r0, J*K), nrow=J, ncol=K) #shape of the Gamma distribution Is not updated through variational inference ! 
  alpha_a0 = 0.5*(v0+1) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  alpha_a1 = 0.5*(v1+1) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  
  res = list()
  
  res[['alpha_lambda']] = alpha_lambda
  res[['omega_lambda']] = init_omega_lambda
  res[['mu_beta']] = init_mu_beta
  res[['sigma_beta']] = init_sigma_beta
  res[["mu_W"]] = init_mu_W
  res[["sigma_W"]] = init_sigma_W
  res[["mu_eta0"]] = init_mu_eta0
  res[["b_Sigma0"]] = init_b_Sigma0
  res[["mu_eta1"]] = init_mu_eta1
  res[["b_Sigma1"]] = init_b_Sigma1
  res[["alpha_a0"]] = alpha_a0
  res[["beta_a0"]] = init_beta_a0
  res[["alpha_a1"]] =alpha_a1
  res[["beta_a1"]] = init_beta_a1
  res[["epsilon"]] = init_epsilon
  res[["ELBO"]] = 0
  
  iter = 1 
  
  ELBO = compute_ELBO(Y, y_ref, X, z, init_epsilon, v0, v1, A, delta0, r0,
                      alpha_lambda, init_omega_lambda,
                      init_mu_beta, init_sigma_beta,
                      init_mu_W, init_sigma_W,
                      init_eta0, init_b_Sigma0,
                      init_eta1, init_b_Sigma1,
                      alpha_a0, init_beta_a0,
                      alpha_a1, init_beta_a1, N,M,J,K)
  
  while(!.has_converged(res[["ELBO"]][iter], ELBO, threshold)){
    
    
    mu_beta_prev = res[['mu_beta']][[iter]]
    sigma_beta_prev = res[['sigma_beta']][[iter]]
    
    mu_W_prev = res[['mu_W']][[iter]]
    sigma_W_prev = res[['sigma_W']][[iter]]
    
    b_Sigma0_prev = res[["b_Sigma0"]][[iter]]
    b_Sigma1_prev = res[["b_Sigma1"]][[iter]]
    
    epsilon_prev =  res[["epsilon"]][[iter]]
    
    omega_lambda = 0.5*(mu_beta_prev^2+sigma_beta_prev)+delta0 #J*K 
    
    mu_beta = sapply(1:J, function(j){
      sapply(1:K, function(k){
        sum(sapply(1:N, function(i) (mu_W_prev[i,j]*X[i,k])/((alpha_lambda[j,k]/omega_lambda[j,k]) + X[i,k]^2) ))
      })
    })
    
    sigma_beta = sapply(1:J, function(j){
      sapply(1:K, function(k){
        sum(sapply(1:N, function(i) ((alpha_lambda[j,k]/omega_lambda[j,k])+X[i,k]^2)*((1-z[i])*(M+v0+J-1)/(res[["b_Sigma0"]][j,j]) + z[i]*(N-M+v1+J-1)/(res[["b_Sigma1"]][j,j])) ))
      })
    })
    
    beta_a0 = sapply(1:J, function(j) v0*(M+v0+J-1)/b_Sigma0_prev[j,j] + 1/A)
    beta_a1 = sapply(1:J, function(j) v1*(N-M+v1+J-1)/b_Sigma1_prev[j,j] + 1/A)
    
    mu_eta0 =  diag(0.5*(v0+1)*beta_a0, J,J)
    mu_eta1 = diag(0.5*(v1+1)*beta_a1, J,J)
    
    b_Sigma0 = 2*v0*mu_eta0  + sum(sapply(1:N, function(i) {
      sum(sapply(1:J, function(j) {
        sum(sapply(1:K, function(k) (mu_W_prev[i,j] + X[i,k]*mu_beta[j,k])^2 + sigma_W_prev[i,j] +(X[i,k]*sigma_beta[j,k])^2))
      }))
    }))
    
    b_Sigma1 = 2*v1*mu_eta1 + sum(sapply(1:N, function(i) {
      sum(sapply(1:J, function(j) {
        sum(sapply(1:K, function(k) (mu_W_prev[i,j] + X[i,k]*mu_beta[j,k])^2 + sigma_W_prev[i,j] +(X[i,k]*sigma_beta[j,k])^2))
      }))
    }))
    
    sigma_W = 2*.lambda_epsilon(epsilon_prev) + 
      sapply(1:N, function(i) {
        sapply(1:J, function(j){
          (K+1)*((1-z[i])*((M+v0+J-1)/b_Sigma0[j,j]) + z[i]*((N-M+v1+J-1)/b_Sigma1[j,j]))
        })
      })
        
    
    
    mu_W = ((X%*%t(mu_beta))*sapply(1:N, function(i){
      sapply(1:J, function(j){
        ((1-z[i])*((M+v0+J-1)/b_Sigma0[j,j])+z[i]*((N-M+v1+J-1)/b_Sigma1[j,j]))*(Y[i,j] - 0.5)
      })
    }))%*%solve(sigma_W)
    
    epsilon = sqrt(mu_W^2 + sigma_W) 
    
    
    res[['alpha_lambda']] = alpha_lambda
    res[['omega_lambda']] = omega_lambda
    res[['mu_beta']] = mu_beta
    res[['sigma_beta']] = sigma_beta
    res[["mu_W"]] = mu_W
    res[["sigma_W"]] = sigma_W
    res[["mu_eta0"]] = mu_eta0
    res[["b_Sigma0"]] = b_Sigma0
    res[["mu_eta1"]] = mu_eta1
    res[["b_Sigma1"]] = b_Sigma1
    res[["alpha_a0"]] = alpha_a0
    res[["beta_a0"]] = beta_a0
    res[["alpha_a1"]] =alpha_a1
    res[["beta_a1"]] = beta_a1
    res[["epsilon"]] = epsilon
    res[["ELBO"]] = c(res[["ELBO"]], ELBO)
    
    rep = rep+1 
    
    ELBO = compute_ELBO(Y, y_ref, X, z, epsilon, v0, v1, A, delta0, r0,
                        alpha_lambda, omega_lambda,
                        mu_beta, sigma_beta,
                        mu_W, sigma_W,
                        eta0, b_Sigma0,
                        eta1, b_Sigma1,
                        alpha_a0, beta_a0,
                        alpha_a1, beta_a1, N,M,J,K)
    
    
  }
  
  
}


