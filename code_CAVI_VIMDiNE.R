.has_converged = function(ELBO_0, ELBO_1, threshold) abs(ELBO1-ELBO0) < threshold

CAVI_MDINE = function(Y, y_ref, X, z, v0, v1, A, delta0, r0,
                      init_omega_lambda, 
                      init_mu_beta, init_sigma_beta, 
                      init_mu_W, init_sigma_W,
                      init_b_Sigma0,init_b_Sigma1,
                      init_beta_a0, init_beta_a1, init_epsilon, threshold=1e-06, niter=10){
  
  N = length(z)
  J = ncol(Y)
  K = ncol(X)
  M = sum(z==0)
  
  res = list()
  
  if(length(z) != nrow(X) | length(z) != nrow(Y) | nrow(Y) != nrow(X)){
    stop("z, Y and X do not have the same number of individuals, please check...")
  }
  
  alpha_lambda = matrix(rep(r0+0.5, J*K), nrow=J, ncol=K) #shape of the Gamma distribution Is not updated through variational inference ! 
  alpha_a0 = 0.5*(v0+J) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  alpha_a1 = 0.5*(v1+J) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  
  res[['alpha_lambda']][[1]] = alpha_lambda
  res[['omega_lambda']][[1]] = init_omega_lambda
  res[['mu_beta']][[1]] = init_mu_beta
  res[['sigma_beta']][[1]] = init_sigma_beta
  res[["mu_W"]][[1]] = init_mu_W
  res[["sigma_W"]][[1]] = init_sigma_W
  res[["b_Sigma0"]][[1]] = init_b_Sigma0
  res[["b_Sigma1"]][[1]] = init_b_Sigma1
  res[["alpha_a0"]][[1]] = alpha_a0
  res[["beta_a0"]][[1]] = init_beta_a0
  res[["alpha_a1"]][[1]] =alpha_a1
  res[["beta_a1"]][[1]] = init_beta_a1
  res[["epsilon"]][[1]] = init_epsilon
  res[["ELBO"]][[1]] = compute_ELBO(Y, y_ref, X, z, init_epsilon, v0, v1, A, delta0, r0,
                      alpha_lambda, init_omega_lambda,
                      init_mu_beta, init_sigma_beta,
                      init_mu_W, init_sigma_W,
                      init_b_Sigma0,init_b_Sigma1,
                      alpha_a0, init_beta_a0,
                      alpha_a1, init_beta_a1, N,M,J,K)
  
  for(iter in 2:niter){
    print(iter)
    
    mu_beta_prev = res[['mu_beta']][[iter-1]]
    sigma_beta_prev = res[['sigma_beta']][[iter-1]]
    
    mu_W_prev = res[['mu_W']][[iter-1]]
    sigma_W_prev = res[['sigma_W']][[iter-1]]
    
    b_Sigma0_prev = res[["b_Sigma0"]][[iter-1]]
    b_Sigma1_prev = res[["b_Sigma1"]][[iter-1]]
    
    epsilon_prev =  res[["epsilon"]][[iter-1]]
    
    omega_lambda = 0.5*(mu_beta_prev^2+sigma_beta_prev)+delta0 #J*K 
    
    mu_beta = sapply(1:K, function(k){
      sapply(1:J, function(j){
        sum(sapply(1:N, function(i) (mu_W_prev[i,j]*X[i,k])/((alpha_lambda[j,k]/omega_lambda[j,k]) + X[i,k]^2) ))
      })
    })
    
    sigma_beta = sapply(1:K, function(k){
      sapply(1:J, function(j){
        sum(sapply(1:N, function(i) ((alpha_lambda[j,k]/omega_lambda[j,k])+X[i,k]^2)*((1-z[i])*(M+v0+J-1)/(b_Sigma0_prev[j,j]) + z[i]*(N-M+v1+J-1)/(b_Sigma1_prev[j,j])) ))
      })
    })
    
    beta_a0 = sapply(1:J, function(j) v0*(M+v0+J-1)/b_Sigma0_prev[j,j] + 1/A[j])
    beta_a1 = sapply(1:J, function(j) v1*(N-M+v1+J-1)/b_Sigma1_prev[j,j] + 1/A[j])
    
    mu_eta0 =  diag(0.5*(v0+J)/beta_a0, J,J)
    mu_eta1 = diag(0.5*(v1+J)/beta_a1, J,J)
    
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
      sapply(1:J, function(j) {
        sapply(1:N, function(i){
          (K+1)*((1-z[i])*((M+v0+J-1)/b_Sigma0[j,j]) + z[i]*((N-M+v1+J-1)/b_Sigma1[j,j]))
        })
      })
        
    print(sigma_W)
    mu_W = ((X%*%t(mu_beta))*sapply(1:J, function(j){
      sapply(1:N, function(i){
        (((1-z[i])*((M+v0+J-1)/b_Sigma0[j,j])+z[i]*((N-M+v1+J-1)/b_Sigma1[j,j]))*(Y[i,j] - 0.5))/sigma_W[i,j]
      })
    }))
    
    
   
    epsilon = sqrt(sapply(1:J, function(j) {
      sapply(1:N, function(i){
        mu_W[i,j]^2*sigma_W[i,j]
      })
    }))
     
    print(epsilon)
    
    ELBO = compute_ELBO(Y, y_ref, X, z, epsilon, v0, v1, A, delta0, r0,
                        alpha_lambda, omega_lambda,
                        mu_beta, sigma_beta,
                        mu_W, sigma_W,
                        b_Sigma0, b_Sigma1,
                        alpha_a0, beta_a0,
                        alpha_a1, beta_a1, N,M,J,K)
    
    res[['alpha_lambda']][[iter]] = alpha_lambda
    res[['omega_lambda']][[iter]] = omega_lambda
    res[['mu_beta']][[iter]]= mu_beta
    res[['sigma_beta']][[iter]] = sigma_beta
    res[["mu_W"]][[iter]] = mu_W
    res[["sigma_W"]][[iter]] = sigma_W
    res[["b_Sigma0"]][[iter]] = b_Sigma0
    res[["b_Sigma1"]][[iter]] = b_Sigma1
    res[["alpha_a0"]][[iter]] = alpha_a0
    res[["beta_a0"]][[iter]] = beta_a0
    res[["alpha_a1"]][[iter]] =alpha_a1
    res[["beta_a1"]][[iter]] = beta_a1
    res[["epsilon"]][[iter]] = epsilon
    res[["ELBO"]][[iter]] = ELBO
  }
  
  
}


