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
  alpha_a0 = 0.5*(v0+J) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  alpha_a1 = 0.5*(v1+J) #shape of the inverse Gamma distribution Is not updated through variational inference ! 
  
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


