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