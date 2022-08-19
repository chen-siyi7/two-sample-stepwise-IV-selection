twosample_mvstepIV = function(p, R, betaZX, betaZY, se_betaZY, n1, n2, sigma2_hat, gamma_hat){
  r = dim(betaZX)[2]
  ZTZ = R
  n = n2
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:p){
    YTY[SNP] = (n2-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  
  testbic = NULL
  for (i in 1:p){
    test11 = diag(x = 0, nrow = p, ncol = p, names = T)
    diag(test11)[i] = 1
    W1 = cbind(test11, gamma_hat)
    solve.W1 = t(W1) %*% ZTZ %*% W1
    beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
    testbic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+log(n2)*sum(diag(test11))
  }
  
  # one by one 
  whichIV = NULL
  whichIV[1] = which.min(testbic)
  BICtest = NULL
  BICtest[1] = testbic[which.min(testbic)]
  for (j in 2:(p-2)){
    testbic = NULL
    for (i in 1:(p-2)){
      test11 = diag(x = 0, nrow = p, ncol = p, names = T)
      diag(test11)[whichIV] = 1
      diag(test11)[i] = 1
      W1 = cbind(test11, gamma_hat)
      solve.W1 = t(W1) %*% ZTZ %*% W1
      W1 = W1[,which(non0>0)]
      beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
      testbic[i] = n2*log(YTY - t(beta1)%*%t(W1)%*%ZTY)+log(n2)*sum(diag(test11))
    }
    whichIV[j] = which.min(testbic)
    BICtest[j] = testbic[which.min(testbic)]
    if(whichIV[j] == whichIV[j-1]) break; 
  }
  invalid.IVs = whichIV[!duplicated(whichIV)]
  # results
  test11 = diag(x = 0, nrow = p, ncol = p, names = T)
  diag(test11)[invalid.IVs] = 1
  W1 = cbind(test11, gamma_hat)
  solve.W1 = t(W1) %*% ZTZ %*% W1
  beta1 = matrix(ginv(solve.W1) %*% t(W1) %*% (ZTY), ncol = 1)
  beta_est = as.numeric(tail(beta1, n = r))
  sigma_u2 = as.numeric(YTY - t(beta1)%*%t(W1)%*%ZTY)
  Varbeta = tail(diag(ginv(solve.W1*n2)), n = r)*sigma_u2+t(beta_est)%*%sigma2_hat%*%beta_est*tail(diag(ginv(solve.W1*n2)), n = r)*n2/n1
  beta_se = as.numeric(sqrt(Varbeta))
  my_list = list('invalidIV' = sort(invalid.IVs), 'beta_est' = beta_est, 'beta_se' = beta_se, 'K' = length(invalid.IVs))
  return(my_list)
}


