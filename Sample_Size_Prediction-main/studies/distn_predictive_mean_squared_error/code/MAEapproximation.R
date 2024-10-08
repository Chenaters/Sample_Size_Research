# simulation to confirm the result of SAWYER
set.seed(3)
 p = 10 # n observations and p predictors
for (n in c(25,50,75,100)){
 N = 1000000  # 

  W = rf(N, df1 = p, df2 = n-p) #empirical distribution of (y_hat-y*)
  X = rnorm(N, mean = 0, sd = 1)
  F_hat = X*sqrt((1+1/n)*(1+W*p/(n-p)))
F_hat_CDF = ecdf(F_hat)  #empirical CDF from 10000 observation

range = 5
list = seq(-range,range,0.01)   # for fitted value
F_hat_fitted = F_hat_CDF(list)
sigma_error = 1  #under assumption
MSE = sigma_error^2*(n+1)*(n-2)/(n*(n-p-2)) # true MSE under assumption
MAE_ep = mean(abs(F_hat))  #empirical MAE
sigma_prime = sqrt(MSE)
MAE_ep;sigma_prime

#MAE approximation from(2a) and (2b)
MAE_1 = sigma_prime*sqrt(2/pi)
MAE_2 = sigma_prime*sqrt(2/pi)-sigma_prime*sqrt(2/pi)*p/(4*(n-2)*(n-p-4))
diff_MAE1 = abs(MAE_1-MAE_ep)
diff_MAE2 = abs(MAE_2-MAE_ep)

#approximation with the 1st term from (1)
F_1_CDF = ecdf(rnorm(1000000, 0, sigma_prime))
F_1_fitted = F_1_CDF(list)

F_1_cdf = function(t) {pnorm(t/sigma_prime)}

#difference between F_hat and F_1
plot(F_1_fitted, x = list, xlim = c(-range,range), col = "red", type = "l");
lines(F_1_cdf(list), x = list, xlim = c(-range,range), col = "blue", type = "l") # add line to plot
lines(F_hat_fitted, x = list, xlim = c(-range,range), col = "green");

diff_F1 = max(abs(F_hat_fitted-F_1_fitted))


#approximation with 1st and 2nd term from (1)
F_2 = function(t){
  pnorm(t/sigma_prime)-(p/(4*sigma_error^3*(n-2)*(n-p-4))*
                          ((t/sigma_prime)^3-3*(t/sigma_prime)))*dnorm(t/sigma_prime)
}
F_2_fitted = F_2(list)#fitted value
F_2_fitted1 = F_2_fitted[which(F_2_fitted<1)]
plot(F_hat_fitted[1:length(F_2_fitted1)], x = list[1:length(F_2_fitted1)], col = "red")
lines(F_2_fitted1, x = list[1:length(F_2_fitted1)], col = 'purple')
diff_F2 = max(abs(F_hat_fitted[1:length(F_2_fitted1)]-F_2_fitted1))


#P_hat and P_1 from (3a)
list_P = seq(0,10,0.001)
P_hat = abs(F_hat)
P_hat_CDF = ecdf(P_hat)
P_hat_fitted = P_hat_CDF(list_P)
P_1_CDF = function(t)
{
  2*pnorm(t, mean = 0, sd = sigma_prime)-1
  }
P_1_fitted = P_1_CDF(list_P)
#comparison between two way of P_1 and P_hat
plot(list_P, P_1_fitted, type = "l", col = "red");
lines(list_P, P_hat_fitted, col = "blue")
diff_P1 = max(abs(P_hat_fitted-P_1_fitted))



##P_2 (3b)
P_2_CDF = function(t)
{
  2*pnorm(t/sigma_prime, mean = 0, sd = 1)-1+(3*(t/sigma_prime)-(t/sigma_prime)^3)*dnorm(t/sigma_prime, mean = 0, sd = 1)*p/(2*(n-2)*(n-p-4))
}
P_2_fitted = P_2_CDF(list_P)
plot(list_P, P_hat_fitted, type = "l", col = "red");
lines(list_P, P_2_fitted, col = "green")
diff_P2 = max(abs(P_hat_fitted-P_2_fitted))   

diff = matrix(c(diff_F1,diff_F2,diff_P1,diff_P2,diff_MAE1,diff_MAE2),nrow=1, ncol=6)

print(diff)

}
