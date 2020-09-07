# rm(list=ls()) 

## load packages
library("tidyverse")
require("MASS")

## true parameters
set.seed(321)
n <- 50        # sample size 
reps <- 10000  # number of simulations 
b0 <- 300      # true intercept for house price in thousands
b1 <- 4        # true parameter for price per thousands of square feet 
b2 <- -7       # true parameter for age of house 
sigma <- 10    # residuals/noise 
e <- rnorm(n,0,sigma) # residuals normally distributed w/ mean 0 and variance 10

# means of individual distributions
mu_sqft <- 3   
mu_age_years <- 10   
mu_age_months <- 120
# variance  
sigma_sqft <- 1  
sigma_age_years <- 5  
sigma_age_months <- 60 
# covariance  
X_sqft <- 0.5 # equivalent to 20% correlation between sqft and age_years
X_age_years <- 17  # equivalent to near perfect correlation between age_years and age_months
X_age_months <- 1.7 # equivalent to 20% correlation between sqft and age_months

## initialize dataframes
# true
beta0_true.hat <- matrix(nrow=reps,ncol=1)
beta1_true.hat <- matrix(nrow=reps,ncol=1)
beta2_true.hat <- matrix(nrow=reps,ncol=1)
beta3_true.hat <- matrix(nrow=reps,ncol=1)
mse_true = matrix(nrow=reps,ncol=1)
pvalue_age_true = matrix(nrow=reps,ncol=1)

# linearity  
beta0_lh.hat <- matrix(nrow=reps,ncol=1)
beta1_lh.hat <- matrix(nrow=reps,ncol=1)
beta2_lh.hat <- matrix(nrow=reps,ncol=1)
mse_lh = matrix(nrow=reps,ncol=1)
pvalue_age_lh = matrix(nrow=reps,ncol=1)

# endogeneity  
beta0_e.hat <- matrix(nrow=reps,ncol=1)
beta1_e.hat <- matrix(nrow=reps,ncol=1)
beta2_e.hat <- matrix(nrow=reps,ncol=1)
mse_e = matrix(nrow=reps,ncol=1)

# multicollinearity  
beta0_m.hat <- matrix(nrow=reps,ncol=1)
beta1_m.hat <- matrix(nrow=reps,ncol=1)
beta2_m.hat <- matrix(nrow=reps,ncol=1)
beta3_m.hat <- matrix(nrow=reps,ncol=1)
beta4_m.hat <- matrix(nrow=reps,ncol=1)
mse_m = matrix(nrow=reps,ncol=1)
pvalue_age_m = matrix(nrow=reps,ncol=1)

## run loop
for(i in 1:reps) {
  
  dat <- mvrnorm(n, mu = c(mu_sqft, mu_age_years, mu_age_months),
                 Sigma = matrix(c(sigma_sqft,    X_sqft,          X_age_months,
                                  X_sqft,        sigma_age_years, X_age_years,
                                  X_age_months,  X_age_years,     sigma_age_months),
                                  ncol = 3, byrow = TRUE),
                 empirical = TRUE)
  
  # fit model with true parameters
  Y <- b0 + b1*dat[,1] + b1*I(dat[,1]^2) + b2*dat[,2] + e; Y

  # true model
  eq1 <- lm(Y~ poly(dat[,1],2,raw=F)+dat[,2]) # notice the use of orthogonal polynomials (raw=F)
  beta0_true.hat[i] <- eq1$coefficients[1]
  beta1_true.hat[i] <- eq1$coefficients[2]
  beta2_true.hat[i] <- eq1$coefficients[3]
  beta3_true.hat[i] <- eq1$coefficients[4]
  mse_true[i] = mean(resid(eq1)^2)
  pvalue_age_true[i] = ifelse(summary(eq1)$coefficients[4,4] > 0.05,1,0)
  
  # linearity  
  eq1 <- lm(Y~dat[,1]+dat[,2])
  beta0_lh.hat[i] <- eq1$coefficients[1]
  beta1_lh.hat[i] <- eq1$coefficients[2]
  beta2_lh.hat[i] <- eq1$coefficients[3]
  mse_lh[i] = mean(resid(eq1)^2)
  pvalue_age_lh[i] = ifelse(summary(eq1)$coefficients[3,4] > 0.05,1,0)
  
  # endogeneity
  eq1 <- lm(Y~dat[,2])
  beta0_e.hat[i] <- eq1$coefficients[1]
  beta1_e.hat[i] <- eq1$coefficients[2]
  mse_e[i] = mean(resid(eq1)^2)
  
  # multicollinearity
  eq1 <- lm(Y~poly(dat[,1],2,raw=F)+dat[,2]+dat[,3])
  beta0_m.hat[i] <- eq1$coefficients[1]
  beta1_m.hat[i] <- eq1$coefficients[2]
  beta2_m.hat[i] <- eq1$coefficients[3]
  beta3_m.hat[i] <- eq1$coefficients[4]
  beta4_m.hat[i] <- eq1$coefficients[5]
  mse_m[i] = mean(resid(eq1)^2)
  pvalue_age_m[i] = ifelse(summary(eq1)$coefficients[4,4] > 0.05,1,0)
  
}

## Store results in df 
df_parameters = data.frame("Intercept_true" = beta0_true.hat,"SquareFeet_true" = beta1_true.hat,"SquareFeet_squared_true" = beta2_true.hat,
                           "Age_true" = beta3_true.hat,"MSE_true" = mse_true,
                           "Intercept_lh" = beta0_lh.hat,"SquareFeet_lh" = beta1_lh.hat,"Age_lh" = beta2_lh.hat,"MSE_lh" = mse_lh,
                           "Intercept_e" = beta0_e.hat,"Age_e" = beta1_e.hat,"MSE_e" = mse_e,
                           "Intercept_m" = beta0_m.hat,"SquareFeet_m" = beta1_m.hat,"SquareFeet_squared_m" = beta2_m.hat,
                           "Age_years_m" = beta3_m.hat,"Age_month_m" = beta4_m.hat,"MSE_m" = mse_m) # store simulation results


## multicollinearity violated ##

# vif
# library(car)
round(vif(lm(Y~poly(dat[,1],2)+dat[,2]+dat[,3])),1) # We can apply the usual VIF rule of thumb if we squared the GVIF^(1/(2*Df)) value

# check MSE (violation)
colors <- c("1" = "red", "2" = "blue")
p1 = ggplot(df_parameters[,c(5,18)],aes(x=MSE_true)) + 
  geom_histogram(aes(fill = "1")) + 
  geom_histogram(aes(x=MSE_m,fill = "2")) + 
  labs(x = "MSE", y = "Frequency", title = "Distribution of 10K MSE",fill="Model") + 
  scale_color_manual(values = colors);p1

p2 = ggplot(df_parameters[,c(4,16)],aes(x=Age_true)) + 
  geom_density(aes(fill = "1")) + 
  geom_density(aes(x=Age_years_m, fill = "2")) + 
  labs(x = expression(hat(beta)[3]), y="Density", title = "Distribution of 10K Coefficient Estimates for age_years", fill="Model") + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks=seq(-20,5,2.5));p2

# compare p-values
mean(pvalue_age_true) # not statistically significant 0% of the time, sigma=10
mean(pvalue_m) # not statistically significant 45% of the time




## linearity violated ##

# linearity violation (violation, fitted vs resid)
resid_vs_fit = data.frame("Residuals" = lm(Y~dat[,1]+dat[,2])$residuals, "Fitted.Values" = lm(Y~dat[,1]+dat[,2])$fitted.values)
p3 = resid_vs_fit %>% ggplot(aes(x=Fitted.Values, y=Residuals)) + geom_point() + geom_smooth(method=loess,se=F) + 
  labs(x = "Fitted Values", title = "Linearity Violation");p3

# check MSE (true model)
p4 = ggplot(df_parameters[,c(5,9)],aes(x=MSE_true)) + 
  geom_histogram(fill = "salmon") + 
  labs(x = "MSE", y = "Frequency", title = "Model 1");p4
mean(mse_true) # 84.03434

# check MSE (violation)
p5 = ggplot(df_parameters[,c(5,9)],aes(x=MSE_lh)) + 
  geom_histogram(fill = "cyan3") + 
  labs(x = "MSE", y = "Frequency", title = "Model 3");p5
mean(mse_lh)   # 112.675

# library(gridExtra)
# library(grid)
grid.arrange(p4, p5, nrow = 1, ncol = 2) 

# Age coefficient
colors <- c("1" = "red", "3" = "blue")
p6 = ggplot(df_parameters[,c(4,8)],aes(x=Age_true)) + 
  geom_density(aes(fill = "1")) + 
  geom_density(aes(x=Age_lh, fill = "3")) + 
  labs(x = expression(hat(beta)), y="Density", title = "Distribution of 10K Coefficient Estimates for age_years", fill="Model") + 
  scale_color_manual(values = colors);p6




## endogeneity violated ##

# mean estimate for age is biased since correlation between sqft and age
mean(beta1_e.hat)    # -4.157
sd(beta1_e.hat)*1.96 # 1.38816

# check MSE (violation)
p7 = ggplot(df_parameters[,c(5,12)],aes(x=MSE_e)) + 
  geom_histogram(fill = "cyan3") + 
  labs(x = "MSE", y = "Frequency", title = "Model 4");p7
mean(mse_true) # 84.03434
mean(mse_e)   # 846.7734

grid.arrange(p4, p7, nrow = 1, ncol = 2) 

# create function to visualize distribution of sample means
sample_means_df = function(col,left_bound,right_bound) {
  dens <- density(df_parameters[,col])
  data <- tibble(x = dens$x, y = dens$y) %>% 
    mutate(variable = case_when(
      (x >= left_bound & x <= right_bound) ~ "On",
      (x >= right_bound & x <= left_bound) ~ "Off",
      TRUE ~ NA_character_))
  data
}

p8 = ggplot(sample_means_df(11,-5.54516,-2.769482), aes(x, y)) + geom_line() +
  geom_area(data = dplyr::filter(sample_means_df(11,-5.54516,-2.769482), variable == 'On'), fill = 'skyblue2') + 
  geom_area(data = dplyr::filter(sample_means_df(11,-5.54516,-2.769482), variable == 'Off'), fill = 'light blue') +
  labs(x=expression(hat(beta)[1]), y="Density", title=expression("Distribution of 10K Coefficient Estimates for age_years")) + 
  annotate(geom="text", x= -4.157642, y=0.15, label="Middle 95% of Data", color="red", size = 8) +
  geom_vline(xintercept = -4.157642, linetype="dotted", color = "blue", size=1) + 
  scale_x_continuous(breaks=seq(-10,0,1));p8