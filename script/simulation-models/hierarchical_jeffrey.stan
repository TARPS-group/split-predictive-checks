
data {
  int<lower=0> I; 
  int<lower=0> J; 
  matrix[I,J]  X; 
}

parameters {
  real mu0;
  real<lower=0> tau0;
  vector[I] theta; 
}

model {
// Priors for random coefficients:
  for (i in 1:I)
    theta[i] ~ normal(mu0, tau0);
    
  for(i in 1:I){
    for (j in 1:J){
      X[i,j] ~ normal(theta[i], 2);
    }
  }
target += -log(tau0);
}

