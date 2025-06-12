require(MASS)
require(MBESS)
require(matrixStats)
require(psych)
################################################################################
set.seed(6202021)

N = 200       # number of examinees
n = 30        # number of items
pe <- 0.20    # proportion of examinees with item preknowledge
pi <- 0.50    # proportion of compromised items

# Generate the binary status of examinee item preknowledge
# 1: examinee has item preknowledge
# 0: examinee has item preknowledge

tmp <- runif(N,0,1)
H  <- ifelse(tmp<=quantile(tmp,pe),1,0)
H
table(H)

# Generate the binary status of item compromise
# 1: item is compromised
# 0: item is not compromised

tmp <- runif(n,0,1)
C  <- ifelse(tmp<=quantile(tmp,pi),1,0)
C
table(C)

################################################################################
#                            RESPONSE TIME DATA
################################################################################

# Generate item parameters

mu_beta        <- 3.5
mu_logalpha    <- 0.5
sigma_beta     <- 0.3
sigma_logalpha <- 0.2
omega_I        <- matrix(c(1,0.25,0.25,1),2,2)

mu_I     <- c(mu_beta,mu_logalpha)
Sigma_I  <- diag(c(sigma_beta,sigma_logalpha))%*%omega_I%*%diag(c(sigma_beta,sigma_logalpha))

item_par <- mvrnorm(n,mu=mu_I,Sigma=Sigma_I)  

beta  <- item_par[,1]    # time intensity parameters
alpha <- exp(item_par[,2])  # time discrimination parameters

# Generate person parameters

mu_taut      <- 0
sigma_taut   <- 0.1
mu_tauc      <- 0.4
sigma_tauc   <- 0.15
omega_P      <- matrix(c(1,0.7,0.7,1),2,2)

mu_P     <- c(mu_taut ,mu_tauc)
Sigma_P  <- diag(c(sigma_taut,sigma_tauc))%*%omega_P%*%diag(c(sigma_taut,sigma_tauc))

tau <- mvrnorm(N,mu_P,Sigma_P)

tau_t <- tau[,1]    # true latent speed parameters 
tau_c <- tau[,2]    # true cheating speed parameters 


# Generate observed response times according to the model

rt <- matrix(nrow=N,ncol=n)

for(i in 1:N){
  for(j in 1:n){
    
    p_t <- beta[j] - tau_t[i]
    p_c <- beta[j] - tau_c[i]
    
    if(H[i] == 1 & C[j] == 1){
      rt[i,j] = exp(rnorm(1,p_c,1/alpha[j]))
    } else {
      rt[i,j] = exp(rnorm(1,p_t,1/alpha[j]))
    }
    
  }
}

# Convert it to data frame

rt           <- as.data.frame(rt)
colnames(rt) <- paste0('RT',1:n)

################################################################################
#                            RESPONSE ACCURACY DATA
################################################################################

# Generate item difficulty parameters

b <- rnorm(n,0,1)
b <- (b-mean(b))/sd(b)
b
describe(b)

# Generate person parameters

mu_t    <- 0      # mean of true latent trait parameters
mu_c    <- 3      # mean of cheating latent trait parameters
sigma_t <- 1      # standard dev. of true latent trait parameters
sigma_c <- 1.25   # standard dev. of cheating latent trait parameters
corr    <- 0.8    # covariance between cheating and true latent trait parameters

1/(1+exp(-mu_t)) # probability of correct for average item without preknowledge
1/(1+exp(-mu_c)) # probability of correct for average item with preknowledge

# some rough idea about effect size, item preknowledge effect
# Odds ratio
exp(mu_c)/exp(mu_t)


th <- mvrnorm(N,
              mu = c(mu_t,mu_c),
              Sigma = matrix(c(sigma_t,corr,corr,sigma_c),2,2))

theta_t <- th[,1] 
theta_c <- th[,2]

describe(theta_t)
describe(theta_c)
describe(theta_c - theta_t)
cor(theta_t,theta_c)

# Generate observed responses

r <- matrix(nrow=N,ncol=n)

for(j in 1:N){
  for(i in 1:n){
    
    p_t <- exp(theta_t[j] - b[i])/(1+exp(theta_t[j] - b[i]))
    p_c <- exp(theta_c[j] - b[i])/(1+exp(theta_c[j] - b[i]))
    
    if(H[j] == 1 & C[i] == 1){
      r[j,i] = rbinom(1,1,p_c)
    } else {
      r[j,i] = rbinom(1,1,p_t)
    }
    
  }
}

# Convert it to data frame

r           <- as.data.frame(r)
colnames(r) <- paste0('R',1:n)

################################################################################
#         Combine Response Time and Response Accuracy Data
################################################################################

d       <- cbind(rt,r)
d$group <- H
d$id    <- 1:nrow(d)

d_long <- reshape(
  data      = d,
  varying   = list(RT = paste0("RT", 1:n), R = paste0("R", 1:n)),
  v.names   = c("RT", "R"),
  timevar   = "item",
  times     = 1:n,
  idvar     = "id",
  direction = "long"
)

# Add item status

d_long$compromised <- NA

for(j in 1:n){
  d_long[d_long$item==j,]$compromised = C[j]
}


describeBy(d_long$RT,list(d_long$group,d_long$compromised),mat=TRUE)
describeBy(d_long$R,list(d_long$group,d_long$compromised),mat=TRUE)

################################################################################

write.csv(d_long,'./data/simdata_long.csv',
          row.names = FALSE)  

write.csv(d,'./data/simdata_wide.csv',
          row.names = FALSE)  
