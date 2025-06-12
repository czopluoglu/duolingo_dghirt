require(cmdstanr)
require(ggplot2)
################################################################################

# Import data (long format)

d_long <- read.csv('./data/simdata_long.csv')
d_wide <- read.csv('./data/simdata_wide.csv')

################################################################################
# Fit DG-HIRT with a single chain
# The goal of this is to obtain initial start values
# We will provide the values obtained from single chain as starting values
# when re-fitting model again later with multiple chains

# Input Data

data_resp <- list(
  I              = length(unique(d_long$item)),
  J              = length(unique(d_long$id)),
  n_obs          = nrow(d_long),
  p_loc          = d_long$id,
  i_loc          = d_long$item,
  RT             = log(d_long$RT),
  Y              = d_long$R
)

# Compile the model syntax

mod <- cmdstan_model('./script/dghirt.stan')

# Fit the model using a single chain with a small number of iterations

fit_init <- mod$sample(
  data            = data_resp,
  seed            = 1234,
  chains          = 1,
  iter_warmup     = 750,
  iter_sampling   = 250,
  refresh         = 10,
  adapt_delta     = 0.99)

fit_init$time()
################################################################################
model_par <- fit_init$summary()


theta <- matrix(c(model_par[grep("^person\\[.*,1\\]$", model_par$variable),]$mean,
                  model_par[grep("^person\\[.*,2\\]$", model_par$variable),]$mean),
                ncol=2,byrow=FALSE)

delta <- matrix(c(model_par[grep("^delta\\[.*,1\\]$", model_par$variable),]$mean,
                  model_par[grep("^delta\\[.*,2\\]$", model_par$variable),]$mean),
                ncol=2,byrow=FALSE)

item <- matrix(c(model_par[grep("^item\\[.*,1\\]$", model_par$variable),]$mean,
                 model_par[grep("^item\\[.*,2\\]$", model_par$variable),]$mean,
                 model_par[grep("^item\\[.*,3\\]$", model_par$variable),]$mean),
               ncol=3,byrow=FALSE)

# A vector of initial P(H=1) parameters, probability of an examinee having item preknowledge
H <- as.vector(model_par[grep("pH", model_par$variable),]$mean)

# A vector of initial P(C=1) parameters, probability of an item being compromised
C <- as.vector(model_par[grep("pC", model_par$variable),]$mean)

# Put the initial estimates together as a list

start <- list(item   = item,
              person = theta,
              delta  = delta,
              pH     = H,
              pC     = C)



fit <- mod$sample(
  data            = data_resp,
  seed            = 1234,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 500,
  refresh         = 10,
  init            = list(start,start,start,start),
  adapt_delta     = 0.99)


# Save the model object with all parameters for future use

fit$save_object(file = "./do_not_upload/model_fit.RDS")

