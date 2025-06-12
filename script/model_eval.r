require(cmdstanr)
require(ggplot2)
################################################################################

d_long <- read.csv('./data/simdata_long.csv')
d_wide <- read.csv('./data/simdata_wide.csv')

fit <- readRDS("./do_not_upload/model_fit.RDS")

# Estimation time
fit$time()
################################################################################
# Model Diagnostics

# Extract model summary

model_summary <- as.data.frame(fit$summary(variables = c("item",
                                                         "person",
                                                         "pC",
                                                         "pH")))

  # Add a column for parameter types extracted from row names
  
  model_summary$type <- gsub("\\[.*$", "", model_summary$variable)
  
  head(model_summary)

# Sampler diagnostics

sampler_params <- posterior::as_draws_df(fit$sampler_diagnostics())

# Effective Sample size

iter <- nrow(sampler_params)   # number of chains x number of sampling iterations

model_summary$n_eff_ratio <- model_summary[,'ess_bulk'] / iter

psych::describeBy(model_summary$n_eff_ratio,
                  model_summary$type,
                  mat=TRUE)[,c('group1','min')]


psych::describeBy(model_summary$ess_bulk,
                  model_summary$type,
                  mat=TRUE)[,c('group1','mean','min','max')]

ggplot(model_summary, aes(x = ess_bulk)) +
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~type, scales = "free_y",nrow=2) +
  labs(
    title = "Distribution of ESS Values by Parameter Type",
    x = "R-hat",
    y = "Frequency"
  ) +
  theme_minimal()



# Split Rhat

ggplot(model_summary, aes(x = rhat)) +
  geom_histogram(binwidth = 0.002, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~type, scales = "free_y",nrow=2) +
  labs(
    title = "Distribution of R-hat Values by Parameter Type",
    x = "R-hat",
    y = "Frequency"
  ) +
  theme_minimal()


psych::describeBy(model_summary$rhat,
                  model_summary$type,
                  mat=TRUE)[,c('group1','mean','min','max')]

psych::describeBy(model_summary$rhat<1.01,
                  model_summary$type,
                  mat=TRUE)[,c('group1','mean','min','max')]

psych::describeBy(model_summary$rhat<1.05,
                  model_summary$type,
                  mat=TRUE)[,c('group1','mean','min','max')]

# Tree depth

table(sampler_params$treedepth__)

# E-BMFI

e_bfmi <- c()
for (n in 1:4) {
  energies = sampler_params[sampler_params$.chain==n,]$energy__
  numer = sum(diff(energies)**2) / length(energies)
  denom = var(energies)
  e_bfmi[n] = numer / denom
  print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
}

# Divergences  

100 * mean(sampler_params$divergent__)


################################################################################

# Person Parameter Estimates

theta_t <- fit$summary(variables ='thetat')
theta_c <- fit$summary(variables ='thetac')

psych::describe(theta_t$mean)
psych::describe(theta_c$mean)

tau_t <- fit$summary(variables ='taut')
tau_c <- fit$summary(variables ='tauc')

psych::describe(tau_t$mean)
psych::describe(tau_c$mean)


# Item Parameter Estimates

ipar <- matrix(fit$summary(variables = c("item"))$mean,
               ncol=3,byrow=FALSE)

psych::describe(ipar)


# pC, probability of item compromise

# Retrive the true item compromise status
# from long format data

C_vec <- c()

for(kk in 1:30){
  C_vec[kk] = unique(d_long[d_long$item==kk,]$compromised)
}

pC <- fit$summary(variables = c("pC"))$mean

psych::describeBy(pC,C_vec,mat=TRUE)[,c('group1','n','mean','sd','min','max')]

plot(density(pC[C_vec==0]),xlim=c(0,1),main="",ylim = c(0,8))
points(density(pC[C_vec==1]),lty=2,type='l')

require(pROC)

auc(C_vec,pC)

roc_analysis <- roc(response = C_vec,
                    predictor = pC)

plot(1-roc_analysis$specificities,
     roc_analysis$sensitivities,
     xlim = c(0,1),ylim=c(0,1),
     xlab = 'False Positive Rate (1-Specificity)',
     ylab = 'True Positive Rate (Sensitivity)',
     type='l')


my_thresholds <- seq(from=0.5,to=0.8,by=0.01)

coords(roc_analysis, 
       my_thresholds, 
       input="threshold", 
       ret=c("threshold","specificity", "sensitivity"))


# pH, probability of item preknowledge


pH <- fit$summary(variables = c("pH"))$mean

psych::describeBy(pH,d_wide$group,mat=TRUE)[,c('group1','n','mean','sd','min','max')]


plot(density(pH[d_wide$group==0]),xlim=c(0,1),main="",ylim = c(0,8))
points(density(pH[d_wide$group==1]),lty=2,type='l')

auc(d_wide$group,pH)


roc_analysis <- roc(response = d_wide$group,
                    predictor = pH)

plot(1-roc_analysis$specificities,
     roc_analysis$sensitivities,
     xlim = c(0,1),ylim=c(0,1),
     xlab = 'False Positive Rate (1-Specificity)',
     ylab = 'True Positive Rate (Sensitivity)',
     type='l')

my_thresholds <- seq(from=0.5,to=0.8,by=0.01)

coords(roc_analysis, 
       my_thresholds, 
       input="threshold", 
       ret=c("threshold","specificity", "sensitivity"))

