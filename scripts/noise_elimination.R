## Noise elimination
library(fitdistrplus)
library(BioNet)
library(ggplot2)

# load data created in "dif_analysis.R"
data_lym <- tT_lym

# BUM parameters
a = 0.15
lambda = 0.73


data_lym_noise <- data.frame(matrix(NA, nrow = 0, ncol = 2))
data_lym_signal <- data.frame(matrix(NA, nrow = 0, ncol = 2))
for (i in 1:nrow(data_lym)){
  h = lambda + (1-lambda)*a*data_lym$adj.P.Val[i]^(a-1)
  p_i = (h-lambda)/h
  x = runif(1,0,1)
  if (x < p_i){
    data_lym_signal <- rbind(data_lym_signal, data_lym_2[i, ])
  }
  else {
    data_lym_noise <- rbind(data_lym_noise, data_lym_2[i, ])
  }
}

# check BUM distribution
fb_noise <- fitBumModel(data_lym_signal$adj.P.Val, plot = FALSE)
hist(fb_noise)
plot(fb_noise)

# data fit to normal distribution
fit_signal <- fitdist(data_lym_signal$logFC, "norm")
fit_noise <- fitdist(data_lym_noise$logFC, "norm")
plot(fit_signal)
plot(fit_noise)

# creating histograms
ggplot(data=data_lym_signal, aes(logFC))+
  geom_density(alpha=1)+
  geom_histogram(color="black", fill="darkblue", alpha=0.5, bins=20)

ggplot(data=data_lym_noise, aes(logFC))+
  geom_density(alpha=1)+
  geom_histogram(color="black", fill="darkgreen", alpha=0.5, bins=20)
