### R code from vignette source 'Tutorial.Rnw'

###################################################
### code chunk number 1: Tutorial.Rnw:50-52
###################################################
options(width=60)
ps.options(family="sans")


###################################################
### code chunk number 2: Tutorial.Rnw:55-59
###################################################
library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)

my_data_2 = read.table("p_v_out_k_const.txt")
#my_data_2 =  table(dataLym)
my_data_2 <- as.numeric(my_data_2)
hist(my_data_2)
fb_2 <- fitBumModel(my_data_2, plot = TRUE)
fb_2
#f <- function(x) pbeta(x, fb_2$a, 1)
#ks_res <- ks.test(my_data_2, f)

dev.new(width = 13, height = 7)
par(mfrow=c(1,2), cex=2.5, lwd = 4)
hist(fb_2)

plot(fb_2)
dev.off()


my_f <- function(arg_1){
  fb_2$lambda*punif(arg_1)+(1-fb_2$lambda)*pbeta(arg_1, fb_2$a, 1)
}
for_test <- c(runif(4000*fb_2$lambda), rbeta((1-fb_2$lambda)*4000,fb_2$a,1))
print(fb_2$a)
print(fb_2$lambda)
ks.test(my_data_2, my_f)
ks.test(for_test, my_f)
hist(for_test)
