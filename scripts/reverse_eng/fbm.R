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
#data(dataLym)
#data(interactome)

my_data_2 = read.table("/home/vladimir/Python_projects/term_project_IB/homework250318/p_v_out.txt")
my_data_2 <- as.numeric(my_data_2)
#hist(my_data_2)
fb_2 <- fitBumModel(my_data_2, plot = FALSE)
#fb_2
#f <- function(x) pbeta(x, fb_2$a, 1)
#ks_res <- ks.test(my_data_2, f)
file_out <- file("/home/vladimir/Python_projects/term_project_IB/homework250318/rout.txt")
write(c(fb_2$a, fb_2$lambda), file_out)
close(file_out)