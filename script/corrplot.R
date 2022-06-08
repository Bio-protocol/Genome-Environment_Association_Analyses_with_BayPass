## Plot population correlation from BayPass covariance matrix
## run as: Rscript corrplot.R
## Kaichi Huang 2022 May

#install.packages("tidyverse")
#install.packages("corrplot")
library(tidyverse)
library(corrplot)
select=dplyr::select

my.cov <- as.matrix(read.table("cache/GSD_RAD_core_mat_omega.out"))
my.pop <- read.table("cache/GSD_RAD.baypass.pop")
dimnames(my.cov)=list(my.pop$V1,my.pop$V1)
pop_order <- as.character(sort(as.numeric(my.pop$V1)))
my.cov <- my.cov[pop_order, pop_order]
my.corr <- cov2cor(my.cov)
png("output/GSD_RAD_core_mat_omega.png", 750, 750)
corrplot(my.corr, mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
dev.off()
