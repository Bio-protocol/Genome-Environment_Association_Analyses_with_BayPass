## Produce POD samples
## run as: Rscript POD.R
## Kaichi Huang 2022 May

#install.packages("mvtnorm")
source("program/baypass_2.3/utils/baypass_utils.R")

nsnp <- 1000 # number of SNPs

# Read in the covariance matrix
omega=as.matrix(read.table("cache/GSD_RAD_core_mat_omega.out"))
# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef=read.table("cache/GSD_RAD_core_summary_beta_params.out",h=T)$Mean
# Original data to obtain total allele count
gt.data<-geno2YN("cache/GSD_RAD.baypass.txt")
# Create the POD
pod.data <- simulate.baypass(omega.mat=omega, nsnp=nsnp, sample.size=gt.data$NN, 
                             beta.pi=pi.beta.coef, pi.maf=0, suffix="GSD_RAD_POD.baypass.txt")

file.rename(from="G.GSD_RAD_POD.baypass.txt", to="cache/GSD_RAD_POD.baypass.txt")
