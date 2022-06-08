## Plot the results of BayPass
## run as: Rscript BayPass_plot.R
## Kaichi Huang 2022 May

library(tidyverse)

# Read in datasets
my.data <- read.table("cache/GSD_RAD_standard_summary_betai_reg.out", h=T) # BayPass output of all SNPs
my.pod <- read.table("cache/GSD_RAD_POD_standard_summary_betai_reg.out", h=T) # BayPass output of the POD data
my.snp <- read.table("cache/GSD_RAD.baypass.snp", h=F) # Locations of SNPs
names(my.snp) <- c("CHR", "BP")
my.snp$CHR <- factor(as.numeric(my.snp$CHR))
my.cov <- read.table("cache/GSD_env.baypass.cov", h=F, sep="\t") # Covariables

# Convert the Positions for plotting
my.snp$BP <- my.snp$BP/1e6
POS <- my.snp$BP[my.snp$CHR == 1]
mid <- 0+max(POS)/2
for (i in 2:length(levels(my.snp$CHR))) {
  max.now <- max(POS)
  mid <- c(mid, max.now+max(my.snp$BP[my.snp$CHR == i])/2)
  POS <- c(POS, my.snp$BP[my.snp$CHR == i] + max.now)
}
my.snp$POS <- POS
my.snp$MRK <- 1:nrow(my.snp)

for (i in unique(my.data$COVARIABLE)) {
  my.plot <- my.data %>% 
    filter(COVARIABLE == i) %>% inner_join(my.snp) %>% select(CHR, POS, BF.dB.)
  my.q <- my.pod %>% 
    filter(COVARIABLE == i) %>% summarize(q99=quantile(BF.dB.,probs=0.99)) %>% pull() # Significance threshold
  p <- ggplot() + theme_classic() +
    geom_point(data=my.plot, mapping=aes(x=POS, y=BF.dB., col=CHR), pch=20, size=1, na.rm=T, show.legend=FALSE) +
    geom_hline(aes(yintercept=my.q), col="red", alpha=0.7, linetype="dotted") +
    xlab("Chromosome") + ylab(expression(BF[is]*" "*(dB))) +
    scale_colour_manual(values=rep(c("chocolate1","grey70"),9)) +
    scale_x_continuous(breaks=mid, labels=levels(my.plot$CHR), expand=c(0.02,0.02)) +
    scale_y_continuous(expand=c(0,0), limits=c(-12,35)) +
    theme(axis.line.x=element_blank()) +
    ggtitle(my.cov$V1[i])
  pdf(paste("output/BayPass_plot",my.cov$V1[i],"pdf",sep="."), width=8, height=3)
  print(p)
  dev.off()
}
