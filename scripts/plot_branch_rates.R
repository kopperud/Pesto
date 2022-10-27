##################################################################
# plot the branch-rate estimates as a function of the event rate #
##################################################################


file = paste0("output/BDS_SCM.log")
samples = read.table(file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)

# discard some burnin (20%)
burnin = 0.20
n_samples = nrow(samples)

# combine the mcmc output
lambda_output   = samples[-c(1:ceiling(n_samples * burnin)),grepl("avg_lambda", colnames(samples))]

# store the parameters
lambda_mean = colMeans(lambda_output)
branch_lambdas_scm = lambda_mean[-length(lambda_mean)]


file = paste0("output/BDS_DA.log")
file = paste0("output/BDS_DBRM.log")
samples = read.table(file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE)

# discard some burnin (20%)
#burnin = 0.20
#n_samples = nrow(samples)

# combine the mcmc output
#lambda_output   = samples[-c(1:ceiling(n_samples * burnin)),grepl("avg_lambda", colnames(samples))]
lambda_output   = samples[,grepl("avg_lambda", colnames(samples))]

# store the parameters
lambda_mean = colMeans(lambda_output)
branch_lambdas_da = lambda_mean[-length(lambda_mean)]

# make the plot
pdf(paste0("figures/branch_rates.pdf"), height=4, width=4)

pch = 4
cex = 0.5
f   = 1.5
m   = 4
range      = range(pretty(branch_lambdas_scm))

plot(branch_lambdas_scm, branch_lambdas_da, xlim=range, ylim=range, pch=pch, cex=cex * f, xaxt="n", yaxt="n", xlab=NA, ylab=NA)
abline(a=0, b=1, lty=2)
axis(1, lwd.tick=1, lwd=0)
#mtext(side=3, text=paste0("k = "),    line=1.2)
mtext(side=1, text="branch-specific speciation rate", line=2.5, cex=0.7)


dev.off()
