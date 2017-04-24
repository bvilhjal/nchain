# Mantel test, introgressed genes
library(ggplot2)
setwd('C:/Users/MariaIzabel/Desktop/MASTER/PHD/nchain')
correl = read.csv('C:/Users/MariaIzabel/Desktop/MASTER/PHD/Methods/Intergenic_LD/corrected_snps/introgressed_gene_stats.csv')

hist(correl$r_scores, breaks = 50, main = 'Kinship vs all', xlab = 'Correlation scores', xlim = c(min(correl$r_scores), max(correl$r_scores)), col = 'grey')
mean(as.numeric(correl$r_scores))

nod_names = c('nodX', 'nodN', 'nodM', 'nodL', 'nodE', 'nodF', 'nodD', 'nodA', 'nodC', 'nodI', 'nodJ', 'nifB', 'nifA','fixX', 'fixC', 'fixB', 'fixA', 'nifH', 'nifD', 'nifK', 'nifE', 'rpoB', 'recA')

# Finding index
for(i in 1:length(nod_names)){
  a[i] = which(as.character(correl$gene) == nod_names[i])
}

nod_cors = correl[a,]
rug(correl[which(correl$origin == "['chrm']"), 5], col = 'blue', ticksize = 0.03, side = 1, lwd = 0.5)
rug(nod_cors$r_scores, ticksize = 0.03, side = 1, lwd = 0.5, col = 'red')
rug(correl[which(correl$origin != "['chrm']" & correl$origin != 0), 5], col = 'green', ticksize = 0.03, side = 1, lwd = 0.5)

abline(v = mean(as.numeric(correl$r_scores)), lty = 3, col = 'black')
abline(v = median(as.numeric(correl$r_scores)), lty = 3, col = 'blue')

write.csv(nod_cors, 'nod_genes_correlations.csv')

# Correlations between R score and number of snps
qplot(correl$r_scores, correl$snps, xlab = 'R scores')
cor(correl$r_scores, correl$snps)

# Correlations between R score and 
qplot(correl$r_scores, correl$members, xlab = 'R scores', ylab = 'Number of members')
cor(correl$r_scores, correl$members)
