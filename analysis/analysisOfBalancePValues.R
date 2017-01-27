library(data.table)
library(xtable)
library(ggplot2)
library(scales)
library(reshape2)

## compute exact distribution of balance p-values for N observations
getExact <- function(N) {
  bins <- 20
  prob <- numeric(bins)
  n <- N/2
  total.p <- 0
  for (K in 0:N) {
    p.binom <- dbinom(K, N, 1/2)        #number of successes
    for (k in 0:n) {                    #number of possible successes observed
      p.hypergeom <- dhyper(k, K, N - K, n)
      p <- K/N
      mu <- p*n
      delta <- abs(k - mu)
      if (mu + delta - 1 < mu - delta) {
        p.value <- 1
      } else {
        p.value <- phyper(mu - delta, K, N - K, n) +
          phyper(mu + delta - 1, K, N - K, n, lower.tail = FALSE)
      }
      idx <- ceiling(p.value*bins)
      prob[idx] <- prob[idx] + p.binom*p.hypergeom
    }
  }
  prob
}

## read in all data

## independent runs
independent.simulation.data <- rbind(
  data.table(read.csv("output_data/50_independent_beta_intercept.tsv", sep = "\t")),
  data.table(read.csv("output_data/100_independent_beta_intercept.tsv", sep = "\t")),
  data.table(read.csv("output_data/200_independent_beta_intercept.tsv", sep = "\t")),
  data.table(read.csv("output_data/400_independent_beta_intercept.tsv", sep = "\t")))
## independent.simulation.data.no.intercept <- rbind(
##   data.table(read.csv("output_data/50_independent_beta.tsv", sep = "\t")),
##   data.table(read.csv("output_data/100_independent_beta.tsv", sep = "\t")),
##   data.table(read.csv("output_data/200_independent_beta.tsv", sep = "\t")),
##   data.table(read.csv("output_data/400_independent_beta.tsv", sep = "\t")))

simulation.data <- data.table(independent.simulation.data)
## simulation.data.no.intercept <- data.table(independent.simulation.data.no.intercept)
setkey(simulation.data, X2N)
## setkey(simulation.data.no.intercept, X2N)

simulation.data[,                
                list(n.simulations=length(set.size),
                     subset.size=mean(subset.size),
                     set.size=mean(set.size),
                     pct.one=sum(subset.size==1)/length(subset.size)),
                by="X2N"]
## simulation.data.no.intercept[,                
##                              list(n.simulations=length(set.size),
##                                   subset.size=mean(subset.size),
##                                   set.size=mean(set.size),
##                                   pct.one=sum(subset.size==1)/length(subset.size)),
##                              by="X2N"]

## Only look at simulations, where the subset size was 1
simulation.data.one <- simulation.data[subset.size==1]
simulation.data.one[,                
                    list(n.simulations=length(set.size)),
                    by="X2N"]
ggplot(simulation.data.one, aes(balance.p.value)) + geom_histogram(binwidth=0.1)

p.value.distribution <- simulation.data.one[,
                                            data.frame(table(cut(balance.p.value, breaks=seq(0,1,by=0.05)))/length(balance.p.value)), 
                                            by=c("X2N")]

## compare the distribution of p-values to the exact distribution
p.value.distribution[order(X2N, Var1),
                     Prob:=c(getExact(50), getExact(100), getExact(200), getExact(400))]
setnames(p.value.distribution, "Var1", "Bin")

## produce a graph that shows that the single covariate is more likely to be unbalanced
## than expected
p.value.distribution.melted <- melt(p.value.distribution, id.vars=c("X2N","Bin"))
pdf("balancePValueDistributionIntercept.pdf", width=16, height=10)
ggplot(p.value.distribution.melted,
       aes(x=Bin, y=value, alpha=variable)) + 
  geom_bar(stat="identity", position="dodge") + 
    facet_grid(X2N ~ ., scales="free") + 
      scale_y_continuous(label=percent) + 
        scale_alpha_manual("Variable", values=c(1,0.4), labels=c("Observed Frequency", "Reference Probability")) + 
          coord_cartesian() +  
            theme_grey(8) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                ggtitle("p-value Distribution")                 
dev.off()

summarize.data <- function(simulation.data) {
  summary.data <- simulation.data[,
                                  list(runs=length(p.value)),
                                  by=c("X2N")]
  setkey(summary.data, X2N)
  summary.data <- summary.data[simulation.data[set.size > 0 & set.size < Inf,
                                               list(mean.subset.size=mean(subset.size),
                                                    sd.subset.size=sd(subset.size),
                                                    mean.set.size=mean(set.size),
                                                    sd.set.size=sd(set.size),
                                                    pct.one=sum(subset.size==1)/length(subset.size)),
                                               by=c("X2N")]]
  summary.data[order(X2N)][,c("X2N","runs","pct.one"),with=FALSE]
  summary.data
}
## print summary statistics for overall data, Table 1
summarize.data(simulation.data)
## print summary statistics when p-value is nearly statistically significant, Table 2
summarize.data(simulation.data[p.value > 0.05 & p.value <= 0.1])
## print summary statistics grouped by p.value group
simulation.data.small.p <- simulation.data[p.value > 0.05 & p.value <= 0.1]
simulation.data.small.p[,p.value.group:=cut(simulation.data.small.p$p.value, breaks=seq(0.05,0.1,by=0.01))]
summary.data <- simulation.data.small.p[,
                                        list(runs=length(p.value)),
                                        by=c("X2N","p.value.group")]
setkey(summary.data, X2N, p.value.group)
summary.data <- summary.data[simulation.data.small.p[set.size > 0 & set.size < Inf,
                                                     list(mean.subset.size=mean(subset.size),
                                                          sd.subset.size=sd(subset.size),
                                                          mean.set.size=mean(set.size),
                                                          sd.set.size=sd(set.size),
                                                          pct.one=sum(subset.size==1)/length(subset.size)),
                                                     by=c("X2N","p.value.group")]]
summary.data[order(X2N, p.value.group)][,c("X2N","p.value.group","runs","pct.one"),with=FALSE]

summary.data <- simulation.data.small.p[,
                                        list(runs=length(p.value)),
                                        by=c("X2N")]
setkey(summary.data, X2N)
summary.data <- summary.data[simulation.data.small.p[set.size > 0 & set.size < Inf,
                                                     list(mean.subset.size=mean(subset.size),
                                                          sd.subset.size=sd(subset.size),
                                                          mean.set.size=mean(set.size),
                                                          sd.set.size=sd(set.size),
                                                          pct.one=sum(subset.size==1)/length(subset.size)),
                                                     by=c("X2N")]]
summary.data

## write data out
## write.csv(simulation.data,
##           "balanced_p_values_intercept_10000.csv",
##           row.names=FALSE)
