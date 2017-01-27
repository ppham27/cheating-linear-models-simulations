library(data.table)
library(xtable)
library(ggplot2)
library(scales)

## read in all data

## independent runs
independent.simulation.data <- rbind(data.table(read.csv("50_independent_beta.tsv", sep = "\t")),
                                     data.table(read.csv("100_independent_beta.tsv", sep = "\t")),
                                     data.table(read.csv("200_independent_beta.tsv", sep = "\t")),
                                     data.table(read.csv("400_independent_beta.tsv", sep = "\t")),
                                     data.table(read.csv("800_independent_beta.tsv", sep = "\t")))

## dependent runs
dependent.simulation.data <- rbind(data.table(read.csv("50_01_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("50_03_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("50_05_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("50_07_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("100_01_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("100_03_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("100_05_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("100_07_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("200_01_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("200_03_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("200_05_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("200_07_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("400_01_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("400_03_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("400_05_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("400_07_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("800_01_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("800_03_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("800_05_R_beta.tsv", sep = "\t")),
                                   data.table(read.csv("800_07_R_beta.tsv", sep = "\t")))

## merge data together
independent.simulation.data[,power:=0.05]
simulation.data <- rbind(dependent.simulation.data, independent.simulation.data)
setkey(simulation.data, X2N, power)

## p value distributions
p.value.distribution <- simulation.data[,
                                        data.frame(table(cut(p.value, breaks=seq(0,1,by=0.05)))/length(p.value)), 
                                        by=c("X2N", "power")]
setnames(p.value.distribution, "Var1", "Bin")
pdf("pValueDistribution.pdf", width=16, height=10)
ggplot(p.value.distribution,
       aes(x=Bin, y=Freq)) + 
  geom_bar(stat="identity") + 
    facet_grid(power ~ X2N, scales="free") + 
      scale_y_continuous(label=percent) + 
        coord_cartesian() +  
          theme_grey(8) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
              ggtitle("p-value Distribution") 
dev.off()

## replace -1 with NA
simulation.data[,set.size := as.numeric(set.size)]
simulation.data[,subset.size := as.numeric(subset.size)]
simulation.data[set.size==-1, set.size := Inf]
simulation.data[subset.size==-1, subset.size := Inf]

## create variables one by one
summary.data <- simulation.data[,
                                list(runs=length(p.value)),
                                by=c("X2N","power")]

## percent where convergence is reached
summary.data <- summary.data[simulation.data[,
                                             list(convergence=gsub("%","\\\\%",percent(sum(set.size < Inf)/length(set.size))),
                                                  initial.significance.pct=sum(p.value <= 0.05)/length(p.value)),
                                             by=c("X2N","power")]]

## quantiles
summary.data <- summary.data[simulation.data[,
                                             list(`25%.set.size`=as.integer(sort(set.size)[ceiling(length(set.size)*0.25)]),
                                                  `80%.set.size`=as.integer(sort(set.size)[ceiling(length(set.size)*0.8)])),
                                             by=c("X2N","power")]]

## mean, average, sd, excluding 0s
summary.data <- summary.data[simulation.data[set.size > 0 & set.size < Inf,
                                             list(mean.subset.size=mean(subset.size),
                                                  median.subset.size=as.integer(median(subset.size)),
                                                  sd.subset.size=sd(subset.size),
                                                  mean.set.size=mean(set.size),
                                                  median.set.size=as.integer(median(set.size)),
                                                  max.set.size=as.integer(max(set.size)),
                                                  sd.set.size=sd(set.size),
                                                  lm.coef.subset.pct=coefficients(lm(subset.size ~ 0 + set.size))["set.size"],
                                                  pct.one=gsub("%","\\\\%",percent(sum(subset.size==1)/length(subset.size)))),
                                             by=c("X2N","power")]]

## power 80 stats
summary.data <- summary.data[simulation.data[summary.data[,c("X2N", "power", "80%.set.size"), with=FALSE]][set.size > 0 & set.size <= `80%.set.size`,
                                                                                                           list(mean.subset.size.80=mean(subset.size),
                                                                                                                sd.subset.size.80=sd(subset.size),
                                                                                                                lm.coef.subset.pct.80=coefficients(lm(subset.size ~ 0 + set.size))["set.size"],
                                                                                                                max.subset.size.80=as.integer(max(subset.size)),
                                                                                                                median.subset.size.80=median(subset.size)),
                                                                                                           by=c("X2N", "power")]]

## effect size stats
summary.data <- summary.data[simulation.data[set.size < Inf,
                                             list(effect.direction.pct=gsub("%","\\\\%",percent(sum(beta > 0)/length(beta))),
                                                  mean.effect.size=mean(beta[beta > 0])),
                                             by=c("X2N","power")]]

## simulation.data[summary.data[,c("X2N", "power", "80%.set.size"), with=FALSE]][set.size > 0 & set.size <= `80%.set.size`,
##                                                                                                          list(median.subset.size.80=median(subset.size),
##                                                                                                               sort(subset.size)[ceiling(length(subset.size)*0.95)],
##                                                                                                               quantile(subset.size)[5]),
##                                                                                                          by=c("X2N", "power")]

summary.data[,c("X2N","power","lm.coef.subset.pct", "lm.coef.subset.pct.80"), with=FALSE][power==0.05]

summary.data[,c("X2N","power","initial.significance.pct"), with=FALSE][power==0.05]


## prep for printing
labeler <- function(text) {
  if (text == "X2N") {
    "$2N$"
  } else if (text == "power") {
    "Power (1 - $\\beta$)"
  } else if (text == "convergence") {
    "Percent Significant"
  } else if (text == "mean.subset.size") {
    "Mean Subset Size"
  } else if (text == "median.subset.size") {
    "Median Subset Size"
  } else if (text == "mean.set.size") {
    "Mean Set Size"
  } else if (text == "max.set.size") {
    "Max Set Size"
  } else if (text == "sd.subset.size") {
    "SD Subset Size"
  } else if (text == "sd.set.size") {
    "SD Set Size"
  } else if (text == "median.set.size") {
    "Median Set Size"
  } else if (text == "80%.set.size") {
    "80\\% Set Size"
  } else if (text == "max.subset.size.80") {
    "80\\% Subset Size"
  } else if (text == "mean.subset.size.80") {
    "Mean Subset Size (80\\%)"
  } else if (text == "lm.coef.subset.pct.80") {
    "Subset Fraction (80\\%)"
  } else if (text == "pct.one") {
    "Subset Size 1 (\\%)"
  } else if (text == "mean.effect.size") {
    "Mean $\\hat{\\beta}_S$ when $\\hat{\\beta}_S > 0$"
  } else if (text == "effect.direction.pct") {
    "Percent $\\hat{\\beta}_S > 0$"
  } else {
    text
  } 
}

print(xtable(summary.data[power==0.05,c("X2N",
                             "convergence"), 
                          with=FALSE],
             caption="Statistical significance is almost always found before we have as many independent covariates as observations. A power of $0.05$ denotes that $\\mathbf{Y}$ and $\\mathbf{X}$ are independent.",
             align=c("r","r","r"),
             label="tab:convergence", 
             digits=3), 
      booktabs=TRUE, 
      include.rownames=FALSE,
      sanitize.text.function=function(columns) { sapply(columns, labeler) })


## general stats
print(xtable(summary.data[,c("X2N",
                             "power",
                             "mean.subset.size", 
                             "pct.one",
                             "mean.set.size", 
                             "sd.set.size"), with=FALSE],
             caption="As the power increases, we need less covariates, and the subset size decreases, too. Power $0.05$ corresponds to the case where $\\mathbf{Y}$ and $\\mathbf{X}$ are independent.",
             label="tab:power_stats",
             align=c("r","r","r","r","r","r","r"),
             digits=3), 
      booktabs=TRUE, 
      include.rownames=FALSE,
      sanitize.text.function=function(columns) { sapply(columns, labeler) })

print(xtable(summary.data[power==0.05,
                          c("X2N",
                            "mean.subset.size", 
                            "pct.one",
                            "mean.set.size", 
                            "sd.set.size"), with=FALSE],
             caption="As $N$ gets larger, we need bigger sets, but not that much bigger. Statistics were calculated after removing the cases of set size 0.",
             align=c("r","r","r","r","r","r"),
             label="tab:independent_stats",
             digits=3), 
      booktabs=TRUE, 
      include.rownames=FALSE,
      sanitize.text.function=function(columns) { sapply(columns, labeler) })

## 80 stats
print(xtable(summary.data[,c("X2N",
                             "power",
                             "80%.set.size", 
                             "max.subset.size.80",
                             "mean.subset.size.80"), 
                          with=FALSE],
             caption="If we only want an 80\\% power test, the set and subset sizes become much smaller.",
             label="tab:power_stats_80",
             digits=3), 
      booktabs=TRUE, 
      include.rownames=FALSE,
      sanitize.text.function=function(columns) { sapply(columns, labeler) })

## effect stats
print(xtable(summary.data[power > 0.05,
                         ,c("X2N",
                             "power",
                             "effect.direction.pct",
                             "mean.effect.size"), 
                          with=FALSE],
             caption="Even at low power, we find the correct effect despite data dredging. However, the size of the effect is overestimated at these low powers.",
             align=rep("r",times=5),
             label="tab:effect_size",
             digits=3), 
      booktabs=TRUE, 
      include.rownames=FALSE,
      sanitize.text.function=function(columns) { sapply(columns, labeler) })



summary.data[,c("X2N","power","median.subset.size.80"),with=FALSE]

write.csv(simulation.data, file="simulation_data.csv", row.names=FALSE)



simulation.data.factor <- data.table(simulation.data)
simulation.data.factor[,X2N:=factor(X2N)]
ggplot(simulation.data.factor[p.value <= 0.1 & p.value > 0.05],
       aes(x=subset.size)) + 
  geom_bar() + 
    facet_wrap( ~ X2N, scales="free")
