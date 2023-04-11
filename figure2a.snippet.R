library(tidyverse)
library(plyr)
library(pairwiseCI)
library(multcompView)
pairwiseLetters <- function (x) {
  df <- data.frame(x$byout)$p.value
  names(df) <- row.names(data.frame(x$byout))
  comp <- multcompLetters(df,compare = "<",threshold = 0.05,Letters = c(letters, LETTERS, "."),reversed = FALSE)
  data.frame(.group = comp$Letters,
             cc_variant = names(comp$Letters))
}
df.MWD <- read.csv("https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/df.MWD.csv", check.names=FALSE)
df.MWD$cc_variant <- factor(df.MWD$cc_variant, levels = c("Fallow", "Mustard","Clover", "Oat" , "Phacelia", "Mix4" , "Mix12"))
pw.MWD <- function(x) {
  pairwiseTest(MWD_cor ~ cc_variant, data = x)
}
pw.list.MWD <- dlply(df.MWD, .(depth), pw.MWD)
(df.pw.MWD <- ldply(pw.list.MWD, .fun=pairwiseLetters))
(df.pw.MWD.pvalues <- ldply(pw.list.MWD, .fun=function (x) {
  data.frame(p.value = x$byout[[1]]["p.value"],
             compnames = x$byout[[1]]["compnames"])
}))
(p1 <- ggplot(df.MWD, aes(x=cc_variant, y=MWD_cor, fill=cc_variant))+
    geom_jitter(shape=21, size=3.5, width = 0.2, alpha=0.3)+
    geom_text(data = df.pw.MWD, aes(x=cc_variant, y=2, label=.group), size=4, vjust = 0.2 )+
    stat_summary(fun.data = "mean_se", geom = "errorbar",width = 0.1, color="black")+
    stat_summary(fun.data = "mean_se",  geom = "point", size=4, shape=21)+
    facet_grid (depth~.,scales = "free")+
    labs(x="", y="MWD (mm)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)