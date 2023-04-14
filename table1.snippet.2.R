library(lme4)
library(lmerTest)
df.MWD  <- read.csv("https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/df.MWD.csv")
df.MWD$cc_variant <- factor(df.MWD$cc_variant, levels = c("Fallow", "Mustard","Clover", "Oat" , "Phacelia", "Mix4" , "Mix12"))
lm.mwd.2 <- lmer(MWD_cor ~ cc_type + (1|depth), data = df.MWD)
df2 <- data.frame(summary(lm.mwd.2)$coefficients, check.names=FALSE)
df2