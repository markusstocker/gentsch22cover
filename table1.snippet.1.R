library(lme4)
library(lmerTest)
df.MWD  <- read.csv("https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/df.MWD.csv")
df.MWD$cc_variant <- factor(df.MWD$cc_variant, levels = c("Fallow", "Mustard","Clover", "Oat" , "Phacelia", "Mix4" , "Mix12"))
lm.mwd.1 <- lmer(MWD_cor ~ cc_variant + (1|depth), data = df.MWD)
df1 <- data.frame(summary(lm.mwd.1)$coefficients, check.names=FALSE)
df1