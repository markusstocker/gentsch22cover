library(lme4)
library(dplyr)
library(lmerTest)
df.20 <- read.csv("https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/CATCHY_aggregate_stability_2020_block2.csv", check.names=FALSE)
for(i in c(3:9, 11:18)) {
  df.20[,i] <- as.factor(df.20[,i])
}
df.20$cc_variant <- recode(df.20$cc_variant, "1" ="Fallow", "2"="Mustard", "3"="Clover", 
                           "4" = "Oat", "5" = "Phacelia","6" ="Mix4", "7"="Mix12")
df.20$cc_type <- ifelse(df.20$cc_variant %in% c("Mix4", "Mix12"), "Mix", "Single")
df.20$cc_type <- as.factor(ifelse(df.20$cc_variant == "Fallow", "Fallow", df.20$cc_type))
df.MWD <- subset(df.20, Fraction=="bulk")
lm.mwd.1 <- lmer(MWD_cor ~ cc_variant + (1|depth), data = df.MWD)
lm.mwd.2 <- lmer(MWD_cor ~ cc_type + (1|depth), data = df.MWD)
df1 <- data.frame(summary(lm.mwd.1)$coefficients, check.names=FALSE)
df2 <- data.frame(summary(lm.mwd.2)$coefficients, check.names=FALSE)