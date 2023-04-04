library(lme4)
library(emmeans)
library(multcomp)
library(tidyverse)
df.MWD <- read.csv("https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/df.MWD.csv", check.names=FALSE)
lm.mwd <- lmer(MWD_cor ~ cc_type + (1|depth), df.MWD)
df.pw.MWD.tot <- cld(emmeans(lm.mwd, list(pairwise ~ cc_type)), Letters=letters, sort=FALSE)
df.MWD$depth2 <- gsub(" cm", "", df.MWD$depth)
(p2 <- ggplot(df.MWD, aes(x=as.numeric(Soil_depth), y = MWD_cor, color=cc_type))+
    geom_smooth(method = "auto", se=T, alpha = 0.2)+
    scale_x_reverse(breaks=as.numeric(df.MWD$Soil_depth), labels = df.MWD$depth2)+
    geom_text(data = df.pw.MWD.tot, aes(y=emmean+SE, x=3, label=.group), size=4, color="black", nudge_x = -0.05)+
    coord_flip() +
    labs(x="Soil depth (cm)", y= expression("MWD (mm)"), color="")+
    theme(legend.position = "bottom")
)