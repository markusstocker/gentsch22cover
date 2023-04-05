library(tidyverse)
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
coord_radar <- function (theta = "x", start = 0, direction = 1, clip = "off") {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),clip = clip,
          is_linear = function(coord) TRUE)
}
df.radar <- read.csv("df.radar.csv")
df.radar$Fraction <- factor(df.radar$Fraction, levels=c("<1","2-1","4-2","8-4","16-8", "bulk"))
df.radar.m <- df.radar %>% 
  filter(cc_variant !="Fallow") %>%
  group_by(cc_variant, depth, Fraction)%>%
  summarise(OC_ratio.m = mean(OC_ratio.pc),
            OC_ratio.se = se(OC_ratio.pc),
            OC_Frac_pc = mean(OC_Frac_pc),
            OC_Frac_pc.fal = mean(OC_Frac_pc.fal),
            OC_ratio.fal = mean(OC_ratio.fal))
df.radar.F <- df.radar %>% 
  filter(cc_variant !="Fallow") %>%
  group_by(cc_variant, depth, Fraction)%>%
  summarise(OC_ratio.m = 100)
ggplot(df.radar.m, aes(x = Fraction, y = OC_ratio.fal, group = cc_variant, color= cc_variant, fill= cc_variant)) +
  geom_polygon(data = df.radar.F, aes(x = Fraction, y = OC_ratio.m), color="black" , fill = NA, size = 0.8, alpha=0.8) +
  geom_polygon(size = 0.8, alpha= 0.5) +
  coord_radar(clip="off") +
  labs(y = "OC (% of Falow)", color="", fill="")+
  facet_grid(depth~cc_variant)+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "gray", size = 0.5,linetype = 2),
        axis.text.x = element_text(color = "gray12", size = 6),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank()
  )