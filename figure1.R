s <- suppressPackageStartupMessages

s(library(orkg))
s(library(tidyverse))

orkg <- ORKG(host="https://incubating.orkg.org/")
orkg$templates$materialize_template(template_id = "R450104")
tp = orkg$templates$list_templates()

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

COL1 <- c("Fallow" = "black", "Mustard" = "red3" , "Clover" = "OliveDrab", "Oat" = "Gold", 
          "Phacelia" ="SteelBlue", "Mix4" = "orchid3", "Mix12"= "orange4")

theme_set(theme_bw())
theme_myBW <- theme(axis.title.x = element_text(size = 10, color = "black"), 
                    axis.title.y = element_text(angle = 90, vjust = 1.5, size = 10, color = "black"),
                    axis.text.x = element_text(size = 10, color = "black"), 
                    axis.text.y = element_text(size = 10, color = "black"), 
                    axis.ticks =element_line(colour="black"),
                    strip.text.x = element_text(size = 10, color = "black"),
                    strip.background = element_blank(),
                    panel.border =element_rect(colour="black", fill=NA), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(size = 12, hjust=0.5),
                    legend.text = element_text(size = 10),
                    legend.text.align=0,
                    legend.title =  element_text(size = 10), 
                    legend.key = element_rect(colour="white", fill = "white"),
                    legend.key.size = unit(5, "mm"),
                    legend.background = element_blank(),
                    legend.position = "bottom")

df.20 <- read.csv("CATCHY_aggregate_stability_2020_block2.csv", check.names=FALSE)

for(i in c(3:9, 11:18)) {
  df.20[,i] <- as.factor(df.20[,i])
}

df.20$cc_variant <- recode(df.20$cc_variant, "1" ="Fallow", "2"="Mustard", "3"="Clover", 
                             "4" = "Oat", "5" = "Phacelia","6" ="Mix4", "7"="Mix12")

df.20$depth <- recode(df.20$depth, "0-10" = "0-10 cm", "20-30"="20-30 cm", "30-40"="30-40 cm")

df.20$cc_type <- ifelse(df.20$cc_variant %in% c("Mix4", "Mix12"), "Mix", "Single")
df.20$cc_type <- as.factor(ifelse(df.20$cc_variant == "Fallow", "Fallow", df.20$cc_type))

coord_radar <- function (theta = "x", start = 0, direction = 1, clip = "off") {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),clip = clip,
          is_linear = function(coord) TRUE)
}

df.radar <- subset(df.20[,c("cc_variant", "depth", "Fraction","Plot", "OC_Frac", "OC_Frac_pc")], Fraction!="bulk")

# subset of fallow, set fallow as 100%
#OC_Frac in mg/g
df.fallow <- df.radar %>% 
  filter(cc_variant == "Fallow" & OC_Frac > 0) %>%
  group_by(cc_variant, depth, Fraction) %>%
  summarise(OC_Frac = mean(OC_Frac),
            OC_Frac_pc = mean(OC_Frac_pc))


df.radar <- merge(df.radar, df.fallow, by =c( "depth", "Fraction"),suffixes = c("",".fal"))

df.radar$OC_ratio.pc <- (df.radar$OC_Frac/df.radar$OC_Frac.fal*100)
# Fallow as 100 %
df.radar$OC_ratio.fal <-  100+df.radar$OC_Frac_pc-df.radar$OC_Frac_pc.fal


#  summarize CC values per treatment
df.radar.m <- df.radar %>% 
  filter(cc_variant !="Fallow") %>%
  group_by(cc_variant, depth, Fraction)%>%
  summarise(OC_ratio.m = mean(OC_ratio.pc),
            OC_ratio.se = se(OC_ratio.pc),
            OC_Frac_pc = mean(OC_Frac_pc),
            OC_Frac_pc.fal = mean(OC_Frac_pc.fal),
            OC_ratio.fal = mean(OC_ratio.fal))

# df only for fallow as polygon in the plot
df.radar.F <- df.radar %>% 
  filter(cc_variant !="Fallow") %>%
  group_by(cc_variant, depth, Fraction)%>%
  summarise(OC_ratio.m = 100)


# set levels into correct order
df.radar.m$Fraction <- factor(df.radar.m$Fraction, levels=c("<1","2-1","4-2","8-4","16-8"))
levels(df.radar.m$Fraction)

ggplot(df.radar.m, aes(x = Fraction, y = OC_ratio.fal, group = cc_variant, color= cc_variant, fill= cc_variant)) +
  geom_polygon(data = df.radar.F, aes(x = Fraction, y = OC_ratio.m), color="black" , fill = NA, size = 0.8, alpha=0.8) +
  geom_polygon(size = 0.8, alpha= 0.5) +
  coord_radar(clip="off") +
  #scale_y_log10()+
  scale_color_manual(values = COL1)+
  scale_fill_manual(values = COL1)+
  labs(y = "OC (% of Falow)", color="", fill="")+
  facet_grid(depth~cc_variant)+
  theme_myBW+
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "gray", size = 0.5,linetype = 2),
        axis.text.x = element_text(color = "gray12", size = 6),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.border = element_blank()
  )

ggsave("Fig.1.png", width=160, height = 110,dpi = 500, units = "mm")

instance <- tp$descriptive_statistical_calculation(
  label="Descriptive statistical calculation for relative proportion of OC in different soil fractions in percentage of fallow level", 
  has_input_dataset="https://github.com/markusstocker/gentsch22cover/blob/main/df.20.csv",
  has_output_dataset=tuple(as.data.frame(df.radar.m), "Relative proportion of OC in different soil fractions in percentage of fallow level"),
  has_output_figure="https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/Fig.1.png",
  has_output_statement="The variability of OC distribution within aggregates was quite large between and within CC treatments."
)
instance$serialize_to_file("article.contribution.3.json", format="json-ld")