s <- suppressPackageStartupMessages

s(library(orkg))
s(library(tidyverse))
s(library(lme4))
s(library(emmeans))
s(library(lmerTest))
s(library(plyr))
s(library(pairwiseCI))
s(library(multcomp))
s(library(multcompView))
s(library(ggpubr))

orkg <- ORKG(host="https://incubating.orkg.org/")
orkg$templates$materialize_template(template_id = "R450109")
tp = orkg$templates$list_templates()

pairwiseLetters <- function (x) {
  df <- data.frame(x$byout)$p.value
  names(df) <- row.names(data.frame(x$byout))
  comp <- multcompLetters(df,compare = "<",threshold = 0.05,Letters = c(letters, LETTERS, "."),reversed = FALSE)
  data.frame(.group = comp$Letters,
             cc_variant = names(comp$Letters))
}

COL <- c("Fallow" = "slategray", "Mustard" = "red3" , "Clover" = "OliveDrab", "Oat" = "Gold", 
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

df.MWD <- subset(df.20, Fraction=="bulk")

# pairwise comparison of CC_variants
pw.MWD <- function(x) {
  # Computes t-tests
  pairwiseTest(MWD_cor ~ cc_variant, data = x)
}

# has_specified_input
df.MWD <- df.MWD[,c("depth","cc_variant","MWD_cor")] 

pw.list.MWD <- dlply(df.MWD, .(depth), pw.MWD)

# apply function on list and produce data frame for plotting
# Small letters denoting significant difference between CC treatments by pairwise comparison 
(df.pw.MWD <- ldply(pw.list.MWD, .fun=pairwiseLetters))

# has_specified_output
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
    scale_fill_manual(values = COL, guide="none")+
    labs(x="", y="MWD (mm)")+
    theme_myBW+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

ggsave("Fig.2a.png", plot = p1, width = 160, height = 100, units = "mm")

instance <- tp$pairwise_test(
  label="Pairwise t-test with MWD response and CC variant predictor", 
  method="http://purl.obolibrary.org/obo/OBI_0000739",
  has_input_dataset=tuple(df.MWD, "Difference of mean weight diameter between the dry and wet sieving method"),
  has_output_dataset=tuple(df.pw.MWD.pvalues, "Pairwise t-test p-values for CC variants at three soil depths"),
  has_output_figure="https://raw.githubusercontent.com/markusstocker/gentsch22cover/main/Fig.2a.png",
  has_input_model=tp$statistical_model(
    label="A pairwise t-test with MWD response and CC variant predictor",
    is_denoted_by=tp$formula(
      label="The formula of the pairwise t-test with MWD response and CC variant predictor",
      has_value_specification=tp$value_specification(
        label="MWD_cor ~ cc_variant",
        has_specified_value="MWD_cor ~ cc_variant"
      )
    )
  ),
)
instance$serialize_to_file("article.contribution.4.json", format="json-ld")