#################### packages and functions ####################
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
s(library(orkg))

orkg <- ORKG(host="https://incubating.orkg.org/")
orkg$templates$materialize_template(template_id = "R450125")
tp = orkg$templates$list_templates()


# function to extract the pairwise letters from the list
pairwiseLetters <- function (x) {
  df <- data.frame(x$byout)$p.value
  names(df) <- row.names(data.frame(x$byout))
  comp <- multcompLetters(df,compare = "<",threshold = 0.05,Letters = c(letters, LETTERS, "."),reversed = FALSE)
  data.frame(.group = comp$Letters,
             cc_variant = names(comp$Letters))
}

# standard error calculation
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# set vector with colors and label
COL <- c("Fallow" = "slategray", "Mustard" = "red3" , "Clover" = "OliveDrab", "Oat" = "Gold", 
         "Phacelia" ="SteelBlue", "Mix4" = "orchid3", "Mix12"= "orange4")

COL1 <- c("Fallow" = "black", "Mustard" = "red3" , "Clover" = "OliveDrab", "Oat" = "Gold", 
          "Phacelia" ="SteelBlue", "Mix4" = "orchid3", "Mix12"= "orange4")

COL.type <- c("Fallow"="black", "Single"= "tomato", "Mix"="darkmagenta")

# customized ggplot theme 
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



###################### read in data ######################
# data from initial soil inventory 2015 block 2
df.15 <- read.csv("CATCHY_soil_data_2015_block2.csv")
# BD = Bulk density in g cm-1
# Clay, silt sand in %
# OC, TN in %

str(df.15)
# change attributes
for(i in c(1:12)) {
  df.15[,i] <- as.factor(df.15[,i])
}

#data from sampling for aggregate fractionation 2020 block 2 
df.20 <- read.csv("CATCHY_aggregate_stability_2020_block2.csv")

str(df.20)
names(df.20)

for(i in c(3:9, 11:18)) {
  df.20[,i] <- as.factor(df.20[,i])
}

# recode cover crop (cc) levels
df.20$cc_variant <- recode(df.20$cc_variant, "1" ="Fallow", "2"="Mustard", "3"="Clover", 
                           "4" = "Oat", "5" = "Phacelia","6" ="Mix4", "7"="Mix12")
levels(df.20$cc_variant)

df.20$depth <- recode(df.20$depth, "0-10" = "0-10 cm", "20-30"="20-30 cm", "30-40"="30-40 cm")

df.20$cc_fac <- ifelse(df.20$cc_variant =="Fallow", "Fallow", "Cover crop")

# produce a new label
df.20$cc_type <- ifelse(df.20$cc_variant %in% c("Mix4", "Mix12"), "Mix", "Single")
df.20$cc_type <- as.factor(ifelse(df.20$cc_variant == "Fallow", "Fallow", df.20$cc_type))
levels(df.20$cc_type)

# order Fractions
df.20$Fraction <- factor(df.20$Fraction, levels=c("<1","2-1","4-2","8-4","16-8", "bulk"))
levels(df.20$Fraction)




########################### mean weight diameter (MWD) in mm ############################
# The higher the MWD as more large scale aggregates are present after water treatment

str(df.20)

df.MWD <- subset(df.20, Fraction=="bulk")
str(df.MWD)
df.MWD$cc_type
levels(df.MWD$cc_type)

# the data represent approximately normally distribution
#ggplot(df.MWD, aes(x=MWD_cor))+
#  geom_histogram(aes(y =..density..))+
#  geom_density(col=2)


# pairwise comparison of CC_variants
pw.MWD <- function(x) {
  pairwiseTest(MWD_cor ~ cc_variant, data = x)
}

pw.list.MWD <- dlply(df.MWD, .(depth), pw.MWD)

# apply function on list and produce data frame for plotting
(df.pw.MWD <-  ldply(pw.list.MWD, .fun=pairwiseLetters))


# comparison of plot B with a mixed model
lm.mwd <- lmer(MWD_cor ~ cc_type + (1|depth), df.MWD)
lm.mwd

df.pw.MWD.tot <- cld(emmeans(lm.mwd, list(pairwise ~ cc_type)), Letters=letters, sort=FALSE)
df.pw.MWD.tot

# plot results
df.MWD$depth2 <- gsub(" cm", "", df.MWD$depth)

df.MWD[c('depth2', 'MWD_cor')]

(p2 <- ggplot(df.MWD, aes(x=as.numeric(Soil_depth), y = MWD_cor, color=cc_type))+
    geom_smooth(method = "auto", se=T, alpha = 0.2)+
    scale_x_reverse(breaks=as.numeric(df.MWD$Soil_depth), labels = df.MWD$depth2)+
    geom_text(data = df.pw.MWD.tot, aes(y=emmean+SE, x=3, label=.group), size=4, color="black", nudge_x = -0.05)+
    coord_flip() +
    scale_color_manual(values = COL.type)+
    labs(x="Soil depth (cm)", y= expression("MWD (mm)"), color="")+
    theme_myBW +
    theme(legend.position = "bottom")
)

ggsave("Fig.2b.svg", plot = p2, scale=0.5)
p2

outTab <- ggplot_build(p2)$data[[1]]
outTab


# overall effects of CC from a LMM (linear mixed effect model)
# Input: df.MWD
# Formula: MWD_cor ~ cc_type + (1|depth)
# Output Figure: Fig.2b.png
# Output Dataset: df.pw.MWD.tot

# Generate Empty Output Dataset
#class(df.pw.MWD.tot) <- "data.frame"

emptDf <- data.frame(matrix(nrow = 1, ncol = 1))

##################### ORKG Output ########################################

# Provisional template (model_fitting_2) to allow for the addition of an output figure.
instance <- tp$model_fitting_2(
  label="Overall effects of CC from a LMM (linear mixed effect model)", 
  
  # Links to Zenodo.org
  has_input_dataset="https://zenodo.org/record/7147566/files/CATCHY_aggregate_stability_2020_block2.csv",
  
  has_input_model=tp$statistical_model(
    label="A linear mixed model with MWD as response and CC type as predictor variable",
    is_denoted_by=tp$formula(
      label="The formula of the linear mixed model with MWD as response and CC type as predictor variable",
      
      # Provisional String for formula/model
      has_value_specification=tp$value_specification(
        label="MWD_cor ~ cc_type + (1|depth)",
        has_specified_value="MWD_cor ~ cc_type + (1|depth)"
      )
    )
  ),
  
  # Empty output dataframe (1 row and 1 column) to account for required field
  has_output_dataset= tuple(emptDf, 'NULL'),
  
  # Data URI with Base64 string (SVG)
  has_output_figure=base64enc::dataURI(file = "Fig.2b.svg", mime = "image/svg+xml"),
  
  has_output_statement= "A comprehensive data evaluation in LMMs indicated, the MWD increased with soil 
  depth and was significantly higher in CC treatments than the fallow."
  
)
instance$serialize_to_file("article.contribution.5.json", format="json-ld")
