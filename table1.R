# written with R version 4.1.3 (2022-03-10) -- "One Push-Up"
# by Norman Gentsch - gentsch@ifbk.uni-hannover.de
# Institute of Soil Science
# Leibniz University Hannover

s <- suppressPackageStartupMessages

s(library(dplyr))
s(library(lme4))
s(library(orkg))

orkg <- ORKG(host="https://sandbox.orkg.org/")
orkg$templates$materialize_template(template_id = "R200727")
tp = orkg$templates$list_templates()

df.20 <- read.csv("CATCHY_aggregate_stability_2020_block2.csv", check.names=FALSE)

for(i in c(3:9, 11:18)) {
  df.20[,i] <- as.factor(df.20[,i])
}

df.20$cc_variant <- recode(df.20$cc_variant, "1" ="Fallow", "2"="Mustard", "3"="Clover", 
                             "4" = "Oat", "5" = "Phacelia","6" ="Mix4", "7"="Mix12")

df.20$cc_type <- ifelse(df.20$cc_variant %in% c("Mix4", "Mix12"), "Mix", "Single")
df.20$cc_type <- as.factor(ifelse(df.20$cc_variant == "Fallow", "Fallow", df.20$cc_type))

# Input data
df.MWD <- subset(df.20, Fraction=="bulk")

# Two Linear Mixed Model (LMM) computations
lm.mwd.1 <- lmer(MWD_cor ~ cc_variant + (1|depth), data = df.MWD)
lm.mwd.2 <- lmer(MWD_cor ~ cc_type + (1|depth), data = df.MWD)

# Output data for the two LMM
df1 <- data.frame(summary(lm.mwd.1)$coefficients, check.names=FALSE)
df2 <- data.frame(summary(lm.mwd.2)$coefficients, check.names=FALSE)

instance <- tp$model_fitting(
  label="Linear mixed model fitting with MWD as response, CC variant as predictor variable, and soil depth as random variable", 
  has_input_dataset="https://github.com/markusstocker/gentsch22cover/blob/main/df.MWD.csv",
  has_input_model=tp$statistical_model(
    label="A linear mixed model with MWD as response and CC variant as predictor variable",
    is_denoted_by=tp$formula(
      label="The formula of the linear mixed model with MWD as response and CC variant as predictor variable",
      has_value_specification=tp$value_specification(
        label="MWD_cor ~ cc_variant + (1|depth)",
        has_specified_value="MWD_cor ~ cc_variant + (1|depth)"
      )
    )
  ),
  has_output_dataset=tuple(df1, "Results of LMM with MWD as response and CC variant as predictor variable")
)
instance$serialize_to_file("article.contribution.1.json", format="json-ld")

instance <- tp$model_fitting(
  label="Linear mixed model fitting with MWD as response, CC type as predictor variable, and soil depth as random variable", 
  has_input_dataset="https://github.com/markusstocker/gentsch22cover/blob/main/df.MWD.csv",
  has_input_model=tp$statistical_model(
    label="A linear mixed model with MWD as response and CC type as predictor variable",
    is_denoted_by=tp$formula(
      label="The formula of the linear mixed model with MWD as response and CC type as predictor variable",
      has_value_specification=tp$value_specification(
        label="MWD_cor ~ cc_type + (1|depth)",
        has_specified_value="MWD_cor ~ cc_type + (1|depth)"
      )
    )
  ),
  has_output_dataset=tuple(df2, "Results of LMM with MWD as response and CC type as predictor variable")
)
instance$serialize_to_file("article.contribution.2.json", format="json-ld")