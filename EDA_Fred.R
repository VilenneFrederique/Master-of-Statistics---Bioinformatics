# Read libraries
library(haven)
library(tidyverse)
library(plotly)

# Read data
renal <- read_sas(data_file = "C:\\Users\\Gebruiker\\Documents\\School\\Master of Statistics - Bioinformatics (2021-2022)\\Longitudinal Data Analysis\\Projects\\Project 1\\renal_transformed.sas7bdat")

# Plot the individual profiles
renal %>% ggplot(aes(x = as.numeric(as.character(time)), y = hc, group = id)) + geom_line() + scale_x_continuous(breaks = seq(0, 10, 1))
# Plot is very busy due to the high amount of measurements
# Now we do see a few things, a lot of lines not being fully completed. 
# So there is missing data.
# There also seems a steep increase after measurement 0. This seems to be quite universal
# So variability between and a lot of variability within
# So things to account for:
# - Within subject variability
# - Missing data, feels a lot like the prostate dataset from Baltimore in class

# Plot in function of gender
renal %>% group_by(male, time) %>% mutate(mean_avg_gender = mean(hc, na.rm = TRUE)) %>% ggplot(aes(x = as.numeric(as.character(time)), y = mean_avg_gender, group = factor(male), col = factor(male))) + geom_line()
# Gender plays quite a big role. So should be in the model at all times

# Plot in function of cardio
renal %>% group_by(cardio, time) %>% mutate(mean_avg_cardio = mean(hc, na.rm = TRUE)) %>% ggplot(aes(x = as.numeric(as.character(time)), y = mean_avg_cardio, group = factor(cardio), col = factor(cardio))) + geom_line()
# Doesn't seem so relevant? Might be ignored later on. Both levels of cardio behave the same.

# Plot in function of reject
renal %>% group_by(reject, time) %>% mutate(mean_avg_reject = mean(hc, na.rm = TRUE)) %>% ggplot(aes(x = as.numeric(as.character(time)), y = mean_avg_reject, group = factor(reject), col = factor(reject))) + geom_line()
# Difference present between both levels

# Plot in function of age
renal %>% ggplot(aes(x = age, y = hc)) + geom_point()
plot_ly(x = as.numeric(as.character(renal$time)), y = renal$hc, z = renal$age, type = "scatter3d", alpha = 0.5)
# Not really as informative but fun to look at lmao
# Anyhow, there should always be a correlation between hc values and age.
# It's described in literature.

# Frequency tables for categorical variables
table(renal$male, renal$cardio)
table(renal$male, renal$reject)
table(renal$cardio, renal$reject)
# Not described yet, need to look back into the definitions

# Some modelling in R as I don't understand a shit of SAS 
library(lme4)
attach(renal)
model1 <- lmer(hc ~ factor(male) + factor(reject) + factor(cardio) + age + (1 | factor(id)))
               