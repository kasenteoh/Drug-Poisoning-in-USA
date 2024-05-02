library(ggplot2)
library(dplyr)
library(tidyverse)
library(nlme)
library(forecast)
library(glmmTMB)
library(lme4)
library(lmtest)


######################################## Cleaning Data ######################################## 
# Import data
df <- read_csv("drug.csv")
df <- df %>%
  select(-c(`State`, `US Crude Rate`)) %>%
  rename(year = names(df)[2]) %>%
  rename(sex = names(df)[3]) %>%
  rename(age = names(df)[4]) %>%
  rename(race = names(df)[5]) %>%
  rename(rate = names(df)[8])

# Adding IDs
df <- transform(df, id = as.numeric(interaction(df$sex, df$age, df$race, drop = TRUE)))

######################################## EDA ######################################## 
################################### Sex ################################### 
# Group the data by year and sex, and calculate the sum of suicides_100k
df_sex <- df %>%
  group_by(year, sex) %>%
  summarise(mean_rate = mean(rate))

# Create the line plot
ggplot(df_sex, aes(x = year, y = mean_rate, color = sex)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", 
       y = "Death Rate per 100k", 
       color = "Sex") +
  ggtitle("Average Drug Poisoning Death Rate - Profile Plot by Sex") +
  theme_minimal()

################################### Age Group ################################### 
# Calculate the average suicide rate per 100k for each age group across each year
avg_by_age_year <- aggregate(rate ~ age + year, data = df, FUN = mean)

# Create the line plot
ggplot(avg_by_age_year, aes(x = year, y = rate, color = age, group = age)) +
  geom_point() + 
  geom_line() +
  labs(x = "Year", y = "Death Rate per 100k", color = "Age Group") +
  ggtitle("Average Drug Poisoning Death Rate - Profile Plot by Age") +
  theme_minimal()

################################### Race ################################### 
# Group the data by year and sex, and calculate the sum of suicides_100k
df_race <- df %>%
  group_by(year, race) %>%
  summarise(mean_rate = mean(rate))

# Create the line plot
ggplot(df_race, aes(x = year, y = mean_rate, color = race)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", 
       y = "Death Rate per 100k", 
       color = "Race") +
  ggtitle("Average Drug Poisoning Death Rate - Profile Plot by Race") +
  theme_minimal()

################################### Empirical Plot ################################### 
df_year <- df %>%
  group_by(year) %>%
  summarise(total_rate = sum(rate), 
            mean_rate = mean(rate), 
            sd_rate = sd(rate), 
            n = n(),
            se_rate = sd_rate / sqrt(n),
            lower_ci = mean_rate - 1.96*se_rate,
            upper_ci = mean_rate + 1.96*se_rate)

ggplot(df_year, aes(x = year, y = mean_rate)) +
  geom_point(size = 3, shape = 21, fill = "white") +
  geom_line() + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  labs(title = "Empirical Summary Plot of Average
       Death Rate with 95% Confidence Interval", 
       x = "Year", 
       y = "Death Rate per 100k") +
  theme_minimal()

######################################## Time Trend ######################################## 
# year
# Random Intercept / Compound Symmetry
RI_l <- lmer(rate ~ year + (1 | id), data = df)

# Random Intercept and Slope
RIAS_l <- lmer(rate ~ year + (1 + year | id), data = df)

# SP(pow)(time)
CAR_l <- lme(rate ~ year, random = ~ 1 | id,
            correlation = corCAR1(form = ~ year | id), data = df)

# year^2
# Random Intercept / Compound Symmetry
RI_q <- lmer(rate ~ I(year ^ 2) + (1 | id), data = df)

# Random Intercept and Slope
RIAS_q <- lmer(rate ~ I(year ^ 2) + (1 + year | id), data = df)

# SP(pow)(time)
CAR_q <- lme(rate ~ I(year ^ 2), random = ~ 1 | id,
           correlation = corCAR1(form = ~ year | id), data = df)

# year^3
# Random Intercept / Compound Symmetry
RI_c <- lmer(rate ~ I(year ^ 3) + (1 | id), data = df)

# Random Intercept and Slope
RIAS_c <- lmer(rate ~ I(year ^ 3) + (1 + year | id), data = df)

# SP(pow)(time)
CAR_c <- lme(rate ~ I(year ^ 3), random = ~ 1 | id,
           correlation = corCAR1(form = ~ year | id), data = df)


AIC(RI_l, RIAS_l, CAR_l, RI_q, RIAS_q, CAR_q, RI_c, RIAS_c, CAR_c)
BIC(RI_l, RIAS_l, CAR_l, RI_q, RIAS_q, CAR_q, RI_c, RIAS_c, CAR_c)

"year is the best time trend, linear time trend."

######################################## Sex Significance ######################################## 
RI_sex <- lme(rate ~ sex, random = ~ 1 | id,
              correlation = corCAR1(form = ~ year | id), data = df)
summary(RI_sex)

"Sex is a significant predictor"

######################################## Age Significance ######################################## 
RI_age <- lme(rate ~ age, random = ~ 1 | id,
              correlation = corCAR1(form = ~ year | id), data = df)
summary(RI_age)

"Age seems to be a significant predictor. However, there are a few groups that 
have a p-value that is greater than 0.05, namely age groups 0-14, 15-24, 65-74, 
and 75+. This result may be a consequence of not many death rates in these four 
age groups. Referring back to plot [age eda], the bottom four lines represent 
these four groups. Hence, some age groups are significant, while others are not."

######################################## Race Significance ######################################## 
RI_race <- lme(rate ~ race, random = ~ 1 | id,
              correlation = corCAR1(form = ~ year | id), data = df)
summary(RI_race)

"Race does not seem to be a significant predictor"

######################################## Sex + Age Significance ######################################## 
RI_sa <- lme(rate ~ sex + age, random = ~ 1 | id,
               correlation = corCAR1(form = ~ year | id), data = df)
summary(RI_sa)

anova(RI_age, RI_sex, RI_sa)
lrtest(RI_age, RI_sa)
lrtest(RI_sex, RI_sa)

"AIC and BIC of RI_sa is the least out of all three models. Using lrtest, 
the model with age and the model with sex both have a log-likelihood value less 
than that of the model with age and sex combined. Additionally, the p-value for 
the Chi-squared tests are minute in both log-ratiotests. Hence, the model with 
sex + age seems to be the best predictors. "

######################################## Sex + Age + Year Significance ########################################
RI_say <- lme(rate ~ sex + age + year, random = ~ 1 | id,
             correlation = corCAR1(form = ~ year | id), data = df)
summary(RI_say)

anova(RI_sa, RI_say)

"From above, the model with sex + age proved to be the model that is able to 
capture the most information. including the year variable as well. The AIC and BIC
with time is far less than that of with only sex and age. Additionally, the 
log-likelihood value is greater for the model including time with a very small 
p-value, all of which support the notion of sex+age+year as the best model. "



######################################## Best Covariance Model ######################################## 
#Using sex + age find best covariance model
# Random Intercept / Compound Symmetry
RI <- lmer(rate ~ sex + age + year + (1 | id), data = df)

# Random Slope
RS  <- lmer(rate ~ sex + age + year + (0 + year | id), data = df)

# Random Intercept and Slope
RIAS <- lmer(rate ~ sex + age + year + (1 + year | id), data = df)

# Independence model  
IN =  glmmTMB(rate ~ sex + age + year, dispformula = ~0, data = df)

# ARMA(1,1)
ARMA11 <- lme(rate ~ sex + age + year, random = ~ 1 | id,
              correlation = corARMA(form = ~ year | id, p = 1, q = 1), 
              data = df)

# CAR
CAR <- lme(rate ~ sex + age + year, random = ~ 1 | id,
            correlation = corCAR1(form = ~ year | id), data = df)


AIC(RI, RS, RIAS, IN, ARMA11, CAR)
BIC(RI, RS, RIAS, IN, ARMA11, CAR)
lrtest(ARMA11, CAR)

"ARMA11 has the lowest AIC and BIC. But CAR is very close. To determine whether
the small difference in AIC and BIC between the two models is negligible, a 
log-ratio test is performed. The log-likelihood score of the ARMA(1, 1) is 
greater than that of the CAR and the Chi-squared further supports this conclusion
with a p-value extremely close to 0. "
