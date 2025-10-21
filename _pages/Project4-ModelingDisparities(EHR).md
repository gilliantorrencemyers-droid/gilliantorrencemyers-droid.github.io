---
title: "Modeling Disparities"
permalink: /modeling-disparities/
---

Download it here!

- [Modeling Disparities (EHR).R](/Project_4_Genomics.R)

```r
## Modeling disparities (EHR)


#======================= Kept from the demo ====================================
  rm(list=ls())
set.seed(13)

#Performance package is missing in tidyr image of PACE.
install.packages("performance")

library("data.table")
library("dplyr")
library("parallel")
library('stringr')
library('ggplot2')
library('performance')

############################################################
###################        Step 1        ###################
###################  Data exploration    ################### 
############################################################


# Load the files for our cohort.
cohort = as.data.frame(fread("C:/Users/gmyer/Downloads/BIOS 4150 GENOMICS & APPL BIOINFORMATICS/ParticipantCohort_4.tsv"))
dim(cohort)
head(cohort)

a1c = as.data.frame(fread("C:/Users/gmyer/Downloads/BIOS 4150 GENOMICS & APPL BIOINFORMATICS/ParticipantA1C.tsv"))
dim(a1c)
head(a1c)

t2d = as.data.frame(fread("C:/Users/gmyer/Downloads/BIOS 4150 GENOMICS & APPL BIOINFORMATICS/ParticipantT2DStatus.tsv"))
dim(t2d)
head(t2d)

############################################################
#####################      Step 2      ##################### 
#####################  Data cleaning   ##################### 
############################################################

# Merge the datasets together.
cohort = merge(cohort, a1c, by = "ParticipantID")
cohort = merge(cohort, t2d, by = "ParticipantID")

# Make SIRE as a categorical variable.
cohort$SIRE = factor(cohort$SIRE, level = c("White", "Black", "Hispanic"))
cohort$T2DStatus = factor(cohort$T2DStatus, level = c(0, 1))

# Check the size and dimensions of the cohort.
dim(cohort)
head(cohort)

# Visualzie the density of age distribution.
ggplot(cohort, aes(x = A1C, fill = SIRE)) + 
  geom_density(alpha = 0.3) + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "bottom")

#=================================================================================



###################################### QUESTION 1 ##########################
# Summarize A1C by SIRE groups.
cohort %>% count(SIRE)
cohort %>%
  group_by(SIRE) %>%
  summarize(
    min = min(A1C),
    max = max(A1C),
    mean = mean(A1C),
    median = median(A1C)
  )

# Summarize T2D by SIRE Groups, needed to convert T2DStatus to a numeric because it was a factor
cohort %>% count(SIRE)
cohort %>%
  group_by(SIRE) %>%
  summarize(
    min = min(as.numeric(as.character(T2DStatus))),
    max = max(as.numeric(as.character(T2DStatus))),
    mean = mean(as.numeric(as.character(T2DStatus))),
    median = median(as.numeric(as.character(T2DStatus)))
  )


#Visualization using box plot for A1C levels
##Needed to create from cohort data frame to also include SIR, so need to use column name A1C
ggplot(cohort, aes(x = SIRE, y = A1C, fill = SIRE)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 1, size = 2, color = "pink")+
    geom_jitter(width = 0.1, alpha = 0.1) +
    labs(title = "A1C levels sorted by SIRE",
          x = "SIRE",
          y = "A1C") +
    theme_minimal() +
    theme(legend.position = "none")+
    coord_flip()
  

#Visualization using Histogram for T2DStatus
#Prevalence of T2D summary table
Prevalence_Summary <- cohort %>%
  group_by(SIRE) %>%
  summarise(
    total_T2D = n(),
    with_T2D = sum(T2DStatus == "1"),
    prevalence = total_T2D/with_T2D
  )

print(Prevalence_Summary)

#bar graph (histogram)
ggplot(Prevalence_Summary, aes(x = SIRE, y = prevalence)) +
       geom_bar(stat = "identity", fill = "lavender") +
       labs(title = "Prevalence by SIRE",
        x = "SIRE",
        y = "Prevalence") +
  theme_minimal()


############################# QUESTION 2 #####################################

# Create a linear model to test for A1C disparities within the outlier-removed dataset.
a1c_mod <- lm(A1C ~ SIRE + Sex + Age, data = cohort)

summary_stats = summary(a1c_mod)$coefficients %>% as.data.frame()
summary(a1c_mod)


############################ QUESTION 3 #####################################
# Create a logistic model to test for T2D disparities.
t2d_mod <- glm(T2DStatus ~ SIRE + Sex + Age, data = cohort, family = "binomial")

summary_stats = summary(t2d_mod)$coefficients %>% as.data.frame()
summary(t2d_mod)



############################ QUESTION 4 #####################################
#Confidence interval for linear regression
confint(a1c_mod)

#all results are statistically significant at 95% confidence interval
#Being black, Hispanic, Male, and aging by a year all increase the likely hood of developing high A1C

#                  2.5 %    97.5 %
#(Intercept)  4.84922338 5.0989804
#SIREBlack    1.04994013 1.1665131
#SIREHispanic 0.87325649 0.9771556
#SexM         0.31052081 0.3953059
#Age          0.01976363 0.0250038


# Get the confidence intervals for the coefficients in logistic regression.
conf_intervals <- confint(t2d_mod) %>% 
  as.data.frame() %>%
  setNames(c("Low", "High")) %>%
  tibble::rownames_to_column("Variable")

summary_stats <- summary_stats %>%
  tibble::rownames_to_column("Variable") %>%
  left_join(conf_intervals, by = "Variable") %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  filter(Variable != "(Intercept)")



#logistic regression odds ratio
exp(confint(t2d_mod))

#odds ratio for t2D
#                  2.5 %     97.5 %
#(Intercept)  0.02162036 0.02952353
#SIREBlack    2.14141408 2.43713693
#SIREHispanic 2.50271961 2.80475481
#SexM         1.13055914 1.25113658
#Age          1.02996258 1.03647394

#Males are more likely to have t2d, as well as black and hispanic people compared to others.

################################### QUESTION 5 ########################################
# Using the forester package.
# Look for collinearity.
check_collinearity(t2d_mod)
#5.1
#Low correlation
#The VIF is very close to 1, ensuring that there is no multicollinearity in the data. 
#This tells us that the independent variables, SIRE, age, and sex are not highly
##correlated with eachother and show that there is some kind of effect between 
##these variables on T2D

#5.2
#A logistic regression is not like a linear regression and does not assume normality/ normally distributed residuals
#I decided to use deviance residuals instead, using a similar method used in  the a1c 5.2 section
dr <- residuals(t2d_mod, type ="deviance")
hist(dr, breaks = 20, main = "Deviance Residuals of T2D", xlab = "Deviance Residuals")

#There are a lot of values around -0.5 and 1.6, but about 90% of deviance residuals are around 0.5
#To look further into this, I looked at the proportion of 0 and 1 to see if there was a large discrepancy in how many 0 and 1 we had
# There is about 83% 0 and 16% 1
table(cohort$T2DStatus)
prop.table(table(cohort$T2DStatus))

#    0     1 
#41891  8109 

# the logistic regression is doing well on 0 which shows in the graph, most people do not have t2d,so it is skewed
#It could be underpredicting the 1, and so the deviance residual analysis here shows that there is a bias in predicting no t2d
#realistically because this is a logistica regression, we should see non constant varianace, however I think due to the extreme skew of 0 and 1 it misrepresents what should be seen

#5.3
#I chose to do boxplots to show if the model fits the group well and to see residuals
y <- cohort$T2DStatus
boxplot(residuals(t2d_mod, type = "deviance") ~ y)

#0 again fits well, its closer to 0 than 1.

#I also chose to do them by SIRE to see if there are any clear patterns
boxplot(dr ~ SIRE, data = cohort)
#This looks interesting, the center of the boxplot is about zero for white and black, but the hispanic has a high box
#The outliers are about the same for white and black, but again hispanic shows a potential misfit to the model
#The spread is similar in white and black, so there may be homosceadasticity however not between all SIRE groups
#It is important to not that for a logistic regression, heomosceadasticity cannot be determined because the variance of residuals is not constant due to the binary outcome

#Overall, it is difficult to determine 5.2 and 5.3 from the nature of the logistical regression,
#However, we can make certain assumptions based off of what we know about the data set (unequal amount of 0 and 1 so likely that 5.2 and 5.3 are inaccurate)
#Additionally, we can say that because this is a logistic regression, heterosceadasticity should be used because of the lack of constant variance
#The lack of homosceadasticity for this does not confirm this and neither does 5.2 lack of constant variance, so mroe tests should be run in the future

boxplot(dr ~ Sex, data = cohort)
#The box plot for sex looks god, the spread and size and height of the boxes are very similar



check_collinearity(a1c_mod)
#5.1
#Low correlation
#The VIF is very close to 1, ensuring that there is no multicollinearity in the data. 
#This tells us that the independent variables, SIRE, age, and sex are not highly
##correlated with eachother and show that there is some kind of effect between 
##these variables on A1C levels

#5.2
#Statistical test for large data sets regarding normality of residuals can apparently be misleading
#due to small deviations skewing the interpretation. I instead opted for a visual assessment, using
#something called a Q-Q plot (Quantile-Quantile plot) to show if a dataset follows
# a specific distribution.

qqnorm(residuals(a1c_mod))
qqline(residuals(a1c_mod), col = "forestgreen")

#it follows the line pretty well, there are a few outliers at the end (~6), but 
#follows the normal distribution well. 

#5.3
#again, graphs are probably more reliable for this large data, it is easier to visualize the differences and compare
boxplot(residuals(a1c_mod) ~ SIRE, data = cohort)

boxplot(residuals(a1c_mod) ~ Sex, data = cohort)

#Levene test to compare variances across categorical data, can be unreliable with larger data sets
install.packages("car")
library(car)
leveneTest(residuals(a1c_mod) ~ SIRE, data = cohort)

#As I suspected, the levene test was unreliable at this scale. It reported an extremley low p-value(0.01961)
#suggesting unequal variance among groups. However, upon further inspection using
#a boxplot, the variances look equal. The size and height of the boxes are about the same, 
#the only concern is the outliers,the amount present across all groups is probably the same, but the spread
#varies across the SIRE groups
#overall, it looks like there is equal residuals across a1c which validates the linear regression results.


#==================================== AI Usage =================================

#I would say ~4-5/10 AI usage was used. There were parts I wasnt sure on for interpreting the data and I wanted to double check to make sure.

```
