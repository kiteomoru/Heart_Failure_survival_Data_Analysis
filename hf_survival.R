# Logistic Regression an Surnival analysis of Heartfailure survival data

#Relative Path
setwd("/Users/kiteomoru/Desktop/Epidata/KAGGLEDATA")

#analysis of medical health insurance data (linear regression)
library(readxl)
library(ggplot2)
library(reshape)
library(psych)
library(car)
library(tidyverse)
library(skimr)
library(janitor)
library(ggpubr)
library(gtsummary)
library(stringr)
library(purrr)
library(tidyr)
library(pacman)
library(survival)
library(survminer)
library(parameters)
theme_set(theme_pubr())
##forestplot
pacman::p_load(easystats)



#load data
#data obtained from Kaggle
HFS_data= read_excel("/Users/kiteomoru/Desktop/Epidata/KAGGLEDATA/survivalheart_failure_clinical_records_dataset.xlsx")
data= as.data.frame(HFS_data)
head(data)

#data exploration
skim(data)

cor_data_m=cor(data)
#visualize

cor_data_m_melt= melt(cor_data_m)
cor_data_m_melt$value= round(cor_data_m_melt$value, digits = 2)
heatmap= ggplot(cor_data_m_melt, aes(X1, X2, fill= value))+
  geom_tile(color = "white", lwd= 1.5, linetype= 1) +
  scale_fill_gradient2(low = "#2C7BB6",
                       mid = "white",
                       high = "#D7191C", breaks=seq(-1,1,0.2), limits= c(-1,1)) +
  coord_fixed() + 
  theme_minimal(base_family="Helvetica")+
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20))

heatmap2=heatmap + 
  geom_text(aes(X1, X2, label = value), color = "black", size = 5) +
  theme(axis.ticks=element_blank())+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
heatmap2

#convert 1 and 0 to 'M' and 'F'
data$sex[data$sex  == 0]<- 'F'
data$sex[data$sex  == 1]<- 'M'


#Convert variables to factors
str(data$diabetes)
str(data$sex)
str(data$high_blood_pressure)
str(data$smoking)
str(data$DEATH_EVENT)

data$diabetes= as.factor(data$diabetes)
data$sex= as.factor(data$sex)
data$smoking = as.factor(data$smoking)
data$high_blood_pressure= as.factor(data$high_blood_pressure)
data$anaemia= as.factor(data$anaemia)
data$DEATH_EVENT= as.factor(data$DEATH_EVENT)


#Statistics
#scatterplot to find relationship among variables
pairs.panels(data)
#visualize Normality distribution
par(mfrow = c(2, 4))
qqPlot(data$age)
qqPlot(data$creatinine_phosphokinase)
qqPlot(data$ejection_fraction)
qqPlot(data$platelets)
qqPlot(data$serum_creatinine)
qqPlot(data$serum_sodium)
qqPlot(data$time)

#Shapiro-Wilk
shapiro.test(data$age)
shapiro.test(data$creatinine_phosphokinase)
shapiro.test(data$ejection_fraction)
shapiro.test(data$platelets)
shapiro.test(data$serum_creatinine)
shapiro.test(data$serum_sodium)
shapiro.test(data$time)


#Wilcoxons test
## testing for significant differences between continuous variables and outcome outcome variables  with a wilcox test
wilcox.test(age ~ DEATH_EVENT, data = data)
wilcox.test(creatinine_phosphokinase ~ DEATH_EVENT, data = data)
wilcox.test(ejection_fraction ~ DEATH_EVENT, data = data)
wilcox.test(platelets ~ DEATH_EVENT, data = data)
wilcox.test(serum_creatinine ~ DEATH_EVENT, data = data)
wilcox.test(serum_sodium ~ DEATH_EVENT, data = data)
wilcox.test(time ~ DEATH_EVENT, data = data)


#Chi-Square test
## testing for significant differences between categorical groups.
chisq.test(data$anaemia, data$DEATH_EVENT)
chisq.test(data$diabetes, data$DEATH_EVENT)
chisq.test(data$high_blood_pressure, data$DEATH_EVENT)
chisq.test(data$sex, data$DEATH_EVENT)
chisq.test(data$smoking, data$DEATH_EVENT)


#Convert age to categories - age_group
data= data %>% 
  mutate(
    # Create categories
    age_group = dplyr::case_when(
     
      age > 39 & age <= 60 ~ "40-60",
      age > 60 & age <= 80 ~ "61-80",
      age > 80            ~ "> 80"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("40-60","61-80", "> 80")
    )
  )
# Death event and Age group
model2 <- glm(DEATH_EVENT ~ age_group, family = "binomial", data = data)
summary(model2)


#Looping multiple univariate models
## define variables of interest 
explanatory_vars1<- c("age_group", "serum_creatinine", "anaemia", "creatinine_phosphokinase", "diabetes","ejection_fraction", 
                      "high_blood_pressure", "platelets", "serum_creatinine", "serum_sodium", "sex", "smoking", "time")
explanatory_vars1 %>% str_c("DEATH_EVENT ~ ", .)

#Multivariate Logistic Regression - glm()
# Death event and all other variables
## run a regression with all variables of interest 
mv_reg <- explanatory_vars1 %>%  ## begin with vector of explanatory column names
  str_c(collapse = "+") %>%     ## combine all names of the variables of interest separated by a plus
  str_c("DEATH_EVENT ~ ", .) %>%    ## combine the names of variables of interest with outcome in formula style
  glm(family = "binomial",      ## define type of glm as logistic,
      data = data)          ## define your dataset

## choose a model using forward selection based on AIC
## you can also do "backward" or "both" by adjusting the direction
final_mv_reg <- mv_reg %>%
  step(direction = "forward", trace = FALSE)

mv_tab_base <- final_mv_reg %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%  ## get a tidy dataframe of estimates 
  mutate(across(where(is.numeric), round, digits = 2))          ## round 

## show results table of final regression -HTML Format
mv_table <- tbl_regression(final_mv_reg, exponentiate = TRUE)
mv_table

## remove the intercept term from your multivariable results- forest plot
final_mv_reg %>% 
  model_parameters(exponentiate = TRUE) %>% 
  plot()


#####Survival plot

data_surv= data
summary(data_surv$time)
# cross tabulate the new event var and the outcome var from which it was created
data_surv$DEATH_EVENT= as.numeric(data_surv$DEATH_EVENT)
data_surv$DEATH_EVENT[data_surv$DEATH_EVENT  == 2]<- 'Alive'
data_surv$DEATH_EVENT[data_surv$DEATH_EVENT  == 1]<- 'Dead'

data_surv= data_surv%>%
  dplyr::mutate(
  # create the event var which is 1 if the patient died and 0 if he was right censored
  outcome = ifelse(is.na(DEATH_EVENT) | DEATH_EVENT == "Alive", 0, 1))
  
data_surv %>% 
  tabyl(DEATH_EVENT, outcome)

summary(data_surv$age_group)

data_surv= data_surv%>%
  select(anaemia, creatinine_phosphokinase, diabetes, ejection_fraction, high_blood_pressure,platelets, serum_creatinine,
         serum_sodium,   serum_sodium, sex, smoking,time,  DEATH_EVENT, age_group, outcome )

##Descriptive table
data_surv %>% 
  tabyl(sex, age_group, show_na = F) %>% 
  adorn_totals(where = "both") %>% 
  adorn_percentages() %>% 
  adorn_pct_formatting() %>% 
  adorn_ns(position = "front")

#Survival Object
# create the survfit object based on gender
data_surv_fit_sex <-  survfit(Surv(time, outcome) ~ sex, data = data_surv)

# set colors
sex <- c("lightgreen", "darkgreen")

# create plot
plot(
  data_surv_fit_sex,
  col = sex,
  xlab = "Days of follow-up",
  ylab = "Survival Probability")

# add legend
legend(
  "topright",
  legend = c("F","M"),
  col = sex,
  lty = 1,
  cex = .9,
  bty = "n")

#Plot
survminer::ggsurvplot(
  data_surv_fit_sex, 
  data = data_surv,          # again specify the data used to fit linelistsurv_fit_sex 
  conf.int = FALSE,              # do not show confidence interval of KM estimates
  surv.scale = "percent",        # present probabilities in the y axis in %
  break.time.by = 10,            # present the time axis with an increment of 10 days
  xlab = "Follow-up days",
  ylab = "Survival Probability",
  pval = T,                      # print p-value of Log-rank test 
  pval.coord = c(40,.91),        # print p-value at these plot coordinates
  risk.table = T,                # print the risk table at bottom 
  legend.title = "Gender",       # legend characteristics
  legend.labs = c("Female","Male"),
  font.legend = 10, 
  palette = "Dark2",             # specify color palette 
  surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
  ggtheme = theme_light()        # simplify plot background
)


# SURVIVAL OBJECT
data_surv_fit_BP <-  survfit(
  Surv(time, outcome) ~ high_blood_pressure,
  data = data_surv
)

# HYPERTENSION plot
ggsurvplot( 
  data_surv_fit_BP,
  data = data_surv,
  size = 1, linetype = "strata",   # line types
  conf.int = T,
  surv.scale = "percent",  
  break.time.by = 10, 
  xlab = "Follow-up days",
  ylab= "Survival Probability",
  pval = T,
  pval.coord = c(40,.91),
  risk.table = T,
  legend.title = "HYPERTENSION",
  legend.labs = c("HIGH BP", "NORMAL BP"),
  font.legend = 10,
  palette = c("#E7B800","#3E606F"),
  surv.median.line = "hv", 
  ggtheme = theme_light()
)

#EF to binary variable
data_surv$ejection_fraction_cat= ifelse(data_surv$ejection_fraction >= 30, 1, 0)
data_surv_fit_EF <-  survfit(
  Surv(time, outcome) ~ ejection_fraction_cat,
  data = data_surv
)

# EJECTION FRACTION plot
ggsurvplot( 
  data_surv_fit_EF,
  data = data_surv,
  size = 1, linetype = "strata",   # line types
  conf.int = T,
  surv.scale = "percent",  
  break.time.by = 10, 
  xlab = "Follow-up days",
  ylab= "Survival Probability",
  pval = T,
  pval.coord = c(40,.91),
  risk.table = T,
  legend.title = "EJECTION FRACTION",
  legend.labs = c("NORMAL", "LOW"),
  font.legend = 10,
  palette = c("#E7B800","#3E606F"),
  surv.median.line = "hv", 
  ggtheme = theme_light()
)





##Survival- Cox Model
#fit the model
datasurv_cox <-  coxph(
  Surv(time, outcome) ~ sex + age_group+ ejection_fraction + serum_creatinine+high_blood_pressure,
  data = data_surv
)


#test the proportional hazard model
datasurv_ph_test <- cox.zph(datasurv_cox)
datasurv_ph_test

#print the summary of the model
summary(datasurv_cox)
survminer::ggcoxzph(datasurv_ph_test)


ggforest(datasurv_cox, data = data_surv)










