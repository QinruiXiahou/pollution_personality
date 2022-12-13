#lead and atmosphere analyses
rm(list=ls())
#report things without scientific notation because it's confusing to my small brain
options(scipen = 999)

memory.limit(size=560000) 

library(nlme)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stargazer)
library(spdep)
library(RANN)
library(spatialreg)
library(tidyr)

##set working directory
setwd('L:/Public/Datasets1/Project-RetrospectiveAnalysis/replication')

##read in data
# lead_density <- read_csv("../../Cross Checks QX/lead_emissions_calc/69-81_PNAS/lead_density_adj_NA_1980.csv") %>%
lead_density <- read_csv("lead_density_NA_newadj.csv") %>%
  mutate(county = str_replace_all(county, pattern = "Newport News City", "Newport News"),
         county = str_replace_all(county, pattern = "Portsmouth City", "Portsmouth"),
         county = str_replace_all(county, pattern = "Richmond City", "Richmond"),
         county = ifelse(geoid==24510, "Baltimore", county)) %>%
  mutate(countystate = paste(county, ", ", state, sep="")) %>%
  select(geoid, Year, countystate, NA1972, isenNA, predNA, lead_den_cnty_km, lead_wden_adj_km) %>%
  mutate(lead_den_cnty_km = lead_den_cnty_km/1000, 
         lead_wden_adj_km = lead_wden_adj_km/1000) # convert the unit from mg/m2 to ug/m2 so that the magnitude is more comparable to lead concentration

usalife <- read.csv("US_b5_plus_lead_lifetimeexposure_2020_all269counties.csv")

##remove the outliers
OUTLIERS <- c("16079", "30049", "48141", "04007", "04021", "32033", "49035", "04003", "04019", "04011", "53053", "30023") #smelters
OUTLIER_NAMES <- unique(filter(lead_density, geoid %in% OUTLIERS)$countystate)

##generate new datasets for each regression
dim(usalife)

usalife_new <- usalife %>%
  pivot_longer(cols = 25:83, names_to = "lead_year", values_to = "lead_exp") %>%
  mutate(lead_year = as.numeric(str_sub(lead_year, 2, 5))) %>%
  mutate(flag_18 = ifelse(2018 - yearborn < 18 & lead_year - yearborn >= 0, 1, ifelse(lead_year - yearborn < 18 & lead_year - yearborn >= 0, 1, 0))) %>%
  group_by(X) %>%
  mutate(avg_leadexp = sum(lead_exp*flag_18, na.rm = TRUE)/18,
         log_avg_leadexp = log(avg_leadexp+1)) %>%
  ungroup() %>%
  select(-flag_18) %>%
  pivot_wider(names_from = lead_year, values_from = c(lead_exp)) %>%
  filter(!name %in% OUTLIER_NAMES) %>%
  mutate(cohort = case_when(yearborn >=1960 & yearborn<1965 ~ "Y6065",
                            yearborn >=1965 & yearborn<1970 ~ "Y6570",
                            yearborn >=1970 & yearborn<1975 ~ "Y7075",
                            yearborn >=1975 & yearborn<1980 ~ "Y7580",
                            yearborn >=1980 & yearborn<1985 ~ "Y8085",
                            yearborn >=1985 & yearborn<1990 ~ "Y8590",
                            yearborn >=1990 & yearborn<1995 ~ "Y9095",
                            yearborn >=1995 & yearborn<2000 ~ "Y9500",
                            yearborn >=2000 & yearborn<2005 ~ "Y0005")) %>% 
  mutate(name2 = name) %>% 
  separate(name2, into = c('county', 'state'), ", ") %>%
  filter(!state %in% c("Puerto Rico", "Hawaii", "Alaska"))


usalife_new_5year <- usalife %>%
  pivot_longer(cols = 25:83, names_to = "lead_year", values_to = "lead_exp") %>%
  mutate(lead_year = as.numeric(str_sub(lead_year, 2, 5))) %>%
  mutate(flag_5 = ifelse(1981 - yearborn < 5 & lead_year - yearborn >= 0, 1, ifelse(lead_year - yearborn < 5 & lead_year - yearborn >= 0, 1, 0))) %>%
  group_by(X) %>%
  mutate(avg_leadexp = sum(lead_exp*flag_5, na.rm = TRUE)/5,
         log_avg_leadexp = log(avg_leadexp+1)) %>%
  ungroup() %>%
  select(-flag_5) %>%
  pivot_wider(names_from = lead_year, values_from = c(lead_exp)) %>%
  filter(!name %in% OUTLIER_NAMES) %>%
  mutate(cohort = case_when(yearborn >=1960 & yearborn<1965 ~ "Y6065",
                            yearborn >=1965 & yearborn<1970 ~ "Y6570",
                            yearborn >=1970 & yearborn<1975 ~ "Y7075",
                            yearborn >=1975 & yearborn<1980 ~ "Y7580",
                            yearborn >=1980 & yearborn<1985 ~ "Y8085",
                            yearborn >=1985 & yearborn<1990 ~ "Y8590",
                            yearborn >=1990 & yearborn<1995 ~ "Y9095",
                            yearborn >=1995 & yearborn<2000 ~ "Y9500",
                            yearborn >=2000 & yearborn<2005 ~ "Y0005")) %>% 
  mutate(name2 = name) %>% 
  separate(name2, into = c('county', 'state'), ", ") %>%
  filter(!state %in% c("Puerto Rico", "Hawaii", "Alaska"))

usalife_new_6981_5year <- usalife %>%
  filter(yearborn >= 1969 & yearborn <= 1977) %>%
  pivot_longer(cols = 25:83, names_to = "lead_year", values_to = "lead_exp") %>%
  mutate(lead_year = as.numeric(str_sub(lead_year, 2, 5))) %>%
  filter(lead_year >= 1969 & lead_year <= 1981) %>%
  mutate(flag_5 = ifelse(1981 - yearborn < 5 & lead_year - yearborn >= 0, 1, ifelse(lead_year - yearborn < 5 & lead_year - yearborn >= 0, 1, 0))) %>%
  left_join(lead_density, by=c("lead_year"="Year", "name"="countystate")) %>%
  group_by(X) %>%
  mutate(avg_leadexp = sum(lead_exp*flag_5, na.rm = TRUE)/5,
         log_avg_leadexp = log(avg_leadexp+1),
         avg_leadden_unadj = sum(lead_den_cnty_km*flag_5, na.rm = TRUE)/5,
         log_avg_leadden_unadj = log(avg_leadden_unadj+1),
         avg_leadden_adj = sum(lead_wden_adj_km*flag_5, na.rm = TRUE)/5,
         log_avg_leadden_adj = log(avg_leadden_adj+1)) %>%
  ungroup() %>%
  rename(lead_den_unadj = lead_den_cnty_km, lead_den_adj = lead_wden_adj_km) %>%
  select(-flag_5) %>%
  pivot_wider(names_from = lead_year, values_from = c(lead_exp, lead_den_unadj, lead_den_adj)) %>%
  filter(!name %in% OUTLIER_NAMES) %>%
  mutate(cohort = case_when(yearborn >=1960 & yearborn<1965 ~ "Y6065",
                            yearborn >=1965 & yearborn<1970 ~ "Y6570",
                            yearborn >=1970 & yearborn<1975 ~ "Y7075",
                            yearborn >=1975 & yearborn<1980 ~ "Y7580",
                            yearborn >=1980 & yearborn<1985 ~ "Y8085",
                            yearborn >=1985 & yearborn<1990 ~ "Y8590",
                            yearborn >=1990 & yearborn<1995 ~ "Y9095",
                            yearborn >=1995 & yearborn<2000 ~ "Y9500",
                            yearborn >=2000 & yearborn<2005 ~ "Y0005")) %>% 
  mutate(name2 = name) %>% 
  separate(name2, into = c('county', 'state'), ", ") %>%
  filter(!state %in% c("Puerto Rico", "Hawaii", "Alaska"))


################################################################################
## 0. Lead concentration, 18 years
################################################################################

# e1 <- lm(se ~ sleadexp + age, data = usalife)
e0 <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new)
cl.e0 <- cluster.vcov(e0, usalife_new$name) # cluster-robust SEs
cl.robust.se.e0 <- sqrt(diag(cl.e0))

# a1 <- lm(sa ~ sleadexp + age, data = usalife)
a0 <- lm(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new)
cl.a0 <- cluster.vcov(a0, usalife_new$name) # cluster-robust SEs
cl.robust.se.a0 <- sqrt(diag(cl.a0))

# c1 <- lm(sc ~ sleadexp + age, data = usalife)
c0 <- lm(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new)
cl.c0 <- cluster.vcov(c0, usalife_new$name) # cluster-robust SEs
cl.robust.se.c0 <- sqrt(diag(cl.c0))

# n1 <- lm(sn ~ sleadexp + age, data = usalife)
n0 <- lm(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new)
cl.n0 <- cluster.vcov(n0, usalife_new$name) # cluster-robust SEs
cl.robust.se.n0 <- sqrt(diag(cl.n0))

# o1 <- lm(so ~ sleadexp + age, data = usalife)
o0 <- lm(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new)
cl.o0 <- cluster.vcov(o0, usalife_new$name) # cluster-robust SEs
cl.robust.se.o0 <- sqrt(diag(cl.o0))

################################################################################
## 1. Lead concentration, 5 years
################################################################################

# e1 <- lm(se ~ s_leadexp + age, data = usalife_new_5year)
e2 <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_5year)
cl.e2 <- cluster.vcov(e2, usalife_new_5year$name) # cluster-robust SEs
cl.robust.se.e2 <- sqrt(diag(cl.e2))

# a1 <- lm(sa ~ s_leadexp + age, data = usalife_new_5year)
a2 <- lm(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_5year)
cl.a2 <- cluster.vcov(a2, usalife_new_5year$name) # cluster-robust SEs
cl.robust.se.a2 <- sqrt(diag(cl.a2))

# c1 <- lm(sc ~ s_leadexp + age, data = usalife_new_5year)
c2 <- lm(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_5year)
cl.c2 <- cluster.vcov(c2, usalife_new_5year$name) # cluster-robust SEs
cl.robust.se.c2 <- sqrt(diag(cl.c2))

# n1 <- lm(sn ~ s_leadexp + age, data = usalife_new_5year)
n2 <- lm(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_5year)
cl.n2 <- cluster.vcov(n2, usalife_new_5year$name) # cluster-robust SEs
cl.robust.se.n2 <- sqrt(diag(cl.n2))

# o1 <- lm(so ~ s_leadexp + age, data = usalife_new_5year)
o2 <- lm(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_5year)
cl.o2 <- cluster.vcov(o2, usalife_new_5year$name) # cluster-robust SEs
cl.robust.se.o2 <- sqrt(diag(cl.o2))

################################################################################
## 2. Lead concentration, 5 years, cohort between 1969 and 1977
################################################################################

# e1 <- lm(se ~ s_leadexp + age, data = usalife_new_6981_5year)
e4 <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.e4 <- cluster.vcov(e4, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.e4 <- sqrt(diag(cl.e4))

# a3 <- lm(sa ~ s_leadexp + age, data = usalife_new_6981_5year)
a4 <- lm(sa  ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.a4 <- cluster.vcov(a4, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.a4 <- sqrt(diag(cl.a4))

# c3 <- lm(sc ~ s_leadexp + age, data = usalife_new_6981_5year)
c4 <- lm(sc  ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.c4 <- cluster.vcov(c4, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.c4 <- sqrt(diag(cl.c4))

# n1 <- lm(sn ~ s_leadexp + age, data = usalife_new_6981_5year)
n4 <- lm(sn  ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.n4 <- cluster.vcov(n4, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.n4 <- sqrt(diag(cl.n4))

# o1 <- lm(so ~ s_leadexp + age, data = usalife_new_6981_5year)
o4 <- lm(so  ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.o4 <- cluster.vcov(o4, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.o4 <- sqrt(diag(cl.o4))

################################################################################
## 3. Lead density, 5 years, cohort between 1969 and 1977
################################################################################

# e1 <- lm(se ~ s_leadden_unadj + age, data = usalife_new_6981_5year)
e6 <- lm(se ~ log_avg_leadden_unadj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.e6 <- cluster.vcov(e6, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.e6 <- sqrt(diag(cl.e6))

e6a <- lm(se ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.e6a <- cluster.vcov(e6a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.e6a <- sqrt(diag(cl.e6a))

# a5 <- lm(sa ~ s_leadden_unadj + age, data = usalife_new_6981_5year)
a6 <- lm(sa ~ log_avg_leadden_unadj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.a6 <- cluster.vcov(a6, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.a6 <- sqrt(diag(cl.a6))

a6a <- lm(sa ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.a6a <- cluster.vcov(a6a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.a6a <- sqrt(diag(cl.a6a))

# c5 <- lm(sc ~ s_leadden_unadj + age, data = usalife_new_6981_5year)
c6 <- lm(sc ~ log_avg_leadden_unadj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.c6 <- cluster.vcov(c6, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.c6 <- sqrt(diag(cl.c6))

c6a <- lm(sc ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.c6a <- cluster.vcov(c6a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.c6a <- sqrt(diag(cl.c6a))

# n1 <- lm(sn ~ s_leadden_unadj + age, data = usalife_new_6981_5year)
n6 <- lm(sn ~ log_avg_leadden_unadj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.n6 <- cluster.vcov(n6, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.n6 <- sqrt(diag(cl.n6))

n6a <- lm(sn ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.n6a <- cluster.vcov(n6a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.n6a <- sqrt(diag(cl.n6a))

# o1 <- lm(so ~ s_leadden_unadj + age, data = usalife_new_6981_5year)
o6 <- lm(so ~ log_avg_leadden_unadj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.o6 <- cluster.vcov(o6, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.o6 <- sqrt(diag(cl.o6))

o6a <- lm(so ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.o6a <- cluster.vcov(o6a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.o6a <- sqrt(diag(cl.o6a))

################################################################################
## Summary statistics
################################################################################
mean_0 <- as.character(round(mean(usalife_new$avg_leadexp, na.rm=TRUE), 3))
mean_2 <- as.character(round(mean(usalife_new_5year$avg_leadexp, na.rm=TRUE), 3))
mean_4 <- as.character(round(mean(usalife_new_6981_5year$avg_leadexp, na.rm=TRUE), 3))
mean_6 <- as.character(round(mean(usalife_new_6981_5year$avg_leadden_unadj, na.rm=TRUE), 3))
mean_6a <- as.character(round(mean(usalife_new_6981_5year$avg_leadden_adj, na.rm=TRUE), 3))

sd_0 <- as.character(round(sd(usalife_new$avg_leadexp, na.rm=TRUE), 3))
sd_2 <- as.character(round(sd(usalife_new_5year$avg_leadexp, na.rm=TRUE), 3))
sd_4 <- as.character(round(sd(usalife_new_6981_5year$avg_leadexp, na.rm=TRUE), 3))
sd_6 <- as.character(round(sd(usalife_new_6981_5year$avg_leadden_unadj, na.rm=TRUE), 3))
sd_6a <- as.character(round(sd(usalife_new_6981_5year$avg_leadden_adj, na.rm=TRUE), 3))


################################################################################
## Export the estimates
################################################################################

## Log

# agreeableness
stargazer(a0,a2,a4,a6,a6a,
          title = "Associations between atmospheric lead concentration/density and agreeableness in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.a0, cl.robust.se.a2, cl.robust.se.a4, cl.robust.se.a6, cl.robust.se.a6a),
          dep.var.labels = c("Agr"),
          covariate.labels = c("log(Lead concentration+1)","log(Lead density (unadj)+1)","log(Lead density (adj)+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Timeframe", "18years", "5years", "5years", "5years", "5years"), 
                           c("Cohort", "all in Schwaba et al.", "all in Schwaba et al.", "1969-1977", "1969-1977", "1969-1977"), 
                           c("Avglead", mean_0, mean_2, mean_4, mean_6, mean_6a),
                           c("Sdlead", sd_0, sd_2, sd_4, sd_6, sd_6a)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the specified timeframe (ug/m3).",
                    "Lead density is the average estimated lead emissions divided by area within the county over the specified timeframe (ug/m2).",
                    "The mean and standard deviation of relevant lead measures (without transformations) are displayed as 'Avglead' and 'Sdlead' in the table."), 
          notes.align = "l",
          type = "html", out = "replication_results/baseresults_agreeableness_appendix.htm")


# consciousness
stargazer(c0,c2,c4,c6,c6a,
          title = "Associations between atmospheric lead concentration/density and conscientiousness in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.c0, cl.robust.se.c2, cl.robust.se.c4, cl.robust.se.c6, cl.robust.se.c6a),
          dep.var.labels = c("Con"),
          covariate.labels = c("log(Lead concentration+1)","log(Lead density (unadj)+1)","log(Lead density (adj)+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Timeframe", "18years", "5years", "5years", "5years", "5years"), 
                           c("Cohort", "all in Schwaba et al.", "all in Schwaba et al.", "1969-1977", "1969-1977", "1969-1977"), 
                           c("Avglead", mean_0, mean_2, mean_4, mean_6, mean_6a),
                           c("Sdlead", sd_0, sd_2, sd_4, sd_6, sd_6a)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the specified timeframe (ug/m3).",
                    "Lead density is the average estimated lead emissions divided by area within the county over the specified timeframe (ug/m2).",
                    "The mean and standard deviation of relevant lead measures (without transformations) are displayed as 'Avglead' and 'Sdlead' in the table."), 
          notes.align = "l",
          type = "html", out = "replication_results/baseresults_conscientiousness_appendix.htm")


# extraversion
stargazer(e0,e2,e4,e6,e6a,
          title = "Associations between atmospheric lead concentration/density and extraversion in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.e0, cl.robust.se.e2, cl.robust.se.e4, cl.robust.se.e6, cl.robust.se.e6a),
          dep.var.labels = c("Ext"),
          covariate.labels = c("log(Lead concentration+1)","log(Lead density (unadj)+1)","log(Lead density (adj)+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Timeframe", "18years", "5years", "5years", "5years", "5years"), 
                           c("Cohort", "all in Schwaba et al.", "all in Schwaba et al.", "1969-1977", "1969-1977", "1969-1977"), 
                           c("Avglead", mean_0, mean_2, mean_4, mean_6, mean_6a),
                           c("Sdlead", sd_0, sd_2, sd_4, sd_6, sd_6a)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the specified timeframe (ug/m3).",
                    "Lead density is the average estimated lead emissions divided by area within the county over the specified timeframe (ug/m2).",
                    "The mean and standard deviation of relevant lead measures (without transformations) are displayed as 'Avglead' and 'Sdlead' in the table."), 
          notes.align = "l",
          type = "html", out = "replication_results/baseresults_extraversion_appendix.htm")


# neuroticism
stargazer(n0,n2,n4,n6,n6a,
          title = "Associations between atmospheric lead concentration/density and neuroticism in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.n0, cl.robust.se.n2, cl.robust.se.n4, cl.robust.se.n6, cl.robust.se.n6a),
          dep.var.labels = c("Neu"),
          covariate.labels = c("log(Lead concentration+1)","log(Lead density (unadj)+1)","log(Lead density (adj)+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Timeframe", "18years", "5years", "5years", "5years", "5years"), 
                           c("Cohort", "all in Schwaba et al.", "all in Schwaba et al.", "1969-1977", "1969-1977", "1969-1977"), 
                           c("Avglead", mean_0, mean_2, mean_4, mean_6, mean_6a),
                           c("Sdlead", sd_0, sd_2, sd_4, sd_6, sd_6a)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the specified timeframe (ug/m3).",
                    "Lead density is the average estimated lead emissions divided by area within the county over the specified timeframe (ug/m2).",
                    "The mean and standard deviation of relevant lead measures (without transformations) are displayed as 'Avglead' and 'Sdlead' in the table."), 
          notes.align = "l",
          type = "html", out = "replication_results/baseresults_neuroticism_appendix.htm")


# openness to experience
stargazer(o0,o2,o4,o6,o6a,
          title = "Associations between atmospheric lead concentration/density and openness in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.o0, cl.robust.se.o2, cl.robust.se.o4, cl.robust.se.o6, cl.robust.se.o6a),
          dep.var.labels = c("Ope"),
          covariate.labels = c("log(Lead concentration+1)","log(Lead density (unadj)+1)","log(Lead density (adj)+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Timeframe", "18years", "5years", "5years", "5years", "5years"), 
                           c("Cohort", "all in Schwaba et al.", "all in Schwaba et al.", "1969-1977", "1969-1977", "1969-1977"), 
                           c("Avglead", mean_0, mean_2, mean_4, mean_6, mean_6a),
                           c("Sdlead", sd_0, sd_2, sd_4, sd_6, sd_6a)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the specified timeframe (ug/m3).",
                    "Lead density is the average estimated lead emissions divided by area within the county over the specified timeframe (ug/m2).",
                    "The mean and standard deviation of relevant lead measures (without transformations) are displayed as 'Avglead' and 'Sdlead' in the table."), 
          notes.align = "l",
          type = "html", out = "replication_results/baseresults_openness_appendix.htm")
