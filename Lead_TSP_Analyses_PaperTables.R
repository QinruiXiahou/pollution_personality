## Working environment setup
rm(list=ls())
options(scipen = 999)

memory.limit(size=56000) 

library(nlme)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stargazer)
library(spdep)
library(RANN)
library(spatialreg)
library(tidyr)
library(haven)
library(ivreg)
library(multiwayvcov)

setwd('L:/Public/Datasets1/Project-RetrospectiveAnalysis/replication')

## Read in data
usalife_new_6981_5year <- readRDS("usalife_new_6981_5year.rds")
usalife_new_6981_1year <- readRDS("usalife_new_6981_1year.rds") %>% filter(yearborn <= 1977)
balancedTSPpanel6977 <- read_dta("balancedTSPpanel6977.dta") %>% 
  mutate(geoid = paste(str_pad(fips_state, 2, side = "left", pad = 0),
                       str_pad(fips_cnty, 3, side = "left", pad = 0),
                       sep = ""))

## Merge datasets
leadtsp_6977 <- usalife_new_6981_5year %>%
  left_join(balancedTSPpanel6977, by=c("geoid", "yearborn"="year")) %>%
  filter(is.na(fips_state)==FALSE) %>%
  mutate(yearFE = as.factor(yearborn))
length(unique(leadtsp_6977$name))

## Summary stats
mean_exp <- as.character(round(mean(usalife_new_6981_5year$avg_leadexp, na.rm=TRUE), 3))
sd_exp <- as.character(round(sd(usalife_new_6981_5year$avg_leadexp, na.rm=TRUE), 3))
mean_den <- as.character(round(mean(usalife_new_6981_5year$avg_leadden_adj, na.rm=TRUE), 3))
sd_den <- as.character(round(sd(usalife_new_6981_5year$avg_leadden_adj, na.rm=TRUE), 3))

mean_exp_1year <- as.character(round(mean(usalife_new_6981_1year$avg_leadexp, na.rm=TRUE), 3))
sd_exp_1year <- as.character(round(sd(usalife_new_6981_1year$avg_leadexp, na.rm=TRUE), 3))
mean_den_1year <- as.character(round(mean(usalife_new_6981_1year$avg_leadden_adj, na.rm=TRUE), 3))
sd_den_1year <- as.character(round(sd(usalife_new_6981_1year$avg_leadden_adj, na.rm=TRUE), 3))

mean_exp_red <- as.character(round(mean(leadtsp_6977$avg_leadexp, na.rm=TRUE), 3))
sd_exp_red <- as.character(round(sd(leadtsp_6977$avg_leadexp, na.rm=TRUE), 3))
mean_den_red <- as.character(round(mean(leadtsp_6977$avg_leadden_adj, na.rm=TRUE), 3))
sd_den_red <- as.character(round(sd(leadtsp_6977$avg_leadden_adj, na.rm=TRUE), 3))
mean_tsp_red <- as.character(round(mean(leadtsp_6977$countymeanTSP, na.rm=TRUE), 3))
sd_tsp_red <- as.character(round(sd(leadtsp_6977$countymeanTSP, na.rm=TRUE), 3))


################################################################################
## Analysis 1.1: lead concentration and personality (log)
################################################################################
log_exp_a <- lm(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_exp_a <- cluster.vcov(log_exp_a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_exp_a <- sqrt(diag(cl.log_exp_a))

log_exp_c <- lm(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_exp_c <- cluster.vcov(log_exp_c, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_exp_c <- sqrt(diag(cl.log_exp_c))

log_exp_e <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_exp_e <- cluster.vcov(log_exp_e, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_exp_e <- sqrt(diag(cl.log_exp_e))

log_exp_n <- lm(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_exp_n <- cluster.vcov(log_exp_n, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_exp_n <- sqrt(diag(cl.log_exp_n))

log_exp_o <- lm(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_exp_o <- cluster.vcov(log_exp_o, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_exp_o <- sqrt(diag(cl.log_exp_o))

stargazer(log_exp_a, log_exp_c, log_exp_e, log_exp_n, log_exp_o,
          title = "Associations between lead concentration and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_exp_a, cl.robust.se.log_exp_c, cl.robust.se.log_exp_e, cl.robust.se.log_exp_n, cl.robust.se.log_exp_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead concentration+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_exp, mean_exp, mean_exp, mean_exp, mean_exp),
                           c("Sdlead", sd_exp, sd_exp, sd_exp, sd_exp, sd_exp)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the first five years of life (ug/m3).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table1_leadexp_log_cohort.htm")


################################################################################
## Analysis 1.2: lead density and personality (log)
################################################################################
log_den_a <- lm(sa ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_den_a <- cluster.vcov(log_den_a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_den_a <- sqrt(diag(cl.log_den_a))

log_den_c <- lm(sc ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_den_c <- cluster.vcov(log_den_c, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_den_c <- sqrt(diag(cl.log_den_c))

log_den_e <- lm(se ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_den_e <- cluster.vcov(log_den_e, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_den_e <- sqrt(diag(cl.log_den_e))

log_den_n <- lm(sn ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_den_n <- cluster.vcov(log_den_n, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_den_n <- sqrt(diag(cl.log_den_n))

log_den_o <- lm(so ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.log_den_o <- cluster.vcov(log_den_o, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.log_den_o <- sqrt(diag(cl.log_den_o))

stargazer(log_den_a, log_den_c, log_den_e, log_den_n, log_den_o,
          title = "Associations between lead density and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_den_a, cl.robust.se.log_den_c, cl.robust.se.log_den_e, cl.robust.se.log_den_n, cl.robust.se.log_den_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead density+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_den, mean_den, mean_den, mean_den, mean_den),
                           c("Sdlead", sd_den, sd_den, sd_den, sd_den, sd_den)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead density is the average estimated lead emissions divided by area within the county over the first five years of life (ug/m2).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table2_leadden_log_cohort.htm")


################################################################################
## Analysis 1.3: lead concentration and personality (log, 1year)
################################################################################
log_exp_a1 <- lm(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_exp_a1 <- cluster.vcov(log_exp_a1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_exp_a1 <- sqrt(diag(cl.log_exp_a1))

log_exp_c1 <- lm(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_exp_c1 <- cluster.vcov(log_exp_c1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_exp_c1 <- sqrt(diag(cl.log_exp_c1))

log_exp_e1 <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_exp_e1 <- cluster.vcov(log_exp_e1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_exp_e1 <- sqrt(diag(cl.log_exp_e1))

log_exp_n1 <- lm(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_exp_n1 <- cluster.vcov(log_exp_n1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_exp_n1 <- sqrt(diag(cl.log_exp_n1))

log_exp_o1 <- lm(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_exp_o1 <- cluster.vcov(log_exp_o1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_exp_o1 <- sqrt(diag(cl.log_exp_o1))

stargazer(log_exp_a1, log_exp_c1, log_exp_e1, log_exp_n1, log_exp_o1,
          title = "Associations between lead concentration and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_exp_a1, cl.robust.se.log_exp_c1, cl.robust.se.log_exp_e1, cl.robust.se.log_exp_n1, cl.robust.se.log_exp_o1),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead concentration+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_exp_1year, mean_exp_1year, mean_exp_1year, mean_exp_1year, mean_exp_1year),
                           c("Sdlead", sd_exp_1year, sd_exp_1year, sd_exp_1year, sd_exp_1year, sd_exp_1year)),
          notes = c("All regressions include county and birth year fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county in the birth year (ug/m3).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table3_leadexp_log_cohort_1year.htm")


################################################################################
## Analysis 1.4: lead density and personality (log, 1year)
################################################################################
log_den_a1 <- lm(sa ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_den_a1 <- cluster.vcov(log_den_a1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_den_a1 <- sqrt(diag(cl.log_den_a1))

log_den_c1 <- lm(sc ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_den_c1 <- cluster.vcov(log_den_c1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_den_c1 <- sqrt(diag(cl.log_den_c1))

log_den_e1 <- lm(se ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_den_e1 <- cluster.vcov(log_den_e1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_den_e1 <- sqrt(diag(cl.log_den_e1))

log_den_n1 <- lm(sn ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_den_n1 <- cluster.vcov(log_den_n1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_den_n1 <- sqrt(diag(cl.log_den_n1))

log_den_o1 <- lm(so ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_1year)
cl.log_den_o1 <- cluster.vcov(log_den_o1, usalife_new_6981_1year$name) # cluster-robust SEs
cl.robust.se.log_den_o1 <- sqrt(diag(cl.log_den_o1))

stargazer(log_den_a1, log_den_c1, log_den_e1, log_den_n1, log_den_o1,
          title = "Associations between lead density and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_den_a1, cl.robust.se.log_den_c1, cl.robust.se.log_den_e1, cl.robust.se.log_den_n1, cl.robust.se.log_den_o1),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead density+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_den_1year, mean_den_1year, mean_den_1year, mean_den_1year, mean_den_1year),
                           c("Sdlead", sd_den_1year, sd_den_1year, sd_den_1year, sd_den_1year, sd_den_1year)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead density is the average estimated lead emissions divided by area within the county in the birth year (ug/m2).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table4_leadden_log_cohort_1year.htm")


################################################################################
## Analysis 1.5: lead concentration and personality (linear)
################################################################################
exp_a <- lm(sa ~ avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.exp_a <- cluster.vcov(exp_a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.exp_a <- sqrt(diag(cl.exp_a))

exp_c <- lm(sc ~ avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.exp_c <- cluster.vcov(exp_c, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.exp_c <- sqrt(diag(cl.exp_c))

exp_e <- lm(se ~ avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.exp_e <- cluster.vcov(exp_e, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.exp_e <- sqrt(diag(cl.exp_e))

exp_n <- lm(sn ~ avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.exp_n <- cluster.vcov(exp_n, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.exp_n <- sqrt(diag(cl.exp_n))

exp_o <- lm(so ~ avg_leadexp + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.exp_o <- cluster.vcov(exp_o, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.exp_o <- sqrt(diag(cl.exp_o))

stargazer(exp_a, exp_c, exp_e, exp_n, exp_o,
          title = "Associations between lead concentration and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.exp_a, cl.robust.se.exp_c, cl.robust.se.exp_e, cl.robust.se.exp_n, cl.robust.se.exp_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("Lead concentration","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_exp, mean_exp, mean_exp, mean_exp, mean_exp),
                           c("Sdlead", sd_exp, sd_exp, sd_exp, sd_exp, sd_exp)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the first five years of life (ug/m3).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table5_leadexp_linear_cohort.htm")


################################################################################
## Analysis 1.6: lead density and personality (linear)
################################################################################
den_a <- lm(sa ~ avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.den_a <- cluster.vcov(den_a, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.den_a <- sqrt(diag(cl.den_a))

den_c <- lm(sc ~ avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.den_c <- cluster.vcov(den_c, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.den_c <- sqrt(diag(cl.den_c))

den_e <- lm(se ~ avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.den_e <- cluster.vcov(den_e, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.den_e <- sqrt(diag(cl.den_e))

den_n <- lm(sn ~ avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.den_n <- cluster.vcov(den_n, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.den_n <- sqrt(diag(cl.den_n))

den_o <- lm(so ~ avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = usalife_new_6981_5year)
cl.den_o <- cluster.vcov(den_o, usalife_new_6981_5year$name) # cluster-robust SEs
cl.robust.se.den_o <- sqrt(diag(cl.den_o))

stargazer(den_a, den_c, den_e, den_n, den_o,
          title = "Associations between lead density and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.den_a, cl.robust.se.den_c, cl.robust.se.den_e, cl.robust.se.den_n, cl.robust.se.den_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("Lead density","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_den, mean_den, mean_den, mean_den, mean_den),
                           c("Sdlead", sd_den, sd_den, sd_den, sd_den, sd_den)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead density is the average estimated lead emissions divided by area within the county over the first five years of life (ug/m2).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table6_leadden_linear_cohort.htm")


################################################################################
## Analysis 1.7: lead concentration and personality (log, reduced sample)
################################################################################
log_exp_a <- lm(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_exp_a <- cluster.vcov(log_exp_a, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_exp_a <- sqrt(diag(cl.log_exp_a))

log_exp_c <- lm(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_exp_c <- cluster.vcov(log_exp_c, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_exp_c <- sqrt(diag(cl.log_exp_c))

log_exp_e <- lm(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_exp_e <- cluster.vcov(log_exp_e, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_exp_e <- sqrt(diag(cl.log_exp_e))

log_exp_n <- lm(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_exp_n <- cluster.vcov(log_exp_n, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_exp_n <- sqrt(diag(cl.log_exp_n))

log_exp_o <- lm(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_exp_o <- cluster.vcov(log_exp_o, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_exp_o <- sqrt(diag(cl.log_exp_o))

stargazer(log_exp_a, log_exp_c, log_exp_e, log_exp_n, log_exp_o,
          title = "Associations between lead concentration and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_exp_a, cl.robust.se.log_exp_c, cl.robust.se.log_exp_e, cl.robust.se.log_exp_n, cl.robust.se.log_exp_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead concentration+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_exp_red, mean_exp_red, mean_exp_red, mean_exp_red, mean_exp_red),
                           c("Sdlead", sd_exp_red, sd_exp_red, sd_exp_red, sd_exp_red, sd_exp_red)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the first five years of life (ug/m3).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table7_leadexp_log_cohort_redsample.htm")


################################################################################
## Analysis 1.8: lead density and personality (log, reduced sample)
################################################################################
log_den_a <- lm(sa ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_den_a <- cluster.vcov(log_den_a, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_den_a <- sqrt(diag(cl.log_den_a))

log_den_c <- lm(sc ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_den_c <- cluster.vcov(log_den_c, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_den_c <- sqrt(diag(cl.log_den_c))

log_den_e <- lm(se ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_den_e <- cluster.vcov(log_den_e, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_den_e <- sqrt(diag(cl.log_den_e))

log_den_n <- lm(sn ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_den_n <- cluster.vcov(log_den_n, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_den_n <- sqrt(diag(cl.log_den_n))

log_den_o <- lm(so ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.log_den_o <- cluster.vcov(log_den_o, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.log_den_o <- sqrt(diag(cl.log_den_o))

stargazer(log_den_a, log_den_c, log_den_e, log_den_n, log_den_o,
          title = "Associations between lead density and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.log_den_a, cl.robust.se.log_den_c, cl.robust.se.log_den_e, cl.robust.se.log_den_n, cl.robust.se.log_den_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("log(Lead density+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avglead", mean_den_red, mean_den_red, mean_den_red, mean_den_red, mean_den_red),
                           c("Sdlead", sd_den_red, sd_den_red, sd_den_red, sd_den_red, sd_den_red)),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "Lead density is the average estimated lead emissions divided by area within the county over the first five years of life (ug/m2).",
                    "The mean and standard deviation are displayed as 'Avglead' and 'Sdlead' in the table."),
          notes.align = "l",
          type = "html", out = "replication_results/table8_leadden_log_cohort_redsample.htm")


################################################################################
## Analysis 2.1: TSP and personality (IV)
################################################################################
tsp_a <- ivreg(sa ~ age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_a <- cluster.vcov(tsp_a, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_a <- sqrt(diag(cl.tsp_a))
summ.tsp_a <- summary(tsp_a, diagnostics=T)

tsp_c <- ivreg(sc ~ age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_c <- cluster.vcov(tsp_c, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_c <- sqrt(diag(cl.tsp_c))
summ.tsp_c <- summary(tsp_c, diagnostics=T)

tsp_e <- ivreg(se ~ age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_e <- cluster.vcov(tsp_e, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_e <- sqrt(diag(cl.tsp_e))
summ.tsp_e <- summary(tsp_e, diagnostics=T)

tsp_n <- ivreg(sn ~ age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_n <- cluster.vcov(tsp_n, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_n <- sqrt(diag(cl.tsp_n))
summ.tsp_n <- summary(tsp_n, diagnostics=T)

tsp_o <- ivreg(so ~ age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_o <- cluster.vcov(tsp_o, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_o <- sqrt(diag(cl.tsp_o))
summ.tsp_o <- summary(tsp_o, diagnostics=T)

stargazer(tsp_a, tsp_c, tsp_e, tsp_n, tsp_o,
          title = "Associations between TSP and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.tsp_a, cl.robust.se.tsp_c, cl.robust.se.tsp_e, cl.robust.se.tsp_n, cl.robust.se.tsp_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("TSP","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avgtsp", mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red),
                           c("Sdtsp", sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red),
                           c(rownames(summ.tsp_a$diagnostics)[1], 
                             round(summ.tsp_a$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_c$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_e$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_n$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_o$diagnostics[1, "p-value"], 3)), 
                           c(rownames(summ.tsp_a$diagnostics)[2], 
                             round(summ.tsp_a$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_c$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_e$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_n$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_o$diagnostics[2, "p-value"], 3))),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "TSP emissions are at the county level in the birth year.",
                    "The mean and standard deviation are displayed as 'Avgtsp' and 'Sdtsp' in the table.",
                    "The p-values of regression diagnostics are presented in 'Weak instruments' and 'Wu-Hausman' rows."),
          notes.align = "l",
          type = "html", out = "replication_results/table9_tsp_redform_cohort.htm")


################################################################################
## Analysis 2.2: TSP and personality (IV first stage)
################################################################################
tsp_stg1 <- lm(countymeanTSP ~ TSPtreat + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.tsp_stg1 <- cluster.vcov(tsp_stg1, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_stg1 <- sqrt(diag(cl.tsp_stg1))

tsp_stg1_exp <- lm(countymeanTSP ~ TSPtreat + log_avg_leadexp + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.tsp_stg1_exp <- cluster.vcov(tsp_stg1_exp, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_stg1_exp <- sqrt(diag(cl.tsp_stg1_exp))

tsp_stg1_den <- lm(countymeanTSP ~ TSPtreat + log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1, data = leadtsp_6977)
cl.tsp_stg1_den <- cluster.vcov(tsp_stg1_den, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_stg1_den <- sqrt(diag(cl.tsp_stg1_den))

stargazer(tsp_stg1, tsp_stg1_exp, tsp_stg1_den,
          title = "First stage results from IV regressions",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.tsp_stg1, cl.robust.se.tsp_stg1_exp, cl.robust.se.tsp_stg1_den),
          dep.var.labels = c("countymeanTSP"),
          covariate.labels = c("TSPtreat", "log(Lead concentration+1)", "log(Lead density+1)", "Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977")),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample."),
          notes.align = "l",
          type = "html", out = "replication_results/table10_tsp_firststage_cohort.htm")


################################################################################
## Analysis 3.1: TSP, lead concentration and personality (IV)
################################################################################
tsp_exp_a <- ivreg(sa ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_exp_a <- cluster.vcov(tsp_exp_a, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_exp_a <- sqrt(diag(cl.tsp_exp_a))
summ.tsp_exp_a <- summary(tsp_exp_a, diagnostics=T)

tsp_exp_c <- ivreg(sc ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_exp_c <- cluster.vcov(tsp_exp_c, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_exp_c <- sqrt(diag(cl.tsp_exp_c))
summ.tsp_exp_c <- summary(tsp_exp_c, diagnostics=T)

tsp_exp_e <- ivreg(se ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_exp_e <- cluster.vcov(tsp_exp_e, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_exp_e <- sqrt(diag(cl.tsp_exp_e))
summ.tsp_exp_e <- summary(tsp_exp_e, diagnostics=T)

tsp_exp_n <- ivreg(sn ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_exp_n <- cluster.vcov(tsp_exp_n, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_exp_n <- sqrt(diag(cl.tsp_exp_n))
summ.tsp_exp_n <- summary(tsp_exp_n, diagnostics=T)

tsp_exp_o <- ivreg(so ~ log_avg_leadexp + age + parentcollege + income + name + cohort - 1 | 
                 countymeanTSP |
                 TSPtreat, data = leadtsp_6977)
cl.tsp_exp_o <- cluster.vcov(tsp_exp_o, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_exp_o <- sqrt(diag(cl.tsp_exp_o))
summ.tsp_exp_o <- summary(tsp_exp_o, diagnostics=T)

stargazer(tsp_exp_a, tsp_exp_c, tsp_exp_e, tsp_exp_n, tsp_exp_o,
          title = "Associations between TSP, lead concentration and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.tsp_exp_a, cl.robust.se.tsp_exp_c, cl.robust.se.tsp_exp_e, cl.robust.se.tsp_exp_n, cl.robust.se.tsp_exp_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("TSP", "log(Lead concentration+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avgtsp", mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red),
                           c("Sdtsp", sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red),
                           c("Avglead", mean_exp_red, mean_exp_red, mean_exp_red, mean_exp_red, mean_exp_red),
                           c("Sdlead", sd_exp_red, sd_exp_red, sd_exp_red, sd_exp_red, sd_exp_red),
                           c(rownames(summ.tsp_exp_a$diagnostics)[1], 
                             round(summ.tsp_exp_a$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_exp_c$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_exp_e$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_exp_n$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_exp_o$diagnostics[1, "p-value"], 3)), 
                           c(rownames(summ.tsp_exp_a$diagnostics)[2], 
                             round(summ.tsp_exp_a$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_exp_c$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_exp_e$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_exp_n$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_exp_o$diagnostics[2, "p-value"], 3))),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "TSP emissions are at the county level in the birth year.",
                    "Lead concentration is the arithmetic mean of non-zero monitor readings within the county over the first five years of life (ug/m3).",
                    "The mean and standard deviation are displayed as 'Avglead/tsp' and 'Sdlead/tsp' in the table.",
                    "The p-values of regression diagnostics are presented in 'Weak instruments' and 'Wu-Hausman' rows."),
          notes.align = "l",
          type = "html", out = "replication_results/table11_tsp_exp_redform_cohort.htm")

################################################################################
## Analysis 3.2: TSP, lead density and personality (IV)
################################################################################
tsp_den_a <- ivreg(sa ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1 | 
                     countymeanTSP |
                     TSPtreat, data = leadtsp_6977)
cl.tsp_den_a <- cluster.vcov(tsp_den_a, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_den_a <- sqrt(diag(cl.tsp_den_a))
summ.tsp_den_a <- summary(tsp_den_a, diagnostics=T)

tsp_den_c <- ivreg(sc ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1 | 
                     countymeanTSP |
                     TSPtreat, data = leadtsp_6977)
cl.tsp_den_c <- cluster.vcov(tsp_den_c, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_den_c <- sqrt(diag(cl.tsp_den_c))
summ.tsp_den_c <- summary(tsp_den_c, diagnostics=T)

tsp_den_e <- ivreg(se ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1 | 
                     countymeanTSP |
                     TSPtreat, data = leadtsp_6977)
cl.tsp_den_e <- cluster.vcov(tsp_den_e, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_den_e <- sqrt(diag(cl.tsp_den_e))
summ.tsp_den_e <- summary(tsp_den_e, diagnostics=T)

tsp_den_n <- ivreg(sn ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1 | 
                     countymeanTSP |
                     TSPtreat, data = leadtsp_6977)
cl.tsp_den_n <- cluster.vcov(tsp_den_n, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_den_n <- sqrt(diag(cl.tsp_den_n))
summ.tsp_den_n <- summary(tsp_den_n, diagnostics=T)

tsp_den_o <- ivreg(so ~ log_avg_leadden_adj + age + parentcollege + income + name + cohort - 1 | 
                     countymeanTSP |
                     TSPtreat, data = leadtsp_6977)
cl.tsp_den_o <- cluster.vcov(tsp_den_o, leadtsp_6977$name) # cluster-robust SEs
cl.robust.se.tsp_den_o <- sqrt(diag(cl.tsp_den_o))
summ.tsp_den_o <- summary(tsp_den_o, diagnostics=T)

stargazer(tsp_den_a, tsp_den_c, tsp_den_e, tsp_den_n, tsp_den_o,
          title = "Associations between TSP, lead density and personality traits in the US",
          no.space = TRUE,
          omit.stat = c("f","ser"),
          digits = 3,
          omit = c("name","cohort"),
          se=list(cl.robust.se.tsp_den_a, cl.robust.se.tsp_den_c, cl.robust.se.tsp_den_e, cl.robust.se.tsp_den_n, cl.robust.se.tsp_den_o),
          dep.var.labels = c("Agr","Con","Ext","Neu","Ope"),
          covariate.labels = c("TSP", "log(Lead density+1)","Age","Parent College","Median County Income"),
          add.lines = list(c("Cohort", "1969-1977", "1969-1977", "1969-1977", "1969-1977", "1969-1977"),
                           c("Avgtsp", mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red, mean_tsp_red),
                           c("Sdtsp", sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red, sd_tsp_red),
                           c("Avglead", mean_den_red, mean_den_red, mean_den_red, mean_den_red, mean_den_red),
                           c("Sdlead", sd_den_red, sd_den_red, sd_den_red, sd_den_red, sd_den_red),
                           c(rownames(summ.tsp_den_a$diagnostics)[1], 
                             round(summ.tsp_den_a$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_den_c$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_den_e$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_den_n$diagnostics[1, "p-value"], 3), 
                             round(summ.tsp_den_o$diagnostics[1, "p-value"], 3)), 
                           c(rownames(summ.tsp_den_a$diagnostics)[2], 
                             round(summ.tsp_den_a$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_den_c$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_den_e$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_den_n$diagnostics[2, "p-value"], 3), 
                             round(summ.tsp_den_o$diagnostics[2, "p-value"], 3))),
          notes = c("All regressions include county and 5-year cohort fixed effects.",
                    "Standard errors are clustered at the county level.",
                    "Counties with smelters or in Alaska, Hawaii, and Puerto Rico are removed from the sample.",
                    "TSP emissions are at the county level in the birth year.",
                    "Lead density is the average estimated lead emissions divided by area within the county over the first five years of life (ug/m2).",
                    "The mean and standard deviation are displayed as 'Avglead/tsp' and 'Sdlead/tsp' in the table.",
                    "The p-values of regression diagnostics are presented in 'Weak instruments' and 'Wu-Hausman' rows."),
          notes.align = "l",
          type = "html", out = "replication_results/table12_tsp_den_redform_cohort.htm")


################################################################################
## Analysis 4.1: summary statistics (usalife_new_6981_5year)
################################################################################
summary(usalife_new_6981_5year)
lapply(usalife_new_6981_5year[, 1:76], sd)

################################################################################
## Analysis 4.2: summary statistics (leadtsp_6977)
################################################################################
summary(leadtsp_6977)
lapply(leadtsp_6977[, 1:76], sd)
lapply(leadtsp_6977[, "countymeanTSP"], sd)

################################################################################
## Analysis 4.3: summary statistics for certain cohort (usalife_new_6981_5year)
################################################################################
cohort69 <- usalife_new_6981_5year %>% filter(yearborn==1969)
cohort77 <- usalife_new_6981_5year %>% filter(yearborn==1977)

summary(cohort69)
summary(cohort77)

lapply(cohort69[, 1:76], sd)
lapply(cohort77[, 1:76], sd)
