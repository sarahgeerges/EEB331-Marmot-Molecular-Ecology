# AUTHOR / EMAIL: Sarah Geerges (sg6324@princeton.edu)
# CREATION DATE: 4/22/26
# MODIFIED DATE: 5/5/26
# PURPOSE: Find correlation between social interactions and relatedness/inbreeding

library(readxl)
library(dplyr) # for data manipulation
library(ggplot2) # for plotting
library(cowplot) # for plotting multiple plots in one
library(lme4)
library(lmerTest)
library(tidyr)
library(glmmTMB)
library(ggeffects)

# set working directory
setwd("/Users/sarahgeerges/Library/CloudStorage/GoogleDrive-sg6324@princeton.edu/My Drive/EEB331")

#### reading in files ####
# metadata 
meta <- read_excel("RADseq_pedigree_metadata.xlsx")
head(meta)
# per comm. marmots whose age was NA were actually adults
meta$age[is.na(meta$age)] <- 'A'
# pairwise interactions
soc <- read.csv("marmot_social_interactions.csv")
head(soc)
# pairwise relatedness 
load("rel.RData")
rel

#### cleaning and merging interactions data with metadata ####
meta_clean <- meta %>%
  select(-dam, -sire, -fur_mark, -cuid) # remove unnecessary columns
head(meta_clean)
# removing marmots w/ no interactions
soc_clean <- soc %>%
  select(-fur_mark_1, -fur_mark_2) %>%
  filter(n_interactions != 0)
head(soc_clean)

# merge by uid (pairwise)
soc_meta <- soc_clean %>%
  left_join(meta_clean, by = c("uid_1" = "uid")) %>%
  rename(yrborn1 = yrborn, sex1 = sex, age1 = age, col_area1 = col_area, year1 = year)
head(soc_meta)
soc_meta <- soc_meta %>%
  left_join(meta_clean, by = c("uid_2" = "uid")) %>%
  rename(yrborn2 = yrborn, sex2 = sex, age2 = age, col_area2 = col_area, year2 = year)
head(soc_meta)

# add a column for interactions by colony (within or between)
soc_meta$colony_type <- ifelse(
  soc_meta$col_area1 == soc_meta$col_area2,
  "Within colony",
  "Between colonies")
head(soc_meta)

# add columns for interactions by sex and age pair
soc_meta <- soc_meta %>%
  mutate(sex_pair = case_when(sex1 == "F" & sex2 == "F" ~ "F-F",
                              sex1 == "M" & sex2 == "M" ~ "M-M",
                              TRUE ~ "F-M")) %>%
  mutate(age_pair = case_when(age1 == "A" & age2 == "A" ~ "A-A",
                              age1 == "Y" & age2 == "Y" ~ "Y-Y",
                              age1 == "P" & age2 == "P" ~ "P-P",
                              (age1 == "A" & age2 == "Y") | (age1 == "Y" & age2 == "A") ~ "A-Y",
                              (age1 == "A" & age2 == "P") | (age1 == "P" & age2 == "A") ~ "A-P",
                              (age1 == "Y" & age2 == "P") | (age1 == "P" & age2 == "Y") ~ "Y-P"))
head(soc_meta)

# converting to long format to reduce pairwise age columns into one column
long_soc_meta <- soc_meta %>%
  pivot_longer(cols = c(age1, age2),
               names_to = "age_of",
               values_to = "age")
View(long_soc_meta)


#### Creating dataset w/ removing pups ####
# pairwise
nopup_soc_meta <- soc_meta %>%
  filter(age1 != 'P') %>%
  filter(age2 != 'P')
head(nopup_soc_meta)
# long format
long_nopup_soc_meta <- long_soc_meta %>%
  filter(age != 'P')
View(long_nopup_soc_meta)

#### Descriptive stats - social int ####
## interactions by colony 
# all marmots
ggplot(soc_meta, aes(x = colony_type)) +
  geom_bar() # within colony interactions most common
# no pups
ggplot(nopup_soc_meta, aes(x = colony_type)) +
  geom_bar() # within colony interactions most common

## interactions by sex
# all marmots
ggplot(soc_meta, aes(x = sex_pair)) +
  geom_bar() # F-F more common; M-M least common
# no pups
ggplot(nopup_soc_meta, aes(x = sex_pair)) +
  geom_bar() # F-F most common; M-M least common 

## interactions by age 
# all marmots
ggplot(long_soc_meta, aes(x = age)) +
  geom_bar() # adults most common
## no pups 
ggplot(long_nopup_soc_meta, aes(x = age)) +
  geom_bar() # adults most common

## interactions by age pair 
# all marmots
ggplot(soc_meta, aes(x = age_pair)) +
  geom_bar() # adult-adult int more common 
# no pups 
ggplot(nopup_soc_meta, aes(x = age_pair)) +
  geom_bar() # adult-adult int more common 

## Histogram with facet of age pair and sex pair -> added to EEB331 google slides
# all marmots
p1 <- ggplot(soc_meta, aes(x = avg_int_per_year)) +
  geom_histogram(bins = 20) + facet_grid(rows = vars(age_pair), cols = vars(sex_pair)) +
  labs(x = "Average Number of Interactions (Year)", y = "Frequency")
p1 
ggsave("plots/histogram_part11.png", plot = p1, width = 10, height = 8, units = "in", dpi=500, bg="white")

## Histogram of Interactions
# all marmots
ggplot(soc_meta, aes(x = avg_int_per_year)) +
  geom_histogram(bins = 20) +
  labs(title = "Distribution of Average Number of Interactions Per Year", 
       x = "Average Number of Interactions", y = "Count")
# no pups
ggplot(nopup_soc_meta, aes(x = avg_int_per_year)) +
  geom_histogram(bins = 20) +
  labs(title = "Distribution of Average Number of Interactions Per Year (excluding pups)", 
       x = "Average Number of Interactions", y = "Count")

#### Extracting relatedness estimators - wang ####
pairwise_rel <- rel$relatedness
wang <- pairwise_rel[,c("ind1.id", "ind2.id", "wang")]
head(wang)

#### Mapping interactions with relatedness ####
# invalid ID names
invalid_ids <- c("cUID_unk_4167_marmot", "RMBL_marmot_05-06",
                 "RMBL_marmot_no_label")
# removing invalid ID names
wang <- wang %>%
  mutate(ind1.id = trimws(ind1.id),
         ind2.id = trimws(ind2.id),
         ind1.id = gsub("^UID_", "", ind1.id),
         ind2.id = gsub("^UID_", "", ind2.id),
         ind1.id = gsub("_so.*$", "", ind1.id),
         ind2.id = gsub("_so.*$", "", ind2.id)) %>%
  filter(!ind1.id %in% invalid_ids,
         !ind2.id %in% invalid_ids)
head(wang)

# standardizing pair order in interactions and wang datasets
soc_meta <- soc_meta %>%
  mutate(uid_min = pmin(uid_1, uid_2),
         uid_max = pmax(uid_1, uid_2)) %>%
  select(uid_min, uid_max, everything(), -uid_1, -uid_2) # keeping new order of ids, and removing original
head(soc_meta) 

wang <- wang %>%
  mutate(uid_min = pmin(ind1.id, ind2.id),
         uid_max = pmax(ind1.id, ind2.id)) %>%
  select(uid_min, uid_max, everything(), -ind1.id, -ind2.id) %>% # keeping new order of ids, and removing original
  mutate(uid_min = gsub("_sort$", "", uid_min),
         uid_max = gsub("_sort$", "", uid_max))
head(wang)

# merging relatedness matrix with interactions
complete_data <- soc_meta %>%
  left_join(wang, by = c("uid_min", "uid_max"))
head(complete_data)

#### Complete dataset excluding pups ####
nopup_complete_data <- complete_data %>%
  filter(age1 != 'P') %>%
  filter(age2 != 'P')
head(nopup_complete_data)

#### Creating separate datasets by sex pair ####
## Female only
# all ages
f_only <- complete_data %>%
  filter(sex_pair == "F-F")
head(f_only)
# no pups
np_f <- f_only %>%
  filter(age1 != 'P') %>%
  filter(age2 != 'P')
head(np_f)

## Female-Male
# all ages
fm <- complete_data %>%
  filter(sex_pair == "F-M")
head(fm)
# no pups
np_fm <- fm %>%
  filter(age1 != 'P') %>%
  filter(age2 != 'P')
head(np_fm)

## Male only
# all ages
m_only <- complete_data %>%
  filter(sex_pair == "M-M")
head(m_only)
# no pups
np_m <- m_only %>%
  filter(age1 != 'P') %>%
  filter(age2 != 'P')
head(np_m)

#### Fitting models ####
#### Negative Binomial ####
# only takes takes integer values (count data) -> will use n_interactions instead (not used in Results Manuscript)
## All individuals
m1_all <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = complete_data)
summary(m1_all)
m2_all <- glmmTMB(n_interactions ~ wang + factor(age1) + factor(age2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = complete_data)
summary(m2_all)
m3_all <- glmmTMB(n_interactions ~ wang + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = complete_data)
summary(m3_all)
m4_all <- glmmTMB(n_interactions ~ wang + factor(age1) + factor(age2) + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = complete_data)
summary(m4_all)
# m4, smallest AIC -> sig pvalue
## repeating w/ no pups
m1_np <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = nopup_complete_data)
summary(m1_np)
m2_np <- glmmTMB(n_interactions ~ wang + factor(age1) + factor(age2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = nopup_complete_data)
summary(m2_np)
m3_np <- glmmTMB(n_interactions ~ wang + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = nopup_complete_data)
summary(m3_np)
m4_np <- glmmTMB(n_interactions ~ wang + factor(age1) + factor(age2) + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = nbinom2, data = nopup_complete_data)
summary(m4_np)
# m2, smallest AIC -> sig pvalue

#### Females only ####
## all ages
m1_f <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = f_only)
summary(m1_f)
# sig pvalue
## no pups
m1_npf <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = np_f)
summary(m1_npf)
# sig pvalue

#### Female-male ####
## all ages
m1_fm <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = fm)
summary(m1_fm)
# non-sig pvalue
## no pups
m1_npfm <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = np_fm)
summary(m1_npfm)
# non-sig pvalue

#### Male only ####
## all ages
m1_m <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = m_only)
summary(m1_m)
# non-sig pvalue
## no pups
m1_npm <- glmmTMB(n_interactions ~ wang + (1 | uid_min) + (1| uid_max), family = nbinom2, data = np_m)
summary(m1_npm)
# non-sig pvalue


#### Gamma GLM ####
# Negative binomial does not take in non-integer values (only count)
# I want to use avg_int_per_year because it reflects the frequency of interations in same year overlap
## All individuals
m1_all <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = complete_data)
summary(m1_all)
m2_all <- glmmTMB(avg_int_per_year ~ wang + factor(age1) + factor(age2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = complete_data)
summary(m2_all)
m3_all <- glmmTMB(avg_int_per_year ~ wang + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = complete_data)
summary(m3_all)
m4_all <- glmmTMB(avg_int_per_year ~ wang + factor(age1) + factor(age2) + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = complete_data)
summary(m4_all)
# m3, smallest AIC -> sig p-value < 0.05
## repeating with no pups
m1_np <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = nopup_complete_data)
summary(m1_np)
m2_np <- glmmTMB(avg_int_per_year ~ wang + factor(age1) + factor(age2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = nopup_complete_data)
summary(m2_np)
m3_np <- glmmTMB(avg_int_per_year ~ wang + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = nopup_complete_data)
summary(m3_np)
m4_np <- glmmTMB(avg_int_per_year ~ wang + factor(age1) + factor(age2) + factor(sex1) + factor(sex2) + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = nopup_complete_data)
summary(m4_np)
# m3, smallest AIC -> sig p-value < 0.05

#### Females only ####
## all ages
m1_f <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = f_only)
summary(m1_f)
# sig p-value < 0.05
## no pups
m1_npf <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = np_f)
summary(m1_npf)
# sig p-value < 0.05

#### Female-male ####
## all ages
m1_fm <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = fm)
summary(m1_fm)
# non-sig p-value > 0.05
## no pups
m1_npfm <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = np_fm)
summary(m1_npfm)
# non-sig p-value > 0.05

#### Male only ####
## all ages
m1_m <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = m_only)
summary(m1_m)
# non-sig p-value
## no pups
m1_npm <- glmmTMB(avg_int_per_year ~ wang + (1 | uid_min) + (1| uid_max), family = Gamma(link = "log"), data = np_m)
summary(m1_npm)
# non-sig p-value

#### Summary ####
# wang relatedness is significantly correlated with interactions among F-F (w/ and w/o pups)
# relatedness is not significantly correlated with interactions among F-M (w and w/o pups)
# relatedness is not significantly correlated with interactions among M-M (w/ and w/o pups)


