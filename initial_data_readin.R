
#Read in and clean-up data

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(udunits2)

current_data <- read_csv("CH4_7.28.csv",
                         trim_ws = TRUE,
                         #0       10        20        30
                         #123456789012345678901234567890123456                     
                         col_types = "ccdddcdccccdddcccddddddccddddddccccc")

current_data %>%
  filter(is.na(`Discard?`)) %>%
  select(-`Discard?`) -> good_data

good_data %>%
  group_by(CH4_flux_unit_V2) %>%
  mutate(count = length(Study_number)) %>%
  select(CH4_flux_unit_V2, count) %>%
  unique(.) %>%
  arrange(desc(count)) -> unit_counts

#still some conversions to do

good_data %>%
  group_by(Manipulation_V2) %>%
  mutate(count = length(Study_number)) %>%
  select(Manipulation_V2, count) %>%
  unique(.) %>%
  arrange(desc(count)) -> manipulation_counts

#"Fertilization + Precipitation as well as Precipitation + Fertilization

#filter out units that we aren't going to use at this time

#select just columns for unit transformation
good_data %>%
  select(Study_number, RN, CH4_annual, CH4_growingseason, CH4_monthly,
         `CH4-C_converted`, Season_converted, N, CH4_flux_unit,
         CH4_flux_unit_V2, SD_CH4_annual, SD_CH4_growingseason,
         SD_CH4_monthly) -> for_transformation

#current csv has zeros for Study # 23035 where NA's should be in Season_converted
for_transformation[for_transformation$Study_number == "23035",]$Season_converted <- NA

for_transformation %>%
  filter(! is.na(CH4_annual)) %>%
  mutate(period = "annual",
         flux = CH4_annual,
         flux_sd = SD_CH4_annual) %>%
  select(- c(CH4_annual,
             SD_CH4_annual,
             CH4_growingseason,
             SD_CH4_growingseason,
             CH4_monthly,
             SD_CH4_monthly)) -> annual

for_transformation %>%
  filter(! is.na(CH4_monthly)) %>%
  mutate(period = "monthly",
         flux = CH4_monthly,
         flux_sd = SD_CH4_monthly) %>%
  select(- c(CH4_annual,
             SD_CH4_annual,
             CH4_growingseason,
             SD_CH4_growingseason,
             CH4_monthly,
             SD_CH4_monthly))-> monthly

for_transformation %>%
  filter(! is.na(CH4_growingseason)) %>%
  mutate(period = "growing season",
         flux = CH4_growingseason,
         flux_sd = SD_CH4_growingseason) %>%
  select(- c(CH4_annual,
             SD_CH4_annual,
             CH4_growingseason,
             SD_CH4_growingseason,
             CH4_monthly,
             SD_CH4_monthly))-> growing_season

annual %>%
  bind_rows(monthly) %>%
  bind_rows(growing_season) %>%
  arrange(Study_number) %>%
  mutate(converted = ifelse(is.na(Season_converted),
                            `CH4-C_converted`,
                            Season_converted),
         converted = coalesce(converted, flux),
         cf = converted/flux,
         flux_sd = flux_sd*cf) %>%
  select(-c(`CH4-C_converted`, Season_converted)) -> combined

#all fluxes are in one column but not yet in the same unit
#line 100-101 will remove any studies with insufficient data for meta analysis
#only run if complete cases are needed
combined %>%
  filter(complete.cases(.)) -> combined

combined %>%
  select(CH4_flux_unit_V2) %>%
  unique(.)

# 1 mg/m2/h         
# 2 kg/ha/yr        
# 3 μg/m2/h         
# 4 mg/m2/d         
# 5 mg/m2/yr        
# 6 g/m2/yr         
# 7 g/ha/d          
# 8 CH4 μmol/mol  #this is the only unit that we can't convert  
# 9 kg/hm2/yr       
# 10 μg/m2/s 

#put all fluxes and their error in mg CH4 per meter squared per hour
#WARNING
#study number 23030 passes through this argument unconverted
combined %>%
  mutate(stnd_flux = ud.convert(converted, CH4_flux_unit_V2, "mg/m2/h"),
         stnd_flux_sd = ud.convert(flux_sd, CH4_flux_unit_V2, "mg/m2/h")) %>%
  select(-c(flux, flux_sd, CH4_flux_unit, CH4_flux_unit_V2,
            converted, cf))-> tidy_standard_unit_data

current_data %>%
  select(Study_number, RN, Study_midyear, Ecosystem_type,
         Manipulation_V2, Manipulation_level, Latitude,
         Longitude, Elevation, Soil_type, Soil_drainage,
         SM_value, SM_sd, SM_depth, SM_unit,
         `same timescale as flux? (T/F)`, `if False, timescale`) %>%
  right_join(tidy_standard_unit_data) -> metadata


metadata %>%
  group_by(Study_number) %>%
  filter(Manipulation_V2 == "None") %>%
  mutate(Control_SM = SM_value,
         Control_SM_sd = SM_sd,
         Control_SM_depth = SM_depth,
         Control_N = N,
         Control = stnd_flux,
         Control_sd = stnd_flux_sd) %>%
  select(-c(RN, Manipulation_V2, Manipulation_level, SM_value, SM_sd, SM_depth, N, stnd_flux, stnd_flux_sd,
            SM_unit, `same timescale as flux? (T/F)`,
            `if False, timescale`)) -> controls

metadata %>%
  group_by(Study_number) %>%
  filter(Manipulation_V2 != "None") %>%
  mutate(Manip_SM = SM_value,
         Manip_SM_sd = SM_sd,
         Manip_SM_depth = SM_depth,
         Manip_N = N,
         Manip = stnd_flux,
         Manip_sd = stnd_flux_sd) %>%
  select(-c(RN, SM_value, SM_sd, SM_depth, N, stnd_flux, stnd_flux_sd)) %>%
  left_join(controls, by = c("Study_number", "Study_midyear",
                             "Ecosystem_type", "Latitude",
                             "Longitude", "Elevation",
                             "Soil_type", "Soil_drainage",
                             "period")) -> tidy_metadata


tidy_metadata %>%
  filter(! is.na(Control)) -> data4meta


library("metafor")
methane_data<- escalc(n1i = Control_N, n2i = Manip_N,
                      m1i = Control, m2i = Manip,
                      sd1i = Control_sd, sd2i = Manip_sd,
                      data = data4meta, measure = "SMD",
                      slab = Study_number, append=TRUE)
ma_model_1 <- rma(yi, vi, data = methane_data)
print(summary(ma_model_1))
forest(ma_model_1, slab = methane_data$Study_number)
title(main="Effect size for methane from all data")
funnel(ma_model_1)

#Forest plot for ecosystem type (Ecosystem_type)
Eco_type<- escalc(n1i = Control_N, n2i = Manip_N,
                      m1i = Control, m2i = Manip,
                      sd1i = Control_sd, sd2i = Manip_sd,
                      data = data4meta, measure = "SMD",
                      slab = Ecosystem_type, append=TRUE)

Eco_model <- rma.mv(yi, vi, data = Eco_type,
random = ~ 1 | Study_number,
mods = ~ Ecosystem_type - 1,
method = "REML")

Coef_meta_Eco <- coef(summary(Eco_model))

print(summary(Eco_model))
forest(Coef_meta_Eco$estimate,ci.lb = Coef_meta_Eco$ci.lb,
       ci.ub =  Coef_meta_Eco$ci.ub ,
       annotate = TRUE,slab=c("Agriculture","Forest","Pasture","Savanna"))
title(main="Ecosystem type")
funnel(Eco_model)


#Forest plot for manipulation effects (Manipulation_V2)
Manip_type<- escalc(n1i = Control_N, n2i = Manip_N,
            m1i = Control, m2i = Manip,
            sd1i = Control_sd, sd2i = Manip_sd,
            data = data4meta, measure = "SMD",
            slab = Manipulation_V2, append=TRUE)

Manip_model <- rma.mv(yi, vi, data = Manip_type,
                     random = ~ 1 | Study_number,
                     mods = ~ Manipulation_V2 - 1,
                     method = "REML")

Coef_meta_Manip <- coef(summary(Manip_model))

print(summary(Manip_model))
forest(Coef_meta_Manip$estimate,ci.lb = Coef_meta_Manip$ci.lb,
       ci.ub =  Coef_meta_Manip$ci.ub ,
       annotate = TRUE,slab=c("Burning","Fertilization","Fumigation","Tilled"))
title(main="Manipulation effects")
funnel(Manip_model)



#Forest plot for manipulation effects in agriculture soil

Ag_mp_data <- filter(data4meta, Ecosystem_type == 'Agriculture', Manipulation_V2 == 'Fertilization')
Ag_type<- escalc(n1i = Control_N, n2i = Manip_N,
            m1i = Control, m2i = Manip,
            sd1i = Control_sd, sd2i = Manip_sd,
            data = Ag_mp_data, measure = "SMD",
            slab = Ecosystem_type, append=TRUE)

ma_model_Ag <- rma.mv(yi, vi, data = Ag_type,
                     random = ~ 1 | Study_number,
                    #mods = ~ Manipulation_V2 - 1,
                     method = "REML")

Coef_meta_Ag <- coef(summary(ma_model_Ag))

print(summary(ma_model_Ag))


forest(Coef_meta_Ag$estimate,ci.lb = Coef_meta_Ag$ci.lb,
       ci.ub = Coef_meta_Ag$ci.ub , annotate = TRUE,slab=c("Agriculture"))


title(main="Manipulation effects in agriculture soil")
funnel(ma_model_Ag)


