set.seed(111)

# Preparing the data ------------------------------------------------------------

# Libraries
library(ggplot2)
library(lavaan)
library(psych)
library(bootnet)
library(qgraph)
library(tidyverse)
library(semTools)
library(responsePatterns)
library(dplyr)
library(lavaan.mi)


# Load data
df <- readxl::read_excel("VP+-+ostr%C3%A1+verze_November+30%2C+2023_12.49.xlsx")
names(df)

# Get PO scenarios content
prompts <- substr(df[1,12:29], 500, nchar(df[1,12:29]))
# Get short names
ids <- substr(prompts, 1, 5)

en_ids <- c("Vote", "Print", "Pray", "Kidney", "Kill", "Speed", "Sales",
            "Litter", "Release", "Hidjab", "Suicide", "Drink", "Dress", "Block",
            "Threat", "Roma", "Trans", "Disinfo")
ids <- en_ids
# Rename
names(df)[12:29] <- ids

# Remove qualtrics column
df <- df[-1,]

# Recode "don't knows" into NA
df[df == "6.0"] <- NA

# Recode likert responses to numeric
df[,12:78] <- lapply(df[,12:78], as.numeric)
df$age <- as.numeric(df$age)


# Remove unnecessary columns
df <- df[,-c(1:10)]

# Keep only people who consented
df <- df[df$consent == "1.0", ]

# Check NA patterns
Amelia::missmap(df, rank.order = FALSE)

# Flag those with < 70 % valid responses per block
df$flag <- ifelse(rowSums(!is.na(df[,2:19]))/ncol(df[,2:19]) < .7 | 
                    rowSums(!is.na(df[,20:31]))/ncol(df[,20:31]) < .7 |
                    rowSums(!is.na(df[,46:58]))/ncol(df[,46:58]) < .7, 1, 0)


table(df$flag)

df <- df[df$flag == 0,]

# Autocorrelation screening
flags <- data.frame(seq(1,nrow(df)))

k <- 5

rps_po <- responsePatterns::rp.acors(df[,2:19], max.lag = k)
responsePatterns::rp.hist(rps_po)
flags$po <- responsePatterns::rp.indices(rps_po)$percentile

rps_rwa <- responsePatterns::rp.acors(df[,20:31], max.lag = k)
responsePatterns::rp.hist(rps_rwa)
flags$rwa <- responsePatterns::rp.indices(rps_rwa)$percentile


rps_sdo <- responsePatterns::rp.acors(df[,32:46], max.lag = k)
responsePatterns::rp.hist(rps_sdo)
flags$sdo <- responsePatterns::rp.indices(rps_sdo)$percentile

rps_pop <- responsePatterns::rp.acors(df[,47:58], max.lag = k)
responsePatterns::rp.hist(rps_pop)
flags$pop <- responsePatterns::rp.indices(rps_pop)$percentile

flags <- flags %>%
  mutate(any = if_any(2:5, ~. > 98))

flags_df <- cbind(flags, df)
flags_df <- flags_df[flags_df$any == TRUE,]

# Straightliners: 104, 642, 222
# View(flags_df[,grep("^SDO", colnames(flags_df))])

# Straightliners: 104, 642
# View(flags_df[,grep("^RWA", colnames(flags_df))])

# Straightliners: 104, 642
# View(flags_df[,grep("^Schulz", colnames(flags_df))])

# View(flags_df[,en_ids])

# Create new df without flagged participants + straightliners
df_clean <- df[-c(104, 222, 642),-1]
Amelia::missmap(df_clean)

# Data description --------------------------------------------------------

psych::describe(df_clean)

# Remove students from lawyers (they get NA instead of 0 because they are halfway between)

df_clean$lawyer[df_clean$job == "9.0"] <- NA

psych::describe(as.numeric(df_clean$age))

round(prop.table(table(df_clean$job)),2)
round(prop.table(table(df_clean$uni)),2)

round(prop.table(table(df_clean$lay_edu)),2)

psych::describeBy(df_clean, df_clean$lawyer)
round(prop.table(table(df_clean$gender, df_clean$lawyer), margin = 2),2)
table(df_clean$gender, df_clean$lawyer)

# Creating indices --------------------------------------------------------

# Scale reliabilities, item analysis
psych::alpha(df_clean[,grep("^RWA", colnames(df_clean), value = TRUE)], check.keys = T)
psych::omega(df_clean[,grep("^RWA", colnames(df_clean), value = TRUE)], check.keys = T)

psych::polychoric(df_clean[,grep("^RWA", colnames(df_clean), value = TRUE)])

psych::alpha(df_clean[,grep("^SDO", colnames(df_clean), value = TRUE)], check.keys = T)
psych::omega(df_clean[,grep("^SDO", colnames(df_clean), value = TRUE)], check.keys = T)

psych::alpha(df_clean[,grep("^Schulz", colnames(df_clean), value = TRUE)], check.keys = T)
psych::omega(df_clean[,grep("^Schulz", colnames(df_clean), value = TRUE)], check.keys = T)

# RWA

rwa_model <- "RWA =~ RWA_1 + RWA_2 + RWA_3 + RWA_4 + RWA_5 + RWA_6 + RWA_7 + RWA_8 + RWA_9 + RWA_10 + RWA_11 + RWA_12
             
              RWA_1 ~~ RWA_2
              RWA_1 ~~ RWA_3
              RWA_2 ~~ RWA_3

              RWA_4 ~~ RWA_5
              RWA_4 ~~ RWA_6
              RWA_5 ~~ RWA_6
              
              RWA_7 ~~ RWA_8
              RWA_7 ~~ RWA_9
              RWA_8 ~~ RWA_9
              
              RWA_10 ~~ RWA_11
              RWA_10 ~~ RWA_12
              RWA_11 ~~ RWA_12
              
"

rwa_fit <- lavaan::cfa(df_clean, 
                       model = rwa_model, estimator = "MLR", ordered = F, 
                       std.lv = T, missing = "fiml")

summary(rwa_fit, fit = T, std = T)
residuals(rwa_fit)

semTools::nullRMSEA(rwa_fit)

## Invariance

rwa_fit_group <- lavaan::cfa(df_clean, 
                             model = rwa_model, estimator = "MLR", 
                             std.lv = T,
                             group = "lawyer", 
                             group.equal = "loadings", missing = "fiml")

summary(rwa_fit_group, fit = T, std = T)

# SDO
sdo_model <- "SDO =~ DOM + ANTI_EG
              DOM =~ SDO_1 + SDO_2 + SDO_3 + SDO_4 + SDO_5 + SDO_6 + SDO_7 + SDO_8  
              ANTI_EG =~ SDO_9 + SDO_10 + SDO_11 + SDO_12 + SDO_13 + SDO_14 + SDO_15

              SDO_6 ~~ SDO_8
              SDO_2 ~~ SDO_3"

sdo_fit <- lavaan::cfa(df_clean[,grep("^SDO", colnames(df_clean), value = TRUE)], 
                       model = sdo_model, estimator = "MLR", std.lv = T, 
                       missing = "fiml", ordered = F)

summary(sdo_fit, fit = T, std = T)
semTools::nullRMSEA(sdo_fit)

## Invariance

sdo_fit_group <- lavaan::cfa(df_clean, 
                             model = sdo_model, estimator = "MLR",
                             group = "lawyer", group.equal = "loadings", 
                             std.lv = T, missing = "fiml")

summary(sdo_fit_group, fit = T, std = T)

# POP
pop_model <- paste(paste0("POP =~", 
                          paste0(grep("^Schulz", colnames(df_clean), value = TRUE), collapse = " + "),
                          sep = " "))

pop_fit <- lavaan::cfa(df_clean, model = pop_model, 
                       estimator = "MLR", std.lv = T,
                       meanstructure = T, missing = "fiml")

summary(pop_fit, fit = T, std = T)

pop_model_3fac <- "
POP =~ 1*POP1 + POP2 + POP3
POP1 =~ Schulz_1 + Schulz_2 + Schulz_3 + Schulz_4
POP2 =~ Schulz_5 + Schulz_6 + Schulz_7 + Schulz_8
POP3 =~ Schulz_9 + Schulz_10 + Schulz_11 + Schulz_12"

pop_fit_3fac <- lavaan::cfa(df_clean, 
                            model = pop_model_3fac, 
                            estimator = "MLR", std.lv = T, missing = "fiml")

anova(pop_fit, pop_fit_3fac)
pop_comp <- semTools::compareFit(pop_fit, pop_fit_3fac)
summary(pop_comp)

summary(pop_fit_3fac, fit = T, std = T)
semTools::nullRMSEA(pop_fit_3fac)


## Invariance

pop_fit_group <- lavaan::cfa(df_clean, 
                             model = pop_model_3fac, estimator = "MLR",
                             group = "lawyer", group.equal = "loadings", 
                             std.lv = T, missing = "fiml")

summary(pop_fit_group, fit = T, std = T)


# Merging factor scores with data

add_scores <- function(model, data){
  idx <- lavInspect(model, "case.idx")
  fscores <- lavPredict(model)
  
  for (fs in colnames(fscores)) {
    data[idx, fs] <- fscores[ , fs]
  }
  return(data)
}

rwa_fit_list <- lavaan::cfa(df_clean, 
                            model = rwa_model, estimator = "MLR", ordered = F, 
                            std.lv = T)

sdo_fit_list <- lavaan::cfa(df_clean, 
                            model = sdo_model, estimator = "MLR", ordered = F, 
                            std.lv = T)

pop_fit_3fac_list <- lavaan::cfa(df_clean, 
                                 model = pop_model_3fac, estimator = "MLR", ordered = F,
                                 std.lv = T)

df_clean <- add_scores(rwa_fit, df_clean)
df_clean <- add_scores(sdo_fit, df_clean)
df_clean <- add_scores(pop_fit_3fac, df_clean)



add_plausible_values <- function(model, data, nDraws = 3) {
  # Step 1: Get row index for each case in the model
  idx <- lavInspect(model, "case.idx")
  
  # Step 2: Get plausible values (factor scores) with the specified number of draws
  fscores <- plausibleValues(model, nDraws = nDraws, method = "EBM")
  
  # Step 3: Add each plausible value draw as separate columns to the data
  for (i in seq_along(fscores)) {
    draw_name_suffix <- paste0("_draw", i)  # Suffix for each draw
    current_draw <- fscores[[i]]
    
    for (fs in colnames(current_draw)) {
      new_column_name <- paste0(fs, draw_name_suffix)  # Create unique column name for each draw
      data[idx, new_column_name] <- current_draw[, fs]  # Assign factor score values to the appropriate rows
    }
  }
  
  # Return the updated data with added plausible values
  return(data)
}




# H1: Variance comparison -------------------------------------------------

# Categorize the PO scenarios
core <- c("Vzdán", "Prode", "Úmysl", "Prodá",
          "Odhaz", "Pití ", "Odesl", "Cílen")

core <- c("Vote", "Kidney", "Kill", "Sales",
          "Litter", "Drink", "Threat", "Disinfo")

periphery <- c("Smlou", "Zaháj", "Překr", "Propu",
               "Veden", "Spách", "Navšt", "Zablo",
               "Nepro", "Změně")

periphery <- c("Print", "Pray", "Speed", "Release",
               "Hidjab", "Suicide", "Dress", "Block",
               "Roma", "Trans")

# Categorize them again for qgraph
colors <- c("core", "peri", "peri", 
            "core", "core", "peri",
            "core", "core", "peri", 
            "peri", "peri", "core",
            "peri", "peri", "core", 
            "peri", "peri", "core")

# Wide to long

psych::alpha(df_clean[df_clean$lawyer == "1.0",ids])

df_long <- df_clean

#df_long[,c("Dress", "Drink")] <- 6 - df_long[,c("Dress", "Drink")] 
df_long[,c("Sales", "Dress", "Drink", "Suicide")] <- 6 - df_long[,c("Sales", "Dress", "Drink", "Suicide")] 


df_long <- df_long[,c(ids, "lawyer")] %>%
  gather(key = "variable", value = "answer", core, periphery) %>%
  mutate(category = ifelse(variable %in% core, "core", "periphery")) %>%
  dplyr::select(lawyer, variable, answer, category)

df_long_law <- df_long[df_long$lawyer == "1.0",]

fligner.test(answer ~ category, df_long_law)

psych::describeBy(df_long_law$answer, df_long_law$category)

# Variance ratio
(1.39^2)/(.97^2)

psych::describeBy(df_long_law$answer, df_long_law$variable)

# H2-H3: Network models ---------------------------------------------------

pol_dims <- c("RWA", "SDO", "POP")

net_law <- bootnet::estimateNetwork(df_clean[df_clean$lawyer == "1.0", c(ids, pol_dims)], 
                                    default = "EBICglasso", corMethod = "spearman",
                                    tuning = .15)

