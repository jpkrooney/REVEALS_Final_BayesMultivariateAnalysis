#### Load necessary packages ####
library(tidyverse)
library(stringr)
library(tableone)
library(patchwork)
library(brms)
library(tidybayes)
library(bayesplot)
library(corrplot)
library(openxlsx)


#### perform set up tasks ####
# Settings for later plots
bayesplot_theme_set(theme_default() + theme(text = element_text(family="Arial")))

# Set brms options for multithreading
options(mc.cores = 8, brms.backend = "cmdstanr")

# Get support functions
source("R/support_fns.R")


##### Load data #####
dflong <- readRDS("../Final_dataset/REVEALS_finallongitudinaldata.rds")
dfbase <- readRDS("../Final_dataset/REVEALS_finalbaselinedata.rds")



##### Data preparation #####
dflong$uin <- as.factor(dflong$uin)
# extract unique IDS
id_list <- unique(dflong$uin)

# Determine number of longitudinal measurements per person
outcome_counts <- dflong %>% count(uin)
dfbase$count_long <- outcome_counts[ match(dfbase$uin, outcome_counts$uin), ]$n

# Exclude any individuals not at stage 2 or 3 at presentation
dfbase <- dfbase[ dfbase$staging %in% c("Kings 2", "Kings 3") | is.na(dfbase$staging), ]
dflong <- dflong[dflong$uin %in% dfbase$uin, ]

# Rename outcome variables with shorter names to make model specification more concise
dflong <- dflong %>% rename(outfvc = fvc_max_l,
                            outsvc = svc_l,
                            outsnip = snip_max_score_cm_h2o,
                            outpeak = peak_cough_flow_max_score)

# Determine number of individual outcomes per person
outfvc_count <- dflong[!is.na(dflong$outfvc), ] %>% count(uin)
outsvc_count <- dflong[!is.na(dflong$outsvc), ] %>% count(uin)
outsnip_count <- dflong[!is.na(dflong$outsnip), ] %>% count(uin)
outpeak_count <- dflong[!is.na(dflong$outpeak), ] %>% count(uin)

dfbase$count_fvc <- outfvc_count[ match(dfbase$uin, outfvc_count$uin), ]$n
dfbase$count_svc <- outsvc_count[ match(dfbase$uin, outsvc_count$uin), ]$n
dfbase$count_snip <- outsnip_count[ match(dfbase$uin, outsnip_count$uin), ]$n
dfbase$count_peak <- outpeak_count[ match(dfbase$uin, outpeak_count$uin), ]$n

#### assign baseline outcome measure to dfbase ####
firstlongs <- dflong %>% 
    group_by(uin) %>% 
    filter(days_from_baseline == first(days_from_baseline)) %>% 
    dplyr::select(uin, outfvc, fvc_max_percent_pred, outsvc, svc_percent_pred,
                  outpeak, outsnip, alsfrs_resp)

dfbase <- dfbase %>%
    left_join(firstlongs)



######## Comparisons of respiratory outcomes using REVEALS dataset #######

#### Make Table 1 of descriptive statistics ####
vars <- c("gender", "site_onset", "staging","height_cm", 
          "weight_kg", "age_dx", "dx_delay_months", "onset2baseline", "onset2censor_months",
          "count_long", "baseline2censor_months", "clinic_followup", "baseline2long_end",
          "outfvc", "fvc_max_percent_pred", "outsvc", "svc_percent_pred",
          "outpeak", "outsnip", "alsfrs_resp")
non_par <- c("onset2censor_months", "count_long", "baseline2censor_months", "dx_delay_months",
             "onset2baseline", "clinic_followup", "baseline2long_end", "alsfrs_resp")

tab1 <- CreateTableOne(vars, strata = "cohort", data = dfbase )
tab1b <- print(tab1, nonnormal = non_par, exact = c("gender", "site_onset", "staging"),
               quote = FALSE, noSpaces = TRUE, includeNA = TRUE,
               printToggle = FALSE, showAllLevels = TRUE)
tab1_total <- CreateTableOne(vars, data = dfbase)
tab1_total <- print(tab1_total, nonnormal = non_par, exact = c("gender", "site_onset", "staging"),
                    quote = FALSE, noSpaces = TRUE,
                    printToggle = FALSE, showAllLevels = TRUE)
table1 <- cbind(tab1b, tab1_total)
write.csv(table1, "Results/Table1_finalREVEALScohort.csv")

# How many in each cohort have only one measurement ?
table(dfbase$count_long, dfbase$cohort)


# Missing rate in outcome variables
dim(dflong)
summary(dflong[, c("outfvc", "outsvc", "outsnip", "outpeak")])

# How many timepoints are missing any outcome variable 
table(rowSums( !is.na(dflong[, c("outfvc", "outsvc", "outsnip",
                                 "outpeak")]) ) == 4)

# Missing values are < 5% in outcome variables. Therefore we can remove these rows without significant bias (only 2 rows are missing the time variable).
# Also will remove any missing site_onset and gender
df2 <- dflong[ rowSums( !is.na(dflong[, c("outfvc", "outsvc", "outsnip",
                                          "outpeak", "site_onset", "gender")]) ) == 6 ,  ]

# How many individuals were excluded?
length(unique(dflong$uin))
length(unique(df2$uin))


#### What is the raw correlation between outcomes ? ####
M <- cor(df2[, c("outfvc", "outsvc", "outpeak", "outsnip",
                 "alsfrs_resp")], use = "pairwise.complete.obs")
corrplot(M, method = "number", type = "full")
# What is the raw correlation between outcomes in spinal patients only ?
M_sp <- cor(df2[ df2$site_onset == "ALS Spinal onset",
                 c("outfvc", "outsvc", "outpeak", "outsnip",
                   "alsfrs_resp")], use = "pairwise.complete.obs")
corrplot(M_sp, method = "number", type = "full")
M_bulb <- cor(df2[ df2$site_onset == "ALS Bulbar onset",
                   c("outfvc", "outsvc", "outpeak", "outsnip",
                     "alsfrs_resp")], use = "pairwise.complete.obs")
corrplot(M_bulb, method = "number", type = "full")

rawCorrstable <- bind_rows(data.frame(cases = "all", outcome = rownames(M), M),
                           data.frame(cases = "spinal", outcome = rownames(M_sp),M_sp),
                           data.frame(cases = "bulbar", outcome = rownames(M_bulb),M_bulb))
rownames(rawCorrstable) <- NULL
write.csv(rawCorrstable, "Results/TableS1_rawCorrelations.csv", row.names = FALSE)

#### Bayesian Multiple outcomes longitudinal models ####
# Updated of models used here: https://www.tandfonline.com/doi/abs/10.1080/21678421.2021.1908362


#### First perform prior predictive check as per Bayesian workflow ####
##### Combined priors
priors_4outs <- c(set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
                  set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
                  set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
                  set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" ))

fit_priors <- brm(
    bf(mvbind(outfvc, outsvc, outsnip, outpeak) ~ days_from_baseline * sexsite +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2, chains=4,
    prior = c(priors_4outs,
              set_prior("normal(0, 0.1)", class = "b",
                        resp = c("outfvc", "outsvc","outsnip", "outpeak"))),
    sample_prior = "only", seed = 1234345,
    iter = 3500, threads = threading(2))

pp1 <- pp_check(fit_priors, resp = "outfvc", ndraws = 200) + ggtitle("FVC") + xlim(-100, 100)
pp2 <- pp_check(fit_priors, resp = "outsvc", ndraws = 200) + ggtitle("SVC") + xlim(-100, 100)
pp3 <- pp_check(fit_priors, resp = "outsnip", ndraws = 200) + ggtitle("SNIP") + xlim(-250, 500)
pp4 <- pp_check(fit_priors, resp = "outpeak", ndraws = 200) + ggtitle("PEAK") + xlim(-1000, 1500)
print(pp1 + pp2 + pp3 + pp4)

tiff("Graphs/ppc_newpriors.tiff", width=1000, height=700, res=108)
    print(pp1 + pp2 + pp3 + pp4)
dev.off()



##### Prior predictive checks are ok - fit model ####
fit_model <- brm(
    bf(mvbind(outfvc, outsvc, outsnip, outpeak) ~ days_from_baseline * sexsite +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_4outs,
    sample_prior = "no", seed = 1234345,
    iter = 3500, threads = threading(2))
saveRDS(fit_model, "Models/brms_finaldata_newmodel_3500.RDS")
#fit_model <- readRDS("Models/brms_finaldata_newmodel_3500.RDS")

# generate trace plots
mcmc_trace(fit_model, pars = c("b_outfvc_days_from_baseline:sexsiteMalebulbaronset",
                               "b_outsvc_days_from_baseline:sexsiteMalebulbaronset",
                               "b_outsnip_days_from_baseline:sexsiteMalebulbaronset",
                               "b_outpeak_days_from_baseline:sexsiteMalebulbaronset") )


# Perform posterior predictive check as per Bayesian workflow
mod_pp1 <- pp_check(fit_model, resp = "outfvc", ndraws = 1000) + ggtitle("FVC")
mod_pp2 <- pp_check(fit_model, resp = "outsvc", ndraws = 1000) + ggtitle("SVC")
mod_pp3 <- pp_check(fit_model, resp = "outsnip", ndraws = 1000) + ggtitle("SNIP")
mod_pp4 <- pp_check(fit_model, resp = "outpeak", ndraws = 1000) + ggtitle("PEAK")
print(mod_pp1 + mod_pp2 + mod_pp3 + mod_pp4)

tiff("Graphs/Final_model_retrodictive.tiff", width=1000, height=700, res=108)
  print(mod_pp1 + mod_pp2 + mod_pp3 + mod_pp4)
dev.off()

#extract model fixed effect coefficients and write to file
coefs_basemodel <- get_model_fixedcoefs(fit_model)
coefs_basemodel <- coefs_basemodel %>% 
    dplyr::select(variable,
                  Estimate_fvc_, Q2.5_fvc_, Q97.5_fvc_, Est.Error_fvc_,
                  Estimate_svc_, Q2.5_svc_, Q97.5_svc_, Est.Error_svc_,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak)
write.csv(coefs_basemodel, "Results/Finaldata_basemodel_fixed_effects.csv", row.names = FALSE)



# Plot of conditional effects
plotdata <- conditional_effects(fit_model, "days_from_baseline:sexsite")


# make new sex variable for plotting facets
plotdata$`outfvc.outfvc_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata$`outfvc.outfvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata$`outsvc.outsvc_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata$`outsvc.outsvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata$`outsnip.outsnip_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata$`outpeak.outpeak_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))


# make new site variable for plotting facets
plotdata$`outfvc.outfvc_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata$`outfvc.outfvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata$`outsvc.outsvc_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata$`outsvc.outsvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata$`outsnip.outsnip_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata$`outpeak.outpeak_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))


# Make plots
post1 <- ggplot(plotdata$`outfvc.outfvc_days_from_baseline:sexsite`,
                aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 4)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("FVC (L)") + facet_wrap(~sex) 
post2 <- ggplot(plotdata$`outsvc.outsvc_days_from_baseline:sexsite`,
                aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 4)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("SVC (L)") + facet_wrap(~sex) 
post3 <- ggplot(plotdata$`outsnip.outsnip_days_from_baseline:sexsite`,
                aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("SNIP (cmH2O )") + facet_wrap(~sex)               
post4 <- ggplot(plotdata$`outpeak.outpeak_days_from_baseline:sexsite`,
                aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("PEAK (L/min)") + facet_wrap(~sex) 

post1 + post2 + post3 + post4 + plot_layout(guides="collect")


tiff("Graphs/FigureS1_cond_effects_mainmodel.tiff", width=1000, height=700, res=108)
    print(post1 + post2 + post3 + post4 + plot_layout(guides="collect"))
dev.off()

# Breakdown by sex and site of onset
df2 %>% 
    group_by(uin) %>% 
    filter(row_number() == 1) %>% 
    group_by(sexsite) %>% 
    summarise(count = n())

#### Extract slopes of fixed effects from fit_model ####
# Make vector of wanted pars
wanted_pars <- c("outfvc_Intercept", "outsvc_Intercept", "outsnip_Intercept", 
                 "outpeak_Intercept", "outfvc_months_from_baseline", "outfvc_site_onsetALSBulbaronset", 
                 "outfvc_months_from_baseline:site_onsetALSBulbaronset", 
                 "outsvc_months_from_baseline", "outsvc_site_onsetALSBulbaronset", 
                 "outsvc_months_from_baseline:site_onsetALSBulbaronset", 
                 "outsnip_months_from_baseline", "outsnip_site_onsetALSBulbaronset", 
                 "outsnip_months_from_baseline:site_onsetALSBulbaronset", "outpeak_months_from_baseline", 
                 "outpeak_site_onsetALSBulbaronset",
                 "outpeak_months_from_baseline:site_onsetALSBulbaronset")

pop_effects_4out <- data.frame( t(fixef(fit_model) ) )

# Rename columns
pop_effects_4out <- pop_effects_4out %>% rename(
    "FVC_spinal_male_Intercept" = "outfvc_Intercept",
    "FVC_bulbar_male_Intercept_Offset" = "outfvc_sexsiteMalebulbaronset",
    "FVC_spinal_female_Intercept_Offset" ="outfvc_sexsiteFemalespinalonset",
    "FVC_bulbar_female_Intercept_Offset" ="outfvc_sexsiteFemalebulbaronset",
    "FVC_spinal_male_Slope" = "outfvc_days_from_baseline",
    "FVC_bulbar_male_Slope_Offset" = "outfvc_days_from_baseline.sexsiteMalebulbaronset",
    "FVC_spinal_female_Slope_Offset" = "outfvc_days_from_baseline.sexsiteFemalespinalonset",
    "FVC_bulbar_female_Slope_Offset" = "outfvc_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SVC_spinal_male_Intercept" = "outsvc_Intercept",
    "SVC_bulbar_male_Intercept_Offset" = "outsvc_sexsiteMalebulbaronset",
    "SVC_spinal_female_Intercept_Offset" = "outsvc_sexsiteFemalespinalonset",                 
    "SVC_bulbar_female_Intercept_Offset" = "outsvc_sexsiteFemalebulbaronset", 
    "SVC_spinal_male_Slope" = "outsvc_days_from_baseline",
    "SVC_bulbar_male_Slope_Offset" = "outsvc_days_from_baseline.sexsiteMalebulbaronset",
    "SVC_spinal_female_Slope_Offset" =  "outsvc_days_from_baseline.sexsiteFemalespinalonset",
    "SVC_bulbar_female_Slope_Offset" =  "outsvc_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SNIP_spinal_male_Intercept" = "outsnip_Intercept",
    "SNIP_bulbar_male_Intercept_Offset" =  "outsnip_sexsiteMalebulbaronset",
    "SNIP_spinal_female_Intercept_Offset"= "outsnip_sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Intercept_Offset"= "outsnip_sexsiteFemalebulbaronset",
    "SNIP_spinal_male_Slope" = "outsnip_days_from_baseline",
    "SNIP_bulbar_male_Slope_Offset" = "outsnip_days_from_baseline.sexsiteMalebulbaronset",
    "SNIP_spinal_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Intercept" = "outpeak_Intercept",
    "PEAK_bulbar_male_Intercept_Offset" = "outpeak_sexsiteMalebulbaronset",
    "PEAK_spinal_female_Intercept_Offset"= "outpeak_sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Intercept_Offset"= "outpeak_sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Slope" = "outpeak_days_from_baseline",
    "PEAK_bulbar_male_Slope_Offset" = "outpeak_days_from_baseline.sexsiteMalebulbaronset",
    "PEAK_spinal_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalebulbaronset")


# Calculated estimates
pop_effects_4out <- pop_effects_4out %>% mutate(
    FVC_bulbar_male_Intercept = FVC_spinal_male_Intercept + FVC_bulbar_male_Intercept_Offset,
    FVC_bulbar_male_Slope = FVC_spinal_male_Slope + FVC_bulbar_male_Slope_Offset,
    FVC_spinal_female_Intercept = FVC_spinal_male_Intercept +
        FVC_spinal_female_Intercept_Offset,
    FVC_bulbar_female_Intercept = FVC_spinal_male_Intercept +
        FVC_bulbar_female_Intercept_Offset,
    FVC_spinal_female_Slope = FVC_spinal_male_Slope + FVC_spinal_female_Slope_Offset,
    FVC_bulbar_female_Slope = FVC_spinal_male_Slope + FVC_bulbar_female_Slope_Offset,
    
    SVC_bulbar_male_Intercept = SVC_spinal_male_Intercept + SVC_bulbar_male_Intercept_Offset,
    SVC_bulbar_male_Slope = SVC_spinal_male_Slope + SVC_bulbar_male_Slope_Offset,
    SVC_spinal_female_Intercept = SVC_spinal_male_Intercept +
        SVC_spinal_female_Intercept_Offset,
    SVC_bulbar_female_Intercept = SVC_spinal_male_Intercept +
        SVC_bulbar_female_Intercept_Offset,
    SVC_spinal_female_Slope = SVC_spinal_male_Slope + SVC_spinal_female_Slope_Offset,
    SVC_bulbar_female_Slope = SVC_spinal_male_Slope + SVC_bulbar_female_Slope_Offset,
    
    SNIP_bulbar_male_Intercept = SNIP_spinal_male_Intercept + SNIP_bulbar_male_Intercept_Offset,
    SNIP_bulbar_male_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_male_Slope_Offset,
    SNIP_spinal_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_spinal_female_Intercept_Offset,
    SNIP_bulbar_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_bulbar_female_Intercept_Offset,
    SNIP_spinal_female_Slope = SNIP_spinal_male_Slope + SNIP_spinal_female_Slope_Offset,
    SNIP_bulbar_female_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_female_Slope_Offset,
    
    PEAK_bulbar_male_Intercept = PEAK_spinal_male_Intercept + PEAK_bulbar_male_Intercept_Offset,
    PEAK_bulbar_male_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_male_Slope_Offset,
    PEAK_spinal_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_spinal_female_Intercept_Offset,
    PEAK_bulbar_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_bulbar_female_Intercept_Offset,
    PEAK_spinal_female_Slope = PEAK_spinal_male_Slope + PEAK_spinal_female_Slope_Offset,
    PEAK_bulbar_female_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_female_Slope_Offset)
    
pop_effects_4out <- data.frame(par = names(pop_effects_4out), t(pop_effects_4out))

# Break reader friendly parameter names into multiple columns
metadata <- pop_effects_4out$par %>% str_split("_", simplify = TRUE) %>% data.frame()
names(metadata) <- names(metadata) <- c("outcome", "onset", "sex", "param", "offsetpar")
pop_effects_4out <- bind_cols(metadata, pop_effects_4out)
pop_effects_4out$par <- NULL
pop_effects_4out <- pop_effects_4out %>% filter( offsetpar!= "Offset")
pop_effects_4out$offsetpar <- NULL

# timescale of slopes currently slopes per day - convert to per average month
pop_effects_4out [ pop_effects_4out$param == "Slope", 5:8] <- (365.25 *
                                                         pop_effects_4out [ pop_effects_4out$param == "Slope", 5:8]) / 12
pop_effects_4out %>% arrange(outcome, param, onset)

# Write to file
write.csv(pop_effects_4out, "Results/TableS3_4outcomes_fixed_effects.csv", row.names = FALSE)



##### Now fit percent of predicted FVC and SVC version of model #####
# Needs new priors
priors_percent4outs <- c(set_prior("normal(75, 5)", class = "Intercept", resp = "fvcmaxpercentpred" ),
                    set_prior("normal(75, 5)", class = "Intercept", resp = "svcpercentpred" ),
                    set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
                    set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" ))
                                

# Prior predictive checks are ok - fit model #
fit_model_pct <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak) ~
           days_from_baseline * sexsite + cohort + (1 | p | uin) +
           (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2 %>% filter(!is.na(fvc_max_percent_pred) & !is.na(svc_percent_pred)), chains=4,
    prior = priors_percent4outs,
    sample_prior = "no", seed = 1234345,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_pct, "Models/brms_finaldata_newmodel_percentfvcsvc_3500.RDS")
#fit_model_pct <- readRDS("Models/brms_finaldata_newmodel_percentfvcsvc_3500.RDS")


mcmc_trace(fit_model_pct, pars = c("b_fvcmaxpercentpred_days_from_baseline:sexsiteMalebulbaronset",
                               "b_svcpercentpred_days_from_baseline:sexsiteMalebulbaronset",
                               "b_outsnip_days_from_baseline:sexsiteMalebulbaronset",
                               "b_outpeak_days_from_baseline:sexsiteMalebulbaronset") )


# Perform posterior predictive check as per Bayesian workflow
mod_pp1 <- pp_check(fit_model_pct, resp = "fvcmaxpercentpred", ndraws = 1000) + ggtitle("FVC %pred")
mod_pp2 <- pp_check(fit_model_pct, resp = "svcpercentpred", ndraws = 1000) + ggtitle("SVC %pred")
mod_pp3 <- pp_check(fit_model_pct, resp = "outsnip", ndraws = 1000) + ggtitle("SNIP %pred")
mod_pp4 <- pp_check(fit_model_pct, resp = "outpeak", ndraws = 1000) + ggtitle("PEAK %pred")
print(mod_pp1 + mod_pp2 + mod_pp3 + mod_pp4)


#extract model fixed effect coefficients and write to file
coefs_basemodel_pct <- get_model_fixedcoefsPCT(fit_model_pct)
write.csv(coefs_basemodel_pct, "Results/Finaldata_basemodel_fixed_effects_percent.csv", row.names = FALSE)


# Plot of conditional effects
plotdata_percent <- conditional_effects(fit_model_pct, "days_from_baseline:sexsite")


# make new sex variable for plotting facets
plotdata_percent$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_percent$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_percent$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_percent$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_percent$`outsnip.outsnip_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_percent$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_percent$`outpeak.outpeak_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_percent$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))


# make new site variable for plotting facets
plotdata_percent$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_percent$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_percent$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_percent$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_percent$`outsnip.outsnip_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_percent$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_percent$`outpeak.outpeak_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_percent$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))


post_pct1 <- ggplot(plotdata_percent$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`,
                aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("FVC (% predicted)") + facet_wrap(~sex)

post_pct2 <- ggplot(plotdata_percent$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`,
                    aes(x=days_from_baseline, y=estimate__, group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite, fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("SVC (% predicted)") + facet_wrap(~sex)

post_pct3 <- ggplot(plotdata_percent$outsnip.outsnip_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("SNIP (cmH2O )") + facet_wrap(~sex)
post_pct4 <- ggplot(plotdata_percent$outpeak.outpeak_days_from_baseline,
                aes(x=days_from_baseline, y=estimate__, col=site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, fill=site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("PEAK (L/min)") + facet_wrap(~sex)

post_pct1 + post_pct2 + post_pct3 + post_pct4 + plot_layout(guides="collect")

tiff("Graphs/Figure1_mainmodel_condeffects_PCT.tiff", width=1000, height=700, res=108)
    print(post_pct1 + post_pct2 + post_pct3 + post_pct4 + plot_layout(guides="collect"))
dev.off()



#### Extract slopes of fixed effects from percent predicted model ####
# Make vector of wanted pars

pop_effects_pct <- data.frame( t(fixef(fit_model_pct) ) )

## Rename columns
pop_effects_pct <- pop_effects_pct %>% rename(
    "FVCPCT_spinal_male_Intercept" = "fvcmaxpercentpred_Intercept",
    "FVCPCT_bulbar_male_Intercept_Offset" = "fvcmaxpercentpred_sexsiteMalebulbaronset",
    "FVCPCT_spinal_female_Intercept_Offset" ="fvcmaxpercentpred_sexsiteFemalespinalonset",
    "FVCPCT_bulbar_female_Intercept_Offset" ="fvcmaxpercentpred_sexsiteFemalebulbaronset",
    "FVCPCT_spinal_male_Slope" = "fvcmaxpercentpred_days_from_baseline",
    "FVCPCT_bulbar_male_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteMalebulbaronset",
    "FVCPCT_spinal_female_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteFemalespinalonset",
    "FVCPCT_bulbar_female_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SVCPCT_spinal_male_Intercept" = "svcpercentpred_Intercept",
    "SVCPCT_bulbar_male_Intercept_Offset" = "svcpercentpred_sexsiteMalebulbaronset",
    "SVCPCT_spinal_female_Intercept_Offset" = "svcpercentpred_sexsiteFemalespinalonset",                 
    "SVCPCT_bulbar_female_Intercept_Offset" = "svcpercentpred_sexsiteFemalebulbaronset", 
    "SVCPCT_spinal_male_Slope" = "svcpercentpred_days_from_baseline",
    "SVCPCT_bulbar_male_Slope_Offset" = "svcpercentpred_days_from_baseline.sexsiteMalebulbaronset",
    "SVCPCT_spinal_female_Slope_Offset" =  "svcpercentpred_days_from_baseline.sexsiteFemalespinalonset",
    "SVCPCT_bulbar_female_Slope_Offset" =  "svcpercentpred_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SNIP_spinal_male_Intercept" = "outsnip_Intercept",
    "SNIP_bulbar_male_Intercept_Offset" =  "outsnip_sexsiteMalebulbaronset",
    "SNIP_spinal_female_Intercept_Offset"= "outsnip_sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Intercept_Offset"= "outsnip_sexsiteFemalebulbaronset",
    "SNIP_spinal_male_Slope" = "outsnip_days_from_baseline",
    "SNIP_bulbar_male_Slope_Offset" = "outsnip_days_from_baseline.sexsiteMalebulbaronset",
    "SNIP_spinal_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Intercept" = "outpeak_Intercept",
    "PEAK_bulbar_male_Intercept_Offset" = "outpeak_sexsiteMalebulbaronset",
    "PEAK_spinal_female_Intercept_Offset"= "outpeak_sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Intercept_Offset"= "outpeak_sexsiteFemalebulbaronset",

    "PEAK_spinal_male_Slope" = "outpeak_days_from_baseline",
    "PEAK_bulbar_male_Slope_Offset" = "outpeak_days_from_baseline.sexsiteMalebulbaronset",
    "PEAK_spinal_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalebulbaronset"
)
    

# Calculated estimates
pop_effects_pct <- pop_effects_pct %>% mutate(
    FVCPCT_bulbar_male_Intercept = FVCPCT_spinal_male_Intercept + FVCPCT_bulbar_male_Intercept_Offset,
    FVCPCT_bulbar_male_Slope = FVCPCT_spinal_male_Slope + FVCPCT_bulbar_male_Slope_Offset,
    FVCPCT_spinal_female_Intercept = FVCPCT_spinal_male_Intercept +
                                            FVCPCT_spinal_female_Intercept_Offset,
    FVCPCT_bulbar_female_Intercept = FVCPCT_spinal_male_Intercept +
                                            FVCPCT_bulbar_female_Intercept_Offset,
    FVCPCT_spinal_female_Slope = FVCPCT_spinal_male_Slope + FVCPCT_spinal_female_Slope_Offset,
    FVCPCT_bulbar_female_Slope = FVCPCT_spinal_male_Slope + FVCPCT_bulbar_female_Slope_Offset,
    
    SVCPCT_bulbar_male_Intercept = SVCPCT_spinal_male_Intercept + SVCPCT_bulbar_male_Intercept_Offset,
    SVCPCT_bulbar_male_Slope = SVCPCT_spinal_male_Slope + SVCPCT_bulbar_male_Slope_Offset,
    SVCPCT_spinal_female_Intercept = SVCPCT_spinal_male_Intercept +
        SVCPCT_spinal_female_Intercept_Offset,
    SVCPCT_bulbar_female_Intercept = SVCPCT_spinal_male_Intercept +
        SVCPCT_bulbar_female_Intercept_Offset,
    SVCPCT_spinal_female_Slope = SVCPCT_spinal_male_Slope + SVCPCT_spinal_female_Slope_Offset,
    SVCPCT_bulbar_female_Slope = SVCPCT_spinal_male_Slope + SVCPCT_bulbar_female_Slope_Offset,
    
    SNIP_bulbar_male_Intercept = SNIP_spinal_male_Intercept + SNIP_bulbar_male_Intercept_Offset,
    SNIP_bulbar_male_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_male_Slope_Offset,
    SNIP_spinal_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_spinal_female_Intercept_Offset,
    SNIP_bulbar_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_bulbar_female_Intercept_Offset,
    SNIP_spinal_female_Slope = SNIP_spinal_male_Slope + SNIP_spinal_female_Slope_Offset,
    SNIP_bulbar_female_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_female_Slope_Offset,
    
    PEAK_bulbar_male_Intercept = PEAK_spinal_male_Intercept + PEAK_bulbar_male_Intercept_Offset,
    PEAK_bulbar_male_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_male_Slope_Offset,
    PEAK_spinal_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_spinal_female_Intercept_Offset,
    PEAK_bulbar_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_bulbar_female_Intercept_Offset,
    PEAK_spinal_female_Slope = PEAK_spinal_male_Slope + PEAK_spinal_female_Slope_Offset,
    PEAK_bulbar_female_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_female_Slope_Offset
)


pop_effects_pct <- data.frame(par = names(pop_effects_pct), t(pop_effects_pct))

## Break reader friendly parameter names into multiple columns
metadata_pct <- pop_effects_pct$par %>% str_split("_", simplify = TRUE) %>% data.frame()
names(metadata_pct) <- c("outcome", "onset", "sex", "param", "offsetpar")
pop_effects_pct <- bind_cols(metadata_pct, pop_effects_pct)
pop_effects_pct$par <- NULL
pop_effects_pct <- pop_effects_pct %>% filter( offsetpar!= "Offset")
pop_effects_pct$offsetpar <- NULL

# timescale of slopes currently slopes per day - convert to per average month
pop_effects_pct[ pop_effects_pct$param == "Slope", 5:8] <- (365.25 *
                            pop_effects_pct[ pop_effects_pct$param == "Slope", 5:8] / 12)

# remove any params not Intercept or slope
pop_effects_pct <- pop_effects_pct %>% filter( param %in% c("Intercept", "Slope"))

# Write to file
write.csv(pop_effects_pct, "Results/Table2_4outcomesPCT_FixefEff.csv", row.names = FALSE)



##### Next will extend the model to include ALSFRS-Resp as a 5th outcome #####
# set priors for 5 outcome model
priors_5outs <- c(
    set_prior("normal(3.2, 0.1)", class = "Intercept", resp = "outfvc" ),
    set_prior("normal(3.1, 0.2)", class = "Intercept", resp = "outsvc" ),
    set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
    set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" ),
    set_prior("normal(10.9, 2)", class = "Intercept", resp = "alsfrsresp" ))

# fit prior model for prior predictive check
fit_priors5outs <- brm(
    bf(mvbind(outfvc, outsvc, outsnip, outpeak, alsfrs_resp) ~ days_from_baseline * sexsite +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2, chains=4,
    prior = c(priors_5outs,
              set_prior("normal(0, 0.1)", class = "b",
                        resp = c("outfvc", "outsvc","outsnip", "outpeak", "alsfrsresp"))),
    sample_prior = "only", seed = 1234345,
    iter = 3500, threads = threading(2))

pp1 <- pp_check(fit_priors5outs, resp = "outfvc", ndraws = 200) + ggtitle("FVC") + xlim(-100, 100)
pp2 <- pp_check(fit_priors5outs, resp = "outsvc", ndraws = 200) + ggtitle("SVC") + xlim(-100, 100)
pp3 <- pp_check(fit_priors5outs, resp = "outsnip", ndraws = 200) + ggtitle("SNIP") + xlim(-250, 500)
pp4 <- pp_check(fit_priors5outs, resp = "outpeak", ndraws = 200) + ggtitle("PEAK") + xlim(-1000, 1500)
pp5 <- pp_check(fit_priors5outs, resp = "alsfrsresp", ndraws = 200) + ggtitle("ALSFRS_resp") + xlim(-15, 15)

print(pp1 + pp2 + pp3 + pp4 + pp5)

tiff("Graphs/ppc_newpriors5outs.tiff", width=1000, height=700, res=108)
    print(pp1 + pp2 + pp3 + pp4 + pp5)
dev.off()



##### Prior predictive checks are ok - fit 5 outcome model #####
fit_model_alsfrsresp_gaussian <- brm(
    bf(mvbind(outfvc, outsvc, outsnip, outpeak, alsfrs_resp) ~ days_from_baseline * sexsite +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outs,
    sample_prior = "no", seed = 1234345,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_alsfrsresp_gaussian, "Models/brms_finaldata_newmodel_3500_5outs_gaussian.RDS")
#fit_model_alsfrsresp_gaussian <- readRDS("Models/brms_finaldata_newmodel_3500_5outs_gaussian.RDS")

# 
retropred1_alsfrsresp_gauss <- pp_check(fit_model_alsfrsresp_gaussian,
                                        resp = "outfvc", ndraws = 1000) + ggtitle("FVC")
retropred2_alsfrsresp_gauss <- pp_check(fit_model_alsfrsresp_gaussian,
                                        resp = "outsvc", ndraws = 1000) + ggtitle("SVC")
retropred3_alsfrsresp_gauss <- pp_check(fit_model_alsfrsresp_gaussian,
                                        resp = "outsnip", ndraws = 1000) + ggtitle("SNIP")
retropred4_alsfrsresp_gauss <- pp_check(fit_model_alsfrsresp_gaussian,
                                        resp = "outpeak", ndraws = 1000) + ggtitle("PEAK")
retropred5_alsfrsresp_gauss <- pp_check(fit_model_alsfrsresp_gaussian, resp = "alsfrsresp",
                                  ndraws = 1000) + ggtitle("ALSFRS-R resp")

retropred1_alsfrsresp_gauss + retropred2_alsfrsresp_gauss + retropred3_alsfrsresp_gauss +
    retropred4_alsfrsresp_gauss + retropred5_alsfrsresp_gauss


# Plot of conditional effects
plotdata_alsfrs_gauss <- conditional_effects(fit_model_alsfrsresp_gaussian, "days_from_baseline:sexsite")


# make new sex variable for plotting facets
plotdata_alsfrs_gauss$`outfvc.outfvc_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gauss$`outfvc.outfvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gauss$`outsvc.outsvc_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gauss$`outsvc.outsvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gauss$`outsnip.outsnip_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gauss$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gauss$`outpeak.outpeak_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gauss$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gauss$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite` $sex <- factor(
    str_split(plotdata_alsfrs_gauss$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))


# make new site variable for plotting facets
plotdata_alsfrs_gauss$`outfvc.outfvc_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gauss$`outfvc.outfvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gauss$`outsvc.outsvc_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gauss$`outsvc.outsvc_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gauss$`outsnip.outsnip_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gauss$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gauss$`outpeak.outpeak_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gauss$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gauss$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gauss$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))



post_alsfrsresp_gauss1 <- ggplot(plotdata_alsfrs_gauss$`outfvc.outfvc_days_from_baseline:sexsite`,
                           aes(x=days_from_baseline, y=estimate__,
                               group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("FVC (L)") + facet_wrap(~sex)
post_alsfrsresp_gauss2 <- ggplot(plotdata_alsfrs_gauss$`outsvc.outsvc_days_from_baseline:sexsite`,
                           aes(x=days_from_baseline, y=estimate__,
                               group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("SVC (L)") + facet_wrap(~sex)
post_alsfrsresp_gauss3 <- ggplot(plotdata_alsfrs_gauss$`outsnip.outsnip_days_from_baseline:sexsite`,
                           aes(x=days_from_baseline, y=estimate__,
                               group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("SNIP (cmH2O )") + facet_wrap(~sex)
post_alsfrsresp_gauss4 <- ggplot(plotdata_alsfrs_gauss$`outpeak.outpeak_days_from_baseline:sexsite`,
                           aes(x=days_from_baseline, y=estimate__,
                               group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("PEAK (L/min)") + facet_wrap(~sex)
post_alsfrsresp_gauss5 <- ggplot(plotdata_alsfrs_gauss$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`,
                           aes(x=days_from_baseline, y=estimate__,
                               group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("ALSFRS-R Resp") + facet_wrap(~sex)

post_alsfrsresp_gauss1 + post_alsfrsresp_gauss2 + post_alsfrsresp_gauss5 +
    post_alsfrsresp_gauss3 + post_alsfrsresp_gauss4  + plot_layout(guides="collect") &
    theme(legend.position='bottom')

tiff("Graphs/FigureS2_extra_outcome_gauss.tiff", width=1000, height=700, res=108)
    post_alsfrsresp_gauss1 + post_alsfrsresp_gauss2 + post_alsfrsresp_gauss5 +
        post_alsfrsresp_gauss3 + post_alsfrsresp_gauss4  + plot_layout(guides="collect") &
    theme(legend.position='bottom')
dev.off()

jpeg("Graphs/alsfrs_gauss.jpg", width=1000, height=700, res=108)
    post_alsfrsresp_gauss5
dev.off()




#### Extract posterior predictions from 5 outcome model at set follow up times ####

# First create new dataset on which to perform predictions
new_data <- dfbase[ , c("uin", "cohort", "sexsite")]
days_from_baseline <- c(0, 183, 365, 548, 730) # 0, 6, 12, 18, 24 months
new_data <- expand_grid(new_data, days_from_baseline)

# Extract posterior predictions
post_samples <- posterior_predict(fit_model_alsfrsresp_gaussian, new_data,
                                  ndraws = 2000, allow_new_levels = TRUE)
n_samples <- dim(post_samples)[1]

# Create a list of datasets for each sample
post_datasets <- lapply(1:dim(post_samples)[1], function(i) {
    data.frame(new_data, post_samples[ i, ,])
})
#dump <- lapply(1:dim(post_samples)[1], function(i) {
#    post_datasets[[i]]$invertalsfrsresp <<- 12 - post_datasets[[i]]$invertalsfrsresp
#})


# calculate correlation of outcomes overall
cor_all <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))
apply(simplify2array(cor_all), 1:2, mean)


# calculate correlation of outcomes for fixed followup times at 0, 6, 12, 18 and 24 months
idx_time0 <- which( post_datasets[[1]]$days_from_baseline ==0)
cor_time0 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time0, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))

idx_time183 <- which( post_datasets[[1]]$days_from_baseline ==183)
cor_time183 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time183, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))

idx_time365 <- which( post_datasets[[1]]$days_from_baseline ==365)
cor_time365 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time365, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))

idx_time548 <- which( post_datasets[[1]]$days_from_baseline ==548)
cor_time548 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time548, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))

idx_time730 <- which( post_datasets[[1]]$days_from_baseline ==730)
cor_time730 <- lapply(1:n_samples, function(i)
    cor(post_datasets[[i]][idx_time730, c("outfvc", "outsvc", "outsnip", "outpeak", "alsfrsresp")]  ))


# Get mean correlations across posterior predictions at each follow-up time
mean_cor_time0 <- data.frame( apply(simplify2array(cor_time0), 1:2, mean) )
mean_cor_time183 <- data.frame( apply(simplify2array(cor_time183), 1:2, mean))
mean_cor_time365 <- data.frame( apply(simplify2array(cor_time365), 1:2, mean))
mean_cor_time548 <- data.frame( apply(simplify2array(cor_time548), 1:2, mean))
mean_cor_time730 <- data.frame( apply(simplify2array(cor_time730), 1:2, mean))

mean_cor_time0 <- data.frame(outcome = row.names(mean_cor_time0),
                             days = 0, mean_cor_time0)
mean_cor_time183 <- data.frame(outcome = row.names(mean_cor_time183),
                               days = 183, mean_cor_time183)
mean_cor_time365 <- data.frame(outcome = row.names(mean_cor_time365),
                               days = 365, mean_cor_time365)
mean_cor_time548 <- data.frame(outcome = row.names(mean_cor_time548),
                               days = 548, mean_cor_time548)
mean_cor_time730 <- data.frame(outcome = row.names(mean_cor_time730),
                               days = 730, mean_cor_time730)

write.xlsx(bind_rows( mean_cor_time0,
                      mean_cor_time183,
                      mean_cor_time365,
                      mean_cor_time548,
                      mean_cor_time730),
           "Results/TableS2_5outcomes_modelled_corrs.xlsx")


#### Extract slopes of fixed effects from 5 outcomes model ####
# Make vector of wanted pars
pop_effects_5out <- data.frame( t(fixef(fit_model_alsfrsresp_gaussian) ) )

## Rename columns
pop_effects_5out <- pop_effects_5out %>% rename(
    "FVC_spinal_male_Intercept" = "outfvc_Intercept",
    "FVC_bulbar_male_Intercept_Offset" = "outfvc_sexsiteMalebulbaronset",
    "FVC_spinal_female_Intercept_Offset" ="outfvc_sexsiteFemalespinalonset",
    "FVC_bulbar_female_Intercept_Offset" ="outfvc_sexsiteFemalebulbaronset",
    "FVC_spinal_male_Slope" = "outfvc_days_from_baseline",
    "FVC_bulbar_male_Slope_Offset" = "outfvc_days_from_baseline.sexsiteMalebulbaronset",
    "FVC_spinal_female_Slope_Offset" = "outfvc_days_from_baseline.sexsiteFemalespinalonset",
    "FVC_bulbar_female_Slope_Offset" = "outfvc_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SVC_spinal_male_Intercept" = "outsvc_Intercept",
    "SVC_bulbar_male_Intercept_Offset" = "outsvc_sexsiteMalebulbaronset",
    "SVC_spinal_female_Intercept_Offset" = "outsvc_sexsiteFemalespinalonset",                 
    "SVC_bulbar_female_Intercept_Offset" = "outsvc_sexsiteFemalebulbaronset", 
    "SVC_spinal_male_Slope" = "outsvc_days_from_baseline",
    "SVC_bulbar_male_Slope_Offset" = "outsvc_days_from_baseline.sexsiteMalebulbaronset",
    "SVC_spinal_female_Slope_Offset" =  "outsvc_days_from_baseline.sexsiteFemalespinalonset",
    "SVC_bulbar_female_Slope_Offset" =  "outsvc_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SNIP_spinal_male_Intercept" = "outsnip_Intercept",
    "SNIP_bulbar_male_Intercept_Offset" =  "outsnip_sexsiteMalebulbaronset",
    "SNIP_spinal_female_Intercept_Offset"= "outsnip_sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Intercept_Offset"= "outsnip_sexsiteFemalebulbaronset",
    "SNIP_spinal_male_Slope" = "outsnip_days_from_baseline",
    "SNIP_bulbar_male_Slope_Offset" = "outsnip_days_from_baseline.sexsiteMalebulbaronset",
    "SNIP_spinal_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Intercept" = "outpeak_Intercept",
    "PEAK_bulbar_male_Intercept_Offset" = "outpeak_sexsiteMalebulbaronset",
    "PEAK_spinal_female_Intercept_Offset"= "outpeak_sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Intercept_Offset"= "outpeak_sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Slope" = "outpeak_days_from_baseline",
    "PEAK_bulbar_male_Slope_Offset" = "outpeak_days_from_baseline.sexsiteMalebulbaronset",
    "PEAK_spinal_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalebulbaronset",
    
    "ALSFRSRESP_spinal_male_Intercept" = "alsfrsresp_Intercept",
    "ALSFRSRESP_bulbar_male_Intercept_Offset" = "alsfrsresp_sexsiteMalebulbaronset",
    "ALSFRSRESP_spinal_female_Intercept_Offset"= "alsfrsresp_sexsiteFemalespinalonset",
    "ALSFRSRESP_bulbar_female_Intercept_Offset"= "alsfrsresp_sexsiteFemalebulbaronset",
    
    "ALSFRSRESP_spinal_male_Slope" = "alsfrsresp_days_from_baseline",
    "ALSFRSRESP_bulbar_male_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteMalebulbaronset",
    "ALSFRSRESP_spinal_female_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteFemalespinalonset",
    "ALSFRSRESP_bulbar_female_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteFemalebulbaronset"
)


# Calculated estimates
pop_effects_5out <- pop_effects_5out %>% mutate(
    FVC_bulbar_male_Intercept = FVC_spinal_male_Intercept + FVC_bulbar_male_Intercept_Offset,
    FVC_bulbar_male_Slope = FVC_spinal_male_Slope + FVC_bulbar_male_Slope_Offset,
    FVC_spinal_female_Intercept = FVC_spinal_male_Intercept +
        FVC_spinal_female_Intercept_Offset,
    FVC_bulbar_female_Intercept = FVC_spinal_male_Intercept +
        FVC_bulbar_female_Intercept_Offset,
    FVC_spinal_female_Slope = FVC_spinal_male_Slope + FVC_spinal_female_Slope_Offset,
    FVC_bulbar_female_Slope = FVC_spinal_male_Slope + FVC_bulbar_female_Slope_Offset,
    
    SVC_bulbar_male_Intercept = SVC_spinal_male_Intercept + SVC_bulbar_male_Intercept_Offset,
    SVC_bulbar_male_Slope = SVC_spinal_male_Slope + SVC_bulbar_male_Slope_Offset,
    SVC_spinal_female_Intercept = SVC_spinal_male_Intercept +
        SVC_spinal_female_Intercept_Offset,
    SVC_bulbar_female_Intercept = SVC_spinal_male_Intercept +
        SVC_bulbar_female_Intercept_Offset,
    SVC_spinal_female_Slope = SVC_spinal_male_Slope + SVC_spinal_female_Slope_Offset,
    SVC_bulbar_female_Slope = SVC_spinal_male_Slope + SVC_bulbar_female_Slope_Offset,
    
    SNIP_bulbar_male_Intercept = SNIP_spinal_male_Intercept + SNIP_bulbar_male_Intercept_Offset,
    SNIP_bulbar_male_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_male_Slope_Offset,
    SNIP_spinal_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_spinal_female_Intercept_Offset,
    SNIP_bulbar_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_bulbar_female_Intercept_Offset,
    SNIP_spinal_female_Slope = SNIP_spinal_male_Slope + SNIP_spinal_female_Slope_Offset,
    SNIP_bulbar_female_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_female_Slope_Offset,
    
    PEAK_bulbar_male_Intercept = PEAK_spinal_male_Intercept + PEAK_bulbar_male_Intercept_Offset,
    PEAK_bulbar_male_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_male_Slope_Offset,
    PEAK_spinal_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_spinal_female_Intercept_Offset,
    PEAK_bulbar_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_bulbar_female_Intercept_Offset,
    PEAK_spinal_female_Slope = PEAK_spinal_male_Slope + PEAK_spinal_female_Slope_Offset,
    PEAK_bulbar_female_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_female_Slope_Offset,
    
    ALSFRSRESP_bulbar_male_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_bulbar_male_Intercept_Offset,
    ALSFRSRESP_bulbar_male_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_bulbar_male_Slope_Offset,
    ALSFRSRESP_spinal_female_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_spinal_female_Intercept_Offset,
    ALSFRSRESP_bulbar_female_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_bulbar_female_Intercept_Offset,
    ALSFRSRESP_spinal_female_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_spinal_female_Slope_Offset,
    ALSFRSRESP_bulbar_female_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_bulbar_female_Slope_Offset
)


pop_effects_5out <- data.frame(par = names(pop_effects_5out), t(pop_effects_5out))

## Break reader friendly parameter names into multiple columns
metadata_5out <- pop_effects_5out$par %>% str_split("_", simplify = TRUE) %>% data.frame()
names(metadata_5out) <- c("outcome", "onset", "sex", "param", "offsetpar")
pop_effects_5out <- bind_cols(metadata_5out, pop_effects_5out)
pop_effects_5out$par <- NULL
pop_effects_5out <- pop_effects_5out %>% filter( offsetpar!= "Offset")
pop_effects_5out$offsetpar <- NULL

# timescale of slopes currently slopes per day - convert to per average month
pop_effects_5out[ pop_effects_5out$param == "Slope", 5:8] <- (365.25 *
                                    pop_effects_5out[ pop_effects_5out$param == "Slope", 5:8] / 12)
pop_effects_5out %>% arrange(outcome, param, onset)

# remove any params not Intercept or slope
pop_effects_5out <- pop_effects_5out %>% filter( param %in% c("Intercept", "Slope"))

# Write to file
write.csv(pop_effects_5out, "Results/TableS4_5outcomeFixedEff.csv", row.names = FALSE)






##### Fit 5 outcomes model using percent predicted of FVC and SVC #####
priors_5outsPCT <- c(
    set_prior("normal(75, 5)", class = "Intercept", resp = "fvcmaxpercentpred" ),
    set_prior("normal(75, 5)", class = "Intercept", resp = "svcpercentpred" ),
    set_prior("normal(65, 20)", class = "Intercept", resp = "outsnip" ),
    set_prior("normal(350, 100)", class = "Intercept", resp = "outpeak" ),
    set_prior("normal(10.9, 2)", class = "Intercept", resp = "alsfrsresp" ))


fit_model_alsfrsresp_gaussianPCT <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~ days_from_baseline * sexsite + cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) + set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 1234345,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_alsfrsresp_gaussianPCT, "Models/brms_finaldata_newmodel_3500_5outs_gaussianPCT.RDS")
#fit_model_alsfrsresp_gaussianPCT <- readRDS("Models/brms_finaldata_newmodel_3500_5outs_gaussianPCT.RDS")



# Plot of conditional effects
plotdata_alsfrs_gaussPCT <- conditional_effects(fit_model_alsfrsresp_gaussianPCT, "days_from_baseline:sexsite")


# make new sex variable for plotting facets
plotdata_alsfrs_gaussPCT$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gaussPCT$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gaussPCT$`outsnip.outsnip_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gaussPCT$`outpeak.outpeak_days_from_baseline:sexsite`$sex <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))
plotdata_alsfrs_gaussPCT$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite` $sex <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[,1], levels = c("Male", "Female"))


# make new site variable for plotting facets
plotdata_alsfrs_gaussPCT$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gaussPCT$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gaussPCT$`outsnip.outsnip_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`outsnip.outsnip_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gaussPCT$`outpeak.outpeak_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`outpeak.outpeak_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))
plotdata_alsfrs_gaussPCT$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$site_onset <- factor(
    str_split(plotdata_alsfrs_gaussPCT$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`$sexsite,
              " ", simplify = TRUE)[ ,2], levels = c("spinal", "bulbar"))



post_alsfrsresp_gauss1PCT <- ggplot(plotdata_alsfrs_gaussPCT$`fvcmaxpercentpred.fvcmaxpercentpred_days_from_baseline:sexsite`,
                                    aes(x=days_from_baseline, y=estimate__,
                                        group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("FVC (% predicted)") + facet_wrap(~sex)


post_alsfrsresp_gauss2PCT <- ggplot(plotdata_alsfrs_gaussPCT$`svcpercentpred.svcpercentpred_days_from_baseline:sexsite`,
                                    aes(x=days_from_baseline, y=estimate__,
                                        group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 100)) + theme_minimal() + theme(panel.spacing = unit(2, "lines")) +
    labs(color='Site of Onset', fill='Site of Onset') + xlab("days from baseline") +
    ylab("SVC (% predicted)") + facet_wrap(~sex)

post_alsfrsresp_gauss3PCT <- ggplot(plotdata_alsfrs_gaussPCT$`outsnip.outsnip_days_from_baseline:sexsite`,
                                    aes(x=days_from_baseline, y=estimate__,
                                        group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 75)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("SNIP (cmH2O )") + facet_wrap(~sex)

post_alsfrsresp_gauss4PCT <- ggplot(plotdata_alsfrs_gaussPCT$`outpeak.outpeak_days_from_baseline:sexsite`,
                                    aes(x=days_from_baseline, y=estimate__,
                                        group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    coord_cartesian(ylim = c(0, 450)) + theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("PEAK (L/min)") + facet_wrap(~sex)

post_alsfrsresp_gauss5PCT <- ggplot(plotdata_alsfrs_gaussPCT$`alsfrsresp.alsfrsresp_days_from_baseline:sexsite`,
                                    aes(x=days_from_baseline, y=estimate__,
                                        group = sexsite, col = site_onset)) + geom_line() +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, linetype=NA, group = sexsite,
                    fill = site_onset), alpha=0.15) +
    theme_minimal() + labs(color='Site of Onset', fill='Site of Onset') +
    xlab("days from baseline") + ylab("ALSFRS-R Resp") + facet_wrap(~sex)

post_alsfrsresp_gauss1PCT + post_alsfrsresp_gauss2PCT + post_alsfrsresp_gauss5PCT +
    post_alsfrsresp_gauss3PCT + post_alsfrsresp_gauss4PCT  + plot_layout(guides="collect") &
    theme(legend.position='bottom')

tiff("Graphs/Figure2_5outs_gaussPCT.tiff", width=1000, height=700, res=108)
    post_alsfrsresp_gauss1PCT + post_alsfrsresp_gauss2PCT + post_alsfrsresp_gauss5PCT +
        post_alsfrsresp_gauss3PCT + post_alsfrsresp_gauss4PCT  + plot_layout(guides="collect") &
        theme(legend.position='bottom')
dev.off()



pop_effects_5outPCT <- data.frame( t(fixef(fit_model_alsfrsresp_gaussianPCT) ) )

## Rename columns
pop_effects_5outPCT <- pop_effects_5outPCT %>% rename(
    "FVCPCT_spinal_male_Intercept" = "fvcmaxpercentpred_Intercept",
    "FVCPCT_bulbar_male_Intercept_Offset" = "fvcmaxpercentpred_sexsiteMalebulbaronset",
    "FVCPCT_spinal_female_Intercept_Offset" ="fvcmaxpercentpred_sexsiteFemalespinalonset",
    "FVCPCT_bulbar_female_Intercept_Offset" ="fvcmaxpercentpred_sexsiteFemalebulbaronset",
    "FVCPCT_spinal_male_Slope" = "fvcmaxpercentpred_days_from_baseline",
    "FVCPCT_bulbar_male_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteMalebulbaronset",
    "FVCPCT_spinal_female_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteFemalespinalonset",
    "FVCPCT_bulbar_female_Slope_Offset" = "fvcmaxpercentpred_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SVCPCT_spinal_male_Intercept" = "svcpercentpred_Intercept",
    "SVCPCT_bulbar_male_Intercept_Offset" = "svcpercentpred_sexsiteMalebulbaronset",
    "SVCPCT_spinal_female_Intercept_Offset" = "svcpercentpred_sexsiteFemalespinalonset",                 
    "SVCPCT_bulbar_female_Intercept_Offset" = "svcpercentpred_sexsiteFemalebulbaronset", 
    "SVCPCT_spinal_male_Slope" = "svcpercentpred_days_from_baseline",
    "SVCPCT_bulbar_male_Slope_Offset" = "svcpercentpred_days_from_baseline.sexsiteMalebulbaronset",
    "SVCPCT_spinal_female_Slope_Offset" =  "svcpercentpred_days_from_baseline.sexsiteFemalespinalonset",
    "SVCPCT_bulbar_female_Slope_Offset" =  "svcpercentpred_days_from_baseline.sexsiteFemalebulbaronset",
    
    "SNIP_spinal_male_Intercept" = "outsnip_Intercept",
    "SNIP_bulbar_male_Intercept_Offset" =  "outsnip_sexsiteMalebulbaronset",
    "SNIP_spinal_female_Intercept_Offset"= "outsnip_sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Intercept_Offset"= "outsnip_sexsiteFemalebulbaronset",
    "SNIP_spinal_male_Slope" = "outsnip_days_from_baseline",
    "SNIP_bulbar_male_Slope_Offset" = "outsnip_days_from_baseline.sexsiteMalebulbaronset",
    "SNIP_spinal_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalespinalonset",
    "SNIP_bulbar_female_Slope_Offset" = "outsnip_days_from_baseline.sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Intercept" = "outpeak_Intercept",
    "PEAK_bulbar_male_Intercept_Offset" = "outpeak_sexsiteMalebulbaronset",
    "PEAK_spinal_female_Intercept_Offset"= "outpeak_sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Intercept_Offset"= "outpeak_sexsiteFemalebulbaronset",
    
    "PEAK_spinal_male_Slope" = "outpeak_days_from_baseline",
    "PEAK_bulbar_male_Slope_Offset" = "outpeak_days_from_baseline.sexsiteMalebulbaronset",
    "PEAK_spinal_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalespinalonset",
    "PEAK_bulbar_female_Slope_Offset" = "outpeak_days_from_baseline.sexsiteFemalebulbaronset",
    
    "ALSFRSRESP_spinal_male_Intercept" = "alsfrsresp_Intercept",
    "ALSFRSRESP_bulbar_male_Intercept_Offset" = "alsfrsresp_sexsiteMalebulbaronset",
    "ALSFRSRESP_spinal_female_Intercept_Offset"= "alsfrsresp_sexsiteFemalespinalonset",
    "ALSFRSRESP_bulbar_female_Intercept_Offset"= "alsfrsresp_sexsiteFemalebulbaronset",
    
    "ALSFRSRESP_spinal_male_Slope" = "alsfrsresp_days_from_baseline",
    "ALSFRSRESP_bulbar_male_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteMalebulbaronset",
    "ALSFRSRESP_spinal_female_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteFemalespinalonset",
    "ALSFRSRESP_bulbar_female_Slope_Offset" = "alsfrsresp_days_from_baseline.sexsiteFemalebulbaronset"
)



pop_effects_5outPCT <- pop_effects_5outPCT %>% mutate(
    FVCPCT_bulbar_male_Intercept = FVCPCT_spinal_male_Intercept + FVCPCT_bulbar_male_Intercept_Offset,
    FVCPCT_bulbar_male_Slope = FVCPCT_spinal_male_Slope + FVCPCT_bulbar_male_Slope_Offset,
    FVCPCT_spinal_female_Intercept = FVCPCT_spinal_male_Intercept +
        FVCPCT_spinal_female_Intercept_Offset,
    FVCPCT_bulbar_female_Intercept = FVCPCT_spinal_male_Intercept +
        FVCPCT_bulbar_female_Intercept_Offset,
    FVCPCT_spinal_female_Slope = FVCPCT_spinal_male_Slope + FVCPCT_spinal_female_Slope_Offset,
    FVCPCT_bulbar_female_Slope = FVCPCT_spinal_male_Slope + FVCPCT_bulbar_female_Slope_Offset,
    
    SVCPCT_bulbar_male_Intercept = SVCPCT_spinal_male_Intercept + SVCPCT_bulbar_male_Intercept_Offset,
    SVCPCT_bulbar_male_Slope = SVCPCT_spinal_male_Slope + SVCPCT_bulbar_male_Slope_Offset,
    SVCPCT_spinal_female_Intercept = SVCPCT_spinal_male_Intercept +
        SVCPCT_spinal_female_Intercept_Offset,
    SVCPCT_bulbar_female_Intercept = SVCPCT_spinal_male_Intercept +
        SVCPCT_bulbar_female_Intercept_Offset,
    SVCPCT_spinal_female_Slope = SVCPCT_spinal_male_Slope + SVCPCT_spinal_female_Slope_Offset,
    SVCPCT_bulbar_female_Slope = SVCPCT_spinal_male_Slope + SVCPCT_bulbar_female_Slope_Offset,
    
    SNIP_bulbar_male_Intercept = SNIP_spinal_male_Intercept + SNIP_bulbar_male_Intercept_Offset,
    SNIP_bulbar_male_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_male_Slope_Offset,
    SNIP_spinal_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_spinal_female_Intercept_Offset,
    SNIP_bulbar_female_Intercept = SNIP_spinal_male_Intercept +
        SNIP_bulbar_female_Intercept_Offset,
    SNIP_spinal_female_Slope = SNIP_spinal_male_Slope + SNIP_spinal_female_Slope_Offset,
    SNIP_bulbar_female_Slope = SNIP_spinal_male_Slope + SNIP_bulbar_female_Slope_Offset,
    
    PEAK_bulbar_male_Intercept = PEAK_spinal_male_Intercept + PEAK_bulbar_male_Intercept_Offset,
    PEAK_bulbar_male_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_male_Slope_Offset,
    PEAK_spinal_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_spinal_female_Intercept_Offset,
    PEAK_bulbar_female_Intercept = PEAK_spinal_male_Intercept +
        PEAK_bulbar_female_Intercept_Offset,
    PEAK_spinal_female_Slope = PEAK_spinal_male_Slope + PEAK_spinal_female_Slope_Offset,
    PEAK_bulbar_female_Slope = PEAK_spinal_male_Slope + PEAK_bulbar_female_Slope_Offset,
    
    ALSFRSRESP_bulbar_male_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_bulbar_male_Intercept_Offset,
    ALSFRSRESP_bulbar_male_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_bulbar_male_Slope_Offset,
    ALSFRSRESP_spinal_female_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_spinal_female_Intercept_Offset,
    ALSFRSRESP_bulbar_female_Intercept = ALSFRSRESP_spinal_male_Intercept +
        ALSFRSRESP_bulbar_female_Intercept_Offset,
    ALSFRSRESP_spinal_female_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_spinal_female_Slope_Offset,
    ALSFRSRESP_bulbar_female_Slope = ALSFRSRESP_spinal_male_Slope + ALSFRSRESP_bulbar_female_Slope_Offset
)



pop_effects_5outPCT <- data.frame(par = names(pop_effects_5outPCT), t(pop_effects_5outPCT))

## Break reader friendly parameter names into multiple columns
metadata_5out <- pop_effects_5outPCT$par %>% str_split("_", simplify = TRUE) %>% data.frame()
names(metadata_5out) <- c("outcome", "onset", "sex", "param", "offsetpar")
pop_effects_5outPCT <- bind_cols(metadata_5out, pop_effects_5outPCT)
pop_effects_5outPCT$par <- NULL
pop_effects_5outPCT <- pop_effects_5outPCT %>% filter( offsetpar!= "Offset")
pop_effects_5outPCT$offsetpar <- NULL

# timescale of slopes currently slopes per day - convert to per average month
pop_effects_5outPCT[ pop_effects_5outPCT$param == "Slope", 5:8] <- (365.25 *
                                                                        pop_effects_5outPCT[ pop_effects_5outPCT$param == "Slope", 5:8] / 12)
pop_effects_5outPCT %>% arrange(outcome, param, onset)

# remove any params not Intercept or slope
pop_effects_5outPCT <- pop_effects_5outPCT %>% filter( param %in% c("Intercept", "Slope"))

# Write to file
write.csv(pop_effects_5outPCT, "Results/Table3_alsfrsresp5outs_MargEff.csv", row.names = FALSE)



##### Fit Models for investigational explanatory Variables #####

# Match baseline vars to df2
df2$baseline_stage <- dfbase[ match(df2$uin, dfbase$uin), ]$staging
df2$baseline_resphx <- dfbase[ match(df2$uin, dfbase$uin), ]$resp_hx
df2$baseline_ever_smoked <- dfbase[ match(df2$uin, dfbase$uin), ]$ever_smoked
df2$baseline_pkyrs <- dfbase[ match(df2$uin, dfbase$uin), ]$baseline_pkyrs

df2$baseline_breathless_lying <- dfbase[ match(df2$uin, dfbase$uin), ]$baseline_breathless_lying
df2$baseline_sob_rest <- dfbase[ match(df2$uin, dfbase$uin), ]$baseline_sob_rest
df2$baseline_sob_active <- dfbase[ match(df2$uin, dfbase$uin), ]$baseline_sob_active
df2$baseline_prod_cough <- dfbase[ match(df2$uin, dfbase$uin), ]$baseline_prod_cough


## fit baseline staging model
fit_model_basestage <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_stage +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_basestage, "Models/brms_finaldata_newmodel_basestage_3500.RDS")
#fit_model_basestage <- readRDS("Models/brms_finaldata_newmodel_basestage_3500.RDS")

#extract model fixed effect coefficients and write to file
coefs_basestage <- get_model_fixedcoefsPCT(fit_model_basestage)
coefs_basestage <- coefs_basestage %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_basestage, "Results/Finaldata_basestagemodel_fixed_effects.xlsx")


## fit baseline resphx model
fit_model_baseresphx <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_resphx +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseresphx, "Models/brms_finaldata_newmodel_baseresphx_3500.RDS")
#fit_model_baseresphx <- readRDS("Models/brms_finaldata_newmodel_baseresphx_3500.RDS")

#extract model fixed effect coefficients and write to file
coefs_baseresphx <- get_model_fixedcoefsPCT(fit_model_baseresphx)
coefs_baseresphx <- coefs_baseresphx %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseresphx, "Results/Finaldata_baseresphxmodel_fixed_effects.xlsx")


# fit model for ever vs never smoked at baseline
fit_model_baseeverSmoke <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_ever_smoked +
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseeverSmoke, "Models/brms_finaldata_newmodel_baseeverSmoke_3500.RDS")
#fit_modelbaseeverSmoke <- readRDS("Models/brms_finaldata_newmodel_baseeverSmoke_3500.RDS")

#extract model fixed effect coefficients and write to file
coefs_baseeverSmoke <- get_model_fixedcoefsPCT(fit_model_baseeverSmoke)
coefs_baseeverSmoke <- coefs_baseeverSmoke %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseeverSmoke, "Results/Finaldata_baselineeverSmokemodel_fixed_effects.xlsx")


# fit model for orthopnea at baseline
fit_model_baseline_breathless_lying <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_breathless_lying + 
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseline_breathless_lying, "Models/brms_finaldata_newmodel_baseline_breathless_lying_3500.RDS")


coefs_baseline_breathlesslying  <- get_model_fixedcoefsPCT(fit_model_baseline_breathless_lying)
coefs_baseline_breathlesslying <- coefs_baseline_breathlesslying %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseline_breathlesslying, "Results/Finaldata_baseline_breathlesslying_fixed_effects.xlsx")


# fit model for SOB at rest at baseline
fit_model_baseline_sob_rest <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_sob_rest + 
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseline_sob_rest, "Models/brms_finaldata_newmodel_baseline_sobrest_3500.RDS")


coefs_baseline_sobrest  <- get_model_fixedcoefsPCT(fit_model_baseline_sob_rest)
coefs_baseline_sobrest <- coefs_baseline_sobrest %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseline_sobrest, "Results/Finaldata_baseline_sobrest_fixed_effects.xlsx")


# Fit model for SOB when active at baseline
fit_model_baseline_sob_active <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_sob_active + 
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseline_sob_active, "Models/brms_finaldata_newmodel_baseline_sobactive_3500.RDS")


coefs_baseline_sobactive  <- get_model_fixedcoefsPCT(fit_model_baseline_sob_active)
coefs_baseline_sobactive <- coefs_baseline_sobactive %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseline_sobactive, "Results/Finaldata_baseline_sobactive_fixed_effects.xlsx")


# Fit productive cough at baseline model
fit_model_baseline_prod_cough <- brm(
    bf(mvbind(fvc_max_percent_pred, svc_percent_pred, outsnip, outpeak, alsfrs_resp) ~
           days_from_baseline * sexsite + baseline_prod_cough + 
           cohort + (1 | p | uin) + (0 + days_from_baseline | uin)) +
        set_rescor(TRUE),
    data = df2, chains=4,
    prior = priors_5outsPCT,
    sample_prior = "no", seed = 343451,
    iter = 3500, threads = threading(2))
saveRDS(fit_model_baseline_prod_cough, "Models/brms_finaldata_newmodel_baseline_prod_cough_3500.RDS")


coefs_baseline_prod_cough  <- get_model_fixedcoefsPCT(fit_model_baseline_prod_cough)
coefs_baseline_prod_cough <- coefs_baseline_prod_cough %>% 
    dplyr::select(variable,
                  Estimate_fvcpct, Q2.5_fvcpct, Q97.5_fvcpct, Est.Error_fvcpct,
                  Estimate_svcpct, Q2.5_svcpct, Q97.5_svcpct, Est.Error_svcpct,
                  Estimate_snip, Q2.5_snip, Q97.5_snip, Est.Error_snip,
                  Estimate_peak, Q2.5_peak, Q97.5_peak, Est.Error_peak,
                  Estimate_alsfrsrr, Q2.5_alsfrsrr, Q97.5_alsfrsrr, Est.Error_alsfrsrr)
write.xlsx(coefs_baseline_prod_cough, "Results/Finaldata_baseline_prod_cough_fixed_effects.xlsx")


#### Assemble Table 4 ####
tab4 <- bind_rows(coefs_basestage %>% filter(variable == "baseline_stageKings3"),
                  coefs_baseresphx %>% filter(variable %in%
                                                  c("baseline_resphxoccasional", "baseline_resphxmodfreq")),
                  coefs_baseeverSmoke %>% filter(variable == "baseline_ever_smokedY"),
                  coefs_baseline_breathlesslying %>% filter(variable == "baseline_breathless_lyingYes"),
                  coefs_baseline_sobrest %>% filter(variable == "baseline_sob_restY"),
                  coefs_baseline_sobactive %>% filter(variable == "baseline_sob_activeY"),
                  coefs_baseline_prod_cough %>% filter(variable == "baseline_prod_coughY"))

write.xlsx(tab4, "Results/Table4_exploratory_model_ParamEsts.xlsx")


