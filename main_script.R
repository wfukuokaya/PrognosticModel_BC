pacman::p_load(lubridate, dcurves, survival, survminer, gtsummary, tidymodels, ggfortify, glmnet, rms, missForest, ggsurvfit, patchwork, broom)
library(tidyverse) # load tidyverse at the end
select <- dplyr::select

df # full analysis set

# define ggplot2 theme...
theme_juog <- function(base_size = 10, 
                       dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.4)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

# risk group
theme_km <- function(base_size = 9, 
                     dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.2)),
          plot.subtitle = element_text(size = rel(1)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

set.seed(234)
imp_df <- missForest(df, maxiter = 1, verbose = TRUE) # random forest imputation
dimp <- imp_df$ximp

# calculate the median follow-up
library(pec)
fu_time <- quantile(prodlim::prodlim(Hist(fu, death_all) ~ 1, data = dimp, reverse = TRUE))

# model development
set.seed(26)
df_split <- initial_split(dimp, prop = 2/3) # split into a 2:1 ratio
train_data <- training(df_split) %>% as_data_frame() # training dataset
test_data <- testing(df_split) %>% as_data_frame() # testing dataset

# time-to-event distribution data
max_os_train <- max(train_data$fu) %>% round(digits = 1)
max_os_test <- max(test_data$fu) %>% round(digits = 1)

# define variable set and survival data
surv <- train_data %>% select(fu, death_cancer) %>% as_data_frame()
fit <- survival::Surv(surv$fu, surv$death_cancer) %>% as.matrix() # survival function
var <- train_data %>% select(age_rc:logcrp_post) %>% mutate_if(is.factor, as.numeric) %>% as.matrix() # covariates

regularization_path <- glmnet(x = var, y = fit, family = "cox", alpha = 1) %>% 
  broom::tidy(return_zeros = TRUE) %>% 
  ggplot2::ggplot(ggplot2::aes(lambda, estimate, group = term, color = term)) + 
  geom_line(size = 0.7) + 
  scale_x_log10() + 
  theme_juog() + 
  scale_color_manual(values = as.vector(pals::glasbey()), 
                     labels = c("Age at surgery", 
                                "Albumin",
                                "Albumin (post-surgery)",
                                "Previous BCG",
                                "Body mass index", 
                                "Carcinoma in situ",
                                "Clinical N stage", 
                                "Clinical T stage",
                                "Urinary diversion",
                                "ECOG performance status",
                                "Hemoglobin",
                                "Hemoglobin (post-surgery)",
                                "Tumor histology",
                                "CRP",
                                "CRP (post-surgery)",
                                "NLR",
                                "NLR (post-surgery)",
                                "LVI",
                                "Tumor multifocality",
                                "Neoadjuvant chemotherapy",
                                "Previous NMIBC",
                                "Pathological N stage",
                                "Pathological T stage",
                                "Resective margin",
                                "Sex", 
                                "Tumor size",
                                "Smoking status",
                                "Procedure",
                                "Previous UTUC")) + 
  theme(legend.key.height = unit(0.25, "cm"),
        legend.key.width = unit(0.24, "cm")) + 
  labs(
    x = "Lambda",
    y = "Coefficient"
  )

# perform adaptive LASSO regression with 10-fold CV
set.seed(3)
lasso.cv <- cv.glmnet(x = var, y = fit, family = "cox", type.measure = "C")

glance_lasso_cv <- lasso.cv %>% broom::glance() # extract data regarding lambda
tidied_lasso_cv <- lasso.cv %>% broom::tidy() # tidied results of 10-fold cross-validation

plot_lasso_cv <- 
  ggplot(tidied_lasso_cv, aes(lambda, estimate)) + 
  geom_line(size = 0.7, color = "#3B4992FF") + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#3B4992FF", alpha = 0.3) + 
  scale_x_log10() + 
  ylim(0.7, 0.85) + 
  theme_juog() + 
  labs(
    x = "Lambda",
    y = "Concordance index (95% CI)"
  ) + 
  geom_vline(xintercept = glance_lasso_cv$lambda.min, size = 0.7) +
  geom_vline(xintercept = glance_lasso_cv$lambda.1se, size = 0.7, lty = 2) 

# solution path with lambda.min
regularization_path_lambda <- regularization_path + 
  geom_vline(xintercept = glance_lasso_cv$lambda.min, 
             color = "grey45", size = 0.7) 
regularization_path_lambda

# extract values of the penalty parameter `lambda` at which predictions are required.
coef_lasso <- (coef(lasso.cv, s = lasso.cv$lambda.min)) 
var_lasso <- as.data.frame(as.matrix(coef_lasso)) %>% 
  add_rownames(var = "rowname") %>% 
  rename(var = rowname, coef = "1") %>% 
  filter(coef != 0)
var_lasso

cox_os_formula <- as.formula(paste0("Surv(fu, death_all) ~ ", paste0(var_lasso$var, collapse = " + "))) # define formula
cox_css_formula <- as.formula(paste0("Surv(fu, death_cancer) ~ ", paste0(var_lasso$var, collapse = " + "))) # define formula
cox_nutrfs_formula <- as.formula(paste0("Surv(nutrfs, nutr) ~ ", paste0(var_lasso$var, collapse = " + "))) # define formula

# LASSO regression
train_reg <- train_data %>% 
  mutate_at(vars(sex, ecog_ps_cat, smoking, utuc, nmibc, bcg, size_tur, 
                 multi_tur, hist_rc, type, div, ct, cn, pt, pn, lvi, cis_rc, rm_other, 
                 nac), 
            funs(as.character)) # numeric to integar

theme_gtsummary_journal(journal = "jama")
lasso_cox_os <- coxph(cox_os_formula, data = train_reg) %>% # LASSO Cox regression model
  tbl_regression(exponentiate = TRUE,
                 label = list(
                   ecog_ps_cat ~ "ECOG performance status",
                   utuc ~ "Previous history of UTUC",
                   pt ~ "Pathoological T stage",
                   pn ~ "Pathological N stage",
                   lvi ~ "Lymphovascular invasion",
                   rm_other ~ "Resective margin",
                   lognlr ~ "Log-transformed NLR",
                   alb ~ "Albumin"
                 ))

lasso_cox_css <- coxph(cox_css_formula, data = train_reg) %>% # LASSO Cox regression model
  tbl_regression(exponentiate = TRUE,
                 label = list(
                   ecog_ps_cat ~ "ECOG performance status",
                   utuc ~ "Previous history of UTUC",
                   pt ~ "Pathoological T stage",
                   pn ~ "Pathological N stage",
                   lvi ~ "Lymphovascular invasion",
                   rm_other ~ "Resective margin",
                   lognlr ~ "Log-transformed NLR",
                   alb ~ "Albumin"
                 ))

lasso_cox_nutrfs <- coxph(cox_nutrfs_formula, data = train_reg) %>% # LASSO Cox regression model
  tbl_regression(exponentiate = TRUE,
                 label = list(
                   ecog_ps_cat ~ "ECOG performance status",
                   utuc ~ "Previous history of UTUC",
                   pt ~ "Pathoological T stage",
                   pn ~ "Pathological N stage",
                   lvi ~ "Lymphovascular invasion",
                   rm_other ~ "Resective margin",
                   lognlr ~ "Log-transformed NLR",
                   alb ~ "Albumin"
                 ))

# multivariable Cox regression with backward variable selection
intercept_only <- coxph(Surv(fu, death_cancer) ~ 1, data = train_reg) # model with intercept only
all <- coxph(Surv(fu, death_cancer) ~ age_rc + sex + bmi + ecog_ps_cat + 
               smoking + utuc + nmibc + bcg + size_tur + multi_tur + hist_rc + type + div + 
               ct + cn + pt + pn + lvi + cis_rc + rm_other + nac + hgb + lognlr + alb + logcrp + hgb_post + lognlr_post + alb_post + logcrp, data = train_reg) # model with full variables
forward <- stats::step(intercept_only, direction = "forward", scope = as.formula(all))

formula_stepwise <- as.formula(Surv(fu, death_cancer) ~ ecog_ps_cat + 
                                 pt + pn + lvi + rm_other + lognlr + alb + utuc)

stepwise_cox_css <- coxph(formula_stepwise, data = train_reg) %>% 
  tbl_regression(
    exponentiate = TRUE,
    label = list(
      ecog_ps_cat ~ "ECOG performance status",
      utuc ~ "Previous history of UTUC",
      pt ~ "Pathoological T stage",
      pn ~ "Pathological N stage",
      lvi ~ "Lymphovascular invasion",
      rm_other ~ "Resective margin",
      lognlr ~ "Log-transformed NLR",
      alb ~ "Albumin"
    ))

data_aic <- tibble(
  num = seq(0, 8), # number of variables
  aic = forward$anova$AIC # changes in AICs
)

# plot variable number and AIC
plot_aic <- ggplot(data = data_aic, 
                   aes(x = num, y = aic, label = sprintf("%.1f",aic))) + 
  geom_point(size = 3, color = "#3B4992FF") + 
  geom_line(size = 0.7, color = "#3B4992FF") + 
  geom_text(vjust = -2, size = 3, show.legend = FALSE, fontface = "bold", color = "#3B4992FF") + 
  theme_juog() + 
  scale_x_continuous(breaks = seq(0, 15, by = 1)) + 
  coord_cartesian(ylim = c(2850, 3150), expand = TRUE) + 
  labs(
    x = "Number of variables",
    y = "Akaike Information Criterion"
  ) 

# table1 for all patients (crude data)
reset_gtsummary_theme()
data_tbl <- df %>% mutate(imp = "Crude")
data_imp <- df_imp %>% mutate(imp = "Imputed")
df_tblimp <- rbind(data_tbl, data_imp)

df_train <- train_data %>% mutate(set = "Training")
df_test <- test_data %>% mutate(set = "Validation")
df_tt <- rbind(df_train, df_test)

# depict survival curves for the total set
# risk group
rfs_fit <- survival::survfit(data = df_tt, Surv(nutrfs, nutr) ~ set)
rfs <- survminer::ggsurvplot(rfs_fit, 
                             risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                             title = "NUTRFS", 
                             legend = "top",  legend.labs = c("Training", "Validation"),
                             legend.title = "Set", 
                             censor = TRUE, censor.shape = "|", censor.size = 1.5,
                             xlab = "Months since radical cystectomy", 
                             ylab = "NUTRFS probability (95% CI)",
                             palette = "aaas", size = 0.5,  break.time.by = 12,
                             ggtheme = theme_grey(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.2,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

rfs$plot <- rfs$plot + 
  theme_km()

rfs$table <- rfs$table + 
  theme_void() + 
  theme(
    text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
css_fit <- survival::survfit(data = df_tt, Surv(fu, death_cancer) ~ set)
css <- survminer::ggsurvplot(css_fit, 
                             risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                             title = "CSS", 
                             legend = "top",  legend.labs = c("Training", "Validation"),
                             legend.title = "Set", 
                             censor = TRUE, censor.shape = "|", censor.size = 1.5,
                             xlab = "Months since radical cystectomy", 
                             ylab = "CSS probability (95% CI)",
                             palette = "aaas", size = 0.5,  break.time.by = 12,
                             ggtheme = theme_grey(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.2,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

css$plot <- css$plot + 
  theme_km() 

css$table <- css$table + 
  theme_void() + 
  theme(
    text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
os_fit <- survival::survfit(data = df_tt, Surv(fu, death_all) ~ set)
os <- survminer::ggsurvplot(os_fit, 
                            risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                            title = "OS", 
                            legend = "top",  legend.labs = c("Training", "Validation"),
                            legend.title = "Set", 
                            censor = TRUE, censor.shape = "|", censor.size = 1.5,
                            xlab = "Months since radical cystectomy", 
                            ylab = "OS probability (95% CI)",
                            palette = "aaas", size = 0.5,  break.time.by = 12,
                            ggtheme = theme_grey(), risk.table.title = "Number at risk",
                            risk.table.col = "strata",
                            tables.height = 0.12, risk.table.fontsize = 3.2,
                            tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os$plot <- os$plot + 
  theme_km() 

os$table <- os$table + 
  theme_void() + 
  theme(
    text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# merge
merge <- survminer::arrange_ggsurvplots(list(rfs, css, os), print = FALSE, nrow = 1, ncol = 3)

tbl_theme <- list(
  "tbl_summary-str:missing_stat" = "{N_miss} ({p_miss})" # display the number of missing values and its percentage
)
set_gtsummary_theme(tbl_theme)

tbl1_imp <- df_tblimp %>%
  select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, size_tur, multi_tur, ct, cn, hist_rc, 
         type, div, pt, pn, lvi, cis_rc, rm_other, nac, hgb, lognlr, alb, logcrp, imp) %>%
  
  tbl_summary(
    by = imp, 
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urothelial carcinoma",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      size_tur ~ "Tumor size",
      multi_tur ~ "Multifocality",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      hist_rc ~ "Tumor histology",
      type ~ "Surgical procedure",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      lvi ~ "Lymphovascular invasion",
      cis_rc ~ "Concomitant carcinoma in situ",
      rm_other ~ "Soft-tissue resective margin",
      nac ~ "Neoadjuvant chemotherapy",
      hgb ~ "Baseline hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      alb ~ "Albumin levels",
      logcrp ~ "C-reactive protein levels (log-transformed)"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) 

# training vs. testing set
tbl1_tt <- df_tt %>%
  select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, size_tur, multi_tur, ct, cn, hist_rc, 
         type, div, pt, pn, lvi, cis_rc, rm_other, nac, hgb, lognlr, alb, logcrp, set) %>%
  tbl_summary(
    by = set,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urothelial carcinoma",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      size_tur ~ "Tumor size",
      multi_tur ~ "Multifocality",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      hist_rc ~ "Tumor histology",
      type ~ "Surgical procedure",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      lvi ~ "Lymphovascular invasion",
      cis_rc ~ "Concomitant carcinoma in situ",
      rm_other ~ "Soft-tissue resective margin",
      nac ~ "Neoadjuvant chemotherapy",
      hgb ~ "Baseline hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      alb ~ "Albumin levels",
      logcrp ~ "C-reactive protein levels (log-transformed)"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_overall() %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Cohort**") %>%
  modify_caption("**Table 1. Patient Characteristics**") %>%
  bold_labels() %>% 
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) 

# create nomogram
library("rms")

data_nomogram <- train_data %>% 
  select(fu, death_cancer, age_rc, utuc, multi_tur, hist_rc, ecog_ps_cat, pt, pn, lvi, rm_other, lognlr, alb) %>% 
  mutate(ecog_ps_cat = case_when(ecog_ps_cat == 0 ~ "0 to 1", TRUE ~ "2 or more"),
         pt = case_when(pt == 0 ~ "pT0 to pT1", pt == 1 ~ "pT2", pt == 2 ~ "pT3", TRUE ~ "pT4") %>% 
           factor(levels = c("pT0 to pT1", "pT2", "pT3", "pT4")),
         pn = case_when(pn == 0 ~ "pN0", TRUE ~ "pN1 to pN3"),
         lvi = case_when(lvi == 1 ~ "Yes", TRUE ~ "No"),
         utuc = case_when(utuc == 0 ~ "no", TRUE ~ "yes"),
         rm_other = case_when(rm_other == 1 ~ "Positive", TRUE ~ "Negative")
  ) %>% as_data_frame()

dist <- datadist(data_nomogram)
options(datadist = "dist")

mod <- cph(cox_css_formula, data = data_nomogram, 
           x = TRUE, y = TRUE, surv = TRUE)
mod <- Newlabels(mod, c(
  ecog_ps_cat = "ECOG PS",
  utuc = "Previous UTUC",
  pt = "Pathological T stage",
  pn = "Pathological N stage",
  lvi = "LVI",
  rm_other = "Resective margin",
  lognlr = "Log-transformed NLR",
  alb = "Albumin"
))
sv <- Survival(mod)
surv60 <- function(x) sv(60, lp = x)

nomogram <- nomogram(mod, fun = list(surv60), 
                     funlabel = c("60-month CSS probability"), lp = FALSE)

print(nomogram) # detail of scoring system

plot(nomogram, 
     total.points.label="Sum of all points",
     tcl = 0.5,
     lmgp = 0,
     ia.space = 0)

# risk score calculation
calculate_score <- function(x) {
  x %>% 
    mutate(
      score = 
        as.numeric(case_when(ecog_ps_cat == 1 ~ 26, TRUE ~ 0)) + 
        as.numeric(case_when(utuc == 1 ~ 24, TRUE ~ 0)) + 
        as.numeric(case_when(pt == 0 ~ 0,
                             pt == 1 ~ 14, 
                             pt == 2 ~ 44,
                             TRUE ~ 57)) + 
        as.numeric(case_when(pn == 1 ~ 22, TRUE ~ 0)) + 
        as.numeric(case_when(lvi == 1 ~ 29, TRUE ~ 0)) + 
        as.numeric(case_when(rm_other == 1 ~ 32, TRUE ~ 0)) + 
        
        as.numeric(case_when(
          lognlr < -2.5 ~ 0,
          lognlr >= -2.5 & lognlr < -2 ~ 5,
          lognlr >= -2 & lognlr < -1.5 ~ 11,
          lognlr >= -1.5 & lognlr < -1 ~ 16,
          lognlr >= -1 & lognlr < -0.5 ~ 21,
          lognlr >= -0.5 & lognlr < 0 ~ 27,
          lognlr >= 0 & lognlr < 0.5 ~ 32,
          lognlr >= 0.5 & lognlr < 1 ~ 37,
          lognlr >= 1 & lognlr < 1.5 ~ 43,
          lognlr >= 1.5 & lognlr < 2 ~ 48,
          lognlr >= 2 & lognlr < 2.5 ~ 53,
          lognlr >= 2.5 & lognlr < 3 ~ 59,
          lognlr >= 3 & lognlr < 3.5 ~ 64,
          lognlr >= 3.5 & lognlr < 4 ~ 69,
          lognlr >= 4 & lognlr < 4.5 ~ 75,
          lognlr >= 4.5 ~ 80
        )) + 
        as.numeric(case_when(
          alb < 1.5 ~ 100,
          alb >= 1.5 & alb < 2 ~ 92,
          alb >= 2 & alb < 2.5 ~ 85,
          alb >= 2.5 & alb < 3 ~ 77,
          alb >= 3 & alb < 3.5 ~ 69,
          alb >= 3.5 & alb < 4 ~ 62,
          alb >= 4 & alb < 4.5 ~ 54,
          alb >= 4.5 & alb < 5 ~ 46,
          alb >= 5 & alb < 5.5 ~ 38,
          alb >= 5.5 & alb < 6 ~ 31,
          alb >= 6 & alb < 6.5 ~ 23,
          alb >= 6.5 & alb < 7 ~ 15,
          alb >= 7 & alb < 7.5 ~ 8,
          alb >= 7.5 ~ 0
        )),
      score_quantile = case_when(
        score <= quantile(score, 0.50) ~ "favorable",
        score >= quantile(score, 0.75) ~ "poor", 
        TRUE ~ "intermediate"
      ) %>% factor(levels = c("favorable", "intermediate", "poor")),
      score_median = case_when(score > median(score) ~ "poor", TRUE ~ "favorable"),
      score_tertile = case_when(score > quantile(score, 0.66) ~ "poor risk", TRUE ~ "good risk") %>% 
        factor(levels = c("good risk", "poor risk"))
    )
}
train_data <- train_data %>% calculate_score()
test_data <- test_data %>% calculate_score()
df_score <- dimp %>% calculate_score()

# survival distribution stratified by risk group...
rfs_fit <- survival::survfit(data = df_score, Surv(nutrfs, nutr) ~ score_tertile)
css_fit <- survival::survfit(data = df_score, Surv(fu, death_cancer) ~ score_tertile)
os_fit <- survival::survfit(data = df_score, Surv(fu, death_all) ~ score_tertile)

cox_rfs_tertile <- coxph(data = test_data, Surv(nutrfs, nutr) ~ score_tertile) %>% tbl_regression(exponentiate = TRUE)
cox_css_tertile <- coxph(data = test_data, Surv(fu, death_cancer) ~ score_tertile) %>% tbl_regression(exponentiate = TRUE)
cox_os_tertile <- coxph(data = test_data, Surv(fu, death_all) ~ score_tertile) %>% tbl_regression(exponentiate = TRUE)

# stratified Kaplan-Meier curves
rfs_model_fit <- survival::survfit(data = df_score, Surv(nutrfs, nutr) ~ score_tertile)
rfs_model <- survminer::ggsurvplot(rfs_model_fit, 
                                   risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                   title = "NUTRFS", 
                                   subtitle = "stratified by the developed model",
                                   legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                   legend.title = "The developed model", 
                                   censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                   xlab = "Months since radical cystectomy", 
                                   ylab = "NUTRFS probability (95% CI)",
                                   palette = "aaas", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.15, risk.table.fontsize = 2.8,
                                   conf.int = TRUE) 

rfs_model$plot <- rfs_model$plot + 
  theme_km() 

rfs_model$table <- rfs_model$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# cancer-specific survival
css_model_fit <- survival::survfit(data = df_score, Surv(fu, death_cancer) ~ score_tertile)
css_model <- survminer::ggsurvplot(css_model_fit, 
                                   risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                   title = "CSS", 
                                   subtitle = "stratified by the developed model",
                                   legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                   legend.title = "The developed model", 
                                   censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                   xlab = "Months since radical cystectomy", 
                                   ylab = "CSS probability (95% CI)",
                                   palette = "aaas", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.15, risk.table.fontsize = 2.8,
                                   conf.int = TRUE) 

css_model$plot <- css_model$plot + 
  theme_km()

css_model$table <- css_model$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
os_model_fit <- survival::survfit(data = df_score, Surv(fu, death_all) ~ score_tertile)
os_model <- survminer::ggsurvplot(os_model_fit, 
                                  risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                  title = "OS", 
                                  subtitle = "stratified by the developed model",
                                  legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                  legend.title = "The developed model", 
                                  censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                  xlab = "Months since radical cystectomy", 
                                  ylab = "OS probability (95% CI)",
                                  palette = "aaas", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.15, risk.table.fontsize = 2.8,
                                  conf.int = TRUE) 

os_model$plot <- os_model$plot + 
  theme_km()

os_model$table <- os_model$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# stratified Kaplan-Meier curves
rfs_cm274_fit <- survival::survfit(data = df_score, Surv(nutrfs, nutr) ~ risk_recurrence)
rfs_cm274 <- survminer::ggsurvplot(rfs_cm274_fit, 
                                   risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                   title = "NUTRFS", 
                                   subtitle = "stratified by CheckMate 274 risk group",
                                   legend = "top", legend.labs = c("Low risk", "High risk"),
                                   legend.title = "CheckMate 274 risk group", 
                                   censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                   xlab = "Months since radical cystectomy", 
                                   ylab = "NUTRFS probability (95% CI)",
                                   palette = "aaas", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.15, risk.table.fontsize = 2.8,
                                   conf.int = TRUE) 

rfs_cm274$plot <- rfs_cm274$plot + 
  theme_km() 

rfs_cm274$table <- rfs_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# cancer-specific survival
css_cm274_fit <- survival::survfit(data = df_score, Surv(fu, death_cancer) ~ risk_recurrence)
css_cm274 <- survminer::ggsurvplot(css_cm274_fit, 
                                   risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                   title = "CSS", 
                                   subtitle = "stratified by CheckMate 274 risk group",
                                   legend = "top", legend.labs = c("Low risk", "High risk"),
                                   legend.title = "CheckMate 274 risk group", 
                                   censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                   xlab = "Months since radical cystectomy", 
                                   ylab = "CSS probability (95% CI)",
                                   palette = "aaas", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.15, risk.table.fontsize = 2.8,
                                   conf.int = TRUE) 

css_cm274$plot <- css_cm274$plot + 
  theme_km()

css_cm274$table <- css_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
os_cm274_fit <- survival::survfit(data = df_score, Surv(fu, death_all) ~ risk_recurrence)
os_cm274 <- survminer::ggsurvplot(os_cm274_fit, 
                                  risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                  title = "OS", 
                                  subtitle = "stratified by CheckMate 274 risk group",
                                  legend = "top", legend.labs = c("Low risk", "High risk"),
                                  legend.title = "CheckMate 274 risk group", 
                                  censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                  xlab = "Months since radical cystectomy", 
                                  ylab = "OS probability (95% CI)",
                                  palette = "aaas", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.15, risk.table.fontsize = 2.8,
                                  conf.int = TRUE) 

os_cm274$plot <- os_cm274$plot + 
  theme_km()

os_cm274$table <- os_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# merge
merge <- survminer::arrange_ggsurvplots(list(rfs_model, rfs_cm274, css_model, css_cm274, os_model, os_cm274), 
                                        print = FALSE, 
                                        nrow = 2, ncol = 3)

nutrfs <- arrange_ggsurvplots(list(rfs_model, rfs_cm274),
                              print = FALSE, nrow = 1, ncol = 2)
css <- arrange_ggsurvplots(list(css_model, css_cm274),
                           print = FALSE, nrow = 1, ncol = 2)
os <- arrange_ggsurvplots(list(os_model, os_cm274),
                          print = FALSE, nrow = 1, ncol = 2)


# predicted survival curves
# predicted survival probability plot
weibull_nutrfs <- flexsurv::flexsurvreg(Surv(nutrfs, nutr) ~ score_quantile, data = df_score, 
                                        dist = "weibullPH")
weibull_nutrfs_predict <- summary(weibull_nutrfs, type = "survival", conf.int = TRUE, tidy = TRUE)
plot_nutrfs_weibull <- ggplot() + 
  geom_line(data = weibull_nutrfs_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_nutrfs_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  labs(
    subtitle = "NUTRFS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) +
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  theme_km() 

# predicted survival probability plot
weibull_css <- flexsurv::flexsurvreg(Surv(fu, death_cancer) ~ score_quantile, data = df_score, 
                                     dist = "weibullPH")
weibull_css_predict <- summary(weibull_css, type = "survival", tidy = TRUE)
plot_css_weibull <- ggplot() + 
  geom_line(data = weibull_css_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_css_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  labs(
    subtitle = "CSS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) +
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  theme_km() 

# predicted survival probability plot
weibull_os <- flexsurv::flexsurvreg(Surv(fu, death_all) ~ score_quantile, data = df_score, 
                                    dist = "weibullPH")
weibull_os_predict <- summary(weibull_os, type = "survival", tidy = TRUE)
plot_os_weibull <- ggplot() + 
  geom_line(data = weibull_os_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_os_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  labs(
    subtitle = "OS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) +
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  theme_km() 

weibull_merge <- (plot_nutrfs_weibull + plot_os_weibull) + plot_annotation(
  title = "Predicted survival curves",
  subtitle = "Estimated by flexible parametric survival analysis according to the developed scores",
) & theme_km() & theme(legend.key.height = unit(0.3, "cm"))

# overlay with Kaplan-Meier curves...

km_data_nutrfs <- survfit(Surv(nutrfs, nutr) ~ score_quantile, data = df_score) %>% broom::tidy(conf.int = TRUE) %>% 
  group_by(strata) %>% mutate(
    strata = case_when(
      strata == "score_quantile=favorable" ~ "favorable",
      strata == "score_quantile=intermediate" ~ "intermediate",
      TRUE ~ "poor"
    ) %>% factor(levels = c("favorable", "intermediate", "poor"))
  )

km_data_css <- survfit(Surv(fu, death_cancer) ~ score_quantile, data = df_score) %>% broom::tidy(conf.int = TRUE) %>% 
  group_by(strata) %>% mutate(
    strata = case_when(
      strata == "score_quantile=favorable" ~ "favorable",
      strata == "score_quantile=intermediate" ~ "intermediate",
      TRUE ~ "poor"
    ) %>% factor(levels = c("favorable", "intermediate", "poor"))
  )

km_data_os <- survfit(Surv(fu, death_all) ~ score_quantile, data = df_score) %>% broom::tidy(conf.int = TRUE) %>% 
  group_by(strata) %>% mutate(
    strata = case_when(
      strata == "score_quantile=favorable" ~ "favorable",
      strata == "score_quantile=intermediate" ~ "intermediate",
      TRUE ~ "poor"
    ) %>% factor(levels = c("favorable", "intermediate", "poor"))
  )

# km + flexible parametric survival curves
km_weibull_nutrfs <- ggplot() +
  # Kaplan-Meier curves
  geom_step(data = km_data_nutrfs, aes(x = time, y = estimate, color = strata, fill = strata), size = 0.5) + 
  geom_ribbon(data = km_data_nutrfs, aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.3) + 
  
  # Flexible parametric survival curves
  geom_line(data = weibull_nutrfs_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_nutrfs_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, by = 0.2)) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) + 
  theme_km() + 
  labs(
    title = "NUTRFS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  ggsci::scale_color_aaas() + 
  ggsci::scale_fill_aaas() 

# km + flexible parametric survival curves
km_weibull_css <- ggplot() +
  # Kaplan-Meier curves
  geom_step(data = km_data_css, aes(x = time, y = estimate, color = strata, fill = strata), size = 0.5) + 
  geom_ribbon(data = km_data_css, aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.3) + 
  
  # Flexible parametric survival curves
  geom_line(data = weibull_css_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_css_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, by = 0.2)) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) + 
  theme_km() + 
  labs(
    title = "CSS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  ggsci::scale_color_aaas() + 
  ggsci::scale_fill_aaas() 

# km + flexible parametric survival curves
km_weibull_os <- ggplot() +
  # Kaplan-Meier curves
  geom_step(data = km_data_os, aes(x = time, y = estimate, color = strata, fill = strata), size = 0.5) + 
  geom_ribbon(data = km_data_os, aes(x = time, ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.3) + 
  
  # Flexible parametric survival curves
  geom_line(data = weibull_os_predict, aes(x = time, y = est, color = score_quantile), size = 0.8) + 
  geom_ribbon(data = weibull_os_predict, aes(x = time, fill = score_quantile, ymin = lcl, ymax = ucl), alpha = 0.3) + 
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, by = 0.2)) + 
  scale_x_continuous(limits = c(0, 132), breaks = seq(0, 132, by = 12)) + 
  theme_km() + 
  labs(
    title = "OS",
    x = "Months since radical cystectomy",
    y = "Survival probability (95% CI)"
  ) + 
  ggsci::scale_color_aaas() + 
  ggsci::scale_fill_aaas() 

km_weibull_merge <- (km_weibull_nutrfs + km_weibull_css + km_weibull_os) + theme_km() 

# harrell"s C
# NUTRFS
# testing dataset
surv_test <- coxph(Surv(nutrfs, nutr) ~ score, data = test_data) %>% summary()
cindex_test <- surv_test$concordance[1] %>% round(digits = 2)

# testing dataset for base model
surv_test_base <- coxph(Surv(nutrfs, nutr) ~ risk_recurrence, data = test_data) %>% summary()
cindex_test_base <- surv_test_base$concordance[1] %>% round(digits = 2)

# testing dataset
surv_test <- coxph(Surv(fu, death_cancer) ~ score, data = test_data) %>% summary()
cindex_test <- surv_test$concordance[1] %>% round(digits = 2)

# testing dataset for base model
surv_test_base <- coxph(Surv(fu, death_cancer) ~ risk_recurrence, data = test_data) %>% summary()
cindex_test_base <- surv_test_base$concordance[1] %>% round(digits = 2)

# testing dataset
surv_test <- coxph(Surv(fu, death_all) ~ score, data = test_data) %>% summary()
cindex_test <- surv_test$concordance[1] %>% round(digits = 2)

# testing dataset for base model
surv_test_base <- coxph(Surv(fu, death_all) ~ risk_recurrence, data = test_data) %>% summary()
cindex_test_base <- surv_test_base$concordance[1] %>% round(digits = 2)

# decision curve analysis
# calculate failure probability at 24 months
surv_train <- coxph(Surv(nutrfs, nutr) ~ score, data = train_data) 
surv_train_base <- coxph(Surv(nutrfs, nutr) ~ risk_recurrence, data = train_data) 

train_data$model = c(1 - (summary(survfit(surv_train, newdata = train_data), times = 24)$surv)) 
train_data$base = c(1 - (summary(survfit(surv_train_base, newdata = train_data), times = 24)$surv)) 

dca_rfs <- dca(Surv(nutrfs, nutr) ~ model + base, 
               data = train_data, 
               time = 24,
               thresholds = seq(0, 1, by = 0.01),
               label = list(base = "CheckMate 274 risk group",
                            model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "Non-urinary tract recurrence at 24 months") +
  guides(colour = guide_legend(nrow = 2)) 

# training dataset
surv_train <- coxph(Surv(fu, death_all) ~ score, data = train_data)
surv_train_base <- coxph(Surv(fu, death_all) ~ risk_recurrence, data = train_data)

# calculate failure probability at 24 months
train_data$model = c(1 - (summary(survfit(surv_train, newdata = train_data), times = 24)$surv)) 
train_data$base = c(1 - (summary(survfit(surv_train_base, newdata = train_data), times = 24)$surv)) 

dca_os <- dca(Surv(fu, death_all) ~ model + base, 
              data = train_data, 
              time = 24,
              thresholds = seq(0, 1, by = 0.01),
              label = list(base = "CheckMate 274 risk group",
                           model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "All-cause mortality at 24 months") +
  guides(colour = guide_legend(nrow = 2)) 

dca_merge24 <- (dca_rfs + dca_os) + theme_juog()

# 60-month event
# decision curve analysis
# calculate failure probability at 60 months
# training dataset
surv_train <- coxph(Surv(nutrfs, nutr) ~ score, data = train_data)
surv_train_base <- coxph(Surv(nutrfs, nutr) ~ risk_recurrence, data = train_data)

train_data$model = c(1 - (summary(survfit(surv_train, newdata = train_data), times = 60)$surv)) 
train_data$base = c(1 - (summary(survfit(surv_train_base, newdata = train_data), times = 60)$surv)) 

dca_rfs <- dca(Surv(nutrfs, nutr) ~ model + base, 
               data = train_data, 
               time = 60,
               thresholds = seq(0, 1, by = 0.01),
               label = list(base = "CheckMate 274 risk group",
                            model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "Non-urinary tract recurrence at 60 months") +
  guides(colour = guide_legend(nrow = 2)) 

# training dataset
surv_train <- coxph(Surv(fu, death_all) ~ score, data = train_data)
surv_train_base <- coxph(Surv(fu, death_all) ~ risk_recurrence, data = train_data)
# calculate failure probability at 60 months
train_data$model = c(1 - (summary(survfit(surv_train, newdata = train_data), times = 60)$surv)) 
train_data$base = c(1 - (summary(survfit(surv_train_base, newdata = train_data), times = 60)$surv)) 

dca_os <- dca(Surv(fu, death_all) ~ model + base, 
              data = train_data, 
              time = 60,
              thresholds = seq(0, 1, by = 0.01),
              label = list(base = "CheckMate 274 risk group",
                           model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "All-cause mortality at 60 months") +
  guides(colour = guide_legend(nrow = 2)) 

dca_merge60 <- (dca_rfs + dca_os) + theme_juog()

dca_merge <- dca_merge24 / dca_merge60

# harrell"s C
surv_test <- coxph(Surv(nutrfs, nutr) ~ score, data = test_data)
surv_test_base <- coxph(Surv(nutrfs, nutr) ~ risk_recurrence, data = test_data)

# decision curve analysis
# calculate failure probability at 24 months
test_data$model = c(1 - (summary(survfit(surv_test, newdata = test_data), times = 24)$surv)) 
test_data$base = c(1 - (summary(survfit(surv_test_base, newdata = test_data), times = 24)$surv)) 

dca_rfs <- dca(Surv(nutrfs, nutr) ~ model + base, 
               data = test_data, 
               time = 24,
               thresholds = seq(0, 1, by = 0.01),
               label = list(base = "CheckMate 274 risk group",
                            model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "Non-urinary tract recurrence at 24 months") +
  guides(colour = guide_legend(nrow = 2)) 

# testing dataset
surv_test <- coxph(Surv(fu, death_all) ~ score, data = test_data)
surv_test_base <- coxph(Surv(fu, death_all) ~ risk_recurrence, data = test_data)

# calculate failure probability at 24 months
test_data$model = c(1 - (summary(survfit(surv_test, newdata = test_data), times = 24)$surv)) 
test_data$base = c(1 - (summary(survfit(surv_test_base, newdata = test_data), times = 24)$surv)) 

dca_os <- dca(Surv(fu, death_all) ~ model + base, 
              data = test_data, 
              time = 24,
              thresholds = seq(0, 1, by = 0.01),
              label = list(base = "CheckMate 274 risk group",
                           model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "All-cause mortality at 24 months") +
  guides(colour = guide_legend(nrow = 2)) 

dca_merge24 <- (dca_rfs + dca_os) + theme_juog()

# 60-month event
# decision curve analysis
# calculate failure probability at 60 months
surv_test <- coxph(Surv(nutrfs, nutr) ~ score, data = test_data)
surv_test_base <- coxph(Surv(nutrfs, nutr) ~ risk_recurrence, data = test_data)

test_data$model = c(1 - (summary(survfit(surv_test, newdata = test_data), times = 60)$surv)) 
test_data$base = c(1 - (summary(survfit(surv_test_base, newdata = test_data), times = 60)$surv)) 

dca_rfs <- dca(Surv(nutrfs, nutr) ~ model + base, 
               data = test_data, 
               time = 60,
               thresholds = seq(0, 1, by = 0.01),
               label = list(base = "CheckMate 274 risk group",
                            model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "Non-urinary tract recurrence at 60 months") +
  guides(colour = guide_legend(nrow = 2)) 

# testing dataset
surv_test <- coxph(Surv(fu, death_all) ~ score, data = test_data)
surv_test_base <- coxph(Surv(fu, death_all) ~ risk_recurrence, data = test_data)

# calculate failure probability at 60 months
test_data$model = c(1 - (summary(survfit(surv_test, newdata = test_data), times = 60)$surv)) 
test_data$base = c(1 - (summary(survfit(surv_test_base, newdata = test_data), times = 60)$surv)) 

dca_os <- dca(Surv(fu, death_all) ~ model + base, 
              data = test_data, 
              time = 60,
              thresholds = seq(0, 1, by = 0.01),
              label = list(base = "CheckMate 274 risk group at 60 months",
                           model = "The developed model")) %>%
  plot(smooth = TRUE) + ggsci::scale_color_aaas() + 
  theme_juog() + 
  labs(subtitle = "All-cause mortality") +
  guides(colour = guide_legend(nrow = 2)) 

dca_merge60 <- (dca_rfs + dca_os) + theme_juog()

dca_merge <- dca_merge24 / dca_merge60

# calibration plot for Cox model
# non-urinary tract recurrence-free survival
# 2-year survival
cox_nutrfs24_training <- cph(Surv(nutrfs, nutr) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
cox_nutrfs24_test <- cph(Surv(nutrfs, nutr) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
calib_train_nutrfs24 <- calibrate(cox_nutrfs24_training, u = 24, cmethod = "KM", B = 1000, m = 440)
calib_test_nutrfs24 <- calibrate(cox_nutrfs24_test, u = 24, cmethod = "KM", B = 1000, m = 220)
df_calib_train_nutrfs24 <- data.frame(predicted = calib_train_nutrfs24[,7],
                                      observed = calib_train_nutrfs24[,9],
                                      se = calib_train_nutrfs24[,10],
                                      group = "Training")
df_calib_test_nutrfs24 <- data.frame(predicted = calib_test_nutrfs24[,7],
                                     observed = calib_test_nutrfs24[,9],
                                     se = calib_test_nutrfs24[,10],
                                     group = "Validation")
df_calib_nutrfs24 <- bind_rows(df_calib_train_nutrfs24, df_calib_test_nutrfs24)

cplot24_nutrfs <- ggplot(data = df_calib_nutrfs24,
                         aes(x = predicted, y = observed,
                             ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                             color = group, fill = group)) +
  geom_smooth(size = 0.8, method = "loess", se = FALSE) + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Two-year NUTRFS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

# 5-year survival
cox_nutrfs60_training <- cph(Surv(nutrfs, nutr) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
cox_nutrfs60_test <- cph(Surv(nutrfs, nutr) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
calib_train_nutrfs60 <- calibrate(cox_nutrfs60_training, u = 60, cmethod = "KM", B = 1000, m = 440)
calib_test_nutrfs60 <- calibrate(cox_nutrfs60_test, u = 60, cmethod = "KM", B = 1000, m = 220)
df_calib_train_nutrfs60 <- data.frame(predicted = calib_train_nutrfs60[,7],
                                      observed = calib_train_nutrfs60[,9],
                                      se = calib_train_nutrfs60[,10],
                                      group = "Training")
df_calib_test_nutrfs60 <- data.frame(predicted = calib_test_nutrfs60[,7],
                                     observed = calib_test_nutrfs60[,9],
                                     se = calib_test_nutrfs60[,10],
                                     group = "Validation")
df_calib_nutrfs60 <- bind_rows(df_calib_train_nutrfs60, df_calib_test_nutrfs60)

cplot60_nutrfs <- ggplot(data = df_calib_nutrfs60,
                         aes(x = predicted, y = observed,
                             ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                             color = group, fill = group)) + 
  geom_smooth(size = 0.8, method = "loess") + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Five-year NUTRFS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

cplot_nutrfs_merge <- (cplot24_nutrfs + cplot60_nutrfs) + theme_juog()

# Cancer-specific survival
# 2-year survival
cox_css24_training <- cph(Surv(fu, death_all) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
cox_css24_test <- cph(Surv(fu, death_all) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
calib_train_css24 <- calibrate(cox_css24_training, u = 24, cmethod = "KM", B = 1000, m = 440)
calib_test_css24 <- calibrate(cox_css24_test, u = 24, cmethod = "KM", B = 1000, m = 220)
df_calib_train_css24 <- data.frame(predicted = calib_train_css24[,7],
                                   observed = calib_train_css24[,9],
                                   se = calib_train_css24[,10],
                                   group = "Training")
df_calib_test_css24 <- data.frame(predicted = calib_test_css24[,7],
                                  observed = calib_test_css24[,9],
                                  se = calib_test_css24[,10],
                                  group = "Validation")
df_calib_css24 <- bind_rows(df_calib_train_css24, df_calib_test_css24)

cplot24_css <- ggplot(data = df_calib_css24,
                      aes(x = predicted, y = observed,
                          ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                          color = group, fill = group)) + 
  geom_smooth(size = 0.8, method = "loess") + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Two-year CSS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

# 5-year survival
cox_css60_training <- cph(Surv(fu, death_all) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
cox_css60_test <- cph(Surv(fu, death_all) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
calib_train_css60 <- calibrate(cox_css60_training, u = 60, cmethod = "KM", B = 1000, m = 440)
calib_test_css60 <- calibrate(cox_css60_test, u = 60, cmethod = "KM", B = 1000, m = 220)
df_calib_train_css60 <- data.frame(predicted = calib_train_css60[,7],
                                   observed = calib_train_css60[,9],
                                   se = calib_train_css60[,10],
                                   group = "Training")
df_calib_test_css60 <- data.frame(predicted = calib_test_css60[,7],
                                  observed = calib_test_css60[,9],
                                  se = calib_test_css60[,10],
                                  group = "Validation")
df_calib_css60 <- bind_rows(df_calib_train_css60, df_calib_test_css60)

cplot60_css <- ggplot(data = df_calib_css60,
                      aes(x = predicted, y = observed,
                          ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                          color = group, fill = group)) + 
  geom_smooth(size = 0.8, method = "loess") + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Five-year CSS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

# overall survival
# 2-year survival
cox_os24_training <- cph(Surv(fu, death_all) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
cox_os24_test <- cph(Surv(fu, death_all) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
calib_train_os24 <- calibrate(cox_os24_training, u = 24, cmethod = "KM", B = 1000, m = 440)
calib_test_os24 <- calibrate(cox_os24_test, u = 24, cmethod = "KM", B = 1000, m = 220)
df_calib_train_os24 <- data.frame(predicted = calib_train_os24[,7],
                                  observed = calib_train_os24[,9],
                                  se = calib_train_os24[,10],
                                  group = "Training")
df_calib_test_os24 <- data.frame(predicted = calib_test_os24[,7],
                                 observed = calib_test_os24[,9],
                                 se = calib_test_os24[,10],
                                 group = "Validation")
df_calib_os24 <- bind_rows(df_calib_train_os24, df_calib_test_os24)

cplot24_os <- ggplot(data = df_calib_os24,
                     aes(x = predicted, y = observed,
                         ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                         color = group, fill = group)) + 
  geom_smooth(size = 0.8, method = "loess") + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Two-year OS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

# 5-year survival
cox_os60_training <- cph(Surv(fu, death_all) ~ score, data = train_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
cox_os60_test <- cph(Surv(fu, death_all) ~ score, data = test_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)
calib_train_os60 <- calibrate(cox_os60_training, u = 60, cmethod = "KM", B = 1000, m = 440)
calib_test_os60 <- calibrate(cox_os60_test, u = 60, cmethod = "KM", B = 1000, m = 220)
df_calib_train_os60 <- data.frame(predicted = calib_train_os60[,7],
                                  observed = calib_train_os60[,9],
                                  se = calib_train_os60[,10],
                                  group = "Training")
df_calib_test_os60 <- data.frame(predicted = calib_test_os60[,7],
                                 observed = calib_test_os60[,9],
                                 se = calib_test_os60[,10],
                                 group = "Validation")
df_calib_os60 <- bind_rows(df_calib_train_os60, df_calib_test_os60)

cplot60_os <- ggplot(data = df_calib_os60,
                     aes(x = predicted, y = observed,
                         ymin = observed - 1.96 * se, ymax = observed + 1.96 * se,
                         color = group, fill = group)) + 
  geom_smooth(size = 0.8, method = "loess") + 
  geom_point(size = 3) + 
  geom_errorbar(width = 0.02) + 
  theme_juog() + 
  geom_abline(linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "Predicted probability",
    y = "Observed probability with 95% CI\n(Kaplan-Meier)",
    subtitle = "Five-year OS probability"
  ) + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas()

cplot_os_merge <- (cplot24_os + cplot60_os) +   
  plot_annotation(
    title = "Calibration plots",
    subtitle = "Actual vs. predicted probability of events at 24 and 60 months"
  ) & 
  theme_juog()

cplot_merge <-  (cplot24_nutrfs + cplot24_css + cplot24_os) / (cplot60_nutrfs + cplot60_css + cplot60_os) + theme_juog()

cplot24_merge <-  (cplot24_nutrfs + cplot24_css + cplot24_os) + theme_juog()
cplot60_merge <-  (cplot60_nutrfs + cplot60_css + cplot60_os) + theme_juog()

# tAUC comparing 2 risk groups...
# time-dependent ROC analysis using timeROC...
library(timeROC)
# overall survival - training set
troc_os_train_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$fu, # time
                    delta = train_data$death_all, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_os_train_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_os_train_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_os_train_cm274)[1] <- "time"
colnames(plot_tdroc_os_train_cm274)[2] <- "tdauc"
colnames(plot_tdroc_os_train_cm274)[3] <- "lcl"
colnames(plot_tdroc_os_train_cm274)[4] <- "ucl"

troc_os_train_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$fu, # time
                    delta = train_data$death_all, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_os_train_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_os_train_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_os_train_model)[1] <- "time"
colnames(plot_tdroc_os_train_model)[2] <- "tdauc"
colnames(plot_tdroc_os_train_model)[3] <- "lcl"
colnames(plot_tdroc_os_train_model)[4] <- "ucl"

plot_tdroc_os_train_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_os_train_model$tdauc, plot_tdroc_os_train_cm274$tdauc),
  lcl = c(plot_tdroc_os_train_model$lcl, plot_tdroc_os_train_cm274$lcl),
  ucl = c(plot_tdroc_os_train_model$ucl, plot_tdroc_os_train_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_train_os_plot <- 
  ggplot(plot_tdroc_os_train_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "OS",
    subtitle = "Training set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

# overall survival - validation set
troc_os_test_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$fu, # time
                    delta = test_data$death_all, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_os_test_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_os_test_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_os_test_cm274)[1] <- "time"
colnames(plot_tdroc_os_test_cm274)[2] <- "tdauc"
colnames(plot_tdroc_os_test_cm274)[3] <- "lcl"
colnames(plot_tdroc_os_test_cm274)[4] <- "ucl"

troc_os_test_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$fu, # time
                    delta = test_data$death_all, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_os_test_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_os_test_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_os_test_model)[1] <- "time"
colnames(plot_tdroc_os_test_model)[2] <- "tdauc"
colnames(plot_tdroc_os_test_model)[3] <- "lcl"
colnames(plot_tdroc_os_test_model)[4] <- "ucl"

plot_tdroc_os_test_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_os_test_model$tdauc, plot_tdroc_os_test_cm274$tdauc),
  lcl = c(plot_tdroc_os_test_model$lcl, plot_tdroc_os_test_cm274$lcl),
  ucl = c(plot_tdroc_os_test_model$ucl, plot_tdroc_os_test_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_test_os_plot <- 
  ggplot(plot_tdroc_os_test_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "OS",
    subtitle = "Validation set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

# cancer-specific survival - training set
troc_css_train_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$fu, # time
                    delta = train_data$death_cancer, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_css_train_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_css_train_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_css_train_cm274)[1] <- "time"
colnames(plot_tdroc_css_train_cm274)[2] <- "tdauc"
colnames(plot_tdroc_css_train_cm274)[3] <- "lcl"
colnames(plot_tdroc_css_train_cm274)[4] <- "ucl"

troc_css_train_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$fu, # time
                    delta = train_data$death_cancer, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_css_train_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_css_train_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_css_train_model)[1] <- "time"
colnames(plot_tdroc_css_train_model)[2] <- "tdauc"
colnames(plot_tdroc_css_train_model)[3] <- "lcl"
colnames(plot_tdroc_css_train_model)[4] <- "ucl"

plot_tdroc_css_train_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_css_train_model$tdauc, plot_tdroc_css_train_cm274$tdauc),
  lcl = c(plot_tdroc_css_train_model$lcl, plot_tdroc_css_train_cm274$lcl),
  ucl = c(plot_tdroc_css_train_model$ucl, plot_tdroc_css_train_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_train_css_plot <- 
  ggplot(plot_tdroc_css_train_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "CSS",
    subtitle = "Training set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

# overall survival - validation set
troc_css_test_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$fu, # time
                    delta = test_data$death_cancer, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_css_test_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_css_test_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_css_test_cm274)[1] <- "time"
colnames(plot_tdroc_css_test_cm274)[2] <- "tdauc"
colnames(plot_tdroc_css_test_cm274)[3] <- "lcl"
colnames(plot_tdroc_css_test_cm274)[4] <- "ucl"

troc_css_test_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$fu, # time
                    delta = test_data$death_cancer, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_css_test_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_css_test_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_css_test_model)[1] <- "time"
colnames(plot_tdroc_css_test_model)[2] <- "tdauc"
colnames(plot_tdroc_css_test_model)[3] <- "lcl"
colnames(plot_tdroc_css_test_model)[4] <- "ucl"

plot_tdroc_css_test_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_css_test_model$tdauc, plot_tdroc_css_test_cm274$tdauc),
  lcl = c(plot_tdroc_css_test_model$lcl, plot_tdroc_css_test_cm274$lcl),
  ucl = c(plot_tdroc_css_test_model$ucl, plot_tdroc_css_test_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_test_css_plot <- 
  ggplot(plot_tdroc_css_test_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "CSS",
    subtitle = "Validation set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

# non-urinary tract recurrence-free survival - training set
troc_nutrfs_train_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$nutrfs, # time
                    delta = train_data$nutr, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_nutrfs_train_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_nutrfs_train_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_nutrfs_train_cm274)[1] <- "time"
colnames(plot_tdroc_nutrfs_train_cm274)[2] <- "tdauc"
colnames(plot_tdroc_nutrfs_train_cm274)[3] <- "lcl"
colnames(plot_tdroc_nutrfs_train_cm274)[4] <- "ucl"

troc_nutrfs_train_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = train_data$nutrfs, # time
                    delta = train_data$nutr, # event
                    cause = 1, # specify the event of interest
                    marker = train_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_nutrfs_train_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_nutrfs_train_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_nutrfs_train_model)[1] <- "time"
colnames(plot_tdroc_nutrfs_train_model)[2] <- "tdauc"
colnames(plot_tdroc_nutrfs_train_model)[3] <- "lcl"
colnames(plot_tdroc_nutrfs_train_model)[4] <- "ucl"

plot_tdroc_nutrfs_train_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_nutrfs_train_model$tdauc, plot_tdroc_nutrfs_train_cm274$tdauc),
  lcl = c(plot_tdroc_nutrfs_train_model$lcl, plot_tdroc_nutrfs_train_cm274$lcl),
  ucl = c(plot_tdroc_nutrfs_train_model$ucl, plot_tdroc_nutrfs_train_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_train_nutrfs_plot <- 
  ggplot(plot_tdroc_nutrfs_train_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "NUTRFS",
    subtitle = "Training set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

# overall survival - validation set
troc_nutrfs_test_cm274 <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$nutrfs, # time
                    delta = test_data$nutr, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$risk_recurrence, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_nutrfs_test_cm274(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_nutrfs_test_cm274 <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_nutrfs_test_cm274)[1] <- "time"
colnames(plot_tdroc_nutrfs_test_cm274)[2] <- "tdauc"
colnames(plot_tdroc_nutrfs_test_cm274)[3] <- "lcl"
colnames(plot_tdroc_nutrfs_test_cm274)[4] <- "ucl"

troc_nutrfs_test_model <- function(cut){
  cutoff <- cut
  
  # compute time-dependent ROC analysis...
  survroc = timeROC(T = test_data$nutrfs, # time
                    delta = test_data$nutr, # event
                    cause = 1, # specify the event of interest
                    marker = test_data$score, 
                    times = cutoff, iid = TRUE)
  
  # estimate confidence intervals...
  confint <- confint(survroc, level = 0.95, n.sim = 1000)
  
  output <- tibble(
    tdauc = 0, 
    lcl = 0, 
    ucl = 0
  )
  
  output$tdauc <- survroc$AUC[2] %>% round(digits = 2)
  output$lcl <- (confint$CI_AUC[1] / 100) %>% round(digits = 2)
  output$ucl <- (confint$CI_AUC[2] / 100) %>% round(digits = 2)
  
  print(output)
}

data_tdroc <- data.frame()
for(i in seq(from = 12, to = 96, by = 6)) {
  data_tdroc = rbind(data_tdroc, troc_nutrfs_test_model(i))
}

tau_vec <- as.data.frame(seq(from = 12, to = 96, by = 6))

plot_tdroc_nutrfs_test_model <- cbind(tau_vec, data_tdroc)
colnames(plot_tdroc_nutrfs_test_model)[1] <- "time"
colnames(plot_tdroc_nutrfs_test_model)[2] <- "tdauc"
colnames(plot_tdroc_nutrfs_test_model)[3] <- "lcl"
colnames(plot_tdroc_nutrfs_test_model)[4] <- "ucl"

plot_tdroc_nutrfs_test_merge <- data.frame(
  month = rep(seq(from = 12, to = 96, by = 6), 2),
  tdauc = c(plot_tdroc_nutrfs_test_model$tdauc, plot_tdroc_nutrfs_test_cm274$tdauc),
  lcl = c(plot_tdroc_nutrfs_test_model$lcl, plot_tdroc_nutrfs_test_cm274$lcl),
  ucl = c(plot_tdroc_nutrfs_test_model$ucl, plot_tdroc_nutrfs_test_cm274$ucl),
  strat = rep(c("The developed model", "CheckMate 274 risk group"), c(15, 15))
)

tdroc_test_nutrfs_plot <- 
  ggplot(plot_tdroc_nutrfs_test_merge, 
         aes(x = month, y = tdauc, ymin = lcl, ymax = ucl,
             group = strat, color = strat, fill = strat)) + 
  scale_x_continuous(breaks = seq(from = 12, to = 96, by = 6)) + 
  geom_line(size = 0.7) + 
  geom_errorbar(width = 1) + 
  geom_point(size = 2) + 
  labs(
    title = "NUTRFS",
    subtitle = "Validation set",
    x = "Months since radical cystectomy",
    y = "Time-dependent AUC (95% CI)"
  ) + 
  theme_juog() + 
  ggsci::scale_color_aaas() + ggsci::scale_fill_aaas() + 
  ylim(0.35, 1)

merge_figure <- (tdroc_train_nutrfs_plot + tdroc_train_css_plot + tdroc_train_os_plot) / (tdroc_test_nutrfs_plot + tdroc_test_css_plot + tdroc_test_os_plot) + theme_juog()
merge_figure_validation <- (tdroc_test_nutrfs_plot + tdroc_test_css_plot + tdroc_test_os_plot) + theme_juog()

# subgroup analysis
# high-risk populations
df_hr <- df_score %>% 
  filter(risk_recurrence == 1) 

# stratified Kaplan-Meier curves
rfs_hr_cm274_fit <- survival::survfit(data = df_hr, Surv(nutrfs, nutr) ~ score_tertile)
rfs_hr_cm274 <- survminer::ggsurvplot(rfs_hr_cm274_fit, 
                                      risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                      title = "NUTRFS", 
                                      subtitle = "CheckMate 274 high-risk population",
                                      legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                      legend.title = "The developed model", 
                                      censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                      xlab = "Months since radical cystectomy", 
                                      ylab = "NUTRFS probability (95% CI)",
                                      palette = "aaas", size = 0.5,  break.time.by = 12,
                                      ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                      risk.table.col = "strata",
                                      tables.height = 0.12, risk.table.fontsize = 2.8,
                                      conf.int = TRUE) 

rfs_hr_cm274$plot <- rfs_hr_cm274$plot + 
  theme_km() 

rfs_hr_cm274$table <- rfs_hr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# cancer-specific survival
css_hr_cm274_fit <- survival::survfit(data = df_hr, Surv(fu, death_cancer) ~ score_tertile)
css_hr_cm274 <- survminer::ggsurvplot(css_hr_cm274_fit, 
                                      risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                      title = "CSS", 
                                      subtitle = "CheckMate 274 high-risk population",
                                      legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                      legend.title = "The developed model", 
                                      censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                      xlab = "Months since radical cystectomy", 
                                      ylab = "CSS probability (95% CI)",
                                      palette = "aaas", size = 0.5,  break.time.by = 12,
                                      ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                      risk.table.col = "strata",
                                      tables.height = 0.12, risk.table.fontsize = 2.8,
                                      conf.int = TRUE) 

css_hr_cm274$plot <- css_hr_cm274$plot + 
  theme_km()

css_hr_cm274$table <- css_hr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
os_hr_cm274_fit <- survival::survfit(data = df_hr, Surv(fu, death_all) ~ score_tertile)
os_hr_cm274 <- survminer::ggsurvplot(os_hr_cm274_fit, 
                                     risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                     title = "OS", 
                                     subtitle = "CheckMate 274 high-risk population",
                                     legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                     legend.title = "The developed model", 
                                     censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                     xlab = "Months since radical cystectomy", 
                                     ylab = "OS probability (95% CI)",
                                     palette = "aaas", size = 0.5,  break.time.by = 12,
                                     ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                     risk.table.col = "strata",
                                     tables.height = 0.12, risk.table.fontsize = 2.8,
                                     conf.int = TRUE) 

os_hr_cm274$plot <- os_hr_cm274$plot + 
  theme_km()

os_hr_cm274$table <- os_hr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# low-risk populations
df_lr <- df_score %>% 
  filter(risk_recurrence == 0) 

# stratified Kaplan-Meier curves
rfs_lr_cm274_fit <- survival::survfit(data = df_lr, Surv(nutrfs, nutr) ~ score_tertile)
rfs_lr_cm274 <- survminer::ggsurvplot(rfs_lr_cm274_fit, 
                                      risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                      title = "NUTRFS", 
                                      subtitle = "CheckMate 274 low-risk population",
                                      legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                      legend.title = "The developed model", 
                                      censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                      xlab = "Months since radical cystectomy", 
                                      ylab = "NUTRFS probability (95% CI)",
                                      palette = "aaas", size = 0.5,  break.time.by = 12,
                                      ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                      risk.table.col = "strata",
                                      tables.height = 0.12, risk.table.fontsize = 2.8,
                                      conf.int = TRUE) 

rfs_lr_cm274$plot <- rfs_lr_cm274$plot + 
  theme_km() 

rfs_lr_cm274$table <- rfs_lr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# cancer-specific survival
css_lr_cm274_fit <- survival::survfit(data = df_lr, Surv(fu, death_all) ~ score_tertile)
css_lr_cm274 <- survminer::ggsurvplot(css_lr_cm274_fit, 
                                      risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                      title = "CSS", 
                                      subtitle = "CheckMate 274 low-risk population",
                                      legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                      legend.title = "The developed model", 
                                      censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                      xlab = "Months since radical cystectomy", 
                                      ylab = "CSS probability (95% CI)",
                                      palette = "aaas", size = 0.5,  break.time.by = 12,
                                      ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                      risk.table.col = "strata",
                                      tables.height = 0.12, risk.table.fontsize = 2.8,
                                      conf.int = TRUE) 

css_lr_cm274$plot <- css_lr_cm274$plot + 
  theme_km()

css_lr_cm274$table <- css_lr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# overall survival
os_lr_cm274_fit <- survival::survfit(data = df_lr, Surv(fu, death_all) ~ score_tertile)
os_lr_cm274 <- survminer::ggsurvplot(os_lr_cm274_fit, 
                                     risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                     title = "OS", 
                                     subtitle = "CheckMate 274 low-risk population",
                                     legend = "top", legend.labs = c("Good risk", "Poor risk"),
                                     legend.title = "The developed model", 
                                     censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                     xlab = "Months since radical cystectomy", 
                                     ylab = "OS probability (95% CI)",
                                     palette = "aaas", size = 0.5,  break.time.by = 12,
                                     ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                     risk.table.col = "strata",
                                     tables.height = 0.12, risk.table.fontsize = 2.8,
                                     conf.int = TRUE) 

os_lr_cm274$plot <- os_lr_cm274$plot + 
  theme_km()

os_lr_cm274$table <- os_lr_cm274$table + 
  theme_void() + 
  theme(
    text = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) # replace table

# summary of subpopulations
merge <- survminer::arrange_ggsurvplots(list(rfs_hr_cm274, rfs_lr_cm274, css_hr_cm274, css_lr_cm274, os_hr_cm274, os_lr_cm274), 
                                        print = FALSE, 
                                        nrow = 2, ncol = 3)

merge_hr <- survminer::arrange_ggsurvplots(list(rfs_hr_cm274, css_hr_cm274, os_hr_cm274), 
                                           print = FALSE, 
                                           nrow = 1, ncol = 3)

# survival probabilities of subgroups
survfit_cm274_good <- os_hr_cm274_fit %>% 
  tbl_survfit(times = c(24, 60), 
              label_header = "**{time} Months**",
              label = "The developed model")

# cross tables
tbl_cross <- df_score %>% 
  mutate(
    risk_recurrence = case_when(risk_recurrence == 1 ~ "high risk", TRUE ~ "low risk") %>% factor(levels = c("low risk", "high risk"))
  ) %>% 
  tbl_cross(row = risk_recurrence, col = score_tertile,
            statistic = "{n} ({p})", 
            percent = "cell",
            label = list(
              score_tertile ~ "The developed model",
              risk_recurrence ~ "CheckMate 274 risk group"
            )
  ) %>% 
  bold_labels()