# =========================================================
# 0) Packages
# =========================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  lme4,
  lmerTest,
  tidyverse,
  lubridate,
  emmeans,
  ggeffects,
  performance,
  htmltools,
  knitr
)

options(
  dplyr.summarise.inform = FALSE,
  contrasts = c("contr.sum", "contr.poly")
)

emmeans::emm_options(lmer.df = "satterthwaite")

# =========================================================
# 1) Global settings
# =========================================================
input_file  <- "temp_for_model_200.csv"
output_html <- "paper_style_model_report_social_context_3.html"

plot_m4_file <- "plot_model4_age_moderation.png"
plot_m5b_social_context_file <- "plot_model5b_social_context_3_moderation.png"
plot_social_context_by_age_file <- "plot_social_context_3_by_age.png"

# =========================================================
# 2) Helpers
# =========================================================
recode_yes_no <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    x_chr %in% c("1", "Yes", "yes", "YES", "Ja", "ja", "TRUE", "true", "T") ~ 1,
    x_chr %in% c("0", "No", "no", "NO", "Nein", "nein", "FALSE", "false", "F") ~ 0,
    TRUE ~ NA_real_
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < .001, "< .001", paste0("= ", formatC(p, digits = 3, format = "f")))
  )
}

yes_no <- function(x) ifelse(isTRUE(x), "Yes", ifelse(isFALSE(x), "No", NA_character_))

safe_r2 <- function(model) {
  out <- try(performance::r2_nakagawa(model), silent = TRUE)
  if (inherits(out, "try-error")) {
    return(data.frame(Marginal_R2 = NA_real_, Conditional_R2 = NA_real_))
  }
  data.frame(
    Marginal_R2 = unname(out$R2_marginal),
    Conditional_R2 = unname(out$R2_conditional)
  )
}

safe_icc <- function(model) {
  out <- try(performance::icc(model), silent = TRUE)
  if (inherits(out, "try-error")) {
    return(data.frame(ICC_adjusted = NA_real_, ICC_unadjusted = NA_real_))
  }
  data.frame(
    ICC_adjusted = unname(out$ICC_adjusted),
    ICC_unadjusted = unname(out$ICC_unadjusted)
  )
}

safe_converged <- function(model) {
  if (inherits(model, "merMod")) {
    return(is.null(model@optinfo$conv$lme4$messages))
  }
  NA
}

safe_singular <- function(model) {
  if (inherits(model, "merMod")) {
    return(tryCatch(lme4::isSingular(model, tol = 1e-4), error = function(e) NA))
  }
  NA
}

fit_stats <- function(model) {
  nobs_value <- tryCatch(
    stats::nobs(model),
    error = function(e) {
      if ("fitted.values" %in% names(model)) {
        return(NROW(model$fitted.values))
      }
      tryCatch(nrow(stats::model.frame(model)), error = function(e2) NA_integer_)
    }
  )
  
  data.frame(
    AIC = tryCatch(AIC(model), error = function(e) NA_real_),
    BIC = tryCatch(BIC(model), error = function(e) NA_real_),
    logLik = tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_),
    nobs = nobs_value
  )
}

diag_stats <- function(model) {
  data.frame(
    Converged = yes_no(safe_converged(model)),
    Singular = yes_no(safe_singular(model))
  )
}

extract_lrt <- function(model_old, model_new, old_name = NULL, new_name = NULL) {
  tab <- as.data.frame(anova(model_old, model_new))
  tab$Model <- rownames(tab)
  rownames(tab) <- NULL
  out <- tab[2, , drop = FALSE]
  if (!is.null(old_name)) out$Compared_with <- old_name
  if (!is.null(new_name)) out$Current_model <- new_name
  out
}

extract_lmm_type3 <- function(model) {
  tab <- as.data.frame(anova(model, type = 3))
  tab$Term <- rownames(tab)
  rownames(tab) <- NULL
  dplyr::select(tab, Term, NumDF, DenDF, `F value`, `Pr(>F)`)
}

extract_ranova <- function(model) {
  tab <- as.data.frame(lmerTest::ranova(model))
  tab$Random_term <- rownames(tab)
  rownames(tab) <- NULL
  needed <- intersect(c("Random_term", "npar", "logLik", "AIC", "LRT", "Df", "Pr(>Chisq)"), names(tab))
  tab[, needed, drop = FALSE]
}

extract_emmeans_lmm <- function(model, factor_name) {
  emm <- emmeans(model, specs = as.formula(paste0("~ ", factor_name)), weights = "proportional")
  as.data.frame(summary(emm, infer = c(TRUE, TRUE)))
}

extract_effect_contrasts_lmm <- function(model, factor_name) {
  emm <- emmeans(model, specs = as.formula(paste0("~ ", factor_name)), weights = "proportional")
  eff <- contrast(emm, "eff")
  as.data.frame(summary(eff, infer = c(TRUE, TRUE)))
}

extract_pairwise_lmm <- function(model, factor_name) {
  emm <- emmeans(model, specs = as.formula(paste0("~ ", factor_name)), weights = "proportional")
  as.data.frame(summary(pairs(emm), infer = c(TRUE, TRUE)))
}

extract_trends <- function(model, specs, var) {
  tr <- emtrends(model, specs = as.formula(paste0("~ ", specs)), var = var)
  as.data.frame(summary(tr, infer = c(TRUE, TRUE)))
}

extract_trend_pairs <- function(model, specs, var) {
  tr <- emtrends(model, specs = as.formula(paste0("~ ", specs)), var = var)
  as.data.frame(summary(pairs(tr), infer = c(TRUE, TRUE)))
}

make_ggeffects_term <- function(var, x, probs = c(.05, .25, .50, .75, .95), digits = 3) {
  vals <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  vals <- unique(round(vals, digits))
  paste0(var, " [", paste(vals, collapse = ","), "]")
}

model_comp_sentence <- function(lrt_df, old_label, new_label) {
  chi <- lrt_df$Chisq[1]
  df  <- lrt_df$Df[1]
  p_col <- names(lrt_df)[grepl("Pr\\(>Chisq\\)", names(lrt_df))]
  p   <- if (length(p_col) == 1) lrt_df[[p_col]][1] else NA_real_
  
  if (is.na(p)) {
    return(paste0(
      "Compared with ", old_label, ", ", new_label,
      " could not be formally compared using a likelihood-ratio test."
    ))
  }
  
  paste0(
    "Compared with ", old_label, ", ", new_label,
    ifelse(p < .05, " improved overall model fit", " did not significantly improve overall model fit"),
    " (\u03C7\u00B2(", df, ") = ", fmt_num(chi, 2), ", p ", fmt_p(p), ")."
  )
}

format_results_df <- function(df) {
  out <- as.data.frame(df)
  if (nrow(out) == 0) return(out)
  
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      if (grepl("Pr\\(|p.value|p_value|^p$", nm, ignore.case = TRUE)) {
        out[[nm]] <- ifelse(
          is.na(out[[nm]]),
          NA_character_,
          ifelse(out[[nm]] < .001, "< .001", formatC(out[[nm]], digits = 3, format = "f"))
        )
      } else {
        out[[nm]] <- round(out[[nm]], 3)
      }
    }
  }
  out
}

table_html <- function(df, caption = NULL) {
  df2 <- format_results_df(df)
  HTML(
    knitr::kable(
      df2,
      format = "html",
      escape = FALSE,
      caption = caption,
      table.attr = 'class="paper-table"'
    )
  )
}

formula_block <- function(x) {
  tags$pre(class = "formula-block", HTML(x))
}

section_block <- function(title, formula_text, paragraphs = NULL, tables = list(), image_file = NULL, image_caption = NULL) {
  children <- list(
    tags$section(
      tags$h2(title),
      tags$h3("Model formula"),
      formula_block(formula_text)
    )
  )
  
  if (!is.null(paragraphs)) {
    for (p in paragraphs) {
      children[[length(children) + 1]] <- tags$p(p)
    }
  }
  
  if (length(tables) > 0) {
    for (tb in tables) {
      children[[length(children) + 1]] <- tags$h3(tb$title)
      children[[length(children) + 1]] <- tags$div(class = "table-wrap", table_html(tb$data, tb$caption))
    }
  }
  
  if (!is.null(image_file)) {
    children[[length(children) + 1]] <- tags$div(
      class = "figure-wrap",
      tags$img(src = image_file, class = "figure-img"),
      if (!is.null(image_caption)) tags$div(class = "figure-caption", image_caption)
    )
  }
  
  tagList(children)
}

# =========================================================
# 3) Read data
# =========================================================
df_raw <- read.csv(input_file, stringsAsFactors = FALSE)

if ("...1" %in% names(df_raw)) {
  df_raw <- df_raw %>% select(-`...1`)
}

# =========================================================
# 4) Preprocess with stepwise filtering log
# =========================================================

# -------------------------
# 4.0 Helper for sample-size logging
# -------------------------

preprocessing_log <- tibble(
  step = character(),
  n_rows = integer(),
  n_participants = integer()
)

log_sample <- function(data, step) {
  n_pid <- if ("PID" %in% names(data)) {
    dplyr::n_distinct(data$PID, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  row <- tibble(
    step = step,
    n_rows = nrow(data),
    n_participants = n_pid
  )
  
  preprocessing_log <<- bind_rows(preprocessing_log, row)
  
  message(
    step, ": ",
    format(nrow(data), big.mark = ","), " rows from ",
    format(n_pid, big.mark = ","), " participants"
  )
  
  invisible(data)
}

filter_step <- function(data, step, ...) {
  out <- data %>% filter(...)
  log_sample(out, step)
  out
}


# -------------------------
# 4.1 Define social-context columns
# -------------------------

all_people_cols_num <- c(
  "pp_partner_num",
  "pp_friends_num",
  "pp_family_num",
  "pp_acquaintances_colleagues_num",
  "pp_other_num",
  "pp_strangers_num"
)

close_tie_cols_num <- c(
  "pp_partner_num",
  "pp_friends_num",
  "pp_family_num"
)


# -------------------------
# 4.2 Initial log
# -------------------------

log_sample(df_raw, "0 Raw data")


# -------------------------
# 4.3 Recode and derive variables
# -------------------------

df <- df_raw %>%
  mutate(
    PID = factor(PID),
    timeStart = as.POSIXct(timeStart, tz = "Europe/Zurich"),
    
    Rating = as.numeric(Rating),
    mean_qi = as.numeric(mean_qi),
    duration_s = as.numeric(duration_s),
    
    AgeClass = trimws(as.character(AgeClass)),
    blsex = trimws(as.character(blsex)),
    bl_home_type = trimws(as.character(bl_home_type)),
    
    pp_other_num = recode_yes_no(pp_other),
    pp_acquaintances_colleagues_num = recode_yes_no(pp_acquaintances_colleagues),
    pp_family_num = recode_yes_no(pp_family),
    pp_strangers_num = recode_yes_no(pp_strangers),
    pp_friends_num = recode_yes_no(pp_friends),
    pp_partner_num = recode_yes_no(pp_partner),
    pp_alone_num = recode_yes_no(pp_alone)
  ) %>%
  mutate(
    time_of_day = case_when(
      hour(timeStart) >= 6  & hour(timeStart) < 12 ~ "morning",
      hour(timeStart) >= 12 & hour(timeStart) < 18 ~ "afternoon",
      !is.na(timeStart) ~ "evening_night",
      TRUE ~ NA_character_
    ),
    weekday_weekend = case_when(
      is.na(timeStart) ~ NA_character_,
      wday(timeStart, week_start = 1) >= 6 ~ "weekend",
      TRUE ~ "weekday"
    ),
    AgeClass = case_when(
      AgeClass %in% c("1", "Jung", "young", "Young", "younger", "Younger") ~ "younger adults",
      AgeClass %in% c("0", "Alt", "older", "Older", "old", "Old") ~ "older adults",
      TRUE ~ NA_character_
    ),
    blsex = case_when(
      blsex %in% c("1", "Männlich", "male", "Male") ~ "male",
      blsex %in% c("2", "Weiblich", "female", "Female") ~ "female",
      blsex %in% c("3", "Divers", "diverse", "Diverse") ~ "diverse",
      TRUE ~ NA_character_
    ),
    bl_home_type = case_when(
      bl_home_type %in% c("1", "Urban", "urban", "Urban (sehr städtisch, zentrumsnah)") ~ "urban",
      bl_home_type %in% c("2", "Sub-urban", "suburban", "Sub-urban (Agglomeration)") ~ "suburban",
      bl_home_type %in% c("3", "Ländlich", "rural", "Rural") ~ "rural",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    n_people_flags = rowSums(across(all_of(all_people_cols_num)), na.rm = TRUE),
    n_close_tie_flags = rowSums(across(all_of(close_tie_cols_num)), na.rm = TRUE),
    
    strict_alone = pp_alone_num == 1 & n_people_flags == 0,
    has_close_tie = n_close_tie_flags > 0,
    valid_social_context = strict_alone | has_close_tie | n_people_flags > 0 | pp_alone_num == 0,
    
    social_context_3 = factor(
      case_when(
        has_close_tie ~ "partner, family and friends",
        strict_alone ~ "alone",
        valid_social_context ~ "other social contacts",
        TRUE ~ NA_character_
      ),
      levels = c("alone", "partner, family and friends", "other social contacts")
    ),
    
    AgeClass = factor(AgeClass, levels = c("younger adults", "older adults")),
    blsex = factor(blsex, levels = c("male", "female", "diverse")),
    bl_home_type = factor(bl_home_type, levels = c("urban", "suburban", "rural")),
    time_of_day = factor(time_of_day, levels = c("morning", "afternoon", "evening_night")),
    weekday_weekend = factor(weekday_weekend, levels = c("weekday", "weekend"))
  )

log_sample(df, "1 After recoding and deriving variables")


# -------------------------
# 4.4 Select modelling variables
# -------------------------

df_model <- df %>%
  select(
    PID, Rating, mean_qi, duration_s,
    AgeClass, blsex, bl_home_type,
    time_of_day, weekday_weekend,
    social_context_3
  )

log_sample(df_model, "2 After selecting modelling variables")


# -------------------------
# 4.5 Stepwise filtering
# -------------------------

df_model <- df_model %>%
  filter_step("3 Keep non-missing PID", !is.na(PID)) %>%
  filter_step("4 Keep non-missing Rating", !is.na(Rating)) %>%
  filter_step("5 Keep non-missing greenspace exposure", !is.na(mean_qi)) %>%
  filter_step("6 Keep non-missing duration", !is.na(duration_s)) %>%
  filter_step("7 Keep non-negative duration", duration_s >= 0) %>%
  filter_step("8 Keep non-missing age group", !is.na(AgeClass)) %>%
  filter_step("9 Keep non-missing sex", !is.na(blsex)) %>%
  filter_step("10 Keep non-missing home type", !is.na(bl_home_type)) %>%
  filter_step("11 Keep non-missing time of day", !is.na(time_of_day)) %>%
  filter_step("12 Keep non-missing weekday/weekend", !is.na(weekday_weekend)) %>%
  filter_step("13 Keep non-missing presence of others", !is.na(social_context_3))


# -------------------------
# 4.6 Keep participants with at least two episodes
# -------------------------

df_model <- df_model %>%
  group_by(PID) %>%
  filter(n() >= 2) %>%
  ungroup()

log_sample(df_model, "14 Keep participants with at least two episodes")


# -------------------------
# 4.7 Duration transformation
# -------------------------

df_model <- df_model %>%
  mutate(
    duration_log = log1p(duration_s)
  )

log_sample(df_model, "15 Final analytic sample after duration-log transformation")


# -------------------------
# 4.8 Print full preprocessing log
# -------------------------

print(preprocessing_log)

# =========================================================
# 5) Within-between decomposition for greenspace exposure
# =========================================================
df_model <- df_model %>%
  group_by(PID) %>%
  mutate(
    gse_person_mean_raw = mean(mean_qi, na.rm = TRUE),
    gseW_raw = mean_qi - gse_person_mean_raw
  ) %>%
  ungroup() %>%
  mutate(
    gseB_raw = gse_person_mean_raw - mean(gse_person_mean_raw, na.rm = TRUE)
  )

within_check <- df_model %>%
  group_by(PID) %>%
  summarise(mean_gseW = mean(gseW_raw, na.rm = TRUE), .groups = "drop")

# =========================================================
# 5.1) Ensure categorical variables are factors before modelling
# =========================================================
df_model <- df_model %>%
  mutate(
    PID = factor(PID),
    AgeClass = factor(AgeClass, levels = c("younger adults", "older adults")),
    blsex = factor(blsex, levels = c("male", "female", "diverse")),
    bl_home_type = factor(bl_home_type, levels = c("urban", "suburban", "rural")),
    time_of_day = factor(time_of_day, levels = c("morning", "afternoon", "evening_night")),
    weekday_weekend = factor(weekday_weekend, levels = c("weekday", "weekend")),
    social_context_3 = factor(social_context_3, levels = c("alone", "partner, family and friends", "other social contacts"))
  )

# =========================================================
# 6) Model controls
# =========================================================
ctrl_lmer <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 2e5)
)

# =========================================================
# 7) Fit models
# =========================================================
# Null model
mnull <- lmer(
  Rating ~ 1 + (1 | PID),
  data = df_model,
  REML = FALSE,
  control = ctrl_lmer
)

# Model 0 comparison model: baseline adjusted model with random intercept only
# This model is used only to test whether adding a random slope improves fit.
m0_random_intercept <- lmer(
  Rating ~ gseW_raw + gseB_raw + duration_s + AgeClass +
    blsex + bl_home_type + time_of_day + weekday_weekend +
    (1 | PID),
  data = df_model,
  REML = FALSE,
  control = ctrl_lmer
)

# Model 1: baseline adjusted model with random intercept and random slope
m1 <- lmer(
  Rating ~ gseW_raw + gseB_raw + duration_s + AgeClass +
    blsex + bl_home_type + time_of_day + weekday_weekend +
    (1 + gseW_raw | PID),
  data = df_model,
  REML = FALSE,
  control = ctrl_lmer
)

# Model 2: age moderation
m2 <- lmer(
  Rating ~ gseW_raw * AgeClass + gseB_raw + duration_s +
    blsex + bl_home_type + time_of_day + weekday_weekend +
    (1 + gseW_raw | PID),
  data = df_model,
  REML = FALSE,
  control = ctrl_lmer
)

# Model 3: presence-of-others moderation
# Presence of others: alone, partner, family and friends, other social contacts.
# No three-way interaction is included.
m3 <- lmer(
  Rating ~ gseW_raw * AgeClass + gseB_raw + duration_s +
    gseW_raw * social_context_3 +
    blsex + bl_home_type + time_of_day + weekday_weekend +
    (1 + gseW_raw | PID),
  data = df_model,
  REML = FALSE,
  control = ctrl_lmer
)

# =========================================================
# 8) Descriptive analytics
# =========================================================

# -----------------------------
# Helper functions
# -----------------------------

count_pct <- function(data, var, denom = NULL) {
  var <- rlang::ensym(var)
  
  out <- data %>%
    count(!!var, name = "n") %>%
    mutate(
      percent = 100 * n / sum(n)
    )
  
  out
}

summarise_cont <- function(data, vars) {
  data %>%
    summarise(
      across(
        all_of(vars),
        list(
          n = ~sum(!is.na(.x)),
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          min = ~min(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = c("variable", ".value"),
      names_pattern = "(.+)_(n|mean|sd|median|min|max)"
    )
}

# =========================================================
# 8) Descriptive analytics
# =========================================================

# -----------------------------
# 8.0 Add reporting variables
# -----------------------------

df_model <- df_model %>%
  mutate(
    duration_min = duration_s / 60,
    AgeClass = factor(AgeClass, levels = c("younger adults", "older adults")),
    blsex = factor(blsex, levels = c("male", "female", "diverse")),
    bl_home_type = factor(bl_home_type, levels = c("urban", "suburban", "rural")),
    social_context_3 = factor(social_context_3, levels = c("alone", "partner, family and friends", "other social contacts")),
    weekday_weekend = factor(weekday_weekend, levels = c("weekday", "weekend")),
    time_of_day = factor(time_of_day, levels = c("morning", "afternoon", "evening_night"))
  )

# Participant-level table:
# age group, gender, and home type should be counted once per participant.
participant_tbl <- df_model %>%
  distinct(PID, AgeClass, blsex, bl_home_type, gse_person_mean_raw, gseB_raw)

# -----------------------------
# 8.1 Helper functions
# -----------------------------

fmt_n <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

fmt_pct1 <- function(x) {
  formatC(x, digits = 1, format = "f")
}

fmt_num2 <- function(x) {
  formatC(x, digits = 2, format = "f")
}

fmt_n_pct <- function(n, denom) {
  if (is.na(n) || is.na(denom) || denom == 0) return("--")
  paste0(fmt_n(n), " (", fmt_pct1(100 * n / denom), ")")
}

fmt_mean_sd <- function(mean_value, sd_value) {
  if (is.na(mean_value) || is.na(sd_value)) return("--")
  paste0(fmt_num2(mean_value), " (", fmt_num2(sd_value), ")")
}

make_participant_cat_row <- function(var_name, value, label) {
  var_sym <- rlang::sym(var_name)
  
  overall_denom <- nrow(participant_tbl)
  young_denom <- participant_tbl %>% filter(AgeClass == "younger adults") %>% nrow()
  older_denom <- participant_tbl %>% filter(AgeClass == "older adults") %>% nrow()
  
  overall_n <- participant_tbl %>%
    filter(!!var_sym == value) %>%
    nrow()
  
  young_n <- participant_tbl %>%
    filter(AgeClass == "younger adults", !!var_sym == value) %>%
    nrow()
  
  older_n <- participant_tbl %>%
    filter(AgeClass == "older adults", !!var_sym == value) %>%
    nrow()
  
  data.frame(
    Characteristic = label,
    Summary = "\\(N\\) (\\%)",
    Overall = fmt_n_pct(overall_n, overall_denom),
    Younger_adults = fmt_n_pct(young_n, young_denom),
    Older_adults = fmt_n_pct(older_n, older_denom),
    stringsAsFactors = FALSE
  )
}

make_episode_cat_row <- function(var_name, value, label) {
  var_sym <- rlang::sym(var_name)
  
  overall_denom <- nrow(df_model)
  young_denom <- df_model %>% filter(AgeClass == "younger adults") %>% nrow()
  older_denom <- df_model %>% filter(AgeClass == "older adults") %>% nrow()
  
  overall_n <- df_model %>%
    filter(!!var_sym == value) %>%
    nrow()
  
  young_n <- df_model %>%
    filter(AgeClass == "younger adults", !!var_sym == value) %>%
    nrow()
  
  older_n <- df_model %>%
    filter(AgeClass == "older adults", !!var_sym == value) %>%
    nrow()
  
  data.frame(
    Characteristic = label,
    Summary = "\\(n\\) (\\%)",
    Overall = fmt_n_pct(overall_n, overall_denom),
    Younger_adults = fmt_n_pct(young_n, young_denom),
    Older_adults = fmt_n_pct(older_n, older_denom),
    stringsAsFactors = FALSE
  )
}

make_episode_cont_row <- function(var_name, label) {
  var_sym <- rlang::sym(var_name)
  
  overall <- df_model %>%
    summarise(
      mean_value = mean(!!var_sym, na.rm = TRUE),
      sd_value = sd(!!var_sym, na.rm = TRUE)
    )
  
  young <- df_model %>%
    filter(AgeClass == "younger adults") %>%
    summarise(
      mean_value = mean(!!var_sym, na.rm = TRUE),
      sd_value = sd(!!var_sym, na.rm = TRUE)
    )
  
  older <- df_model %>%
    filter(AgeClass == "older adults") %>%
    summarise(
      mean_value = mean(!!var_sym, na.rm = TRUE),
      sd_value = sd(!!var_sym, na.rm = TRUE)
    )
  
  data.frame(
    Characteristic = label,
    Summary = "Mean (SD)",
    Overall = fmt_mean_sd(overall$mean_value, overall$sd_value),
    Younger_adults = fmt_mean_sd(young$mean_value, young$sd_value),
    Older_adults = fmt_mean_sd(older$mean_value, older$sd_value),
    stringsAsFactors = FALSE
  )
}

section_row <- function(label) {
  data.frame(
    Characteristic = paste0("\\textit{", label, "}"),
    Summary = "",
    Overall = "",
    Younger_adults = "",
    Older_adults = "",
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 8.2 Raw descriptive tables
# -----------------------------

sample_tbl <- data.frame(
  N_episodes = nrow(df_model),
  N_participants = n_distinct(df_model$PID)
)

sample_by_age_tbl <- df_model %>%
  group_by(AgeClass) %>%
  summarise(
    n_episodes = n(),
    n_participants = n_distinct(PID),
    .groups = "drop"
  )

participant_age_tbl <- participant_tbl %>%
  count(AgeClass, name = "n_participants") %>%
  mutate(percent = 100 * n_participants / sum(n_participants))

participant_sex_overall_tbl <- participant_tbl %>%
  count(blsex, name = "n_participants") %>%
  mutate(percent = 100 * n_participants / sum(n_participants))

participant_sex_by_age_tbl <- participant_tbl %>%
  count(AgeClass, blsex, name = "n_participants") %>%
  group_by(AgeClass) %>%
  mutate(percent_within_age = 100 * n_participants / sum(n_participants)) %>%
  ungroup()

participant_home_overall_tbl <- participant_tbl %>%
  count(bl_home_type, name = "n_participants") %>%
  mutate(percent = 100 * n_participants / sum(n_participants))

participant_home_by_age_tbl <- participant_tbl %>%
  count(AgeClass, bl_home_type, name = "n_participants") %>%
  group_by(AgeClass) %>%
  mutate(percent_within_age = 100 * n_participants / sum(n_participants)) %>%
  ungroup()

episode_cont_overall_tbl <- df_model %>%
  summarise(
    Rating_mean = mean(Rating, na.rm = TRUE),
    Rating_sd = sd(Rating, na.rm = TRUE),
    mean_qi_mean = mean(mean_qi, na.rm = TRUE),
    mean_qi_sd = sd(mean_qi, na.rm = TRUE),
    duration_min_mean = mean(duration_min, na.rm = TRUE),
    duration_min_sd = sd(duration_min, na.rm = TRUE)
  )

episode_cont_by_age_tbl <- df_model %>%
  group_by(AgeClass) %>%
  summarise(
    n_episodes = n(),
    n_participants = n_distinct(PID),
    Rating_mean = mean(Rating, na.rm = TRUE),
    Rating_sd = sd(Rating, na.rm = TRUE),
    mean_qi_mean = mean(mean_qi, na.rm = TRUE),
    mean_qi_sd = sd(mean_qi, na.rm = TRUE),
    duration_min_mean = mean(duration_min, na.rm = TRUE),
    duration_min_sd = sd(duration_min, na.rm = TRUE),
    .groups = "drop"
  )

social_overall_tbl <- df_model %>%
  count(social_context_3, name = "n_episodes") %>%
  mutate(percent = 100 * n_episodes / sum(n_episodes))

social_by_age_tbl <- df_model %>%
  count(AgeClass, social_context_3, name = "n_episodes") %>%
  group_by(AgeClass) %>%
  mutate(percent_within_age = 100 * n_episodes / sum(n_episodes)) %>%
  ungroup()

weekday_overall_tbl <- df_model %>%
  count(weekday_weekend, name = "n_episodes") %>%
  mutate(percent = 100 * n_episodes / sum(n_episodes))

weekday_by_age_tbl <- df_model %>%
  count(AgeClass, weekday_weekend, name = "n_episodes") %>%
  group_by(AgeClass) %>%
  mutate(percent_within_age = 100 * n_episodes / sum(n_episodes)) %>%
  ungroup()

time_of_day_overall_tbl <- df_model %>%
  count(time_of_day, name = "n_episodes") %>%
  mutate(percent = 100 * n_episodes / sum(n_episodes))

time_of_day_by_age_tbl <- df_model %>%
  count(AgeClass, time_of_day, name = "n_episodes") %>%
  group_by(AgeClass) %>%
  mutate(percent_within_age = 100 * n_episodes / sum(n_episodes)) %>%
  ungroup()

participant_cont_tbl <- participant_tbl %>%
  summarise(
    gse_person_mean_raw_n = sum(!is.na(gse_person_mean_raw)),
    gse_person_mean_raw_mean = mean(gse_person_mean_raw, na.rm = TRUE),
    gse_person_mean_raw_sd = sd(gse_person_mean_raw, na.rm = TRUE),
    gseB_raw_n = sum(!is.na(gseB_raw)),
    gseB_raw_mean = mean(gseB_raw, na.rm = TRUE),
    gseB_raw_sd = sd(gseB_raw, na.rm = TRUE)
  )

# -----------------------------
# 8.3 Within-person centering check
# -----------------------------

within_check <- df_model %>%
  group_by(PID) %>%
  summarise(
    mean_gseW = mean(gseW_raw, na.rm = TRUE),
    .groups = "drop"
  )

within_tbl <- data.frame(
  Mean_of_person_specific_centered_gse = mean(within_check$mean_gseW, na.rm = TRUE),
  SD_of_person_specific_centered_gse = sd(within_check$mean_gseW, na.rm = TRUE),
  Min = min(within_check$mean_gseW, na.rm = TRUE),
  Max = max(within_check$mean_gseW, na.rm = TRUE)
)

# -----------------------------
# 8.4 Social-context coding check
# -----------------------------

social_coding_check_tbl <- data.frame(
  coding_rule = c(
    "alone: pp_alone == 1 and all other person indicators == 0",
    "partner, family and friends: partner OR friends OR family present",
    "other social contacts: remaining valid contexts"
  ),
  n_episodes = c(
    sum(df_model$social_context_3 == "alone", na.rm = TRUE),
    sum(df_model$social_context_3 == "partner, family and friends", na.rm = TRUE),
    sum(df_model$social_context_3 == "other social contacts", na.rm = TRUE)
  )
) %>%
  mutate(percent = 100 * n_episodes / sum(n_episodes))

# -----------------------------
# 8.5 Paper-ready descriptive table
# -----------------------------

n_participants_overall <- n_distinct(df_model$PID)
n_participants_young <- participant_tbl %>% filter(AgeClass == "younger adults") %>% nrow()
n_participants_older <- participant_tbl %>% filter(AgeClass == "older adults") %>% nrow()

n_episodes_overall <- nrow(df_model)
n_episodes_young <- df_model %>% filter(AgeClass == "younger adults") %>% nrow()
n_episodes_older <- df_model %>% filter(AgeClass == "older adults") %>% nrow()

descriptive_table_paper <- bind_rows(
  section_row("Sample size"),
  data.frame(
    Characteristic = "Participants",
    Summary = "\\(N\\)",
    Overall = fmt_n(n_participants_overall),
    Younger_adults = fmt_n(n_participants_young),
    Older_adults = fmt_n(n_participants_older),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Characteristic = "Walking episodes",
    Summary = "\\(n\\)",
    Overall = fmt_n(n_episodes_overall),
    Younger_adults = fmt_n(n_episodes_young),
    Older_adults = fmt_n(n_episodes_older),
    stringsAsFactors = FALSE
  ),
  
  section_row("Participant-level characteristics"),
  make_participant_cat_row("AgeClass", "younger adults", "Age group: younger adults"),
  make_participant_cat_row("AgeClass", "older adults", "Age group: older adults"),
  make_participant_cat_row("blsex", "male", "Gender: male"),
  make_participant_cat_row("blsex", "female", "Gender: female"),
  make_participant_cat_row("blsex", "diverse", "Gender: diverse"),
  make_participant_cat_row("bl_home_type", "urban", "Home type: urban"),
  make_participant_cat_row("bl_home_type", "suburban", "Home type: suburban"),
  make_participant_cat_row("bl_home_type", "rural", "Home type: rural"),
  
  section_row("Episode-level continuous variables"),
  make_episode_cont_row("Rating", "Perceived restoration, ROS"),
  make_episode_cont_row("mean_qi", "Greenspace exposure, NDVI"),
  make_episode_cont_row("duration_min", "Episode duration, minutes"),
  
  section_row("Episode-level categorical variables"),
  make_episode_cat_row("social_context_3", "alone", "Presence of others: alone"),
  make_episode_cat_row("social_context_3", "partner, family and friends", "Presence of others: partner, family and friends"),
  make_episode_cat_row("social_context_3", "other social contacts", "Presence of others: other social contacts"),
  make_episode_cat_row("weekday_weekend", "weekday", "Weekday"),
  make_episode_cat_row("weekday_weekend", "weekend", "Weekend"),
  make_episode_cat_row("time_of_day", "morning", "Morning"),
  make_episode_cat_row("time_of_day", "afternoon", "Afternoon"),
  make_episode_cat_row("time_of_day", "evening_night", "Evening/night")
)

# -----------------------------
# 8.6 Print key outputs
# -----------------------------

print(sample_tbl)
print(sample_by_age_tbl)

print(participant_age_tbl)
print(participant_sex_overall_tbl)
print(participant_sex_by_age_tbl)
print(participant_home_overall_tbl)
print(participant_home_by_age_tbl)

print(episode_cont_overall_tbl)
print(episode_cont_by_age_tbl)

print(social_overall_tbl)
print(social_by_age_tbl)
print(weekday_overall_tbl)
print(weekday_by_age_tbl)
print(time_of_day_overall_tbl)
print(time_of_day_by_age_tbl)

print(within_tbl)
print(social_coding_check_tbl)

print(descriptive_table_paper)

# -----------------------------
# 8.7 Export tables
# -----------------------------

descriptive_dir <- "descriptive_outputs"
dir.create(descriptive_dir, showWarnings = FALSE)

write.csv(sample_tbl, file.path(descriptive_dir, "sample_size_overall.csv"), row.names = FALSE)
write.csv(sample_by_age_tbl, file.path(descriptive_dir, "sample_size_by_age.csv"), row.names = FALSE)

write.csv(participant_age_tbl, file.path(descriptive_dir, "participant_age_group_overall.csv"), row.names = FALSE)
write.csv(participant_sex_overall_tbl, file.path(descriptive_dir, "participant_gender_overall.csv"), row.names = FALSE)
write.csv(participant_sex_by_age_tbl, file.path(descriptive_dir, "participant_gender_by_age.csv"), row.names = FALSE)
write.csv(participant_home_overall_tbl, file.path(descriptive_dir, "participant_home_type_overall.csv"), row.names = FALSE)
write.csv(participant_home_by_age_tbl, file.path(descriptive_dir, "participant_home_type_by_age.csv"), row.names = FALSE)

write.csv(episode_cont_overall_tbl, file.path(descriptive_dir, "episode_continuous_overall.csv"), row.names = FALSE)
write.csv(episode_cont_by_age_tbl, file.path(descriptive_dir, "episode_continuous_by_age.csv"), row.names = FALSE)

write.csv(social_overall_tbl, file.path(descriptive_dir, "episode_social_context_overall.csv"), row.names = FALSE)
write.csv(social_by_age_tbl, file.path(descriptive_dir, "episode_social_context_by_age.csv"), row.names = FALSE)
write.csv(weekday_overall_tbl, file.path(descriptive_dir, "episode_weekday_weekend_overall.csv"), row.names = FALSE)
write.csv(weekday_by_age_tbl, file.path(descriptive_dir, "episode_weekday_weekend_by_age.csv"), row.names = FALSE)
write.csv(time_of_day_overall_tbl, file.path(descriptive_dir, "episode_time_of_day_overall.csv"), row.names = FALSE)
write.csv(time_of_day_by_age_tbl, file.path(descriptive_dir, "episode_time_of_day_by_age.csv"), row.names = FALSE)

write.csv(participant_cont_tbl, file.path(descriptive_dir, "participant_continuous_variables.csv"), row.names = FALSE)
write.csv(within_tbl, file.path(descriptive_dir, "within_person_centering_check.csv"), row.names = FALSE)
write.csv(social_coding_check_tbl, file.path(descriptive_dir, "social_context_coding_check.csv"), row.names = FALSE)

write.csv(descriptive_table_paper, file.path(descriptive_dir, "descriptive_table_paper_ready.csv"), row.names = FALSE)
# =========================================================
# 9) Extract model results
# =========================================================
# Null model
mnull_fit  <- fit_stats(mnull)
mnull_r2   <- safe_r2(mnull)
mnull_icc  <- safe_icc(mnull)
mnull_diag <- diag_stats(mnull)

# Model 0: baseline adjusted random-intercept model
m0_fit   <- fit_stats(m0_random_intercept)
m0_comp  <- extract_lrt(mnull, m0_random_intercept, "Null model", "Model 0")
m0_type3 <- extract_lmm_type3(m0_random_intercept)
m0_r2    <- safe_r2(m0_random_intercept)
m0_icc   <- safe_icc(m0_random_intercept)
m0_diag  <- diag_stats(m0_random_intercept)

# Model 1: baseline adjusted random-slope model
m1_fit   <- fit_stats(m1)
m1_comp  <- extract_lrt(m0_random_intercept, m1, "random-intercept baseline model", "Model 1")
m1_type3 <- extract_lmm_type3(m1)
m1_r2    <- safe_r2(m1)
m1_icc   <- safe_icc(m1)
m1_diag  <- diag_stats(m1)

# Model 2: age moderation
m2_fit         <- fit_stats(m2)
m2_comp        <- extract_lrt(m1, m2, "Model 1", "Model 2")
m2_type3       <- extract_lmm_type3(m2)
m2_r2          <- safe_r2(m2)
m2_icc         <- safe_icc(m2)
m2_diag        <- diag_stats(m2)
m2_slopes      <- extract_trends(m2, "AgeClass", "gseW_raw")
m2_slope_pairs <- extract_trend_pairs(m2, "AgeClass", "gseW_raw")

# Model 3: presence-of-others moderation
m3_fit                        <- fit_stats(m3)
m3_comp                       <- extract_lrt(m2, m3, "Model 2", "Model 3")
m3_type3                      <- extract_lmm_type3(m3)
m3_r2                         <- safe_r2(m3)
m3_icc                        <- safe_icc(m3)
m3_diag                       <- diag_stats(m3)
m3_social_context_slopes      <- extract_trends(m3, "social_context_3", "gseW_raw")
m3_social_context_slope_pairs <- extract_trend_pairs(m3, "social_context_3", "gseW_raw")

# =========================================================
# 10) Figures
# =========================================================
gseW_plot_term <- make_ggeffects_term("gseW_raw", df_model$gseW_raw)

# Model 2: age moderation
pred_m2 <- ggpredict(m2, terms = c(gseW_plot_term, "AgeClass"))

p_m2 <- ggplot(pred_m2, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  theme_minimal(base_size = 12) +
  labs(
    # title = "Predicted restorative rating across age groups and greenspace exposure",
    x = "Within-person greenspace exposure",
    y = "Predicted restorative rating",
    color = "Age group",
    fill = "Age group"
  )

ggsave(plot_m4_file, p_m2, width = 7, height = 5, dpi = 300)

# Descriptive Presence of others by age group
p_social_context_by_age <- ggplot(
  social_by_age_tbl,
  aes(x = social_context_3, y = percent_within_age, fill = AgeClass)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Three-category presence of others by age group",
    x = "Presence of others",
    y = "Percentage within age group (%)",
    fill = "Age group"
  )

ggsave(plot_social_context_by_age_file, p_social_context_by_age, width = 7, height = 5, dpi = 300)

# Model 3: social-context moderation of the within-person greenspace effect
pred_m3_social_context <- ggpredict(
  m3,
  terms = c(gseW_plot_term, "social_context_3")
)

p_m3_social_context <- ggplot(
  pred_m3_social_context,
  aes(x = x, y = predicted, color = group, fill = group)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  theme_minimal(base_size = 12) +
  labs(
    # title = "Predicted restorative ratings across presence of others and greenspace exposure",
    x = "Within-person greenspace exposure",
    y = "Predicted restorative rating",
    color = "Presence of others",
    fill = "Presence of others"
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(
      hjust = 0.5,
      size = 12
    ),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(plot_m5b_social_context_file, p_m3_social_context, width = 7, height = 5, dpi = 300)

# =========================================================
# 11) Narrative text
# =========================================================
intro_text <- paste0(
  "This report summarises a refined mixed-effects modelling strategy for episode-level restorative ratings. ",
  "The model sequence first evaluates within-person greenspace exposure and age moderation. ",
  "It then introduces a three-category social-context variable distinguishing episodes spent alone, episodes with partner, family and friendss, and episodes with other people. ",
  "partner, family and friendss are defined as episodes in which a partner, friend, or family member was present. ",
  "Alone episodes are defined strictly as pp_alone = 1 with all other person indicators equal to 0. ",
  "The remaining valid presence of others are coded as other social contacts. ",
  "All adjusted models control for log-transformed episode duration. ",
  "Greenspace exposure is decomposed into within-person and between-person components but is not standardized. ",
  "No three-way interaction between greenspace exposure, age group, and presence of others is estimated. ",
  "Categorical predictors are interpreted using omnibus Type III tests, estimated marginal means, and pairwise contrasts."
)

m1_text <- paste0(
  model_comp_sentence(m1_comp, "the corresponding random-intercept baseline model", "Model 1"),
  " This comparison evaluates whether allowing the within-person greenspace slope to vary across participants improves model fit."
)
m2_text <- model_comp_sentence(m2_comp, "Model 1", "Model 2")
m3_text <- model_comp_sentence(m3_comp, "Model 2", "Model 3")

# =========================================================
# 12) Build HTML report
# =========================================================
css <- "
body { font-family: Georgia, 'Times New Roman', serif; background: #ffffff; color: #222; line-height: 1.6; margin: 0; padding: 0; }
.container { max-width: 1080px; margin: 0 auto; padding: 40px 48px 60px 48px; }
h1 { font-size: 30px; margin-bottom: 10px; }
h2 { font-size: 24px; margin-top: 36px; margin-bottom: 8px; border-bottom: 1px solid #d9d9d9; padding-bottom: 6px; }
h3 { font-size: 18px; margin-top: 24px; margin-bottom: 8px; }
p { font-size: 15px; margin-top: 8px; margin-bottom: 12px; }
.formula-block { background: #f7f7f7; border: 1px solid #e5e5e5; padding: 12px 14px; font-family: Menlo, Consolas, monospace; font-size: 13px; white-space: pre-wrap; }
.paper-table { width: 100%; border-collapse: collapse; font-size: 14px; margin-top: 6px; margin-bottom: 18px; }
.paper-table th, .paper-table td { border: 1px solid #d8d8d8; padding: 6px 8px; text-align: left; }
.paper-table th { background: #f2f2f2; }
.table-wrap { overflow-x: auto; }
.figure-wrap { margin-top: 18px; margin-bottom: 26px; }
.figure-img { max-width: 100%; border: 1px solid #e0e0e0; }
.figure-caption { font-size: 13px; color: #555; margin-top: 6px; }
"

html_doc <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$title("Paper-style mixed model report"),
    tags$style(HTML(css))
  ),
  tags$body(
    tags$div(
      class = "container",
      
      tags$h1("Paper-style mixed model report"),
      tags$p(intro_text),
      
      tags$h2("Data overview"),
      tags$p("The analytic sample retained episodes with non-missing values on the outcome, greenspace exposure, duration, age group, biological sex, home type, time-of-day, weekday/weekend status, and the three-category social-context variable. Duration was log-transformed before being included as a control variable. Greenspace exposure was decomposed into within-person and between-person components without standardization."),
      tags$h3("Preprocessing and filtering steps"),
      tags$div(class = "table-wrap", table_html(preprocessing_log)),
      tags$h3("Sample size"),
      tags$div(class = "table-wrap", table_html(sample_tbl)),
      tags$h3("Within-person centering check"),
      tags$p("The person-mean centred within-person greenspace variable should have participant-specific means close to zero."),
      tags$div(class = "table-wrap", table_html(within_tbl)),
      tags$h3("Three-category social-context coding"),
      tags$div(class = "table-wrap", table_html(social_coding_check_tbl)),
      tags$h3("Three-category social-context distribution"),
      tags$div(class = "table-wrap", table_html(social_overall_tbl)),
      tags$h3("Three-category social-context distribution by age group"),
      tags$div(class = "table-wrap", table_html(social_by_age_tbl)),
      tags$div(
        class = "figure-wrap",
        tags$img(src = plot_social_context_by_age_file, class = "figure-img"),
        tags$div(class = "figure-caption", "Figure 1. Distribution of the three-category social-context variable by age group.")
      ),
      
      section_block(
        title = "Null model",
        formula_text = "Rating ~ 1 + (1 | PID)",
        tables = list(
          list(title = "Fit statistics", data = mnull_fit),
          list(title = "R-squared", data = mnull_r2),
          list(title = "Intraclass correlation", data = mnull_icc),
          list(title = "Model diagnostics", data = mnull_diag)
        )
      ),
      
      section_block(
        title = "Model 0. Baseline adjusted random-intercept model",
        formula_text = "Rating ~ gseW_raw + gseB_raw + duration_s + AgeClass + blsex + bl_home_type + time_of_day + weekday_weekend + (1 | PID)",
        tables = list(
          list(title = "Fit statistics", data = m0_fit),
          list(title = "Model comparison with Null model", data = m0_comp),
          list(title = "Type III tests", data = m0_type3),
          list(title = "R-squared", data = m0_r2),
          list(title = "Intraclass correlation", data = m0_icc),
          list(title = "Model diagnostics", data = m0_diag)
        )
      ),
      
      section_block(
        title = "Model 1. Baseline adjusted model",
        formula_text = "Rating ~ gseW_raw + gseB_raw + duration_s + AgeClass + blsex + bl_home_type + time_of_day + weekday_weekend + (1 + gseW_raw | PID)",
        paragraphs = list(
          "This baseline adjusted model includes within-person and between-person greenspace exposure, episode duration, age group, gender, living environment, time of day, and weekday/weekend status, with random intercepts and random slopes for within-person greenspace exposure.",
          m1_text
        ),
        tables = list(
          list(title = "Fit statistics", data = m1_fit),
          list(title = "Model comparison with random-intercept baseline model", data = m1_comp),
          list(title = "Type III tests", data = m1_type3),
          list(title = "R-squared", data = m1_r2),
          list(title = "Intraclass correlation", data = m1_icc),
          list(title = "Model diagnostics", data = m1_diag)
        )
      ),
      
      section_block(
        title = "Model 2. Age moderation model",
        formula_text = "Rating ~ gseW_raw * AgeClass + gseB_raw + duration_s + blsex + bl_home_type + time_of_day + weekday_weekend + (1 + gseW_raw | PID)",
        paragraphs = list(
          "This model tests whether the within-person greenspace slope differs by age group.",
          m2_text,
          "The key follow-up quantities are the group-specific within-person greenspace slopes estimated separately for each age group."
        ),
        tables = list(
          list(title = "Fit statistics", data = m2_fit),
          list(title = "Model comparison with Model 1", data = m2_comp),
          list(title = "Type III tests", data = m2_type3),
          list(title = "R-squared", data = m2_r2),
          list(title = "Intraclass correlation", data = m2_icc),
          list(title = "Model diagnostics", data = m2_diag),
          list(title = "Within-person greenspace slopes by age group", data = m2_slopes),
          list(title = "Slope contrast between age groups", data = m2_slope_pairs)
        ),
        image_file = plot_m4_file,
        image_caption = "Figure 2. Predicted restorative ratings across within-person greenspace exposure in the age moderation model."
      ),
      
      section_block(
        title = "Model 3. Presence-of-others moderation model",
        formula_text = "Rating ~ gseW_raw * AgeClass + gseB_raw + duration_s + gseW_raw * social_context_3 + blsex + bl_home_type + time_of_day + weekday_weekend + (1 + gseW_raw | PID)",
        paragraphs = list(
          "This model tests whether the within-person greenspace slope differs among alone, close-tie, and other-people episodes.",
          "No three-way interaction with age group is included.",
          m3_text
        ),
        tables = list(
          list(title = "Fit statistics", data = m3_fit),
          list(title = "Model comparison with Model 2", data = m3_comp),
          list(title = "Type III tests", data = m3_type3),
          list(title = "R-squared", data = m3_r2),
          list(title = "Intraclass correlation", data = m3_icc),
          list(title = "Model diagnostics", data = m3_diag),
          list(title = "Within-person greenspace slopes by presence of others", data = m3_social_context_slopes),
          list(title = "Slope contrasts by presence of others", data = m3_social_context_slope_pairs)
        ),
        image_file = plot_m5b_social_context_file,
        image_caption = "Figure 4. Predicted restorative ratings across within-person greenspace exposure by three-category presence of others."
      )
    )
  )
)

# =========================================================
# 13) Save HTML
# =========================================================
save_html(html_doc, file = output_html)

message("Analysis complete.")
message("HTML report saved to: ", output_html)
message("Figure saved to: ", plot_m4_file)
message("Figure saved to: ", plot_m5b_social_context_file)
message("Figure saved to: ", plot_social_context_by_age_file)