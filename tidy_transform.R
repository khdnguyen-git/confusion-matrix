# Load
if (!require("pacman")) {
  install.packages("pacman")
  library(pacman)
}

pacman::p_load(tidyverse,
               rio,
               here,
               janitor,
               psych)


raw_data <- import(here("input", "PublicationFileFinal 35.csv")) %>% 
  clean_names() %>% 
  rename("sens" = "sensitivity",
         "spec" = "specificity",
         "n_tot" = "total_n",
         "n_positive_miss" = "n_for", # Screened positive but no follow-up 
         "asd_prev_natl" = "asd_prev_n",
         "asd_prev_cdc" = "asd_prev_c",
         "asd_prev_rep" = "report_prev",
         "asd_prev_world" = "glob_prev")
  

# 1.0 - Clean
redo_tidyverse <-  raw_data %>% 
  mutate(
    sens_calc_1 = tp / (tp + fn),
    spec_calc_1 = tn / (tn + fp),
    ppv_calc_1 = tp / (tp + fp),
    npv_calc_1 = tn / (tn + fn)) %>% 
  
  mutate(
    tp_calc_sens_fn_1 = (sens * fn) / (1 - sens),
    tp_calc_ppv_fp_1 = (ppv * fp) / (1 - ppv),
    fn_calc_sens_tp = (tp * (1 - sens))/sens,
    fn_calc_npv_tn = (tn * (1 - npv)/npv),
    
    fp_calc_spec_tn = (tn * (1 - spec)) / spec,
    fp_calc_ppv_tp = (tp * (1 - ppv)) / ppv,
    tn_calc_spec_fp = (spec * fp) / (1 - spec),
    tn_calc_npv_fn = (npv * fn) / (1 - npv)) %>% 
  
  mutate(
    tn_calc_sub_1 = n - (fn + fp + tp),
    tn_calc_sub_2 = n - (fn_calc_sens_tp + fp_calc_ppv_tp + tp),
    fn_calc_sub_1 = n - (tn + fp + tp), # Only one entry, but add _1 for naming consistency
    fp_calc_sub_1 = n - (tn + fn + tp), 
    tp_calc_sens_fn_2 = (sens * fn_calc_sens_tp) / (1 - sens),
    tp_calc_ppv_fp_2 = (ppv * fp_calc_ppv_tp) / (1 - ppv)) %>% 
  
  mutate(n_analyze = n, .after = n,
         tp_rep = tp) %>%                    # Reused output 
  
  mutate(tn_rep = coalesce(                  # Reused output 
    ifelse(!is.na(tn), tn, tn_calc_sub_1),
    tn_calc_spec_fp,
    tn_calc_npv_fn,
    tn_calc_sub_2,
    0
  )) %>% 
  
  mutate(fp_calc_sub_2 = n - (tn_rep + fn + tp)) %>% 
  
  mutate(fp_rep = coalesce(                  # Reused output 
    ifelse(!is.na(fp), fp, fp_calc_sub_1),
    fp_calc_spec_tn,
    fp_calc_ppv_tp,
    fp_calc_sub_2,
    0
  )) %>% 
  
  mutate(fn_rep = coalesce(                  # Reused output 
    ifelse(!is.na(fn), fn, fn_calc_sub_1),
    fn_calc_sens_tp,
    fn_calc_npv_tn,
    ((tn_rep + tp_rep + fp_rep) / n),
    0
  )) %>%
  
  mutate(sens_calc_2 = tp / (tp + fn_rep),
         spec_calc_2 = tn_rep / (tn_rep + fp_rep),
         ppv_calc_2 = tp / (tp + fp_rep),
         npv_calc_2 = tn / (tn_rep + fn_rep),
         asd_rep = tp + fn_rep) %>% 
  
# 1.2 - Reclaim
  mutate(tp_reclaim_1p = tp_rep + 1 * n_positive_miss,                      # Reused output 
         tp_reclaim_0p = tp_rep + 0 * n_positive_miss,                      # Reused output
         tp_reclaim_ppv = tp_rep + ppv_calc_2 * n_positive_miss,            # Reused output
         fp_reclaim_1p = fp_rep + (1 - 1) * n_positive_miss,                # Reused output
         fp_reclaim_0p = fp_rep + (1 - 0) * n_positive_miss,                # Reused output
         fp_reclaim_ppv = fp_rep + (1 - ppv_calc_2) * n_positive_miss) %>%  # Reused output
  
# 2.0 - Missingness calc
  
  # National prevalence
  mutate(asd_exp_natl = as.integer(n * asd_prev_natl)) %>% 
  mutate(asd_diff_natl = asd_exp_natl - (tp + fn_rep),
         asd_diff_natl_cat = ifelse(asd_diff_natl > 0, 1, 
                                    ifelse(asd_diff_natl < 0, 0, asd_diff_natl))) %>% 
  
  mutate(fn_calc_natl = fn_rep + ifelse(asd_diff_natl > 0, asd_diff_natl, 0),
         tn_calc_natl = tn_rep - (asd_diff_natl * asd_diff_natl_cat)) %>% 
  
  # CDC prevalence
  mutate(asd_exp_cdc = as.integer(n * asd_prev_cdc)) %>%
  mutate(asd_diff_cdc = asd_exp_cdc - (tp + fn_rep),
         asd_diff_cdc_cat = ifelse(asd_diff_cdc > 0, 1, 
                                   ifelse(asd_diff_cdc < 0, 0, asd_diff_cdc))) %>% 
  
  mutate(fn_calc_cdc = fn_rep + ifelse(asd_diff_cdc > 0, asd_diff_cdc, 0),
         tn_calc_cdc = tn_rep - (asd_diff_cdc * asd_diff_cdc_cat)) %>% 
  
  # World prevalence
  mutate(asd_exp_world = as.integer(n * asd_prev_world)) %>%
  mutate(asd_diff_world = asd_exp_world - (tp + fn_rep),
         asd_diff_world_cat = ifelse(asd_diff_world > 0, 1, 
                                   ifelse(asd_diff_world < 0, 0, asd_diff_world))) %>% 
  
  mutate(fn_calc_world = fn_rep + ifelse(asd_diff_world > 0, asd_diff_world, 0),
         tn_calc_world = tn_rep - (asd_diff_world * asd_diff_world_cat)) %>% 
  
  # Diagnostic accuracy - National
  mutate(sens_calc_natl = tp_rep / (tp_rep + fn_calc_natl),
         spec_calc_natl = tn_calc_natl / (tn_calc_natl + fp_rep),
         ppv_calc_natl = tp_rep / (tp_rep + fp_rep),
         npv_calc_natl = tn_calc_natl / (tn_calc_natl + fn_calc_natl)) %>% 
  
  # Diagnostic accuracy - CDC
  mutate(sens_calc_cdc = tp_rep / (tp_rep + fn_calc_cdc),
         spec_calc_cdc = tn_calc_cdc / (tn_calc_cdc + fp_rep),
         ppv_calc_cdc = tp_rep / (tp_rep + fp_rep),
         npv_calc_cdc = tn_calc_cdc / (tn_calc_cdc + fn_calc_cdc)) %>% 
  
  # Diagnostic accuracy - World
  mutate(sens_calc_world = tp_rep / (tp_rep + fn_calc_world),
         spec_calc_world = tn_calc_world / (tn_calc_world + fp_rep),
         ppv_calc_world = tp_rep / (tp_rep + fp_rep),
         npv_calc_world = tn_calc_world / (tn_calc_world + fn_calc_world)) %>% 
  
# 2.0.3 - Missingness calc; with 100% missing data
  
  # National prevalence
  mutate(asd_exp_natl_100 = as.integer((n + n_positive_miss) * asd_prev_natl)) %>% 
  mutate(asd_diff_natl_100 = asd_exp_natl_100 - (tp_reclaim_1p + fn_rep),
         asd_diff_natl_cat_100 = ifelse(asd_diff_natl_100 > 0, 1, 
                                    ifelse(asd_diff_natl_100 < 0, 0, asd_diff_natl_100))) %>% 
  
  mutate(fn_calc_natl_100 = fn_rep + ifelse(asd_diff_natl_100 > 0, asd_diff_natl_100, 0),
         tn_calc_natl_100 = tn_rep - (asd_diff_natl_100 * asd_diff_natl_cat_100)) %>% 
  
  # CDC prevalence
  mutate(asd_exp_cdc_100 = as.integer((n + n_positive_miss) * asd_prev_cdc)) %>%
  mutate(asd_diff_cdc_100 = asd_exp_cdc_100 - (tp_reclaim_1p + fn_rep),
         asd_diff_cdc_cat_100 = ifelse(asd_diff_cdc_100 > 0, 1, 
                                        ifelse(asd_diff_cdc_100 < 0, 0, asd_diff_cdc_100))) %>% 
  
  mutate(fn_calc_cdc_100 = fn_rep + ifelse(asd_diff_cdc_100 > 0, asd_diff_cdc_100, 0),
         tn_calc_cdc_100 = tn_rep - (asd_diff_cdc_100 * asd_diff_cdc_cat_100)) %>% 
  
  # World prevalence
  mutate(asd_exp_world_100 = as.integer((n + n_positive_miss) * asd_prev_world)) %>%
  mutate(asd_diff_world_100 = asd_exp_world_100 - (tp_reclaim_1p + fn_rep),
         asd_diff_world_cat_100 = ifelse(asd_diff_world_100 > 0, 1, 
                                       ifelse(asd_diff_world_100 < 0, 0, asd_diff_world_100))) %>% 
  
  mutate(fn_calc_world_100 = fn_rep + ifelse(asd_diff_world_100 > 0, asd_diff_world_100, 0),
         tn_calc_world_100 = tn_rep - (asd_diff_world_100 * asd_diff_world_cat_100)) %>%  
  
  mutate(sens_calc_og_100 = tp_reclaim_1p / (tp_reclaim_1p + fn_rep),
         spec_calc_og_100 = tn_rep / (tn_rep + fp_reclaim_1p),
         ppv_calc_og_100 = tp_reclaim_1p / (tp_reclaim_1p + fp_reclaim_1p),
         npv_calc_og_100 = tn_rep / (tn_rep + fn_rep)) %>% 
  
  # Diagnostic accuracy - National
  mutate(sens_calc_natl_100 = tp_reclaim_1p / (tp_reclaim_1p + fn_calc_natl_100),
         spec_calc_natl_100 = tn_calc_natl_100 / (tn_calc_natl_100 + fp_reclaim_1p),
         ppv_calc_natl_100 = tp_reclaim_1p / (tp_reclaim_1p + fp_reclaim_1p),
         npv_calc_natl_100 = tn_calc_natl_100 / (tn_calc_natl_100 + fn_calc_natl_100)) %>% 
  
  # Diagnostic accuracy - CDC
  mutate(sens_calc_cdc_100 = tp_reclaim_1p / (tp_reclaim_1p + fn_calc_cdc_100),
         spec_calc_cdc_100 = tn_calc_cdc_100 / (tn_calc_cdc_100 + fp_reclaim_1p),
         ppv_calc_cdc_100 = tp_reclaim_1p / (tp_reclaim_1p + fp_reclaim_1p),
         npv_calc_cdc_100 = tn_calc_cdc_100 / (tn_calc_cdc_100 + fn_calc_cdc_100)) %>% 
  
  # Diagnostic accuracy - World
  mutate(sens_calc_world_100 = tp_reclaim_1p / (tp_reclaim_1p + fn_calc_world_100),
         spec_calc_world_100 = tn_calc_world_100 / (tn_calc_world_100 + fp_reclaim_1p),
         ppv_calc_world_100 = tp_reclaim_1p / (tp_reclaim_1p + fp_reclaim_1p),
         npv_calc_world_100 = tn_calc_world_100 / (tn_calc_world_100 + fn_calc_world_100)) %>% 
  
# 2.0.3 - Missingness calc; with 0% missing data
  
  # National prevalence
  mutate(asd_exp_natl_000 = as.integer((n + n_positive_miss) * asd_prev_natl)) %>% 
  mutate(asd_diff_natl_000 = asd_exp_natl_000 - (tp_reclaim_0p + fn_rep),
         asd_diff_natl_cat_000 = ifelse(asd_diff_natl_000 > 0, 1, 
                                        ifelse(asd_diff_natl_000 < 0, 0, asd_diff_natl_000))) %>% 
  
  mutate(fn_calc_natl_000 = fn_rep + ifelse(asd_diff_natl_000 > 0, asd_diff_natl_000, 0),
         tn_calc_natl_000 = tn_rep - (asd_diff_natl_000 * asd_diff_natl_cat_000)) %>% 
  
  # CDC prevalence
  mutate(asd_exp_cdc_000 = as.integer((n + n_positive_miss) * asd_prev_cdc)) %>%
  mutate(asd_diff_cdc_000 = asd_exp_cdc_000 - (tp_reclaim_0p + fn_rep),
         asd_diff_cdc_cat_000 = ifelse(asd_diff_cdc_000 > 0, 1, 
                                       ifelse(asd_diff_cdc_000 < 0, 0, asd_diff_cdc_000))) %>% 
  
  mutate(fn_calc_cdc_000 = fn_rep + ifelse(asd_diff_cdc_000 > 0, asd_diff_cdc_000, 0),
         tn_calc_cdc_000 = tn_rep - (asd_diff_cdc_000 * asd_diff_cdc_cat_000)) %>% 
  
  # World prevalence
  mutate(asd_exp_world_000 = as.integer((n + n_positive_miss) * asd_prev_world)) %>%
  mutate(asd_diff_world_000 = asd_exp_world_000 - (tp_reclaim_0p + fn_rep),
         asd_diff_world_cat_000 = ifelse(asd_diff_world_000 > 0, 1, 
                                         ifelse(asd_diff_world_000 < 0, 0, asd_diff_world_000))) %>% 
  
  mutate(fn_calc_world_000 = fn_rep + ifelse(asd_diff_world_000 > 0, asd_diff_world_000, 0),
         tn_calc_world_000 = tn_rep - (asd_diff_world_000 * asd_diff_world_cat_000)) %>%  
  
  mutate(sens_calc_og_000 = tp_reclaim_0p / (tp_reclaim_0p + fn_rep),
         spec_calc_og_000 = tn_rep / (tn_rep + fp_reclaim_0p),
         ppv_calc_og_000 = tp_reclaim_0p / (tp_reclaim_0p + fp_reclaim_0p),
         npv_calc_og_000 = tn_rep / (tn_rep + fn_rep)) %>% 
  
  # Diagnostic accuracy - National
  mutate(sens_calc_natl_000 = tp_reclaim_0p / (tp_reclaim_0p + fn_calc_natl_000),
         spec_calc_natl_000 = tn_calc_natl_000 / (tn_calc_natl_000 + fp_reclaim_0p),
         ppv_calc_natl_000 = tp_reclaim_0p / (tp_reclaim_0p + fp_reclaim_0p),
         npv_calc_natl_000 = tn_calc_natl_000 / (tn_calc_natl_000 + fn_calc_natl_000)) %>% 
  
  # Diagnostic accuracy - CDC
  mutate(sens_calc_cdc_000 = tp_reclaim_0p / (tp_reclaim_0p + fn_calc_cdc_000),
         spec_calc_cdc_000 = tn_calc_cdc_000 / (tn_calc_cdc_000 + fp_reclaim_0p),
         ppv_calc_cdc_000 = tp_reclaim_0p / (tp_reclaim_0p + fp_reclaim_0p),
         npv_calc_cdc_000 = tn_calc_cdc_000 / (tn_calc_cdc_000 + fn_calc_cdc_000)) %>% 
  
  # Diagnostic accuracy - World
  mutate(sens_calc_world_000 = tp_reclaim_0p / (tp_reclaim_0p + fn_calc_world_000),
         spec_calc_world_000 = tn_calc_world_000 / (tn_calc_world_000 + fp_reclaim_0p),
         ppv_calc_world_000 = tp_reclaim_0p / (tp_reclaim_0p + fp_reclaim_0p),
         npv_calc_world_000 = tn_calc_world_000 / (tn_calc_world_000 + fn_calc_world_000)) %>% 
                                    
# 2.0.3 - Missingness calc using PPV data
  
  # National prevalence
  mutate(asd_exp_natl_ppv = as.integer((n + n_positive_miss) * asd_prev_natl)) %>% 
  mutate(asd_diff_natl_ppv = asd_exp_natl_ppv - (tp_reclaim_ppv + fn_rep),
         asd_diff_natl_cat_ppv = ifelse(asd_diff_natl_ppv > 0, 1, 
                                        ifelse(asd_diff_natl_ppv < 0, 0, asd_diff_natl_ppv))) %>% 
  
  mutate(fn_calc_natl_ppv = fn_rep + ifelse(asd_diff_natl_ppv > 0, asd_diff_natl_ppv, 0),
         tn_calc_natl_ppv = tn_rep - (asd_diff_natl_ppv * asd_diff_natl_cat_ppv)) %>% 
  
  # CDC prevalence
  mutate(asd_exp_cdc_ppv = as.integer((n + n_positive_miss) * asd_prev_cdc)) %>%
  mutate(asd_diff_cdc_ppv = asd_exp_cdc_ppv - (tp_reclaim_ppv + fn_rep),
         asd_diff_cdc_cat_ppv = ifelse(asd_diff_cdc_ppv > 0, 1, 
                                       ifelse(asd_diff_cdc_ppv < 0, 0, asd_diff_cdc_ppv))) %>% 
  
  mutate(fn_calc_cdc_ppv = fn_rep + ifelse(asd_diff_cdc_ppv > 0, asd_diff_cdc_ppv, 0),
         tn_calc_cdc_ppv = tn_rep - (asd_diff_cdc_ppv * asd_diff_cdc_cat_ppv)) %>% 
  
  # World prevalence
  mutate(asd_exp_world_ppv = as.integer((n + n_positive_miss) * asd_prev_world)) %>%
  mutate(asd_diff_world_ppv = asd_exp_world_ppv - (tp_reclaim_ppv + fn_rep),
         asd_diff_world_cat_ppv = ifelse(asd_diff_world_ppv > 0, 1, 
                                         ifelse(asd_diff_world_ppv < 0, 0, asd_diff_world_ppv))) %>% 
  
  mutate(fn_calc_world_ppv = fn_rep + ifelse(asd_diff_world_ppv > 0, asd_diff_world_ppv, 0),
         tn_calc_world_ppv = tn_rep - (asd_diff_world_ppv * asd_diff_world_cat_ppv)) %>%  
  
  mutate(sens_calc_og_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fn_rep),
         spec_calc_og_ppv = tn_rep / (tn_rep + fp_reclaim_ppv),
         ppv_calc_og_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fp_reclaim_ppv),
         npv_calc_og_ppv = tn_rep / (tn_rep + fn_rep)) %>% 
  
  # Diagnostic accuracy - National
  mutate(sens_calc_natl_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fn_calc_natl_ppv),
         spec_calc_natl_ppv = tn_calc_natl_ppv / (tn_calc_natl_ppv + fp_reclaim_ppv),
         ppv_calc_natl_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fp_reclaim_ppv),
         npv_calc_natl_ppv = tn_calc_natl_ppv / (tn_calc_natl_ppv + fn_calc_natl_ppv)) %>% 
  
  # Diagnostic accuracy - CDC
  mutate(sens_calc_cdc_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fn_calc_cdc_ppv),
         spec_calc_cdc_ppv = tn_calc_cdc_ppv / (tn_calc_cdc_ppv + fp_reclaim_ppv),
         ppv_calc_cdc_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fp_reclaim_ppv),
         npv_calc_cdc_ppv = tn_calc_cdc_ppv / (tn_calc_cdc_ppv + fn_calc_cdc_ppv)) %>% 
  
  # Diagnostic accuracy - World
  mutate(sens_calc_world_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fn_calc_world_ppv),
         spec_calc_world_ppv = tn_calc_world_ppv / (tn_calc_world_ppv + fp_reclaim_ppv),
         ppv_calc_world_ppv = tp_reclaim_ppv / (tp_reclaim_ppv + fp_reclaim_ppv),
         npv_calc_world_ppv = tn_calc_world_ppv / (tn_calc_world_ppv + fn_calc_world_ppv))


