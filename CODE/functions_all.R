# Load libraries
if (require(pacman)){
  install.packages("pacman")
  library(pacman)
}

pacman::p_load(tidyverse, 
               rio, 
               here, 
               janitor, 
               stringr, 
               rlang)


# Import raw data and rename for clarity ####
raw_data <- import(here("input", "PublicationFileFinal 35.csv")) %>% 
  clean_names() %>% 
  rename("sens" = "sensitivity",
         "spec" = "specificity",
         "n_tot" = "total_n",
         "n_positive_miss" = "n_for",
         "asd_prev_natl" = "asd_prev_n",
         "asd_prev_cdc" = "asd_prev_c",
         "asd_prev_rep" = "report_prev",
         "asd_prev_world" = "glob_prev")

# Architecture ####

# Function should do:
  # Calculate using any missing p% 
  # Naming for missing p%
  # Includes any: tp tn fp fn spec sens ppv npv, n n_positive_miss (if not missing, run reclaimn)
  # Extra: include any prev, any conditions
  # Round...

# Proposed steps
  # Step 1: Clean, calculate from existing
  # Step 2a: If not reclaim -> Skip to step 3
  # Step 2b: If reclaim -> Calculate from reclaimn p
  # Step 3: Epi adjustment

# Function codes ####

# Function 1: Calculate new Confusion Matrix values from original data ####
# Input: original dataframe
# Output: tp_rep, fn_rep, fp_rep, tn_rep

recalculated_missing_data <-  function(df){
  zero_cols <- c("tp","tn","fp","fn","n","n_positive_miss")
  
  calculated_df <- df %>% 
    mutate(across(all_of(zero_cols), ~ifelse(is.na(.), 0, .))) %>% 
    mutate(
      sens_calc1 = tp / (tp + fn),
      spec_calc1 = tn / (tn + fp),
      ppv_calc1 = tp / (tp + fp),
      npv_calc1 = tn / (tn + fn)) %>% 
    
    mutate(
      tp_calc_sens_fn1 = (sens * fn) / (1 - sens),
      tp_calc_ppv_fp1 = (ppv * fp) / (1 - ppv),
      
      fn_calc_sens_tp = (tp * (1 - sens))/sens,
      fn_calc_npv_tn = (tn * (1 - npv)/npv),
      
      fp_calc_spec_tn = (tn * (1 - spec)) / spec,
      fp_calc_ppv_tp = (tp * (1 - ppv)) / ppv,
      
      tn_calc_spec_fp = (spec * fp) / (1 - spec),
      tn_calc_npv_fn = (npv * fn) / (1 - npv)) %>% 
    
    mutate(
      tn_calc_sub1 = n - (fn + fp + tp),
      tn_calc_sub2 = n - (fn_calc_sens_tp + fp_calc_ppv_tp + tp)) %>% 
    
    mutate(fn_calc_sub1 = n - (tn + fp + tp), .after = fn_calc_npv_tn) %>%                     # Only one entry, but add 1 for naming consistency
  
    mutate(fp_calc_sub1 = n - (tn + fn + tp), .after = fp_calc_ppv_tp) %>% 
      
    mutate(tp_calc_sens_fn2 = (sens * fn_calc_sens_tp) / (1 - sens), .after = tp_calc_ppv_fp1,
           tp_calc_ppv_fp2 = (ppv * fp_calc_ppv_tp) / (1 - ppv)) %>% 
    
    mutate(tp_rep = tp, .after = tp_calc_ppv_fp2) %>% 
    
    mutate(tn_rep = coalesce(
      ifelse(!is.na(tn), tn, tn_calc_sub1),
      tn_calc_spec_fp,
      tn_calc_npv_fn,
      tn_calc_sub2,
      0
    ),
    .after = tn_calc_sub2) %>% 
    
    mutate(fp_calc_sub2 = n - (tn_rep + fn + tp), .after = fp_calc_sub1) %>%
    
    mutate(fp_rep = coalesce(
      ifelse(!is.na(fp), fp, fp_calc_sub1),
      fp_calc_spec_tn,
      fp_calc_ppv_tp,
      fp_calc_sub2,
      0
    ),
    .after = fp_calc_sub2) %>% 
    
    mutate(fn_rep = coalesce(
      ifelse(!is.na(fn), fn, fn_calc_sub1),
      fn_calc_sens_tp,
      fn_calc_npv_tn,
      ((tn_rep + tp_rep + fp_rep) / n),      
      0
    ),
    .after = fn_calc_sub1) %>%
    
    mutate(sens_calc2 = tp / (tp + fn_rep), .after = sens_calc1) %>% 
    mutate(spec_calc2 = tn_rep / (tn_rep + fp_rep), .after = spec_calc1) %>% 
    mutate(ppv_calc2 = tp / (tp + fp_rep), .after = ppv_calc1) %>% 
    mutate(npv_calc2 = tn / (tn_rep + fn_rep), .after = npv_calc1) %>% 
    mutate(asd_rep = tp + fn_rep) 
}  

test_df <- recalculated_missing_data(raw_data)
export(test_df, "test_df.xlsx")

# Function 1x: Pairwise subtraction for validation ####
# Input: recalculated df from Function 1
# Output: .txt file with console messages about which studies to check
# Note: not usable for reclaimn input

thresholds <- list(
  tp = 10,
  tn = 50,
  fp = 1,
  fn = 1,
  sens = 0.02,
  spec = 0,
  ppv = 0,
  npv = 0
)

validation <- function(df, thresholds, console_output) { # Threshold list needs to be made separately
  
  sink(console_output, append = FALSE) # Starting now, every output messages will be exported to a .txt file
  
  original_vars <- c("tp", "tn", "fp", "fn", "sens", "spec", "ppv", "npv")
  
  # Loop 1.1: over the list of vars above
  for (var in original_vars) {
    
    all_vars <- names(df)[str_starts(names(df), var)]  # Find all vars that starts with original vars (to loop over)
    
    # Start subtraction
    for (i in seq_along(all_vars)) {   # Loop 1.2: get first column col1
      for (j in seq(i + 1, length(all_vars))) {   # Loop 1.3: get second column col2
        
        col1 <- all_vars[i]
        col2 <- all_vars[j]
        
        # These are checks for NAs 
        
        if (!is.na(col1) && is.na(col2)) {
          cat("------------------------- ", paste0(var), "'s loop is finished. Moving to next var -------------------------\n\n", sep = "")
          next
          
        } else if (col1 == col2) {
          next
          
        } else if (is.na(col1) || is.na(col2)) {
          cat(col1, "or", col2, "is missing. Skipped validation.\n\n", sep = "")
          next
        }
        
        comparison_df <- df %>% # Store results in a new df 
          filter(!is.na(.data[[col1]]) & !is.na(.data[[col2]])) %>%
          mutate(
            difference = abs(.data[[col1]] - .data[[col2]]),   # Pairwise subtraction
            flag = ifelse(difference > thresholds[[var]], TRUE, FALSE)  # Flag if above defined thresholds
          )
        
        # Get corresponding studies where flag is TRUE
        flagged_studies <- comparison_df$study[comparison_df$flag == TRUE]
        
        if (length(flagged_studies) > 0){ # If this space isn't empty
          
          cat("Check these studies because of differences > ", thresholds[[var]], " between ", col1, " and ", col2, ":\n", sep = "", # Output this message
              
              paste(flagged_studies, collapse = ", "), "\n\n") 
          
        } else {                         # If there are no flagged studies (0 differences)
          
          if (sum(comparison_df$difference) == 0) {  # If combined sum is 0 (for this loop)
            
            cat("No difference between ", col1, " and ", col2, ".\n\n", sep = "")  # Output this message
            
          }
        }
      }
    }
  }
  cat("------------------------ Out of variables to loop. Validation finished ------------------------")
  sink() # Close output
  
  return(comparison_df) # Return result df 
}

output_filename <- "test.txt"
validation_result <- validation(test_df, thresholds, output_filename)
export(validation_result, "result.xlsx")


# Function 2a_a: Epi adjustment for when no cases are missing from screened (n_positive_miss = 0) ####
# For asd, fixed prevalence values
# Input: Function 1's resulted df
# Output: sens, spec, ppv, npv using national, cdc, world prevalence
   
adjust_no_missed_asd <-  function(df){
    df %>% 
    mutate(asd_exp_natl = as.integer(n * asd_prev_natl),
           asd_exp_cdc = as.integer(n * asd_prev_cdc),
           asd_exp_world = as.integer(n * asd_prev_world)) %>% 

    mutate(asd_diff_natl = asd_exp_natl - (tp + fn_rep),
           asd_diff_cdc = asd_exp_cdc - (tp + fn_rep),
           asd_diff_world = asd_exp_world - (tp + fn_rep)) %>% 
    
    mutate(asd_diff_natl_cat = ifelse(asd_diff_natl > 0, 1, 
                                      ifelse(asd_diff_natl < 0, 0, asd_diff_natl)),
           asd_diff_cdc_cat = ifelse(asd_diff_cdc > 0, 1, 
                                     ifelse(asd_diff_cdc < 0, 0, asd_diff_cdc)),
           asd_diff_world_cat = ifelse(asd_diff_world > 0, 1, 
                                       ifelse(asd_diff_world < 0, 0, asd_diff_world))) %>% 
    
    mutate(fn_calc_natl = fn_rep + ifelse(asd_diff_natl > 0, asd_diff_natl, 0),
           tn_calc_natl = tn_rep - (asd_diff_natl * asd_diff_natl_cat)) %>% 
    
    mutate(fn_calc_cdc = fn_rep + ifelse(asd_diff_cdc > 0, asd_diff_cdc, 0),
           tn_calc_cdc = tn_rep - (asd_diff_cdc * asd_diff_cdc_cat)) %>% 
    
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
           npv_calc_world = tn_calc_world / (tn_calc_world + fn_calc_world))
}
# Function 2b_a: Epi adjustment for when there are kids missed from screened (n_positive_miss != 0) ####
# For asd, fixed prevalence values
# Input: Function 1's resulted df
# Output: sens, spec, ppv, npv using national, cdc, world prevalence

adjust_missed_asd <- function(df, p){
  df %>% 
  mutate(!!paste0("tp_reclaim_", p, "%") := tp_rep + (p/20) * n_positive_miss,             
         !!paste0("fp_reclaim_", p, "%") := fp_rep + (1- (p/20)) * n_positive_miss) %>% 
    
  mutate(!!paste0("asd_exp_natl_", p, "%") := as.integer((n + n_positive_miss) * asd_prev_natl),
         !!paste0("asd_exp_cdc_", p, "%") := as.integer((n + n_positive_miss) * asd_prev_cdc),
         !!paste0("asd_exp_world_", p, "%") := as.integer((n + n_positive_miss) * asd_prev_world)) %>% 
                                              
  mutate(!!paste0("asd_diff_natl_", p, "%") := !!sym(paste0("asd_exp_natl_", p, "%")) - !!sym(paste0("tp_reclaim_", p, "%")) + fn_rep,
         !!paste0("asd_diff_cdc_", p, "%") := !!sym(paste0("asd_exp_cdc_", p, "%")) - !!sym(paste0("tp_reclaim_", p, "%")) + fn_rep,
         !!paste0("asd_diff_world_", p, "%") := !!sym(paste0("asd_exp_world_", p, "%")) - !!sym(paste0("tp_reclaim_", p, "%")) + fn_rep) %>% 
  
  mutate(!!paste0("asd_diff_natl_cat_", p, "%") := ifelse(!!sym(paste0("asd_diff_natl_", p, "%")) > 0, 1,
                                                             ifelse(!!sym(paste0("asd_diff_natl_", p, "%")) < 0, 0, !!sym(paste0("asd_diff_natl_", p, "%")))),
         !!paste0("asd_diff_cdc_cat_", p, "%") := ifelse(!!sym(paste0("asd_diff_cdc_", p, "%")) > 0, 1,
                                                         ifelse(!!sym(paste0("asd_diff_cdc_", p, "%")) < 0, 0, !!sym(paste0("asd_diff_cdc_", p, "%")))),
         !!paste0("asd_diff_world_cat_", p, "%") := ifelse(!!sym(paste0("asd_diff_world_", p, "%")) > 0, 1,
                                                         ifelse(!!sym(paste0("asd_diff_world_", p, "%")) < 0, 0, !!sym(paste0("asd_diff_world_", p, "%"))))) %>% 

  mutate(!!paste0("fn_calc_natl_", p, "%") := fn_rep + ifelse(!!sym(paste0("asd_diff_natl_", p, "%")) > 0, !!sym(paste0("asd_diff_natl_", p, "%")), 0),
         !!paste0("tn_calc_natl_", p, "%") := tn_rep + !!sym(paste0("asd_diff_natl_", p, "%")) * !!sym(paste0("asd_diff_natl_cat_", p, "%")),
         
         !!paste0("fn_calc_cdc_", p, "%") := fn_rep + ifelse(!!sym(paste0("asd_diff_cdc_", p, "%")) > 0, !!sym(paste0("asd_diff_cdc_", p, "%")), 0),
         !!paste0("tn_calc_cdc_", p, "%") := tn_rep + !!sym(paste0("asd_diff_cdc_", p, "%")) * !!sym(paste0("asd_diff_cdc_cat_", p, "%")),
         
         !!paste0("fn_calc_world_", p, "%") := fn_rep + ifelse(!!sym(paste0("asd_diff_world_", p, "%")) > 0, !!sym(paste0("asd_diff_world_", p, "%")), 0),
         !!paste0("tn_calc_world_", p, "%") := tn_rep + !!sym(paste0("asd_diff_world_", p, "%")) * !!sym(paste0("asd_diff_world_cat_", p, "%"))) %>% 
    
  # Diagnostic accuracy - original data 
  mutate(!!paste0("sens_calc_og_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + fn_rep),
         !!paste0("spec_calc_og_", p, "%") := tn_rep / (tn_rep + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("ppv_calc_og_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("npv_calc_og_", p, "%") := tn_rep / (tn_rep + fn_rep)) %>% 
    
  # Diagnostic accuracy - National 
  mutate(!!paste0("sens_calc_natl_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fn_calc_natl_", p, "%"))), 
         !!paste0("spec_calc_natl_", p, "%") := !!sym(paste0("tn_calc_natl_", p, "%")) / (!!sym(paste0("tn_calc_natl_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("ppv_calc_natl_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("npv_calc_natl_", p, "%") := !!sym(paste0("tn_calc_natl_", p, "%")) / (!!sym(paste0("tn_calc_natl_", p, "%")) + !!sym(paste0("fn_calc_natl_", p, "%")))) %>% 
  
  # Diagnostic accuracy - CDC  
  mutate(!!paste0("sens_calc_cdc_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fn_calc_cdc_", p, "%"))), 
         !!paste0("spec_calc_cdc_", p, "%") := !!sym(paste0("tn_calc_cdc_", p, "%")) / (!!sym(paste0("tn_calc_cdc_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("ppv_calc_cdc_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("npv_calc_cdc_", p, "%") := !!sym(paste0("tn_calc_cdc_", p, "%")) / (!!sym(paste0("tn_calc_cdc_", p, "%")) + !!sym(paste0("fn_calc_cdc_", p, "%")))) %>% 
  
  # Diagnostic accuracy - World
  mutate(!!paste0("sens_calc_world_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fn_calc_world_", p, "%"))), 
         !!paste0("spec_calc_world_", p, "%") := !!sym(paste0("tn_calc_world_", p, "%")) / (!!sym(paste0("tn_calc_world_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("ppv_calc_world_", p, "%") := !!sym(paste0("tp_reclaim_", p, "%")) / (!!sym(paste0("tp_reclaim_", p, "%")) + !!sym(paste0("fp_reclaim_", p, "%"))),
         !!paste0("npv_calc_world_", p, "%") := !!sym(paste0("tn_calc_world_", p, "%")) / (!!sym(paste0("tn_calc_world_", p, "%")) + !!sym(paste0("fn_calc_world_", p, "%")))) 
}

# Wrapper (master) function for asd dataset ####
# It will output a dataset, so just input datasets and p% and run
master_function_asd <- function(df) {
  cleaned_data <- recalculated_missing_data(df)
  
  # Divide into 2 datasets based on n_positive_miss
  no_miss_data <- subset(cleaned_data, n_positive_miss == 0)
  miss_data <- subset(cleaned_data, n_positive_miss != 0)
  
  # Function 2a: Apply to rows with n_positive_miss == 0
  if (nrow(no_miss_data) > 0) {
    cat("Found rows with no missed cases -> Running adjustment without reclaim... \n")
    
    adjusted_result_no_miss <- adjust_no_missed_asd(no_miss_data)
    adjusted_result_no_miss_name <- "adjusted_result_no_miss_asd"
    assign(adjusted_result_no_miss_name, adjusted_result_no_miss, envir = .GlobalEnv)
  }
  
  # Function 2b: Apply to rows with n_positive_miss != 0
  if (nrow(miss_data) > 0) {
    cat("Found rows with missed cases -> Running adjustment with reclaim... \n")
    if (!"p" %in% colnames(df)) {
      repeat {
        p_input <- as.numeric(readline(prompt = "Enter the missing percentage (0 to 100): "))
        if (p_input >= 0 && p_input <= 100 && !is.na(p_input)) {
          df$p_missed <- p_input
          p <- p_input
          break
        } else {
          cat("Invalid percentage. Please enter a number between 0 to 100 \n")
        }
      }  
    } else {
      p <- df$p 
    } 
    adjusted_result_miss <- adjust_missed_asd(miss_data, p = p)
    adjusted_result_name_miss <- paste0("adjusted_result_miss_asd_", p)
    assign(adjusted_result_name_miss, adjusted_result_miss, envir = .GlobalEnv)
  }
  
  # Merge back results from 2a and 2b
  if (nrow(no_miss_data) > 0 && nrow(miss_data) > 0) {
    adjusted_result <- bind_rows(adjusted_result_no_miss, adjusted_result_miss)
    adjusted_result_name <- paste0("adjusted_result_asd_", p)
    assign(adjusted_result_name, adjusted_result, envir = .GlobalEnv)
    cat("The output dataset is: ", adjusted_result_name, "\n")
  } else if (nrow(no_miss_data) > 0) {
    cat("Only dataset without missed cases was found. \n")
  } else if (nrow(miss_data) > 0) {
    cat("Only dataset with missed cases was found. \n")
  } else {
    cat("No data to process.\n")
  }
}
# Optional - Excel export for validation
test_df <- master_function_asd(raw_data)
export(test_df, "test_df.xlsx")



# Function 2a_b: Epi adjustment for when no cases are missing from screened (n_positive_miss = 0) ####
# For any condition, any prevalence
# Input: Function 1's resulted df
# Output: sens, spec, ppv, npv using user-defined prevalence

adjust_no_missed_any_prev <-  function(df, prev){
  df %>% 
    mutate(expected = as.integer(n * prev)) %>% 
    
    mutate(difference = expected - (tp + fn_rep),
           difference_multiplier = ifelse(difference > 0, 1, 
                                          ifelse(difference < 0, 0, difference))) %>% 
    
    mutate(fn_calc = fn_rep + ifelse(difference > 0, difference, 0),
           tn_calc = tn_rep - (difference * difference_multiplier)) %>% 
    
    
    # Diagnostic accuracy - any prevalence
    mutate(sens_calc = tp_rep / (tp_rep + fn_calc),
           spec_calc = tn_calc / (tn_calc + fp_rep),
           ppv_calc = tp_rep / (tp_rep + fp_rep),
           npv_calc = tn_calc / (tn_calc + fn_calc)) %>% 
    
    mutate(across(contains("calc"), ~ round(., digits = 4))) 
  
}

# Function 2b_b: Epi adjustment for when there are kids missed from screened (n_positive_miss != 0) ####
# For any condition, any prevalence
# Input: Function 1's resulted df
# Output: sens, spec, ppv, npv using user-defined prevalence

adjust_missed_any_prev <- function(df, p, prev){
  df %>% 
    mutate(!!paste0("tp_reclaim_", p) := tp_rep + (p/100) * n_positive_miss,             
           !!paste0("fp_reclaim_", p) := fp_rep + ((1- (p/100)) * n_positive_miss)) %>% 
    
    mutate(!!paste0("expected_at_", p) := as.integer((n + n_positive_miss) * prev),
           !!paste0("difference_at_", p) := !!sym(paste0("expected_at_", p)) - (!!sym(paste0("tp_reclaim_", p)) + fn_rep),
           !!paste0("difference_multiplier_", p) := ifelse(!!sym(paste0("difference_at_", p)) > 0, 1, 0)) %>%
    
    mutate(!!paste0("fn_calc_", p) := fn_rep + ifelse(!!sym(paste0("difference_at_", p)) > 0, !!sym(paste0("difference_at_", p)), 0),
           !!paste0("tn_calc_", p) := tn_rep + !!sym(paste0("difference_at_", p)) * !!sym(paste0("difference_multiplier_", p))) %>%    
    
    
    # Diagnostic accuracy - original data 
    mutate(!!paste0("sens_calc_og_", p) := !!sym(paste0("tp_reclaim_", p)) / (!!sym(paste0("tp_reclaim_", p)) + fn_rep),
           !!paste0("spec_calc_og_", p) := tn_rep / (tn_rep + !!sym(paste0("fp_reclaim_", p))),
           !!paste0("ppv_calc_og_", p) := !!sym(paste0("tp_reclaim_", p)) / (!!sym(paste0("tp_reclaim_", p)) + !!sym(paste0("fp_reclaim_", p))),
           !!paste0("npv_calc_og_", p) := tn_rep / (tn_rep + fn_rep)) %>% 
    
    # Diagnostic accuracy - new prevalence
    mutate(!!paste0("sens_calc_", p) := !!sym(paste0("tp_reclaim_", p)) / (!!sym(paste0("tp_reclaim_", p)) + !!sym(paste0("fn_calc_", p))), 
           !!paste0("spec_calc_", p) := !!sym(paste0("tn_calc_", p)) / (!!sym(paste0("tn_calc_", p)) + !!sym(paste0("fp_reclaim_", p))),
           !!paste0("ppv_calc_", p) := !!sym(paste0("tp_reclaim_", p)) / (!!sym(paste0("tp_reclaim_", p)) + !!sym(paste0("fp_reclaim_", p))),
           !!paste0("npv_calc_", p) := !!sym(paste0("tn_calc_", p)) / (!!sym(paste0("tn_calc_", p)) + !!sym(paste0("fn_calc_", p)))) %>% 
    
    mutate(across(contains("calc"), ~ round(., digits = 4))) 
}

# Wrapper function for any condition, any prevalence ####
# It will output a dataset, so just input dataset and p% and run
# Condition list can be edited

master_function_any <- function(df) {
  
  conditions <- c(
    "AD - Adjustment Disorder",
    "ADHD - Attention Deficit Hyperactivity Disorder",
    "ADD - Attention Deficit Disorder",
    "ASD - Autism Spectrum Disorder",
    "ASPD - Antisocial Personality Disorder",
    "AVPD - Avoidant Personality Disorder",
    "BED - Binge Eating Disorder",
    "BDD - Body Dysmorphic Disorder",
    "BPD - Borderline Personality Disorder",
    "CD - Conduct Disorder",
    "CRSD - Circadian Rhythm Sleep Disorders",
    "EDNOS - Eating Disorder Not Otherwise Specified",
    "GAD - Generalized Anxiety Disorder",
    "MDD - Major Depressive Disorder",
    "NES - Night Eating Syndrome",
    "OCD - Obsessive-Compulsive Disorder",
    "ODD - Oppositional Defiant Disorder",
    "ON - Orthorexia",
    "PPD - Paranoid Personality Disorder",
    "PTED - Post-Traumatic Embitterment Disorder",
    "PTSD - Post-Traumatic Stress Disorder",
    "SZA - Schizoaffective Disorder",
    "SAFD - Seasonal Affective Disorder",
    "SEAD - Separation Anxiety Disorder",
    "SOAD - Social Anxiety Disorder",
    "SPD - Schizoid Personality Disorder"
    # Abbreviated name - Condition name"
  )
  # Start ####
  choice <- menu(conditions, title = "Select a condition: ")
  if(choice == 0) {
    cat("Please select a valid input")
  }
  else{
    cat("Condition for analysis: ", conditions[choice] , "\n")
    cond <- word(conditions[choice], 1)
  }
  # Check for prev in df
  if (!"prev" %in% colnames(df)) {
    repeat { 
      prev_input <- as.numeric(readline(prompt = "Enter the condition's prevalence (0 to 1): "))
      
      if(prev_input >= 0 && prev_input <= 1 && !is.na(prev_input)) {
        df$prev <-  prev_input
        prev <- prev_input 
        break
      }
      else{
        cat("Invalid prevalence. Please enter a number between 0 to 1 \n")
      }
    }
  }
  else{
    prev <- df$prev 
  }
  # Function 1
  cleaned_data <- recalculated_missing_data(df)
  
  # Divide into 2 datasets based on n_positive_miss
  no_miss_data <- subset(cleaned_data, n_positive_miss == 0)
  miss_data <- subset(cleaned_data, n_positive_miss != 0)
  
  # Function 2a: Apply to rows with n_positive_miss == 0
  if (nrow(no_miss_data) > 0) {
    cat("Found rows with no missed cases -> Running step 2a... \n")
    
    adjusted_result_no_miss <- adjust_no_missed_any_prev(no_miss_data, prev = prev)
    adjusted_result_no_miss_name <- paste0("adjusted_result_no_miss_", cond)
    assign(adjusted_result_no_miss_name, adjusted_result_no_miss, envir = .GlobalEnv)
  }
  
  # Function 2b: Apply to rows with n_positive_miss != 0
  if (nrow(miss_data) > 0) {
    cat("Found rows with missed cases -> Running step 2b... \n")
    if (!"p" %in% colnames(df)) {
      repeat{
        p_input <- as.numeric(readline(prompt = "Enter the missing percentage (0 to 100): "))
        if (p_input >= 0 && p_input <= 100 && !is.na(p_input)) {
          df$p_missed <-  p_input
          p <- p_input
          break
        }
        else{
          cat("Invalid percentage. Please enter a number between 0 to 100 \n")
        }
      }  
    }
    
    else{
      p <- df$p 
    } 
    adjusted_result_miss <- adjust_missed_any_prev(miss_data, p = p, prev = prev)
    adjusted_result_name_miss <- paste0("adjusted_result_miss_", cond, "_", p)
    assign(adjusted_result_name_miss, adjusted_result_miss, envir = .GlobalEnv)
  }
  # Merge back results from 2a and 2b
  if (nrow(no_miss_data) > 0 && nrow(miss_data) > 0) {
    adjusted_result <- bind_rows(adjusted_result_no_miss, adjusted_result_miss)
    adjusted_result_name <- paste0("adjusted_result_", cond, "_", p)
    assign(adjusted_result_name, adjusted_result, envir = .GlobalEnv)
    cat("The output dataset is: ", adjusted_result_name, "\n")
  } 
  else if (nrow(no_miss_data) > 0) {
    cat("Only dataset without missed cases was found. \n")
  } 
  else if (nrow(miss_data) > 0) {
    cat("Only dataset with missed cases was found. \n")
  } 
  else {
    cat("No data to process.\n")
  }
}


# Making testing datasets ####
raw_slice <- slice_head(raw_data, n = 10)
raw_slice_0miss <- raw_slice %>% 
  mutate(n_positive_miss = 0)


# Testing ####
master_function_any(raw_slice)
master_function_any(raw_slice_0miss)







