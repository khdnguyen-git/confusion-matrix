# Confusion-matrix

## 1. Input Requirements

The input dataset should have the following variable names:

- `tp`, `tn`, `fp`, `fn`: confusion matrix values
- `n`: total sample size
- `n_positive_miss`: number of missed positive cases
- `sens`, `spec`, `ppv`, `npv`: sensitivity, specificity, positive predictive value, negative predictive value
- `asd_prev_natl`, `asd_prev_cdc`, `asd_prev_world`: prevalence values for ASD (national, CDC, world)
  - Not required if analyzing other conditions.
  - Can be renamed to other conditions by using the `Find/Replace All` option RStudio or any text editors. 

For non-ASD conditions:
- Include a `prev` column with the condition's prevalence
  - Can be skipped as the user will be prompted to provide a `prev` value if missing from the input dataset.

If recalculation is needed:
- Include a `p` column with the percentage of missed cases
  - Can be skipped as the user will be prompted to provide a `p` value if missing from the input dataset.


## 2. Functions
- **Function 1** - `recalculated_basic()`: Calculates new confusion matrix values from original data
- **Function 1x** - `validation()`: Performs pairwise column subtraction on the output of **Function 1** for validation 
- **Function 2a_a** - `adjust_no_missed_asd()`: Performs epidemiological adjustment when no cases are missing (ASD prevalences)
- **Function 2a_b** - `adjust_no_missed_any_prev()`: Performs epidemiological adjustment for any condition when no cases are missing
- **Function 2b_a** - `adjust_with_missed_asd()`: Performs epidemiological adjustment when cases are missing (ASD prevalences)
- **Function 2b_b** - `adjust_missed_any_prev()`: Performs epidemiological adjustment for any condition when cases are missing
- **Function 3a** - `master_function_asd()`: Wrapper function for ASD dataset analysis
- **Function 3b** - `master_function_any()`: Wrapper function for any condition analysis


## 3. Output 

### 3.a Naming Scheme for Output Datasets:

1. For ASD only (from **Function 3a** `master_function_asd`):
   - Without missed cases: `adjusted_result_no_miss_asd`
   - With missed cases: `adjusted_result_miss_asd_{p}`
       - `{p}` is the percentage of missed cases

2. For other conditions (from **Function 3b** `master_function_any`):
   1) Without missed cases (observations where `n_positive_miss` = 0) -> `adjusted_result_no_miss_{cond}`
   2) With missed cases (observations where `n_positive_miss` != 0) -> `adjusted_result_miss_{cond}_{p}`
   3) Combined results (if both scenarios exist) -> `adjusted_result_{cond}_{p}`
   4) Validation results: `validation_result` (from the `validation` **Function 1x**)

Where:
- `{cond}` is the name of the condition (e.g., "ADHD")
- `{p}` is the percentage of missed positive cases

Output variable names include suffixes `calc_`, `reclaim_`, followed by prevalence and/or missing percentage.

### 3.b Purpose:

1) Datasets without missed cases: These contain recalculations for studies where all cases were initially captured. They use the original data without adjustments for missed cases.

2) Datasets with missed cases: These contain recalculations for studies where some cases were initially missed. They incorporate adjustments based on the percentage of missed cases.

3) Combined dataset: When both types of studies (with and without missed cases) are present in the input data, this dataset merges the results for a comprehensive view.

4) Validation dataset: This helps identify potential discrepancies in the recalculations from **Function 1**


## 4. Usage
- For ASD analysis: `master_function_asd(R dataset)`
- For other conditions: `master_function_any(R dataset)`


## 5. Psychological Conditions

- The code includes a list of psychological conditions that can be selected for analysis. 
- This list can be modified to any condition by editing the `conditions` vector in the `master_function_any()` function.
    - For readability and simplicity sake, use the format: "Abbreviated condition name - Condition name".

# Note
- **Function 1x** run on the results of **Function 1**. It is not included in the master functions `master_function_asd` or `master_function_any`
- The combined dataset `adjusted_result_{cond}_{p}`, when available, should be considered the "final" dataset as it includes all combined adjusted values.
- The code uses the `tidyverse` package for data manipulation
