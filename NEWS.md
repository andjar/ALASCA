# ALASCA 1.0.19

* Bug fix: The prediction plots should now work for models with separated effects too

# ALASCA 1.0.18

* Replacing `aes_string()` due to changes in the `ggplot2` package
* Moved `finalize()` from public to private re changes in the `R6` package
* Bug fix: Plots all variables instead of subset

# ALASCA 1.0.17

* Minor bug fix related to 1.0.16

# ALASCA 1.0.16

* Bug fix: Plots show wrong number of variables when compared to `n_limit` for odd numbers
  * Solves issue [https://github.com/andjar/ALASCA/issues/12](Plot function prnts output as showing x number of variables, but not actually)

# ALASCA 1.0.14

* New function: `predict_scores()`. It accepts a data table with columns `variable` and `value`, and returns a new data table with a score column

# ALASCA 1.0.11

* New features:
  * `optimize_PCs = TRUE` (default: `FALSE`) will check if significant principal components have to be re-ordered during bootstrapping. This may happen if PCs are close in explanatory value so that they are shuffled during bootstrap. If this happens, it will trigger a warning
  * `waterfall = TRUE` (default: `FALSE`) in `plot()` will replace points with bars for loadings. This can be very nice in combination with `loading_group_column`

# ALASCA 1.0.10

* New feature: `plot()` now accepts the argument `sort_loadings` to control the order of the loading variables
  * `sort_loadings = "loading"` sorts the variables by loading (default)
  * `sort_loadings = "alpha"` sorts the variables alphabetically
  * `sort_loadings = c(...)` sorts the variables in the same order as `c(...)`, where `...` is the variables of interest. Note that it may be required to increase `n_limit` (or use `n_limit = 0`) to ensure that all variables are shown. It is recommended to use `n_limit = 0`, i.e., `plot(..., sort_loadings = c(...), n_limit = 0)`

# ALASCA 1.0.9

* Fix: Error when running LMMs without scaling

# ALASCA 1.0.8

* Fix: Error when running linear models with Rfast (`Error in crossprod(x, y) : requires numeric/complex matrix/vector arguments`)

# ALASCA 1.0.7

* Fix: Error for small datasets
* Fix: Error for custom stratification columns

# ALASCA 1.0.5

* New feature: Black-and-white mode for more plot types
* Fix: Crash when combining `use_Rfast = FALSE` and another random intercept than `ID`
* Fix: Crash when trying to use only a three-way interaction as effect. Still a but unstable

# ALASCA 1.0.4

* New feature: Plot effects in gray scale/black-and-white with symbols instead of colors. Can be tested with `plot(..., bw = TRUE)` (or `grayscale = TRUE` or `greyscale = TRUE`) or similarly with `ALASCA(..., plot.bw = TRUE)`
* Fix: Prediction plot with single-variable effect (e.g., `time`) did not color the groups correctly

# ALASCA 1.0.3

* Fix: Crash when some participants are missing certain measurements

# ALASCA 1.0.2

* New feature: Permutation testing (`validation_method = "permutation"`)
  * Simple permutation testing where data labels are shuffled at two levels: either the participant is re-assigned (e.g., a participant is randomly moved to a new group (or not)), or labels are shifted *within* participant (e.g. the time labels for a participant are shuffled)
  * By default, the first effect is assumed to be shuffled *within* participant and the others shuffled across *participant*. The default can be overwritten by specifying `permutation_within_participants` and `permutation_across_participants`, e.g. `permutation_within_participants = c("time")`
  * The participants should only belong to *one* group for each of the variables in `permutation_across_participants` and samples will be be reassigned as a block for `permutation_within_participants` (i.e., if a participant has two samples in group A and one sample in group B, then the two former samples will be reassigned *together* and not individually)
* Improved performance: `df["value"][ rowNumers ]` is somewhat faster than `df[ rowNumers, value ]`

# ALASCA 1.0.1

* Fix: Error when using another column name than `ID` for ID

# ALASCA 1.0.0

* First release with new framework