# ALASCA 0.0.0.98

## Improvements

* Nicer default plots
* Minor improvements in performance

## New features

* `plotresiduals()` will plot residuals as a QQ-plot

## Bug fixes

* ALASCA failes on large datasets when the number of PCs is different between variables

# ALASCA 0.0.0.98

## New features

* `plot.ALASCA()` will now plot components if `component` is a vector with more than one value
* `plotComponents()` plots two components against each other

# ALASCA 0.0.0.97

## Improvements

* Better plotting with the `viridis()` palette, different line types, and inherited plot settings from the ALASCA object

## Bug fixes

* `assessGroupDifferences()` should now be working

# ALASCA 0.0.0.96

## New features

* `assessGroupDifferences()` now supports `lm()`
* One can now use Cox regression to take limits of detection into account. To do this, specify `method = "KM"` (using `survival::coxph()`) for ASCA or `method = "KMM"` (using `coxme::coxme()`) for RM-ASCA. You should be careful with scaling; `scaleFun <- function(x){return(x)}` is an alternative. Next, you need to make a data frame and provide it to the function `lowerLimit = lowerLimit`. `lowerLimit` should have a column `value` with the lower cutoff values, and `variable` for the variable names (must correspond to the variables in `df`!)
  * Note: This function is in alpha and has not been thoroughly tested!

## Breaking changes

* The choice to use `Rfast` is now moved to a separate parameter: `useRfast` which defaults to `TRUE`. The goal is to implement `Rfast` for linear models as well.

# ALASCA 0.0.0.95

## New features

* `assessGroupDifferences()` applies `emmeans::emmeans()` to calculate pairwise group differences

# ALASCA 0.0.0.94

## New features

* `flipIt()` can now selectively flip either time or group by specifying `effect = "time"` or `effect = "group"` (default: `effect = "both"`)

# ALASCA 0.0.0.93

## New features

* Set `save = TRUE` when running `ALASCA()` to automatically save the model and the plots you make. (Only tested on Windows)

## Bux fix

* Fixed: `ALASCA::summary()` is now working

# ALASCA 0.0.0.92

## Bux fix

* Fixed: `flipIt()` did not flip the upper limit for score of the group effect correctly

# ALASCA 0.0.0.91

## New features

* `validationMethod = "permutation"` now produces actual P values.
  * P values are kept in `$pvals`
  * The lowest P values will be set to 1/nValRuns, so it may be wise to use eg. `nValRuns=1001` to be able to detect P<.001
  * Plotting this model adds asterisks indicating significance
* `validationMethod = "bootstrap"` is now functional

# ALASCA 0.0.0.9

## Known bugs

* Reference level cannot be set when using multiple interactions
* `plotPart()` sometimes fails to provide more than a blank plot