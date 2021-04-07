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