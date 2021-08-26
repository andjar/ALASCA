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