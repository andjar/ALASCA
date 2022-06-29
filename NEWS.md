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