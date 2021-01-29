# Installation

```{r}
devtools::install_github("andjar/RMASCA",ref="main")
```

# Notice

This is work in practice. The validation is not fully working yet.

# Usage

For a workng example, see the introductory vignette.

RMASCA expects a data frame, eg. `df`, with at least the following columns (with these exact names)

* `time` Either factor, string or integer. Defines when a sample is taken.
* `variable` Either a factor or a string. The measured variable.
* `group` Either a factor, string or integer. Defines the group of a participant

In addition you need to define your model, including at least one random effect (usually participant), for example `mod <- value ~ time*group + (1|partid)`. The participant column has to be specified with `participantColumn`. Usually, you also want to test the robustness of your model.

```{r}
mod.RMASCA <- RMASCA(df = df, formula = mod, participantColumn = "partid", validate = TRUE)
```