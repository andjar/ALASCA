# Assorted, linear ASCA (ALASCA)
## Installation

```{r}
devtools::install_github("andjar/ALASCA",ref="main")
```

## Notice

<div class="alert alert-danger" role="alert">
  This is work in process and should not be relied upon at the moment.
</div>

## Usage

For a working example, see [Get Started](articles/ALASCA.html). For different case examples, see

* [Repeated-measures observational data](articles/pregnancy.html)
* [Repeated-measures intervention data](articles/metabolomics.html)
* [Single-measures observational data](articles/personality.html)

ALASCA expects a data frame, eg. `df`, with at least the following columns (with these exact names)

* `time` Either factor, string or integer. Defines when a sample is taken.
* `variable` Either a factor or a string. The measured variable.
* `group` Either a factor, string or integer. Defines the group of a participant

In addition you need to define your model, including at least one random effect (usually participant), for example `mod <- value ~ time*group + (1|partid)`. The participant column has to be specified with `participantColumn`. Usually, you also want to test the robustness of your model.

```{r}
mod.ALASCA <- ALASCA(df = df, formula = mod, participantColumn = "partid", validate = TRUE)
```

## Background
This implementation of ALASCA is based on the MATLAB version presented in [Repeated measures ASCA+ for analysis of longitudinal intervention studies with multivariate outcome data](https://www.medrxiv.org/content/10.1101/2020.12.03.20243097v1) by Torfinn S. Madssen, Guro F. Giskeødegård, Age K. Smilde and Johan A. Westerhuis. They are, however, not involved in the development of this R package.