# Assorted, linear ASCA (ALASCA) Functions
## Installation

```{r}
devtools::install_github("andjar/ALASCA", ref="main")
```

## Notice

<div class="alert alert-danger" role="alert">
  The documentation is under revision
</div>

## Usage

Let `df` be a long-format data frame with repeated measurements and the columns `ID`, `time`, `group`, `variable` and `value`;

```{r}
mod.ALASCA <- ALASCA(df = df, formula = value ~ time * group + (1|ID), validate = TRUE)
```

## Background
This implementation of ALASCA is based on the MATLAB version presented in [Repeated measures ASCA+ for analysis of longitudinal intervention studies with multivariate outcome data](https://www.medrxiv.org/content/10.1101/2020.12.03.20243097v1) by Torfinn S. Madssen, Guro F. Giskeødegård, Age K. Smilde and Johan A. Westerhuis. They are, however, not involved in the development of this R package.
