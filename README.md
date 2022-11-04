# Assorted, linear ASCA (ALASCA) Functions
The [ALASCA package](https://andjar.github.io/ALASCA) is described in the paper [ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods](https://www.frontiersin.org/articles/10.3389/fmolb.2022.962431/full). The paper contains several examples of how the package can be used.

## Installation

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools:install_github(“andjar/ALASCA”, ref = “main”)
```

## Citation
If you have utilized the ALASCA package, please consider citing:

Jarmund AH, Madssen TS and Giskeødegård GF (2022) ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods. *Front. Mol. Biosci.* 9:962431. doi: 10.3389/fmolb.2022.962431

```
@ARTICLE{10.3389/fmolb.2022.962431,
  AUTHOR={Jarmund, Anders Hagen and Madssen, Torfinn Støve and Giskeødegård, Guro F.},
  TITLE={ALASCA: An R package for longitudinal and cross-sectional analysis of multivariate data by ASCA-based methods},
  JOURNAL={Frontiers in Molecular Biosciences},
  VOLUME={9},
  YEAR={2022},
  URL={https://www.frontiersin.org/articles/10.3389/fmolb.2022.962431},       
  DOI={10.3389/fmolb.2022.962431},      
  ISSN={2296-889X}
}
```