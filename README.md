# Simulation and Application Data of the Paper: A Regularized Cox’s Hierarchical Model for Incorporating Annotation Information in Predictive Omic Studies

## To implement the algorithm described in the paper, install **xrnet** package with survival module

``` r
# Survival module only available under Development branch
devtools::install_github("USCbiostats/xrnet", ref = "development")
```

OS-specific prerequisites
    -   *Windows*: Install
        [RTools](https://cran.r-project.org/bin/windows/Rtools/) (not an
        R package)
        
    -   *Mac*: If using R version &gt;= 4.3.0 and above, verify your GNU Fortran
        version is &gt;= 12.2. If you have an older version, go
        [here](https://cran.r-project.org/bin/macosx/tools/) to install
        the required version
        
    -   *Mac*: If using R version &gt;= 3.6.0 and &lt; 4.3.0, verify your GNU Fortran
        version is &gt;= 6.1. If you have an older version, go
        [here](https://cran.r-project.org/bin/macosx/tools/) to install
        the required version    

## Simulation data generated in the paper, section 3 

The simulation data are generated, and analyzed using the function **simCoxExt()** in the script, *fast_regularized_cox.R*. The simulation results in section 3 of the paper can be reproduced with the script, *sim_Cox_V2.R*.

## Application in the paper, section 4.2, Anti-PD1 Immunotherapy Predictive Biomarker for Melanoma Survival

Data used in section 4.2 of the paper is in the folder, anti-PD1. 

