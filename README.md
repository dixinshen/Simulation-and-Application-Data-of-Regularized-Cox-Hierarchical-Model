# Simulation and Application Data of the Paper: A Regularized Cox’s Hierarchical Model for Incorporating Annotation Information in Predictive Omic Studies

## 1.  To implement the algorithm described in the paper, install **xrnet** package with survival module

1.  OS-specific prerequisites
    -   *Windows*: Install
        [RTools](https://cran.r-project.org/bin/windows/Rtools/) (not an
        R package)
    -   *Mac*: If using R version &gt;= 4.3.0, verify your GNU Fortran
        version is &gt;= 12.2. If you have an older version, go
        [here](https://cran.r-project.org/bin/macosx/tools/) to install
        the required version
2.  Install the **xrnet** package with the *install\_github()* function
    
<!-- end list -->

``` r
# Survival module only available under Development branch
devtools::install_github("USCbiostats/xrnet", ref = "development")
```

## 2.  Simulation data generated in the paper, section 3 

The simulation data are generated, and analyzed using the function **simCoxExt()** in the script, *simulation_function.R*. The simulation results in section 3 of the paper can be reproduced with the script, *sim_xrnetcox.R*.

## 3.  Application in the paper, section 4.2, Anti-PD1 Immunotherapy Predictive Biomarker for Melanoma Survival

Data used in section 4.2 of the paper is in the folder, *anti-PD1 Data*.

## 4.  C++ Implementation of Regularized Cox’s Hierarchical Model

The C++ implement of Regularized Cox's Hierarchical Model, which can be imported to R by Rcpp, are in the folder, *xrnet cox C++ implementation*.

