# Predicting causal relationships from biological data: Comparison of computational methods for automated causal discovery applied to mass cytometry data of human immune cells

## backShift applied to mass cytometry data - R Code

### 0) Dependencies
Install R packages `backShift`, `CompareCausalNetworks`, and `analyzeMC`. `backShift` and `CompareCausalNetworks` are available from [https://github.com/christinaheinze](). `analyzeMC` is provided in this submission. The latter can be installed as follows.

```
install.packages("/PATH/TO/analyzeMC_0.1.1.tar.gz", type = "source", repos = NULL)
```

### 1) Create .rds files from .mat files
Adjust file paths in script `scripts/readMC.R` and run the script to obtain .rds files from .mat files. 

### 2) Run backShift on inhibitor data
Set files paths and adjust options in `scripts/configCytoEuler_inhibitor.R`. Run `scripts/massCyto_runs_inhibitor.R` which uses the config file `scripts/configCytoEuler_inhibitor.R`. 

### 3) Aggregate results
Adjust file paths in script `scripts/aggregate.R` and run the script to aggregate the results as described in the manuscript.
