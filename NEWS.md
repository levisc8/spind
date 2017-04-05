## `spind` v2.0.1
-Updated vignette to newest version. Old version is not well developed.
-added citation() to package 



## `spind` v2.0.0 
introduces functions for implementing GEEs and WRMs in the context of spatial models in R. We are also including multiple model selection options for you to find the best model possible for your data. We have retained all of the original functions to calculate spatially corrected measures of model fit.

## GEEs
We've added wrappers for the `gee` and `geese` functions found in packages `gee` and `geepack` so that they can be implemented for spatial models. We've also added S3 methods for `summary` and `predict` to make them easy to interact with. 

## WRMs
In addition to wrappers for GEEs, we've added wrappers for functions found in the `waveslim` package that allow you to implement wavelet-revised models on spatial data sets.  As with GEEs, we've also added S3 methods for `predict` and `summary`. 

## Multi-model inference and model selection
We've added functions to conduct multi-model inference and stepwise model selection (backwards only for the time being), as well as a slew of helper and utility functions to examine your model in greater detail. 

## Previous functionality
The spatial indices included in the original release should still work as they did before.
