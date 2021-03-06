## `spind` v2.2.1

- Removed unneeded functions from Namespace

- Makes WRM a lot noisier when the model breaks.

## `spind` v2.2.0
- deprecated `customize_plot` and `plot` arguments in favor of new `plot` methods
for `GEE` and `WRM`. 
    + `GEE` and `WRM` now return the `ggplot` object in their respective classes.
    Users can modify them from there rather than use the clunky interface to 
    modify them within the model function call.
    
- `rvi.plot` and `th.indep` now return
`ggplot` as list items along with the values that they returned before. 
`customize_plot` arguments have been deprecated in every function it appears in.
`plot.ROC` in `th.indep` is also deprecated.

- `covar.plot`'s `customize_plot` argument is also deprecated, but the `plot` 
argument is not, as this controls which statistic is plotted (variance or covariance)
rather than whether or not a plot is generated at all. This function no longer
prints the `ggplot` object by default, but returns it as part of a list.

- deprecated `color.maps` arguments to `upscale` and `adjusted.actuals`.

- Some documentation clarification and typo fixes.

- Updated a couple error messages for consistency across functions.

- Behind the scenes, moved all plotting functions that use `ggplot2` to tidy 
evaluation equivalents. If you're interested, more on that [here](https://github.com/tidyverse/ggplot2/blob/master/NEWS.md).


## `spind` v2.1.3
-Updated citation information.

## `spind` v2.1.2
-Updated plot outputs from `GEE` and `WRM` so that y-axes aren't absurdly packed when autocorrelation values are very large
-Updated documentation for some modeling functions to reflect that factors as expanatory variables are not supported by `spind`
-Added some continuous integration functionality to the development branch to ensure that bugs are caught faster. No package functionality is affected by this

## `spind` v2.1.1
-fix bug in `predict.WRM` that prevented calculating smooth components
-fix bug in `scaleWMRR` that prevented acfft from being called internally
-added `trace` argument to `mmiGEE`
-updated vignette to teach users how to customize plots

## `spind` v2.1.0
-Updated `step.spind` model hierarchy recognition algorithm for improved efficiency and generality. Additionally reformatted code so no error message is produced when the initial model is the best model.
-Fixed examples in `wavevar` and `wavecovar` which called a non-existant function.
-Updated the package vignette for clarity and formatting. Additionally, removed pointless WRM example (padding with mean values) because the particular case being demonstrated was identical to first example.
-Added prettier plots from `GEE`, `WRM`, `covar.plot`, `rvi.plot`, and `th.indep` using `ggplot2` graphics. These functions should produce publication quality graphics, rather than the bare bones, minimalist approach in the previous versions of the package. 
-Provide optional color palettes for graphical outputs of `adjusted.actuals` and `upscale` so they are more aesthetically pleasing.
-Updated argument names for a couple functions to ensure consistency throughout the package.


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
