[![Build Status](https://travis-ci.org/levisc8/spind.svg?branch=master)](https://travis-ci.org/levisc8/spind)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/levisc8/spind?branch=master&svg=true)](https://ci.appveyor.com/project/levisc8/spind)
[![cran checks](https://badges.cranchecks.info/worst/spind.svg)](https://badges.cranchecks.info/pkgs/spind)
[![CRAN_Status_Badge](https://badges.cranchecks.info/summary/spind.svg)](http://cran.r-project.org/package=spind)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/spind)](https://cran.r-project.org/package=spind)
[![codecov](https://codecov.io/gh/levisc8/spind/branch/master/graph/badge.svg)](https://codecov.io/gh/levisc8/spind)
[![DOI](https://zenodo.org/badge/DOI/10.3897/BDJ.6.e20760.svg)](https://doi.org/10.3897/BDJ.6.e20760)
[![DOI](https://zenodo.org/badge/DOI/10.1111/ecog.02593.svg)](https://doi.org/10.1111/ecog.02593)


# spind 

`spind` is a package for dealing with spatially autocorrelated data in instances of 2-D gridded data sets. We provide you with a suite of functions that can compute spatial models based on Generalized Estimating Equations (GEEs) and Wavelet-Revised Models (WRMs) based on the work of Gudrun Carl, Ingolf Kuehn, and Sam Levin. We also provide an array of options to allow for scale-specific investigations, conduct model selection procedures, and assess goodness of fit for your model using spatially corrected statistics.  

A more complete introduction to the package and its functionality is in the vignette (`vignette('spind_vignette', 'spind'`)) The most recent version will be in the /vignettes folder in this repo, as updates are not sent to CRAN as frequently. Additionally, you can find papers introducing [`spind 1.0` in Ecography](https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.02593) and [`spind > 2.0`](https://bdj.pensoft.net/articles.php?id=20760) in Biodiversity Data Journal. 

This is the development version of `spind` and may not be totally stable/bug free. The most recent stable version of this package is available on CRAN. If you encounter any bugs, please create an issue in this repo and we will try to resolve them as soon as we can. 

**Install stable version from CRAN or development version from Github**

*from CRAN*

`install.packages('spind')`

*from Github*

`devtools::install_github('levisc8/spind')`
