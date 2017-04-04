## Resubmission

This is a resubmission of Spind_2.0.0. In this version,
we have corrected the code for WRM.R which caused
the example to pass testing on x64, but not i386.
We have retested with winbuilder using r-devel,
r-release, and r-oldrelease. 

All passed with no ERRORS OR WARNINGS.

R-oldrelease had the following NOTE:
* checking dependencies in R code ... NOTE
Namespace in Imports field not imported from: 'sp'
  All declared Imports should be used.

'sp' is listed in IMPORTS because it is listed
DEPENDS for 'splancs', which our package does use.

This note was not present in checks on more recent 
versions of R, so we're guessing it's OK. If not, 
please advise. 

ALL builds had the following NOTE: 
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sam Levin <levisc8@gmail.com>'
New submission

Based on information I found online, this is not
actually an issue with the package. Please advise
if this is incorrect.

Also tested on:
* ubuntu 12.04 (on travis-ci), R 3.3.2
* macOS Sierra 10.12.3

*Built cleanly on these two platforms with no 
ERRORS, WARNINGS, or NOTES


## Downstream dependencies
There are currently no downstream dependencies for this package.
