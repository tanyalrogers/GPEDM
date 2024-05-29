# GPEDM 0.0.0.9008
* Add `fitSmap` and related functions. 

# GPEDM 0.0.0.9007
* Enable use of `update.default` with `fitGP` by storing function call in output.
* Fix bug so that parameter values from prior model can be used in `initpars`. Initial values for `ve` and `sigma2` previously had to be constrained to (0,1), but can now be in (0,5).
* Create `predict_seq` function.
* Add `fixedpars` option to fix phi, ve, sigma2
* Add `returnGPgrad` option to `predict.GP`

# GPEDM 0.0.0.9006

* Added 2 new vignettes, 2 empirical datasets, and `posdef` function to find nearest `rhomatrix` that is positive definite.  
* Covariance matrix is sometimes non-symmetric due to numerical issues. Force to be symmetric by averaging: `(C+t(C))/2`.  

# GPEDM 0.0.0.9005

* Fixed bugs discovered by Vadim, related to prediction with `newdata`. This now works with (1) rhomatrix, (2) a single prediction row, (3) column numbers specified instead of names.

# GPEDM 0.0.0.9004

* Added `predictmethod="lto"` option to `predict.GP`, which leaves out points with the same time index (relevant for multi-population models).

# GPEDM 0.0.0.9003

* Added `getrhomatrix` function, and `rhomatrix` argument to `fitGP`. This allows you to specify a fixed pairwise rho matrix for a hierarchical model.
* Added vignette for VS-EDM and new example dataset.

# GPEDM 0.0.0.9002

* First logged development version.
* VS-EDM with augmentation added but not fully tested.
