# GPEDM 0.0.0.9004

* Added `predictmethod="lto"` option to `predict.GP`, which leaves out points with the same time index (relevant for multi-population models).

# GPEDM 0.0.0.9003

* Added `getrhomatrix` function, and `rhomatrix` argument to `fitGP`. This allows you to specify a fixed pairwise rho matrix for a hierarchical model.
* Added vignette for VS-EDM and new example dataset.

# GPEDM 0.0.0.9002

* First logged development version.
* VS-EDM with augmentation added but not fully tested.
