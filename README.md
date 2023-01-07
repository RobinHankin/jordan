Jordan algebras in R
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/jordan.png" width = "150" align="right" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/jordan?color=green)](https://cran.r-project.org/package=jordan)
![](https://cranlogs.r-pkg.org/badges/grand-total/jordan?color=green)
![](https://cranlogs.r-pkg.org/badges/jordan?color=green)
![](https://cranlogs.r-pkg.org/badges/last-week/jordan?color=green)

A *Jordan algebra* is a non-associative algebra over the reals with a
bilinear multiplication that satisfies the following identities:

![xy=yx](https://latex.codecogs.com/png.latex?xy%3Dyx "xy=yx")

![(xy)(xx)=x(y(xx))](https://latex.codecogs.com/png.latex?%28xy%29%28xx%29%3Dx%28y%28xx%29%29 "(xy)(xx)=x(y(xx))")

(the second identity is known as the Jordan identity). In literature one
usually indicates multiplication by juxtaposition but one sometimes sees
![x\circ y](https://latex.codecogs.com/png.latex?x%5Ccirc%20y "x\circ y").
Package idiom is to use an asterisk, as in `x*y`. There are five types
of Jordan algebras, all of which have their own representation in the
package, which uses S4 methods.

# Further information

For more detail, see the package vignette

`vignette("jordan")`
