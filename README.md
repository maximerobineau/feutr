# feutr
A Package for Evaluating Forecast Efficiency and Unbiasedness in R

This package currently operates as functions. The intention is to eventually convert them into methods and publish it on CRAN.

You can clone or download this repository and reference it in a .R file using `source("path_to_feutr.R")`.

The minimal structure should be `feutr(Data, h)`, where `Data` is a data.frame object of dimensions `n` by `n+1` (with the first column containing dates) and `h` is the forecast horizon (a number or a vector) to study.

Additionnal arguments may be passed to perform different tests and control for various factors. Users may refer to the source code until more comprehensive documentation is available.
