# Install 

devtools::install_github("bitmask/labnetmet")
library(labnetmet)

`


# Labelled network metrics

This R package provides a number of metrics that defined distances between small, labelled networks.  The distances and networks can also be visualized.  See the examples for usage. 

## Transitive metric

Distance between two graphs is defined as the number of pairs of nodes that have a path between them (through any number of nodes) in one graph but not in the other. 
