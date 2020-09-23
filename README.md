# PRSflex
In Progress


## Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("caret") # if not already installed 
install.packages("stats") # if not already installed 
devtools::install_github("gabrielrvsc/HDeconometrics") # if not already installed 


devtools::install_github("SharonLutz/PRSflex")
```
Note: caret uses ‘ggplot2’, which was built under R version 4.0.2. You may need to update your version of R.
## Example
```
library(PRSflex)
?PRSflex
PRSflex(nSim=50)
```
