# singR

singR package is built on SING method
<https://github.com/thebrisklab/SING>. SING is used to extract joint and
individual non-gaussian components from different datasets. This is a
tutorial example supporting the paper **Simultaneous Non-Gaussian
Component Analysis (SING) for Data Integration in Neuroimaging Benjamin
Risk, Irina Gaynanova** <https://arxiv.org/abs/2005.00597v1>

## Installation

You can install singR from github with:

``` r
library(devtools)
install_github("thebrisklab/singR")
```

## Quick start guide

#### Load package and data

##### for the quick start, we use a compressed version to accelerate the computation, which is a subpart of correlation matrix and 2k-resolution dtseries data.

``` r
# Load the package
library(singR)

# Read and visualize data
load(file = "Small_Simulated_data.Rdata")
# It contains dX, dY, mj, new_sIx,new_sIy,new_sjx,new_sjy

## True Data and signchange
Sxtrue = t(new_sjx) #dim(Sxtrue) px x n
Sytrue = t(new_sjy)

Sxtrue = signchange(Sxtrue) #sign degree amplification
Sytrue = signchange(Sytrue)


# plot for the true component
out_true1 = plotNetwork_change(Sytrue[,1], title='Truth',qmin=0.005, qmax=0.995, path = 'new_mmp.csv') 
out_true2 = plotNetwork_change(Sytrue[,2], title='Truth',qmin=0.005, qmax=0.995, path = 'new_mmp.csv') 

# function plotNetwork_change is tailored for this compressed version data.
# The original function is called plotNetwork, which is set for the standard version of correlation matrix.
p1=out_true1$netmatfig
p2=out_true1$loadingsfig
p3=out_true2$netmatfig
p4=out_true2$loadingsfig
```

![](fig/Truth_Comp1_matrix.png)![](fig/Truth_Comp1_Hist.png)![](fig/Truth_Comp2_matrix.png)![](fig/Truth_Comp2_Hist.png)
