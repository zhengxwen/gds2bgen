gds2bgen: Format Conversion Between GDS and BGEN
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Description

This package provides functions for format conversion between [bgen](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) files to [SeqArray](https://github.com/zhengxwen/SeqArray) GDS files.


## Package Maintainer

Dr. Xiuwen Zheng ([zhengxwen@gmail.com](zhengxwen@gmail.com))


## Installation

* Installation from Github:
```R
library("devtools")
install_github("zhengxwen/gds2bgen")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.


## Copyright Notice

This package includes the sources of the bgen library written by Gavin Band and
Jonathan Marchini (https://bitbucket.org/gavinband/bgen), Boost (the C++
libraries, https://www.boost.org) and Zstandard (https://zstd.net).


## Citations for GDS

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Examples

```R
library(gds2bgen)
```


## Also See

[gds2bcf](https://github.com/zhengxwen/gds2bcf): Format Conversion Between GDS and BCF