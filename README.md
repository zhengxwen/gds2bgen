gds2bgen: Format Conversion from BGEN to GDS
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Description

This package provides functions for format conversion from [bgen](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) files to [SeqArray GDS](https://github.com/zhengxwen/SeqArray) files.


## Version

v0.9.0


## Package Maintainer

Dr. Xiuwen Zheng ([zhengxwen@gmail.com](zhengxwen@gmail.com))


## Installation

* Installation from Github:
```R
library("devtools")
install_github("zhengxwen/gds2bgen")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.

Or manually intall the package
```sh
git clone https://github.com/zhengxwen/gds2bgen
cd gds2bgen/src
tar -vxzf gavinband-bgen-0b7a2803adb5.tar.gz
cd gavinband-bgen-0b7a2803adb5
./waf configure
./waf
cp build/libbgen.a ..
cp build/3rd_party/zstd-1.1.0/libzstd.a ..
rm -rf build
cd ../../..
R CMD INSTALL gds2bgen
```


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

bgen_fn <- system.file("extdata", "example.8bits.bgen", package="gds2bgen")
seqBGEN_Info(bgen_fn)

## bgen file: gds2bgen/extdata/example.8bits.bgen
## # of samples: 500
## # of variants: 199
## compression method: zlib
## layout version: v1.2
## sample id: sample_001, sample_002, sample_003, sample_004, ...


# example.8bits.bgen ==> example.gds, using 4 cores
seqBGEN2GDS(bgen_fn, "example.gds", float.type="packed8",
    geno=TRUE, dosage=TRUE, prob=TRUE, parallel=4)


library(SeqArray)
(f <- seqOpen("example.gds"))
seqClose(f)

## File: example.gds (137.7K)
## +    [  ] *
## |--+ description   [  ] *
## |--+ sample.id   { Str8 500 LZMA_ra(7.02%), 393B } *
## |--+ variant.id   { Int32 199 LZMA_ra(33.9%), 277B } *
## |--+ position   { Int32 199 LZMA_ra(60.6%), 489B } *
## |--+ chromosome   { Str8 199 LZMA_ra(15.7%), 101B } *
## |--+ allele   { Str8 199 LZMA_ra(11.8%), 101B } *
## |--+ genotype   [  ] *
## |  |--+ data   { Bit2 2x500x0 LZMA_ra, 18B } *
## |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
## |  \--+ extra   { Int16 0 LZMA_ra, 18B }
## |--+ phase   [  ]
## |  |--+ data   { Bit1 500x0 LZMA_ra, 18B } *
## |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
## |  \--+ extra   { Bit1 0 LZMA_ra, 18B }
## |--+ annotation   [  ]
## |  |--+ id   { Str8 199 LZMA_ra(18.6%), 321B } *
## |  |--+ qual   { Float32 199 LZMA_ra(11.8%), 101B } *
## |  |--+ filter   { Int32 199 LZMA_ra(11.3%), 97B } *
## |  |--+ info   [  ]
## |  \--+ format   [  ]
## |     |--+ DS   [  ] *
## |     |  \--+ data   { PackedReal8U 500x199 LZMA_ra(55.6%), 54.0K } *
## |     \--+ GP   [  ] *
## |        \--+ data   { PackedReal8U 500x398 LZMA_ra(38.8%), 75.3K } *
## \--+ sample.annotation   [  ]
```


## Also See

