# dbnorm
![dbnorm](https://user-images.githubusercontent.com/37698532/69883192-4f890280-12d3-11ea-8b3b-13f110949e98.jpg)



*dbnorm* package, drift across batch normalization, contains R functions which allow visualization and 
removal of technical heterogeneity from large metabolomics dataset. By including advanced statistical tools, the *dbnorm* package allows user to inspect the structure and quality of large metabolomics datasets both in macroscopic and microscopic scales at the sample batch level and metabolic feature level, respectively.
It allows users to efficiently correct drift across batch and to adjust large metabolomics datasets for technical variation which helps improving the estimation of the biological mechanisms underlying disease condition or medical state.
*dbnorm* includes 11 distinct functions for pre-processing of data and estimation of missing values, 
conventional functions for batch effect correction based on statistical models, as well as functions using advanced statistical 
tools to generate several diagnosis plots to help users to choose the statistical model which better fits to their data 
structure. *dbnorm* implements several statistical models, including ComBat(parametric and non-parametric)-model[PMID: 16632515]  from *svs* package [PMID: 22257669]
, that was already in use for metabolomics data normalization, and ber function [DOI 10.1007/s12561-013-9081-1], priorly developed for microarray gene expression data,
, that we propose here as a new approach for correction of drift across batch in metabolomics datasets. 

## Getting started
### Step1: installation
Install package dependencies in CRAN and Bioconductor via cole bellow:

> R installer
```
install.packages(c("ggplot2","parallel","reshape2","plyr",
"mixOmics","knitr","tibble","installr","fs","Rcp", "rmarkdown","processx","backports",
"bootstrap","boot","caret","dplyr","stringr","ggfortify","factoextra","NormalizeMets","MASS","ber",
"RColorBrewer","RCurl","lattice","data.table","igraph","tidyr","scales",
"e1071","fpc"))
```
> Bioconductor installer
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("pcaMethods","limma","impute","sva","BiocParallel","genefilter","Biobase"))


```
### Step2: install the *dbnorm*
*dbnorm* is freely available from GitHub. 
The package documentation, including  user manual is available within the downloaded R package file. 
If all package dependencies were installed, you will be able to install the *dbnorm*. Users can either manually download the .tar.gz file or
clon the github.


## Usage

TODO: Write usage instructions

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History

TODO: Write history

## Credits

TODO: Write credits

## License

TODO: Write license
]]></content>
  <tabTrigger>readme</tabTrigger>
</snippet>
