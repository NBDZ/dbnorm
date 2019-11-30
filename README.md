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
If all package dependencies were installed, you will be able to install the *dbnorm*. Users can either manually download the .tar.gz file or clon the github.

#### Manual downloading
```
cd ~/Downloads
R CMD INSTALL dbnorm.tar.gz
```
> Rstudio users
```
First, download the tar.gz file from GitHun then in Rstudio from Tools icon, select Install.packages and 
then Browse to the tar.gz file.
```
#### clon the github
```
git clone 
git@github.com:NBDZ/dbnorm.git
R CMD build Mdbnorm
R CMD INSTALL dbnorm.tar.gz
```
> Rstudio users
```
Download the tar.gz file from GitHun then in Rstudio from Tools icon, select Install.packages and 
then Browse to the tar.gz file.
```
## Functions
In this section, we breifly introduce the functions implemented in the *dbnorm* and explaine the expected outcome. 
##### Upload your data
data<-read.csv(" path/to/directory/folde/mydata.csv",sep = ",", header = T, row.names = 1)
> data must be normalized and scaled on the log2 to account for the high abun dance features (variables) by which technical heterogeneity might be overlooked. The input data must be in csv format with the independent exper-
iments in the rows and the features (variables) in the columns. The **batch**
levels must be frames in the first column. 


>emvd

```
Missing values in the forms of Zero values and NA values are estimated by the lowest detected value in the entire experiment






TODO: Write license
]]></content>
  <tabTrigger>readme</tabTrigger>
</snippet>
