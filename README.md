# dbnorm (V-0.2.2)
A package for drift across batches normalization and visualization

![image](https://user-images.githubusercontent.com/37698532/69905457-f93fc080-13b3-11ea-9824-5e7fe8ee183a.png)
*dbnorm* contains R functions which allow visualization and 
removal of technical heterogeneity from large metabolomics dataset. By including advanced statistical tools, the *dbnorm* package allows user to inspect the structure and quality of large metabolomics datasets both in macroscopic and microscopic scales at the sample batch level and metabolic feature level, respectively.
It allows users to efficiently correct drift across batch and to adjust large metabolomics datasets for technical variation which helps improving the estimation of the biological mechanisms underlying disease condition or medical state.
*dbnorm* includes 11 distinct functions for pre-processing of data and estimation of missing values, 
conventional functions for batch effect correction based on statistical models, as well as functions using advanced statistical 
tools to generate several diagnosis plots to help users to choose the statistical model which better fits to their data 
structure. *dbnorm* includes several statistical models such as, ComBat(parametric and non-parametric)-model [PMID:16632515]  from sva package [PMID:22257669] ,that was already in use for metabolomics data normalization, and ber function [DOI:10.1007/s12561-013-9081-1], priorly developed for microarray gene expression data, that we propose here as a new approach for correction of drift across batch in metabolomics datasets. 
## A glimpse into the "dbnorm"
![image](https://user-images.githubusercontent.com/37698532/99926565-00638b00-2cf7-11eb-8101-2436a8d17c40.png)
![link to dbnorm](https://www.nature.com/articles/s41598-021-84824-3)
## Getting started
### Step1: installation
Install package dependencies in CRAN and Bioconductor via codes bellow:

> R installer
```
install.packages(c("ggplot2","parallel","reshape2","plyr",
"knitr","tibble","installr","fs","rmarkdown","processx","backports",
"bootstrap","boot","caret","dplyr","stringr","ggfortify","factoextra","MASS",
"RColorBrewer","RCurl","lattice","data.table","igraph","tidyr","scales",
"e1071","fpc","rlang","glue","digest"))

Call all those package by library:
e.g. library (MASS)
```
> Bioconductor installer
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


BiocManager::install(c("pcaMethods","limma","impute","sva","BiocParallel","genefilter","Biobase","mixOmics","statTarget", "multtest"))
Call all those package by library:
e.g. library (sva)

```
> R install package from archive
```
require(usethis)
require(devtools)
URL: https://cran.r-project.org/src/contrib/Archive/
devtools::install_version("ber", version = "4.0", repos = "http://cran.us.r-project.org")
devtools::install_version("NormalizeMets", version = "0.25", repos = "http://cran.us.r-project.org")
devtools::install_version("metabolomics", version = "0.1.4", repos = "http://cran.us.r-project.org")

library(ber)
library(NormalizeMets)
library(Rcpp)
library(metabolomics)
```
### Step2: install the *dbnorm*
*dbnorm* is freely available from GitHub. 
The package documentation, including  user manual is available within the downloaded R package file. 
If all package dependencies were installed, you will be able to install the *dbnorm*. Users can either manually download the **tar.gz** file or clon the GitHub.

#### Manual downloading
```
cd ~/Downloads
R CMD INSTALL dbnorm tar.gz
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
R CMD INSTALL dbnorm tar.gz
```
> Rstudio users
```
Download the tar.gz file from GitHun then in Rstudio from Tools icon, select Install.packages and 
then Browse to the tar.gz file.
```
## Instructions
In this section, we briefly introduce and explain the functions implemented in the *dbnorm* with the expected outcome. Data to be uploaded must be normalized and scaled on the log2 to account for the high abundance features (variables) by which technical heterogeneity might be overlooked. The input data must be in **.csv** format with the independent experiments in the rows and the features (variables) in the columns, with the `batch` levels considered in the first column.

##### *Upload your data ana call the package;*

- Example:

> data<-read.csv(" path/to/directory/folde/mydata.csv",sep = ",", header = T, row.names = 1)

> library(dbnorm)

### Functions

> - emvd

This function allows you to estimate missing values (Zero or/and NA values) by the lowest detected value in the entire experiment.

```
df<-data[-1] #keep data matrix by removing batch level in the first column
>emvd [df]
f<- emvd[df] # save the imputed data matrix
```
> - emvf

This function allows you to estimate missing values (Zero or/and NA values) for each feature (variable) by the lowest value detected for the corresponding feature (variable), applied on the column. 

```
df<-data[-1] #keep data matrix by removing batch level in the first column
> emvf [df]
f <- emvf[df] # save the imputed data matrix
```
> - Visodbnorm ; 
*Visualization of drift across batch normalization*

This function performs batch effect adjustment via three statistical models implemented in the *dbnorm*, namely two-stage procedure as described by Giordan (2013)[DOI:10.1007/s12561-013-9081-1] and/or empirical Bayes methods in two setting of parametric and non-parametric as described by Johnson et al.(2007) [PMID: 16632515] and in sva package by Leek et al.(2012)[PMID: 22257669]. Meanwhile, the graphical inferences in the context of unsupervised learning algorithms create visual inspection to inform users about the spatial separation of the sample sets analyzed in the different analytical runs alongside the distribution of features (variables) in the raw and treated datasets. This function is suggested for less than 2000 features (variables).

- Value

Graphical check such as *PCA* plot and *Scree* plot compiled into a **PDF** (saved in the working directory) and three **.csv** files (saved in a folder, intiate with *Rtmpe*, in Users's Temporary directory: "C:\Users\ “%USERNAME”\AppData\Local\Tem") for adjusted data based on the models implemented in the package.The *RLA* plots are visualized in the **Viewer** panel in the **rstudio** console. 
```

Visodbnorm(data)
```
>- dbnormSCORE ; 
*Adjusted coefficient of determination for a data normalized for across batch signal drift*

This function gives a quick notification about the performance of the statistical models, two-stage procedure [DOI:10.1007/s12561-013-9081-1] and/or empirical Bayes methods in two setting of parametric and non-parametric as described in [PMID: 16632515] and by sva package [PMID:22257669], implemented in the dbnorm package, in accommodating technical variability. Subsequently, the  adjusted coefficient of determination or Adjusted R-Squared is calculated for each variable estimated in a regression model for its dependency to the batch level in the raw data and treated data via either of those models. Immediately, the performance of applied models are presented by a score calculated based on the maximum variability explained by the batch level, notify the consistency of model performance for all detected features (variables), facilitating quick comparison of the models for selecting one of those models, which is more appropriate to the data structure. This function is suggested for less than 2000 features (variables) for better computational speed.

- Value

Graphical check such as *Correlation* plot and *Score* plot compiled into a **PDF** file (saved in the working directory) and **.csv** files (saved in a folder, intiate with *Rtmpe*, in Users's Temporary directory: "C:\Users\ “%USERNAME”\AppData\Local\Tem") in the two vector data matrix of Adjusted R-squared for each model and a *Table* of score for the maximum Adjusted R-squared detected for the applied models.

```
dbnormSCORE (data)
```


>- ProfPlotraw
>- ProfPlotBer
>- ProfPlotBagging
>- ProfPlotComPara
>- ProfPlotComPara
>- ProfPlotComNPara;
*Visualization of analytical heterogeneity on the profile of features (variables) in raw data and after correction via ber-, ber-bagging, parametric ComBat and non-parametric ComBat*

These functions allow users to adjust the data for batch effect using either of models implemented in the package described earlier, and inform about the presence of across batch signal drift or batch effect in the raw and treated data represented via the shifted probability density function (PDF) plots of the features (variables) detected in an experiment.

- Value

Graphical check such as the plots compiled into a **PDF** file (saved in the working directory) and a **.csv** file (saved in a folder, intiate with *Rtmpe*, in Users's Temporary directory: "C:\Users\ “%USERNAME”\AppData\Local\Tem") of corrected dataset via either of applied function.

```
profplotraw (data)
ProfPlotBer (data)
ProfPlotBagging (data)
ProfPlotComPara (data)
ProfPlotComPara (data)
ProfPlotComNPara (data)
```

> - dbnormBer
> - dbnormBagging
> - dbnormPcom
> - dbnormNPcom;
*Data normalization for across batches signal drift using either of  ber-,ber-bagging, parametric ComBat- and non parametric ComBat- models and unsupervised clustering and regression analysis of corrected data*


To increase computational processing of big data, in these functions, statistical models and  graphical checks implemented in “Visodborm” decomposed in to several separated functions each of these performing a unique batch effect correction with respective result and graphical checks.

- Value

Graphical check such as *PCA* plot, *Scree* plot and *Correlation* plot compiled into a **PDF** (saved in the working directory) and the **.csv** (saved in a folder, intiate with *Rtmpe*, in Users's Temporary directory: "C:\Users\ “%USERNAME”\AppData\Local\Tem") for corrected dataset based on either of applied model and the two column matrix of Adjusted-R squared. The *RLA* plots are visualized in the **Viewer** panel in the **rstudio** console. 

```
dbnormBer(data)
dbnormBagging (data)
dbnormPcom(data)
dbnormNPcom(data)

```

> - hclustdbnorm;
*Hierarchical clustering analysis of original data and corrected data for batch effect*
 
This function allows users to evaluate dissimilarity between identical samples (quality control replicates or analytical replicates) analyzed in different batches, prior and after correction using, ber, ber_bagging and parametric and non-parametric ComBat . Pearson distance and average method for clustering were considered.Bagging model is performed using partial bagging with n=150 bootstrap samples

```
hclustdbnorm (data)
```
# License
Distributed under the GPL license. See LICENSE for details.

