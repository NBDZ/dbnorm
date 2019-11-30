# dbnorm  
- a package for drift across batch normalization and visualization

![image](https://user-images.githubusercontent.com/37698532/69902291-11034e80-138c-11ea-8cad-d3b8dadd1493.png)
*dbnorm* contains R functions which allow visualization and 
removal of technical heterogeneity from large metabolomics dataset. By including advanced statistical tools, the *dbnorm* package allows user to inspect the structure and quality of large metabolomics datasets both in macroscopic and microscopic scales at the sample batch level and metabolic feature level, respectively.
It allows users to efficiently correct drift across batch and to adjust large metabolomics datasets for technical variation which helps improving the estimation of the biological mechanisms underlying disease condition or medical state.
*dbnorm* includes 11 distinct functions for pre-processing of data and estimation of missing values, 
conventional functions for batch effect correction based on statistical models, as well as functions using advanced statistical 
tools to generate several diagnosis plots to help users to choose the statistical model which better fits to their data 
structure. *dbnorm* implements several statistical models, including ComBat(parametric and non-parametric)-model[PMID: 16632515]  from *svs* package [PMID: 22257669]
, that was already in use for metabolomics data normalization, and ber function [DOI 10.1007/s12561-013-9081-1], priorly developed for microarray gene expression data,
, that we propose here as a new approach for correction of drift across batch in metabolomics datasets. 
## A glimpse into the "dbnorm"
![dbnorm](https://user-images.githubusercontent.com/37698532/69902657-8e30c280-1390-11ea-9770-811e4a18eb85.jpg)
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

Call all those package by library:
e.g. library (ber)
```
> Bioconductor installer
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("pcaMethods","limma","impute","sva","BiocParallel","genefilter","Biobase"))
Call all those package by library:
e.g. library (sva)

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
##### *Upload your data*
###### Example:
> data<-read.csv(" path/to/directory/folde/mydata.csv",sep = ",", header = T, row.names = 1)
> library(dbnorm)

data must be normalized and scaled on the log2 to account for the high abun dance features (variables) by which technical heterogeneity might be overlooked. The input data must be in csv format with the independent experiments in the rows and the features (variables) in the columns. The `batch` levels must be frames in the first column. 

>**emvd**

This function allows you to estimate missing values in the forms of Zero values and/or NA values by the lowest detected value in the entire experiment.

```
df<-data[-1] #keep data matrix by removing batch level in the first column
>emvd [df]
f<- emvd[df] # save the imputed data matrix
```
>**emvf**

This function allows you to estimate missing value for each feature (variable) by the lowest
value detected for the corresponding feature (variable). This function is applicable in all sort of high-throughput experiment. Both Zero values and NA values are imputed.

```
df<-data[-1] #keep data matrix by removing batch level in the first column
> emvf [df]
f <- emvf[df] # save the imputed data matrix
```
>**Visdbnorm**; 
*Visualization of drift across batch normalization*

This function performs batch effect adjustment via three statistical models implemented in the
dbnorm, namely two-stage procedure as described by Giordan (2013)[DOI 10.1007/s12561-013-9081-1]
and/or empirical Bayes methods in two setting of parametric and non-parametric as described by
Johnson et al.(2007) [PMID: 16632515] and in sva package by Leek et al.(2012)[PMID: 22257669]. Meanwhile, the graphical inferences
in the context of unsupervised learning algorithms create visual inspection to inform users about
the spatial separation of the sample sets analyzed in the different analytical runs alongside the distribution
of variables (features) in the raw and treated datasets. This function is suggested for less
than 2k variables (features).

- Value

Graphical check such as **PCA** plot, **RLA** plot and **Scree** plot together with data frame of corrected data in *csv* format 
```
ff<- data.frame(data[1],f]# frame the batch level with imputed matrix
Visdbnorm(ff)
```
>**ACDdbnorm**; 
*Adjusted coefficient of determination for a data normalized for across batch signal drift*

This function gives a quick notification about the performance of the statistical models implemented
in the dbnorm package such as two-stage procedure [DOI 10.1007/s12561-013-9081-1] and/or empirical Bayes methods in two setting of
parametric and non-parametric as described in [PMID: 16632515] and by sva package [PMID: 22257669]. Then adjusted coefficient of determination or Adjusted R-Squared is calculated for each variable
estimated in a regression model for its dependency to the batch level in the raw data and treated
data via either of those models. Immediately, the performance of applied models are presented by
a score calculated based on the maximum variability explained
by the batch level, notify the consistency of model performance for all detected variables (features), facilitating quick
comparison of the models for selecting one of those models, which is more appropriate to the data
structure. This function is suggested for less than 2k features.

- Value

Graphical check such as **Correlation** plot and **Score** plot together with the two vector data matrix of Adjusted R-squared for each model and a *Table* of score for the maximum Adjusted R-squared for each model

```
ff<- data.frame(data[1],f]# frame the batch level with imputed matrix
ACDdbnorm(ff)
```
>  **profplotraw**, **profplotber**, **profplotpcomr** and **profplotnpcomr**; 
*Visualization of analytical heterogeneity on the profile of variables
(features) in raw, ber- parametric ComBat and non-parametric ComBat corrected data*

These functions allow users to adjust the data for batch effect base on either of models implemented in the package decribed earlier, and inform about the presence of across batch signal drift or batch effect in the raw and treated data represented via the shifted probability density function plots (pdf plots) of variables (features) detected in an experiment.

- Value

Graphical check such as **pdf** plot together with the data frame of corrected data in *csv* format

```
ff<- data.frame(data[1],f]# frame the batch level with imputed matrix
profplotraw(ff)
profplotber(ff)
profplotpcomr(ff)
profplotnpcomr(ff)
```
> **dbnormBer**,**dbnormPcom** and **dbnormNPcom**; *unsupervised clustering and regression analysis of data corrected via ber-, parametric ComBat- and non parametric ComBat- model*

These functions allow users to adjust the data for across batch signal drift or batch effect using either of
models implemented in the package decribed earlier. These functions includ advanced statistical tools to inspect the structure and quality of high throughput experiment both in macroscopic and microscopic scale at the sample batch level and metabolic
feature level, respectively. Notably, using this function users perform unsupervised clustering analysis
of the raw data and the treated dataset. In parallel, Adjusted-R squared value for each variable
(feature) estimated by regression model is calculated, which demonstrate the dependency of variable
(feature) to the batch level in either of those datasets. In addition, for quick notification about
the performance of the applied model we considered a score, which is calculated based on the maximum variability. This score notifies  the consistency of model performance for the detected variables (features).

- Value

Graphical check such as **PCA** plot, **RLA** plot, **Scree** plot and **Correlation** plot together with data frame of corrected data in *csv* format and two column matrix of Adjusted-R squared.

```
ff<- data.frame(data[1],f]# frame the batch level with imputed matrix
dbnormBer(ff)
dbnormPcom(ff)
dbnormNPcom(ff)

```

TODO: Write license
]]></content>
  <tabTrigger>readme</tabTrigger>
</snippet>
