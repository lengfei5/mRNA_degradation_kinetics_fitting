mRNA kinetic model fitting 
====================
This R package is fitting circadian (temproal) profiles of pre-mRNA and mRNA with a kineticc model,
dissecting the contributions of rhythmic transcription and mRNA degradation 
and finally inferring mRNA kineitc parameters, e.g. half-life, rhythmic amplitude and phase of rhythmic degradation.

## TO-DO
Here is the reminder for the improvement:

* Now the code is designed just for fitting one gene. 
Ideally the code can easily fit all genes in the data in parallel.
Thus the parallization should be taken into consideration now. 

* Since we are also planning to add the option for "Gaussian model", we should think how to design the fucntions in such way 
that they can be easily to be adapted to do it. 
And keep in mind that JW has done it, at least partially, 
in the inital effort (in the folder`origin/`).  

* Not sure we should change S3 class to S4 (more strict in the definition and less error-prone in usage)

* Headers in all scripts should probably removed or modified

* connect the general parameter boundaries (modifiable by used) and gene-specific boundaries (refine the boundaries by the gene data) 

## Directory structure
* **run_modelFitting_forAll.R** (#run_modelFitting_forAll.R)-- the script showing how to run the main function and to specify the parameters
* **data/** -- data example (the read count table used in the PNAS paper) 
* **R/** -- scripts for the main function
* **origin/** -- origin scripts based on which we implement this pacakge and also the scripts used for the PNAS paper (before cleaning)

## Code structure
Here is the structure of inital code:

#### run_modelFitting_forAll.R

The script mainly show how to use the package :

- install required R packages in case they are absent
  - `R/configure.R`

- Important data table, time points, pre-mRNA and mRNA lengths...  

- Estimate scaling factors for each sample (e.g. 48 samples in our case), 
  dispersion parameter for each time point (e.g. 12 time points in our case) and store all those parameters in a so-called 
  MDfitDataSet object (S3 class)
  
  `MDfitDataSet(P, M, length.P, length.M, zt, fitType.dispersion = "local")`    
  in `R/preprocess_prepare_for_fitting.R`  
  - `print.MDfitDataSet` -- define simply print function for object MDfitDataSet
  - `calculate.SizeFactors.DESeq2` -- size factor from DESeq2
  - `calculate.scaling.factors.DESeq2` -- size factor * constant close to libary size (here is arbitrarily defined, just a constant) 
  - `calculate.dispersions.for.each.time.point.DESeq2` -- dispersion estiamtion for each time point

- Run the main function after specifying the parameters
  
  **This is the main function for the pacakge**  
  `make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index, debug, outliers.removal, identifiablity.analysis.gamma)`  
  in `R/fitting_degradation_do_stepbystep.R`  
  `mds` -- the MDfitDataSet object <br />
  `gene.index` -- the index of gene to fit (e.g.  gene.index = 1 (first), 2 (second), 4 (4th))<br />
  `debug` -- TRUE or FALSE, print the results for each step<br />
  `outlier.remove` -- TRUE or FALSE, remove the outliers or not<br />
  `identifiability.anlaysis.gamma` -- TRUE or FLASE, perform identifiability analysis for gamma (degradatio rate) or not <br />
  
  - `R/utilities_generalFunctions.R` -- general utility function in this script
    - `set.scaling.factors(mds$scaling.factors)` -- set scaling factors  
    - `set.time.points(mds$zt) ` -- set time points  
    - `set.nb.data.param()` -- set number of data points  
    - `norm.RPKM() ` -- convert read counts to RPKM using scaling factors 
    - `convert.nb.reads() ` -- convert the RPKM (the output of model) to read counts
    - `set.bounds.general()` -- set general parameter boundaries (which can be modified by the users)   
      - `set.general.bounds.int()` -- general param boundaries for pre-mRNAs
      - `set.general.bounds.degr.splicing()` -- general param boundaries for mRAN degradation and splicing
  
      
  - extracting data and required parameters for one gene and wrap them into a list called `GeneDataSet`  
  
  - `make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, debug = debug)`  
    function to fit the model and optimize parameter  
    which is in the `R/optimization_params.R`
    
    
  - `detect.ouliters.loglike(param.fits.results, GeneDataSet);`   
    function to detect outliers using the output of optimization function as input arguments  
    which is in `R/outliers_detection.R`  
    
  
  - `Identifiablity.analysis.gamma.all.models(param.fits.results, GeneDataSet)`   
    function to analyze the identifiability for parameter gamma  
    which is in `R/identifiability_analysis.R`
      - `set.bounds.gene(M, S, model)` -- set gene-specific parameter boundaries  
      - `Identifiablity.analysis.gamma.each.model() ` -- identifibility analysis for each model 
        - `f2min.profile() ` -- log profile-likelihood function
            
  - `my.model.selection.one.gene.loglike(param.fits.results, method = 'BIC', outlier.m = outlier.m, outlier.s = outlier.s)`  
    function to do the model selection   
    in `R/model_selection.R`  
    
    
  - `transform.parameter.combinations.cleaning(param.fits.results, res.model.sel, res.nonident.analysis.gamma.all.models)`  
    function to transform the estimated parameters to more biology-relevant parameters (e.g. degradation rate to half-life)  
    and also converted the parameter combination used in the model fitting  
    and clean the parameter estimation and calculate the error bars  
    in `R/params_transformation_cleaning.R`
      - `transform.parameter.combinations()` -- convert the parameter combinations 
      - `parameter.cleaning()` -- filter and clean the estimated parameter and also test parameter if identifiable 
    
    
    
- Compare the output with the origin code


## Installation
#### Prerequisites
* R 3.4.1 (currently tested by JW), >= R 3.0.0 should work as well (but to check) 

```
library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
library(numDeriv)
library(Matrix)
library(DESeq2)
library(voom)

```

#### Cloning the git repository
```
cd dir_to_place_repository
git clone https://github.com/lengfei5/mRNA_degradation_kinetics_fitting

```

#### Installing the R package
A step by step series of examples that tell you how to get a development env running
Say what the step will be
```
Give the example
```

## Getting Started
 
 Running the test for one gene (see `run_modelFitting_for_all.R`)

#### install dependencies step 
```
source("R/configure.R")
```

#### import data example from data/
```
data.version = "data_example_readCount"
dataDir = "data/"
load(file = paste0(dataDir, "fitting_degrdation_all_", data.version, ".Rdata"))

ZT.int = grep('.count.premRNA', colnames(T))
ZT.ex = grep('.count.mRNA', colnames(T))
length.int = which(colnames(T) == "length.premRNA")
length.ex = which(colnames(T) == "length.mRNA")
zt = seq(0,94,by = 2)

```
#### creat a MDfitDataSet object (a S3 class) to store the table and estimated scaling factors and dispersion parameters
```
source("R/preprocess_prepare_for_fitting.R")
mds = MDfitDataSet(P = T[, ZT.int], M = T[, ZT.ex], length.P = T[, length.int], length.M = T[, length.ex], zt=zt, fitType.dispersion = "local")
#save(mds, file = "data/MDfitDataSet_example.Rdata")
```
#### Or you can import a saved MDfitDataSet object (a S3 class) 
```
load(file = "data/MDfitDataSet_example.Rdata")
```

#### parameter required to specify and test the function 
```
outliers.removal = TRUE;
debug = TRUE;
absolute.signal = TRUE
parametrization = 'cosine.beta'
identifiablity.analysis.gamma = TRUE

gg = 'Per3'
gene.index = which(T$gene==gg)

source("R/fitting_degradation_do_stepbystep.R")

ptm <- proc.time()
res.fit = make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index = gene.index, debug = debug,
                                                                            outliers.removal = outliers.removal,
                                                                            identifiablity.analysis.gamma = identifiablity.analysis.gamma);
proc.time() - ptm

```


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors


See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

