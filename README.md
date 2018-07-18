mRNA kinetic model fitting 
====================
This R package is fitting circadian (temproal) profiles of pre-mRNA and mRNA with a kineticc model,
dissecting the contributions of rhythmic transcription and mRNA degradation 
and finally inferring mRNA kineitc parameters, e.g. half-life, rhythmic amplitude and phase of rhythmic degradation.

## TO-DO
Here is the reminder for the improvement:

* Since we are also planning to add the option for "Gaussian noise", we should think how to design the fucntions in such way 
that they can be easily to be adapted to do it. 
And also keep in mind that JW has done it, at least partially, 
in the inital effort (in the folder`origin/`).   
  - The empirical Bayes for the variance in gaussian noise have been implemented in limma pacakge for microarray; however, it can not be  
    directly borrowed, because the limma estimate first the gene-wide variance by fitting a GLM, which is applicable in our case. 
    we need to understand how it works in some detailed steps and to ajust it for our case.
  
* Now the code is designed just for fitting one gene. 
Ideally the code can easily fit all genes in the data in parallel.
Thus the parallization should be taken into consideration now. 

* Not sure we should change S3 class to S4 (more strict in the definition and less error-prone in usage)

* Parameter cleaning in the last step is not clear how to integrate from the origin code

* Headers in all scripts should probably removed or modified

* Connect the general parameter boundaries (modifiable by used) and gene-specific boundaries (refine the boundaries by the gene data) 

## Directory structure
* **[run_modelFitting_forAll.R]** -- the script showing how to run the main function and to specify the parameters
* **data/** -- data example (the read count table used in the PNAS paper) 
* **R/** -- scripts for the main function
* **origin/** -- origin scripts based on which we implement this pacakge and also the scripts used for the PNAS paper (before cleaning)

## Code structure
Here is the structure of inital codes:

### Script showing usage
#### run_modelFitting_forAll.R
**The script mainly show the usage of the package**

- install required R packages in case they are absent
  - `R/configure.R`

- Important data table, time points, pre-mRNA and mRNA lengths...  

- Estimate scaling factors for each sample (e.g. 48 samples in our case), 
  dispersion parameter for each time point (e.g. 12 time points in our case) and store all those parameters in a so-called 
  MDfitDataSet object (S3 class)
  
  `MDfitDataSet(P, M, length.P, length.M, zt, fitType.dispersion = "local")`    

- Run the main function after specifying the parameters
  
  `make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index, debug, outliers.removal, identifiablity.analysis.gamma)`  
  in `R/fitting_degradation_do_stepbystep.R`  
    
- Compare the output with the origin code


### Data processing function
#### MDfitDataSet(P, M, length.P, length.M, zt, fitType.dispersion = "local")
**Generate a MDfitDataSet object (S3 class) in which the data were processed and prepared for the main function**

in `R/preprocess_prepare_for_fitting.R`  

- `print.MDfitDataSet()` -- define simply print function for object MDfitDataSet

- `calculate.SizeFactors.DESeq2()` -- size factor from DESeq2

- `calculate.scaling.factors.DESeq2()` -- size factor * constant close to libary size (here is arbitrarily defined, just a constant) 

- `calculate.dispersions.for.each.time.point.DESeq2()` -- dispersion estiamtion for each time point


### Main function 
#### make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index, debug, outliers.removal, identifiablity.analysis.gamma)

in `R/fitting_degradation_do_stepbystep.R`  
`mds` -- the MDfitDataSet object <br />
`gene.index` -- the index of gene to fit (e.g.  gene.index = 1 (first), 2 (second), 4 (4th))<br />
`debug` -- TRUE or FALSE, print the results for each step<br />
`outlier.remove` -- TRUE or FALSE, remove the outliers or not<br />
`identifiability.anlaysis.gamma` -- TRUE or FLASE, perform identifiability analysis for gamma (degradatio rate) or not <br />
  
- **`R/utilities_generalFunctions.R`** -- general utility function in this script
  - `set.scaling.factors(mds$scaling.factors)` -- set scaling factors  
  - `set.time.points(mds$zt) ` -- set time points  
  - `set.nb.data.param()` -- set number of data points  
  - `norm.RPKM() ` -- convert read counts to RPKM using scaling factors 
  - `convert.nb.reads() ` -- convert the RPKM (the output of model) to read counts
  - `set.bounds.general()` -- set general parameter boundaries (which can be modified by the users)   
    - `set.general.bounds.int()` -- general param boundaries for pre-mRNAs
    - `set.general.bounds.degr.splicing()` -- general param boundaries for mRAN degradation and splicing
  
  - **`kinetic_model.R`**     
    as general utility function becasue called by many steps:   
    parameter optimizaiton in `R/error_function.R`;  
    outlier detection in `R/outliers_detection.R`;  
    identifiability analysis in `R/identifiability_analysis.R`     
    - `compute.s.beta()` -- pre-mRNA concentration from the kinetic model
    - `compute.s.beta.v1()` -- NOT USED pre-mRNA but based on slightly different formular using mean instead of minimun in the model
    - `compute.m.beta()` -- main function to calculte mRNA concentration
      - `integrate.m()` -- intergrate method (not working always)  
        - `Gamma()`
        - `f2integrate()`
      - `simulate.m()` -- by simulation
        - `dmdt()`
      
- extracting data and required parameters for one gene and wrap them into a list called `GeneDataSet`  
  
  
- **`make.fits.with.all.models.for.one.gene(GeneDataSet = GeneDataSet, debug = debug)`**  
  Function to fit the model and optimize parameter  
  in the `R/optimization_params.R`
  
  - `make.fit.spec.model()` -- fit data for specific model (M1-M4)
    - `calculate.error.for.flat.model()` -- for M1
      - `NB.error()` -- calcualte -2loglikelihood for NB (in `R/error_functions.R`)
    - `make.optimization()` -- fit M2-M4
      - fit only pre-mRNA for M2 and M4 
        - `Sampling.Initial.Values.for.fitting.S()` -- gene-specific initial values for pre-mRNA parameters, in `set_bounds_initialize_value.R` 
        - `set.bounds.gene.s()` -- gene-specific parameter boundaries, in `set_bounds_initialize_value.R`
        - `f2min.int()` -- -2loglikelihood for pre-mRNA, in `error_functions.R`
      
      - fit only mRNA for M4 by fixing the pre-mRNA parameters from previous step
        - `Sampling.Initial.Values.for.fitting.M()` -- gene-specific initial values for mRNA degradation, in `set_bounds_initialize_value.R` 
        - `set.bounds.gene.m()` -- gene-specific parameter boundaries, in `set_bounds_initialize_value.R`
        - `f2min.mrna()` -- -2loglikelihood for mRNA, in `error_functions.R`
      
      - fit both pre-mRNA and mRNA
        - `Sampling.Initial.Values.for.fitting.M.S()` -- gene-specific initial values, in `set_bounds_initialize_value.R` 
        - `set.bounds.gene()` -- gene-specific parameter boundaries, in `set_bounds_initialize_value.R`
        - `f2min()` -- -2loglikelihood, in `error_functions.R`
          - `sigmoid.bound.contraint`   
          use a signmoid function to impose smooth constrain for relative amplitude of rhythmic degradation
        
      - compute the Standard Error of estimates using hessian function
    
- **`detect.ouliters.loglike(param.fits.results, GeneDataSet);`**   
  function to detect outliers using the output of optimization function as input arguments  
  in `R/outliers_detection.R`  
    
  
- **`Identifiablity.analysis.gamma.all.models(param.fits.results, GeneDataSet)`**   
  function to analyze the identifiability for parameter gamma  
  which is in `R/identifiability_analysis.R`
    - `set.bounds.gene(M, S, model)` -- set gene-specific parameter boundaries  
    - `Identifiablity.analysis.gamma.each.model() ` -- identifibility analysis for each model 
      - `f2min.profile() ` -- log profile-likelihood function
            
- **`my.model.selection.one.gene.loglike(param.fits.results, method = 'BIC', outlier.m = outlier.m, outlier.s = outlier.s)`**  
  function to do the model selection   
  in `R/model_selection.R`  
    
    
- **`transform.parameter.combinations.cleaning(param.fits.results, res.model.sel, res.nonident.analysis.gamma.all.models)`**  
  function to transform the estimated parameters to more biology-relevant parameters (e.g. degradation rate to half-life)  
  and also converted the parameter combination used in the model fitting  
  and clean the parameter estimation and calculate the error bars  
  in `R/params_transformation_cleaning.R`
    - `transform.parameter.combinations()` -- convert the parameter combinations 
    - `parameter.cleaning()` -- filter and clean the estimated parameter and also test parameter if identifiable 
    

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
 
Running the test for one gene (see [run_modelFitting_for_all.R](#run_modelFitting_for_all.R))


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors


See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

