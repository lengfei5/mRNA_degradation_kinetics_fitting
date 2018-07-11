mRNA kinetic model fitting 
====================
This R package is fitting circadian (temproal) profiles of pre-mRNA and mRNA with a kineticc model,
dissecting the contributions of rhythmic transcription and mRNA degradation 
and finally inferring mRNA kineitc parameters, e.g. half-life, rhythmic amplitude and phase of rhythmic degradation.

## TO-DO

1. now the code is designed just for fitting one gene. 
Ideally the code can easily fit all genes in the data in parallel.
Thus the parallization should be taken into consideration now. 
2. Since we are also planning to add the option for "Gaussian model", we should think how to design the fucntions in such way 
that they can be easily to be adapted to do so 

## Directory structure
* **run_modelFitting_forAll.R** -- the script showing how to run the main function and to specify the parameters
* **data/** -- data example (the read count table used in the PNAS paper) 
* **R/** -- scripts for the main function
* **origin/** -- origin scripts based on which we implement this pacakge and also the scripts used for the PNAS paper (before cleaning)

## Code structure

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

* **Laura Symul ** - *Initial work* - [PurpleBooth]
* **Jingkui Wang ** - *Initial work* - [PurpleBooth]
* **Felix Naef ** - *Initial work* - [PurpleBooth]
* **Pal Westernmark ** - *Initial work* - [PurpleBooth]

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

