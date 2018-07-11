mRNA kinetic model fitting 
====================
This R package is fitting circadian (temproal) profiles of pre-mRNA and mRNA with a kineticc model,
dissecting the contributions of rhythmic transcription and mRNA degradation 
and finally inferring mRNA kineitc parameters, e.g. half-life, rhythmic amplitude and phase of rhythmic degradation.

## TODO

In general, I found DESeq2 R package is good example of implementing complicated functions in R, which we can inspire very much from it.

1. now the code is design just for fitting one gene. Ideally the code can easily fit all genes in the data in parallel.
Thus the parallization should be taken into consideration now. 
2. Since we are also planning to add the option for "Gaussian model", we should think how to design the fucntions in such way 
that they can be easily to be adapted to do so 

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
And repeat

```
until finished
```
End with an example of getting some data out of the system or using it for a little demo

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.


#### Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
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

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

