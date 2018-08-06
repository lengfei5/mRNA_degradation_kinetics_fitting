##########################################################################
##########################################################################
## Project: 
## Script purpose: In additional to the global variables, here are some global and general utility functions
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Fri Jun 15 14:53:53 2018
##########################################################################
##########################################################################
###############
## utility functions for model fitting
## these functions are also mainly for the read count table
## they are not used in gaussian model
###############

### extracting the data for one specific gene
ExtractGeneData = function(mds, gene.index)
{
  outlier.m = c();
  outlier.s = c();
  zt = mds$zt; 
  mode = mds$mode;
  if(mds$mode == 'NB')
  {
    R.m = unlist(mds$M[gene.index, ])
    R.s = unlist(mds$P[gene.index, ])
    L.m = T$length.mRNA[gene.index];
    L.s = T$length.premRNA[gene.index];
    
    alpha.m = unlist(mds$dispersions.M[gene.index, ])
    alpha.s = unlist(mds$dispersions.P[gene.index, ])
    
    M = norm.RPKM(R.m, L.m)
    S = norm.RPKM(R.s, L.s)
    
    GeneDataSet = list(R.m = R.m, R.s = R.s, 
                       L.m = L.m, L.s = L.s,
                       Norm.m = M, Norm.s = S, 
                       alpha.m = alpha.m, alpha.s = alpha.s,
                       outlier.m = outlier.m, outlier.s = outlier.s,
                       zt = zt,
                       mode = mode)
  }else{
    
    GeneDataSet = list(Norm.m = unlist(mds$M[gene.index, ]),
                       Norm.s = unlist(mds$P[gene.index, ]), 
                       var.m = unlist(mds$var.M[gene.index, ]),
                       var.s = unlist(mds$var.P[gene.index, ]),
                       outlier.m = c(),
                       outlier.s = c(),
                       zt = mds$zt,
                       mode = mds$mode)
  }
  
  return(GeneDataSet)
}

## time points, data, parameters
set.nb.param = function(absolute.signal=TRUE) ## requires modification 
{
  nb.param <<- c(2,6,5,8);
  #n.param <<- c(0,6,5,8); #the one before 
  #n.param = c(0,5,3,7) from Laura's function
}

set.time.points = function(time.points)
{
  zt <<- time.points;
}

set.nb.data = function(GeneDataSet)
{
  nb.data <<- length(GeneDataSet$Norm.m) + length(GeneDataSet$Norm.s) - length(GeneDataSet$outlier.m) - length(GeneDataSet$outlier.s); 
}

###############################
# import kinetic model and set general parameter boundaries 
###############################
source("R/kinetic_model.R", local = TRUE)

set.general.bounds.int = function(lower.user = NULL, 
                                  upper.user = NULL,
                                  parametrization =c('cosine.beta'), 
                                  absolute.signal = TRUE)
{
  Min.int.max = 2^20; ## intron signal mean 2^9;
  Min.int.min = 2^(-20);
  
  Amp.int.min = 0;
  Amp.int.max = 200;
  
  phase.int.min = 0;
  phase.int.max = 24; 
  
  beta.int.min = 1;
  beta.int.max = 5;
  
  param.synthesis.upper = c(Min.int.max, Amp.int.max, phase.int.max, beta.int.max)
  param.synthesis.lower = c(Min.int.min, Amp.int.min, phase.int.min, beta.int.min)
  names(param.synthesis.upper) = c("Min.int", "Amp.int", "phase.int", "beta.int")
  names(param.synthesis.lower) = names(param.synthesis.upper);
  
  return(list(lower = param.synthesis.lower, upper = param.synthesis.upper))
}

set.general.bounds.degr.splicing = function(lower.user = NULL, 
                                            upper.user = NULL,
                                            parametrization =c('cosine.beta'), 
                                            absolute.signal = TRUE)
{
  ### minimum and maximum of gamma
  gamma.min = log(2)/24
  gamma.max = log(2)/(10/60);  
  
  ### relative amplitudes of gamma (normalized)
  eps.gamma.min = 0;
  eps.gamma.max = 1;
  
  ### phases of gamma (normalized)
  phase.gamma.min = 0;
  phase.gamma.max = 24; 
  
  ### ratio between splicing rate and degradation rate   
  splicing.k.max = 10^5 ## 1s of splicing time;
  splicing.k.min = 10^(-1) ## 30 min of splicing time;
  
  degr.splicing.upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max)
  degr.splicing.lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min)
  
  names(degr.splicing.upper) = c("gamma", "eps.gamma", "phase.gamma", "splicing.k")
  names(degr.splicing.lower) =  names(degr.splicing.upper)
  
  return(list(lower = degr.splicing.lower, upper = degr.splicing.upper))
  
}

set.bounds.general = function(model = 4, lower.user = NULL, upper.user = NULL,
                              parametrization =c('cosine.beta'), absolute.signal = TRUE)
{
  bounds.int = set.general.bounds.int();
  bounds.degr.splicing = set.general.bounds.degr.splicing();
  
  upper = c(bounds.degr.splicing$upper, bounds.int$upper);
  lower = c(bounds.degr.splicing$lower, bounds.int$lower);
  #upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max,  param.synthesis.upper)
  #lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min,  param.synthesis.lower)
  #names(upper) = c("gamma", "eps.gamma", "phase.gamma", "splicing.k", 
  #                         "Min.int", "Amp.int", "phase.int", "beta.int")
  #names(lower) =  names(upper)
  
  if(absolute.signal){
    if(model==1){ upper = upper[c(1, 4, 5)];  lower = lower[c(1, 4, 5)]; }
    if(model==2){ upper = upper[c(1, 4, 5:8)];  lower = lower[c(1, 4, 5:8)];}
    if(model==3){ upper = upper[c(1:5)]; lower = lower[c(1:5)]; }
  }
  
  return(list(lower = lower, upper = upper))
}

###############################
# general functions for NB mode 
###############################
set.scaling.factors = function(sfs)
{
  scaling.factors <<- unlist(sfs);
  #scaling.factors <<- c(39920608, 42250245, 38121270, 45609244, 41511752, 45781196, 43722568, 39638552, 30496638, 30573333, 54950572, 47158379,
  #                      31722765, 39931646, 36317783, 35382708, 47293167, 42408985, 39842283, 40230336, 43691685, 39237518, 51051196, 44778546,
  #                      43858841, 42791401, 42357301, 49782402, 44628140, 44561463, 43485553, 47853067, 43318817, 45055723, 30180984, 46825671,
  #                      43270558, 37496344, 40971385, 45828360, 37065376, 35776330, 45025514, 43026714, 43116633, 35173387, 28538212, 36707156);
  
}

norm.RPKM = function(nb.reads, length)
{
  #set.scaling.factors();
  return(nb.reads/length/scaling.factors*10^9);
}

convert.nb.reads = function(rpkm, length)
{
  #set.scaling.factors();
  return(rpkm*length*scaling.factors/10^9);
}
