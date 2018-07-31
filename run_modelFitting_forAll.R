##########################################################################
##########################################################################
## Project: fitting temporal profiles of pre-mRNA and mRNA with a kinetic model to infer the rhythmic transcriptional 
## and post-transcriptional regulation
## Script purpose: this script is the function of the highest order for R package
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Fri Jun  1 12:01:00 2018
##########################################################################
##########################################################################

######################################
######################################
## Section: test fitting for read count table
######################################
######################################
rm(list=ls())

####################
## import data example and create an object
## prepare the table, geneNames, geneLengths, sizeFactors, dispersion estiamtion and variance estimation 
####################
dataDir = "data/"
load(file = paste0(dataDir, "fitting_degradation_all_data_example_readCount_rpkm.Rdata"))

TEST.readCount.NB = FALSE

source("R/preprocess_prepare_for_fitting.R")

if(TEST.readCount.NB){
  # zt = seq(0,94,by = 2)
  #ZT.int = grep('.count.premRNA', colnames(T))
  #ZT.ex = grep('.count.mRNA', colnames(T))
  #length.int = which(colnames(T) == "length.premRNA")
  #length.ex = which(colnames(T) == "length.mRNA")
  
  ####################
  ## creat a MDfitDataSet object (a S3 class)
  ####################
  #mds = MDfitDataSet(P = T[, ZT.int], M = T[, ZT.ex], length.P = T[, length.int], length.M = T[, length.ex], zt=zt,
  #                   mode = "NB", fitType.dispersion = "local")
  #save(mds, file = "data/MDfitDataSet_example.Rdata")
  load(file = "data/MDfitDataSet_example.Rdata")
  
}else{
  #zt = seq(0,94,by = 2)
  #ZT.int = intersect(grep('.rpkm.premRNA', colnames(T)), grep("ZT", colnames(T)))
  #ZT.ex = intersect(grep('.rpkm.mRNA', colnames(T)), grep("ZT", colnames(T)))
  #mds = MDfitDataSet(P = T[, ZT.int], M = T[, ZT.ex], zt=zt, mode = "logNormal")
  #save(mds, file = "data/MDfitDataSet_example_logNormal.Rdata")
  
  load(file = "data/MDfitDataSet_example_logNormal.Rdata")
}

####################
## parameter required to specify
####################
outliers.removal = FALSE;
debug = TRUE;
identifiablity.analysis.gamma = FALSE

gg = 'Per3'
gene.index = which(T$gene==gg)

####################
## test the current functions 
####################
source("R/fitting_degradation_do_stepbystep.R")

ptm <- proc.time()
res.fit = make.fits.with.all.models.for.one.gene.remove.outliers(mds, gene.index = gene.index, debug = debug,
                                                                            outliers.removal = outliers.removal,
                                                                            identifiablity.analysis.gamma = identifiablity.analysis.gamma);
proc.time() - ptm

###########################
## compare with origine function
###########################
load(file = paste0(dataDir, "fitting_results_for_examples_2compare_with_new_functions.Rdata"))
keep2compare[which(keep2compare$gene == gg), ]

rm(list = lsf.str())
source('origin/functions_origin.R')
ptm <- proc.time()
param.fits.results.v0 = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = gene.index, debug = TRUE,
                                                                            zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = outliers.removal,
                                                                            Identifiablity.Analysis.by.Profile.Likelihood.gamma 
                                                                            = Identifiablity.Analysis);
proc.time() - ptm

index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
res.fit = as.numeric(param.fits.results[-index])
names(res.fit) = names(param.fits.results)[-index]

Check.dispersion.parameters = FALSE
if(Check.dispersion.parameters){
  load(file = "archives/fitting_degradation_all_data_example_readCount.Rdata")
}

####################
## test a list of examples of circadian genes
####################
Test.circadian.gene.examples = FALSE
if(Test.circadian.gene.examples){
  examples = c("Rorc", "Cp",  "Fus", "Per3", "Cbs", "Per1", "Per3", "Nr1d1", "Tef", "Dbp", "Hlf", "Bhlhe40", "Nfil3", "Per2", "Clock", "Wee1")
  circ.examples = c('Per2', 'Cry1', 'Per3','Per1', 'Npas2', 'Arntl', 'Clock', 'Cry2', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Rorc','Rora',    
               'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Gsk3a','Csnk1d', 
               'Por', 'Tymp', 'Ppara', 'Hnrnpc', 'Elmo3', 'Camta1', 'Abcb11', 'Tfrc', 'Loxl4', 'Rcan1', 'Cdkn1a', 'Tubb2a', 'Mfsd2', 'Ppard',
               'Nedd4l', 'Cbs', 'Gsta3', 'Eps8l2', 'Insig2', 'Rbm3')
  
  mm = match(unique(examples), T[,1])
  mm = mm[which(!is.na(mm)==TRUE)]
  T[mm, c(1:3)]
  
  #source('functions.R')
  source('origin/functions_origin.R')
  
  keep2compare = c()
  #for(j in match(controls, T[,1]))
  for(j in mm)
  {
    cat(T[j, 1], '\n');
    
    # gg = 'Per2'
    #j = which(T$gene==gg)
    #source('functions.R')
    #T$mRNA.outlier[j] = paste(outlier.m, sep='', collapse = ';')
    #T$premRNA.outlier[j] = paste(outlier.s, sep='', collapse = ';')
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
    proc.time() - ptm
    
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    
    outers = c(unlist(strsplit(param.fits.results[index[1]], ';')), unlist(strsplit(param.fits.results[index[2]], ';')))
    outers = outers[which(outers != "NA")]
    nb.outliers = length(outers)
    
    test1 = c(res.fit, my.model.selection.one.gene.loglike(res.fit, nb.data = (96-nb.outliers), method = 'BIC'))
    print(test1);
    m = test1[which(names(test1)=='BIC.best.model')]
    name1 = paste('gamma.m', m, sep='');
    name2 = paste('gamma.stderr.m', m, sep='')
    eval(parse(text=paste('print(c(half.lfie=log(2)/test1[which(names(test1)==name1)], 
                          test1[which(names(test1)==name2)]/test1[which(names(test1)==name1)]))',
                          sep='')));
   
    keep2compare = rbind(keep2compare, test1) 
  }
  
  keep2compare = data.frame(T[mm[1:7], 1],  keep2compare, stringsAsFactors = FALSE)
  colnames(keep2compare)[1] = 'gene'
  
  save(keep2compare, file = paste0(dataDir, "fitting_results_for_examples_2compare_with_new_functions.Rdata"))
  
}


