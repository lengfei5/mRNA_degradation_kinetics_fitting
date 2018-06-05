##########################################################################
##########################################################################
## Project: fitting temporal profiles of pre-mRNA and mRNA with a kinetic model to infer the rhythmic transcriptional and post-transcriptional regulation
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

## import data example
data.version = "data_example_readCount"
dataDir = "data/"
load(file = paste0(dataDir, "fitting_degrdation_all_", data.version, ".Rdata"))

ZT.int = grep('.count.premRNA', colnames(T))
ZT.ex = grep('.count.mRNA', colnames(T))
zt = seq(0,94,by = 2)
outliers = FALSE;

debug = TRUE;

gg = 'Rorc'
j = which(T$gene==gg)
#source('functions.R')

ptm <- proc.time()
param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = debug,
                                                                            zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = outliers);
proc.time() - ptm

###########################
## compare with origine function
###########################
source('functions_origin.R')
ptm <- proc.time()
param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                            zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
proc.time() - ptm

index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
res.fit = as.numeric(param.fits.results[-index])
names(res.fit) = names(param.fits.results)[-index]


####################
## test a list of examples of circadian genes
####################
Test.circadian.gene.examples = FALSE
if(Test.circadian.gene.examples){
  examples = c('Per2', 'Cry1', 'Per3','Per1', 'Npas2', 'Arntl', 'Clock', 'Cry2', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Rorc','Rora',    
               'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Gsk3a','Csnk1d', 
               'Por', 'Tymp', 'Ppara', 'Hnrnpc', 'Elmo3', 'Camta1', 'Abcb11', 'Tfrc', 'Loxl4', 'Rcan1', 'Cdkn1a', 'Tubb2a', 'Mfsd2', 'Ppard',
               'Nedd4l', 'Cbs', 'Gsta3', 'Eps8l2', 'Insig2', 'Rbm3')
  mm = match(unique(examples), T[,1])
  mm = mm[which(!is.na(mm)==TRUE)]
  T[mm, c(1:3)]
  
  
  source('functions.R')
  
  keep5 = c()
  #for(j in match(controls, T[,1]))
  for(j in mm)
  {
    cat(T[j, 1], '\n');
    
    gg = 'Per2'
    j = which(T$gene==gg)
    source('functions.R')
    #T$mRNA.outlier[j] = paste(outlier.m, sep='', collapse = ';')
    #T$premRNA.outlier[j] = paste(outlier.s, sep='', collapse = ';')
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene.remove.outliers(T = T, gene.index = j, debug = TRUE,
                                                                                zt = zt, i.ex = ZT.ex, i.int = ZT.int, outliers = TRUE);
    proc.time() - ptm
    
    index = match(c('outlier.m', 'outlier.s'), names(param.fits.results))
    res.fit = as.numeric(param.fits.results[-index])
    names(res.fit) = names(param.fits.results)[-index]
    
    missing.data = length(unlist(strsplit(param.fits.results[index[1]], ';'))) + length(unlist(strsplit(param.fits.results[index[2]], ';')))
    test1 = c(res.fit, my.model.selection.one.gene.loglike(res.fit, nb.data = (96-missing.data), method = 'BIC'))
    print(test1);
    m = test1[which(names(test1)=='BIC.best.model')]
    name1 = paste('gamma.m', m, sep='');name2 = paste('gamma.stderr.m', m, sep='')
    eval(parse(text=paste('print(c(half.lfie=log(2)/test1[which(names(test1)==name1)], test1[which(names(test1)==name2)]/test1[which(names(test1)==name1)]))', sep='')));
    
  }
  
  keep5 = data.frame(T[mm[c(1:14)], 1], keep5, stringsAsFactors = FALSE)
  keep5[, c(1, grep('best.model', colnames(keep5)))]
}


