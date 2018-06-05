######################################
######################################
## Section: should be moved into the function scripts
######################################
######################################
### install some packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# usage
packages <- c("fdrtool", "circular", "preprocessCore", "gtools", "biomaRt", "numDeriv", "Matrix")
ipak(packages)

library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
#library(plotrix)
library(numDeriv)
library('Matrix')
#library('matrixcalc')
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

#####################
#### We work only on the absolute value of RNA-seq and cosine.beta parametrization
#####################

###### splicing demonstration and intron selection
selection = function(vect, cutoff=-4)
{
    length(which(vect>cutoff))>=6
}

packing = function(vect)
{
    return(paste(vect, sep='', collapse=','))
}
unpacking = function(vect)
{
    return(as.numeric(unlist(strsplit(as.character(vect), ','))))
}
unpacking.plots = function(vect)
{
    gene.name = vect[1]
    gene.length = vect[2]
    nb.introns = as.numeric(vect[3])
    #strand = as.numeric(vect[29])
    
    rainbow = rainbow(24,s = 0.85, v = 0.85)
    for(kk in c(5, 8, 11, 14))
    {
        test = unpacking(vect[kk]);
        #if(strand==-1) test = test[c(nb.introns:1)]
        if(kk ==5)
        {
            xlim=c(1, nb.introns)
            ylim = c(-6,6)
            plot(c(1:nb.introns), test, type='l', xlab='Index of Introns', ylab='RPKM expression', xlim=xlim, ylim=ylim, main=paste(gene.name, ', gene.length:', gene.length, sep=''), col=rainbow[kk-4], lty = (kk-1)%/%24+1)
            legend(x = 0.8*xlim[2],y = ylim[2]*(kk-4)/24, legend = paste('ZT',2*(kk-5)), col =  rainbow[kk-4], lty = (kk-1)%/%12+1, bty = 'n' )
            abline(h=-6, col='darkgray', lwd=2.0)
            
        }else{
            
            points(c(1:nb.introns), test, type='l', col=rainbow[kk-4], lty = (kk-1)%/%24+1)
            legend(x = 0.8*xlim[2],y = ylim[2]*(kk-4)/24, legend = paste('ZT',2*(kk-5)), col =  rainbow[kk-4], lty = (kk-1)%/%12+1, bty = 'n' )
        }
    }
    #return(as.numeric(unlist(strsplit(as.character(vect), ','))))
}


######## FUNCTIONS FOR SIMULATED DATA and for the optimization and model selection
set.scaling.factors = function()
{
  #ptm = proc.time()
  #load(file='Scaling_factors_48_samples.Rdata')
  #proc.time() - ptm;
  #xx = scaling.factors;
  #scaling.factors <<- xx;
  #proc.time() - ptm;
  
  #ptm = proc.time()
  scaling.factors <<- c(39920608, 42250245, 38121270, 45609244, 41511752, 45781196, 43722568, 39638552, 30496638, 30573333, 54950572, 47158379,
                        31722765, 39931646, 36317783, 35382708, 47293167, 42408985, 39842283, 40230336, 43691685, 39237518, 51051196, 44778546,
                        43858841, 42791401, 42357301, 49782402, 44628140, 44561463, 43485553, 47853067, 43318817, 45055723, 30180984, 46825671,
                        43270558, 37496344, 40971385, 45828360, 37065376, 35776330, 45025514, 43026714, 43116633, 35173387, 28538212, 36707156);
  #proc.time() - ptm;
  
}

mean.substract = function(vect)
{
	return(vect-mean(vect[which(!is.na(vect))]))
}

index.outliers = function(data.xx)
{
	c = 1.5
	#data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
	Q1 = quantile(data.xx, 0.25,type=5)
	Q3 = quantile(data.xx, 0.75, type=5)
	IQD = Q3 - Q1
	lower = Q1 - c*IQD
	upper = Q3 + c*IQD
	index = which(data.xx<lower|data.xx>upper)
	#boxplot(data.xx);abline(h=Q1);abline(h=Q3);
}

variance.mean.counts = function(x, normalization=TRUE)
{
  if(normalization)
  {
    load(file='size_factors_48_samples.Rdata')
    #set.scaling.factors();
    x = x/size.factors;
  }
  # x = exon
  # var.mean = matrix(NA, ncol=3, nrow=(nrow(x)*12))
  var.mean = c()
  for(n in 1:12)
  {
    index = c(n, (n+12), (n+24), (n+36))
    mean = apply(x[, index], 1, mean)
    var = apply(x[, index], 1, var)
    var.mean = rbind(var.mean, cbind(rep(n, nrow(x)), mean, var))
  }
  #sd = mean(apply(xx, 2, sd))
  return(var.mean)
}

variance.mean.rpkm.gene = function(x)
{
  # var.mean = matrix(NA, ncol=3, nrow=(nrow(x)*12))
  x = as.numeric(x);
  var.mean = c()
  for(n in 1:12)
  {
    index = c(n, (n+12), (n+24), (n+36))
    #mean = apply(x[index], 1, mean)
    #var = apply(x[, index], 1, var)
    var.mean = rbind(var.mean, c(mean(x[index]), var(x[index])))
  }
  #sd = mean(apply(xx, 2, sd))
  return(var.mean)
}

variance.estimate.replicates = function(x, nb.replicates=4, log=TRUE)
{
  if(!log)
  {
    x = log(x)
  }
  x = as.numeric(x)
  #xx = matrix(x, ncol=4)
  if(nb.replicates==4) xx = rbind(x[c(1:12)], x[c(13:24)], x[25:36], x[37:48]);
  if(nb.replicates==2) xx = rbind(x[c(1:12)], x[c(13:24)]);
  sd = mean(apply(xx, 2, sd))
  
  return(sd)
}

index.outliers.loglike = function(data.xx, c=1.5)
{
  #c = 3
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
  #data.xx = data.xx[which(!is.na(data.xx)==TRUE)]
  Q1 = quantile(data.xx, 0.25,type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx>upper)
  #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
  return(index)
}


plot.dispersion.examples = function(T, index)
{
  ##index = mm; T=F;
  pdfname = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/dispersion_examples_clock_genes_logscale_clock_genes.pdf';
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  
  pdf(pdfname, width=6, height=6)
  for(n in index)
  {
    # n = index[8]
    vm.ex = variance.mean.counts(T[n, ZT.ex], normalization = TRUE)
    vm.int = variance.mean.counts(T[n, ZT.int], normalization = TRUE)
    
    lims = range(vm.ex[,2], vm.int[,2])
    lims = c(0.1, lims[2])
    xx = lseq(lims[1], lims[2], 100)
    par(mfrow=c(1,1))
    
    plot(vm.ex[,2], vm.ex[,3], log='xy', xlim=lims, ylim=c(0.01, range(vm.ex[,3], vm.int[,3])[2]), col='blue',
         main=paste(T$gene[n], '  alpha  mRNA ', signif(T$alpha.mRNA[n], d=2), '  ; premRNA ', signif(T$alpha.premRNA[n], d=2),  sep=''), 
         xlab='mean of normalized counts', ylab='variance of normalized counts')
    points(vm.int[,2], vm.int[,3], col='green')
    points(xx, (xx+xx^2*T$alpha.mRNA[n]), type='l', col='blue', lwd=2.0)
    points(xx, (xx+xx^2*T$alpha.premRNA[n]), type='l', col='green', lwd=2.0)
    #points(xx, (xx+xx^2*0.03), type='l', col='darkgray', lwd=2.0, lty=2)
    abline(0, 1, col='red', lwd=2.0, lty=2)
    
  }
  dev.off()
  
}

make.fitting.premrna = function(T = T, i = 1, debug = FALSE, parametrization = c('cosine.beta'),zt = seq(0,46,by = 2),i.ex = ZT.ex, i.int = ZT.int, absolute.signal=TRUE)
{
	# parametrization = 'cosine'; i.int = ZT.int;zt = seq(0,46,by = 2);i = n;
	S = as.numeric(T[i, i.int])
	#source('model_modify_new_RNA_seq/functions.R')
	parametrization = 'cosine.beta';
	#x = seq(0, 46,by=2)
	#fc = 10
	#beta = 1
	#phase = 15
	#mean = 10
	#yy = 10*1/(1+(fc-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)*(1+(fc-1)*((1+cos(2*pi/24*(x-phase)))/2)^beta)
	#plot(x, yy, type='l', col=1)
	#S = yy
	fitting.factor = 3;
	Nfit.S = 4*fitting.factor;
	index = i
	
	### model 1 flat
	s = rep(exp(mean(log(S))),length(S))
	
	err = sum((log(S)-log(s))^2)
	param.fit = err;
	names(param.fit) = 'error.m1';
	
	#cat(err, '\n');
	#gamma.init = c(rep(log(2)/lseq(log(2)/lower[1],log(2)/upper[1], length = 5), Nfit%/%5), rep(log(2)/5,Nfit%%5))
	if(!is.null(T$rel.amp.premRNA[index]))
	{
		rel.ampl.int = T$rel.amp.premRNA[index]; 
		if(rel.ampl.int>=1.0|rel.ampl.int<=0) rel.ampl.int = 0.5;
		phase.int = T$phase.premRNA[index]
	}else{
		rel.ampl.int = (max(S)-min(S))/2;
		phase.int = zt[which.max(S)]
	}
	mean.int.init = rep(mean(S), Nfit.S)
	
	fold.change.init = rep(min((1+ rel.ampl.int)/(1-rel.ampl.int),1000),Nfit.S);
	fold.change.init = fold.change.init*sample(seq(1/max(fold.change.init),2-(1/max(fold.change.init)),len = 1000),Nfit.S)
	fold.change.init = c(fold.change.init[1:ceiling(Nfit.S/2)], lseq(1.5, 50, length=(Nfit.S-ceiling(Nfit.S/2))))
	
	phase.init = (rep(phase.int,Nfit.S)+rnorm(Nfit.S,sd = 6))%%24
	
	beta.min = 1.0;
	beta.max = 5;
	beta.init = lseq(beta.min,beta.max, length = Nfit.S)
	PAR.INIT.S = cbind(mean.int.init,fold.change.init, phase.init, beta.init)
	
	## fitting S 
	errors.fit.s = rep(NA, Nfit.S)
	min.mean = min(S);
	max.mean = max(S)
	for(fit.nb.s in 1:Nfit.S)
	{
		par.init.s = PAR.INIT.S[fit.nb.s,]
		#cat(par.init.s, '\n')
		opt.s = optim(par.init.s, f2min.int, S=S, parametrization=parametrization, absolute.signal = absolute.signal, zt = zt, method = 'L-BFGS-B', lower = c(min.mean, 1,0,1), upper = c(max.mean, 1000, 24, 5))
		res.fit.s = opt.s$par
		errors.fit.s[fit.nb.s] = opt.s$value
		#cat(opt.s$value, '\n');
		#cat(opt.s$par, '\n');
		#cat('..........\n')
		eval(parse(text = paste('res.fit.s.', fit.nb.s, ' = res.fit.s', sep = '')))
		
	}
	
	## choose the best-fitting parameters for S and select the best model with BIC
	imin.s = which.min(errors.fit.s); 
	eval(parse(text = paste('res.fit.s = res.fit.s.', imin.s, sep = '')))
	
	param.fit2 = c(errors.fit.s[imin.s], res.fit.s);
	names(param.fit2) = paste(c('error', colnames(PAR.INIT.S)),'.m2',sep = '')
	
	param.fits.results = c(param.fit, param.fit2)
	
	set.nb.param(absolute.signal=absolute.signal);
	index = match(c('error.m1', 'error.m2'), names(param.fits.results))
	error.m1 = param.fits.results[index[1]]
	error.m2 = param.fits.results[index[2]]
	
	## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
	BIC.m1 = log(48)*n.param[1] + 48*log(error.m1)
	BIC.m2 = log(48)*n.param[2] + 48*log(error.m2)
		
	BIC = c(BIC.m1, BIC.m2)
	BIC = c(BIC, which.min(BIC))
	names(BIC) = c('BIC.m1', 'BIC.m2', 'BIC.best.model')
	
	return(c(param.fits.results, BIC))
	#s = compute.s.beta(x, res.fit.s[2], res.fit.s[3], res.fit.s[4])*res.fit.s[1]
	#ss = compute.s.beta(x, fc, phase, beta)*mean
	#plot(x, S, type='l', col='blue')
	#points(x,s,type='p',col='red' )
	#points(x,ss,type='l',col='orange' )
	#cat(res.fit.s, '\n')
	#cat(c(mean, fc, phase, beta), '\n')
	
}


################################
################################
### Exploration and Preparation for model selection and fitting
###         estimate alpha and detect outliers
################################
################################
test.NB.distribution.read.counts = function(T = T)
{
  #### alpha for each feature (from Deseq)
  cat('TEST THE FIT on THE REAL DATA\n')
  data.version = '_total_counts_v2'
  load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
  source("/Users/jiwang/Degradation_Liver/Main_Code_Degradation/f24_modified_1.0.r")
  T = R;
  
  examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
  mm = match(examples, T[,1])
  #mm = mm[order(-T[mm, 5])]
  mm = mm[which(!is.na(mm)==TRUE)]
  T[mm, c(1:5)]
  
  #### alpha for each feature and each condition (from Desq)
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  zt = seq(0,94,by = 2)
  zt.p = seq(0, 94, by=2)
  #zt = seq(0, 22, by=2)
  #source('functions.R')
 
  alpha.per.condition = TRUE
  aa = cbind(alphas[c(1:nrow(T)), ], alphas[-c(1:nrow(T)), ])
  
  pdfname = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/dispersion_examples_clock_genes_logscale_clock_fitting_NB_alpha_gene_conditions_specific_all_timepoints_likelihood_contribution_outliers.pdf';
  outlier.s1 = c()
  outlier.m1 = c()
  pdf(pdfname, width=20, height=12)
  for(n in mm)
  {
    #i = n;
    
    keep = keep5
    i = which(T[,1]=='Per2')
    n = i;
    param.fits.results =  unlist(keep[which(keep[,1]==T[i, 1]), -c(1)])  
    aa = alphas;
    
    R.m = T[n, ZT.ex]
    R.s = T[n, ZT.int]
    a.m = as.numeric(rep(aa[n, c(1:12)], 4))
    a.s = as.numeric(rep(aa[n, c(13:24)], 4))
    
    #model = 4;
    gamma = param.fits.results[27];eps.gamma = param.fits.results[28];phase.gamma = param.fits.results[29];splicing.k =  param.fits.results[30];
    Min = param.fits.results[31]; Amp=param.fits.results[32]; phase=param.fits.results[33]; beta=param.fits.results[34];
    m4 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s4 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    read.m4 = convert.nb.reads(m4, T$length.mRNA[i])
    read.s4 = convert.nb.reads(s4, T$length.premRNA[i])
    #read.s4 = s1;
    #read.m4 = m1;
    
    #model = 3;
    gamma = param.fits.results[16];eps.gamma = param.fits.results[17];phase.gamma = param.fits.results[18];splicing.k =  param.fits.results[19];
    Min = param.fits.results[20]; Amp=0; phase=12; beta=1;
    m3 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta);
    s3 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    read.m3 = convert.nb.reads(m3, T$length.mRNA[i])
    read.s3 = convert.nb.reads(s3, T$length.premRNA[i])
    
    #model = 2;
    gamma = param.fits.results[3];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
    Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
    m2 = compute.m.beta(zt.p, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    s2 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    read.m2 = convert.nb.reads(m2, T$length.mRNA[i])
    read.s2 = convert.nb.reads(s2, T$length.premRNA[i])
    #read.s2 = s2;
    #read.m2 = m2;
    
    ## find outliers
    loglike.m = matrix(NA, nrow=3, ncol=48)
    loglike.s = matrix(NA, nrow=3, ncol=48)
    for(kk in c(2:4))
    {
      #kk = 4
      eval(parse(text=paste('averg = read.m', kk, sep=''))); 
      err1 = qnbinom(0.025, size = 1/a.m, mu=averg)
      err2 = qnbinom(0.975, size = 1/a.m, mu=averg)
      cat(which(R.s<err1 | R.s>err2), '\n')
      err1 = err1 + 1;
      loglike = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(averg), log = TRUE)
      cat(index.outliers.loglike(loglike), '\n')
      loglike.m[(kk-1), ] = loglike
      
      eval(parse(text=paste('averg = read.s', kk, sep=''))); 
      err1 = qnbinom(0.025, size = 1/a.s, mu=averg)
      err2 = qnbinom(0.975, size = 1/a.s, mu=averg)
      cat(which(R.s<err1 | R.s>err2), '\n')
      err1 = err1 + 1;
      loglike = -2*dnbinom(as.numeric(R.s), size=1/a.s, mu=as.numeric(averg), log = TRUE)
      cat(index.outliers.loglike(loglike), '\n')
      loglike.s[(kk-1), ] = loglike
    }
    
    outlier.m = intersect(index.outliers.loglike(loglike.m[3,]), intersect(index.outliers.loglike(loglike.m[1,]), index.outliers.loglike(loglike.m[2,])))
    outlier.s = intersect(index.outliers.loglike(loglike.s[3,]), intersect(index.outliers.loglike(loglike.s[1,]), index.outliers.loglike(loglike.s[2,])))
    
    apply(cbind(loglike.s, loglike.m), 1, sum)
    apply(cbind(loglike.s[, setdiff(c(1:48), outlier.s)], loglike.m[, setdiff(c(1:48), outlier.m)]), 1, sum)
    
    par(mfrow=c(2,3))
    
    plot(zt, (R.s+1), ylim=range(c(err1, err2)), log='y', main=paste(T[n, 1], ' : premRNA M2'))
    points(zt.p, read.s2, type='b', col='green', pch=10, cex=0.6)
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='green', lwd=1.2, cex=0.7)
    
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    points(zt.p[index.out.s], R.s[index.out.s], col='red', cex=1.0)
    if(length(index.out.s)>0)
    {
      test = cbind(index = index.out.s, mean.M2 = averg[index.out.s], loglike.M2 =loglike[index.out.s])
    }
    
    averg = read.s4
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    err1 = qnbinom(0.05, size = 1/a.s, mu=averg)+1
    err2 = qnbinom(0.95, size = 1/a.s, mu=averg)
    
    plot(zt, (R.s+1), ylim=range(err1, err2), log='y', main=paste(T[n, 1], ' : premRNA M4'))
    points(zt.p, read.s4, type='b', col='green', pch=10, cex=0.6)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='green', lwd=1.2, cex=0.7)
    loglike = -2*dnbinom(as.numeric(R.s), size=1/a.s, mu=as.numeric(averg), log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    points(zt.p[index.out.s], R.s[index.out.s], col='red', cex=1.0)
    if(length(index.out.s)>0)
    {
      test = cbind(test, mean.M4 = averg[index.out.s], loglike.M4 =loglike[index.out.s])
    }
    outlier.s1 = rbind(outlier.s1, test)
    
    averg = read.m2
    err1 = qnbinom(0.05, size = 1/a.m, mu=averg)+1
    err2 = qnbinom(0.95, size = 1/a.m, mu=averg)
    index.out.m = which(R.m<err1 | R.m>err2) 
    plot(zt, (R.m+1), ylim=range(c(err1, err2)), log='y', main=paste(T[n, 1], ' : mRNA M2'))
    points(zt.p, read.m2, type='b', col='blue', pch=10, cex=0.6)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='blue', lwd=1.2, cex=0.7)
    loglike = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(averg), log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    points(zt.p[index.out.m], R.m[index.out.m], col='red', cex=1.0)
    if(length(index.out.m)>0)
    {
      test = cbind(index = index.out.m, mean.M2 = averg[index.out.m], loglike.M2 =loglike[index.out.m])
    }
    
    averg = read.m4
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    err1 = qnbinom(0.05, size = 1/a.m, mu=averg)+1
    err2 = qnbinom(0.95, size = 1/a.m, mu=averg)
    
    plot(zt, (R.m+1), ylim=range(err1, err2), log='y', main=paste(T[n, 1], ' : mRNA M4'))
    points(zt.p, read.m4, type='b', col='blue', cex=0.6, pch=10)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='blue', lwd=1.2, cex=0.7)
    loglike = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(averg), log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    points(zt.p[index.out.m], R.m[index.out.m], col='red', cex=1.0)
    if(length(index.out.m)>0)
    {
      test = cbind(test, mean.M4 = averg[index.out.m], loglike.M4 =loglike[index.out.m])
    }
    outlier.m1 = rbind(outlier.m1, test)
    
  }
  dev.off()
  
  ##### NB vs lognormal distribution
  
  res.version = '_total_all_v1';
  load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
  BIC = my.BIC.loglike(T = T)
  T = data.frame(T, BIC, stringsAsFactors = FALSE);
  T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
  
  examples = c('Npas2', 'Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Bhlhe41', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
  mm = match(examples, T[,1])
  #mm = mm[order(-T[mm, 5])]
  mm = mm[which(!is.na(mm)==TRUE)]
  T[mm, c(1:5)]
  
  
  outlier.s2 = c();
  outlier.m2 = c();
  #Tt = T
  
  ZT.int = grep('.abs.premRNA', colnames(T))
  ZT.ex = grep('.abs.mRNA', colnames(T))
  zt = seq(0,94,by = 2)
  
  #zt = seq(0, 22, by=2)
  #source('functions.R')
  zt.p = seq(0, 94, by=2)
  #alpha.per.condition = TRUE
  #aa = cbind(alphas[c(1:nrow(T)), ], alphas[-c(1:nrow(T)), ])
  pdfname = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/Examples_clock_genes_logscale_clock_fitting_lognormal_outliers.pdf';
  
  pdf(pdfname, width=20, height=12)
  for(n in mm)
  {
    i = n;
    
    param.fits.results =  unlist(T[n, c(112:159)])  
    model = 4;
    gamma = param.fits.results[27];eps.gamma = param.fits.results[28];phase.gamma = param.fits.results[29];splicing.k =  param.fits.results[30];
    Min = param.fits.results[31]; Amp=param.fits.results[32]; phase=param.fits.results[33]; beta=param.fits.results[34];
    m4 = compute.m.beta(zt.p, (gamma+0.5*eps.gamma), eps.gamma/(2*gamma+eps.gamma), phase.gamma, splicing.k, Min, Amp, phase, beta)
    s4 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    #read.m4 = convert.nb.reads(m4, T$length.mRNA[i])
    #read.s4 = convert.nb.reads(s4, T$length.premRNA[i])
    #read.s4 = s4;
    #read.m4 = m4;
    model = 2;
    gamma = param.fits.results[3];eps.gamma = 0;phase.gamma = 12;splicing.k =  param.fits.results[4];
    Min = param.fits.results[5]; Amp=param.fits.results[6]; phase=param.fits.results[7]; beta=param.fits.results[8];
    m2 = compute.m.beta(zt.p, (gamma+0.5*eps.gamma), eps.gamma/(2*gamma+eps.gamma), phase.gamma, splicing.k, Min, Amp, phase, beta)
    s2 = compute.s.beta(zt.p, Min, Amp, phase, beta);
    #read.m2 = convert.nb.reads(m2, T$length.mRNA[i])
    #read.s2 = convert.nb.reads(s2, T$length.premRNA[i])
    read.s2 = s2;
    read.m2 = m2;
    
    R.m = T[n, ZT.ex]
    R.s = T[n, ZT.int]
    sigma.ss = variance.estimate.replicates(log(R.s), nb.replicates=4)
    sigma.mm = variance.estimate.replicates(log(R.m), nb.replicates=4)
    
    par(mfrow=c(2,2))
    
    averg = read.s2
    #err1 = qnbinom(0.05, size = 1/a.s, mu=averg)
    err1 = qlnorm(0.05, meanlog = log(averg), sdlog = sigma.ss)
    err2 = qlnorm(0.95, meanlog = log(averg), sdlog = sigma.ss)
    index.out.s = which(R.s<err1 | R.s>err2) 
    
    plot(zt, (R.s), ylim=range(c(err1, err2)), log='y', main=paste(T[n, 1], ' : premRNA M2'))
    #plot(zt, R.s, log='y', main=paste(T[n, 1], ' : premRNA M2'))
    points(zt.p, read.s2, type='b', col='green', pch=10, cex=0.6)
    points(zt.p[index.out.s], R.s[index.out.s], col='red', cex=1.0)
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='green', lwd=1.2, cex=0.7)
    #dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
    loglike = -2*dlnorm(as.numeric(R.s), meanlog=log(averg), sdlog=sigma.ss, log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    
    if(length(index.out.s)>0)
    {
      test = cbind(index = index.out.s, mean.M2 = averg[index.out.s], loglike.M2 = loglike[index.out.s])
    }
    
    averg = read.s4
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    err1 = qlnorm(0.05, meanlog = log(averg), sdlog = sigma.ss)
    err2 = qlnorm(0.95, meanlog = log(averg), sdlog = sigma.ss)
    plot(zt, (R.s), ylim=range(err1, err2), log='y', main=paste(T[n, 1], ' : premRNA M4'))
    points(zt.p, read.s4, type='b', col='green', pch=10, cex=0.6)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='green', lwd=1.2, cex=0.7)
    loglike = -2*dlnorm(as.numeric(R.s), meanlog=log(averg), sdlog=sigma.ss, log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    points(zt.p[index.out.s], R.s[index.out.s], col='red', cex=1.0)
    if(length(index.out.s)>0)
    {
      test = cbind(test, mean.M4 = averg[index.out.s], loglike.M4 =loglike[index.out.s])
    }
    outlier.s2 = rbind(outlier.s2, test)
    
    averg = read.m2
    err1 = qlnorm(0.05, meanlog = log(averg), sdlog = sigma.mm)
    err2 = qlnorm(0.95, meanlog = log(averg), sdlog = sigma.mm)
    
    index.out.m = which(R.m<err1 | R.m>err2) 
    
    plot(zt, (R.m), ylim=range(c(err1, err2)), log='y', main=paste(T[n, 1], ' : mRNA M2'))
    points(zt.p, read.m2, type='b', col='blue', pch=10, cex=0.6)
    points(zt.p[index.out.m], R.m[index.out.m], col='red', cex=1.0)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='blue', lwd=1.2, cex=0.7)
    loglike = -2*dlnorm(as.numeric(R.m), meanlog=log(averg), sdlog=sigma.mm, log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    
    if(length(index.out.m)>0)
    {
      test = cbind(index = index.out.m, mean.M2 = averg[index.out.m], loglike.M2 = loglike[index.out.m])
    }
    
    averg = read.m4
    #err = sqrt(T$alpha.premRNA[n]*averg^2 + averg)
    err1 = qlnorm(0.05, meanlog = log(averg), sdlog = sigma.mm)
    err2 = qlnorm(0.95, meanlog = log(averg), sdlog = sigma.mm)
    plot(zt, (R.m), ylim=range(err1, err2), log='y', main=paste(T[n, 1], ' : mRNA M4'))
    points(zt.p, read.m4, type='b', col='blue', cex=0.6, pch=10)
    points(zt.p[index.out.m], R.m[index.out.m], col='red', cex=1.0)
    arrows(zt.p, err1, zt.p, err2, length=0.05, angle=90, code=3, col='blue', lwd=1.2, cex=0.7)
    loglike = -2*dlnorm(as.numeric(R.m), meanlog=log(averg), sdlog=sigma.mm, log = TRUE)
    text(zt, averg, signif(loglike, d=3), cex=0.7, pos=3, offset=0.5, adj=0.5)
    
    if(length(index.out.m)>0)
    {
      test = cbind(test, mean.M4 = averg[index.out.m], loglike.M4 = loglike[index.out.m])
    }
    outlier.m2 = rbind(outlier.m2, test)
  }
  
  dev.off()
  
  par(mfrow=c(2,2))
  plot(outlier.s1[, 5], (outlier.s1[, 3] - outlier.s1[,5]), cex=0.5, ylim = c(-4.5, 4.5*5))
  abline(h=4.5, col='red', lwd=2.0)
  plot(outlier.m1[, 5], (outlier.m1[, 3] - outlier.m1[,5]), cex=0.5,  ylim = c(-4.5, 4.5*5))
  abline(h=4.5, col='red', lwd=2.0)
  plot(outlier.s2[, 5], (outlier.s2[, 3] - outlier.s2[,5]), cex=0.5,  ylim = c(-4.5, 4.5*5))
  abline(h=4.5, col='red', lwd=2.0)
  plot(outlier.m2[, 5], (outlier.m2[, 3] - outlier.m2[,5]), cex=0.5,  ylim = c(-4.5, 4.5*5))
  abline(h=4.5, col='red', lwd=2.0)
  
  par(mfrow=c(2,2))
  hist((outlier.m2[, 3] - outlier.m2[,5]), xlim=c(-5, 30), breaks=seq(-5, 30, by=1))
  hist((outlier.s2[, 3] - outlier.s2[,5]), xlim=c(-5, 30), breaks=seq(-5, 30, by=1))
  hist((outlier.m1[, 3] - outlier.m1[,5]), xlim=c(-5, 30), breaks=seq(-5, 30, by=1))
  hist((outlier.s1[, 3] - outlier.s1[,5]), xlim=c(-5, 30), breaks=seq(-10, 40, by=1))
  
  ### compare NB and lognormal distribution
  xx = rnbinom(n=10000, size=1/0.05, mu=1000)
  outs = seq(1000, 4000, by=1000)
  yy = rlnorm(n=10000, meanlog=log(1000), sdlog = sqrt(0.05))
  yy1 = rlnorm(n=10000, meanlog=log(500), sdlog = sqrt(0.05))
  yy2 = rlnorm(n=10000, meanlog=log(2000), sdlog = sqrt(0.05))
  
  xlims = range(xx, yy, yy1, yy2)
  ylims = range(density(xx)$y, density(yy)$y, density(yy1)$y, density(yy2)$y)
  
  plot(density(xx), xlim=xlims, ylim=ylims, col='blue', lwd=2.0)
  points(density(yy), col='red', lwd=2.0, type='l')
  points(density(yy1), col='red', lwd=2.0, type='l')
  points(density(yy2), col='red', lwd=2.0, type='l')
  
  for(n in 1:100)
  {
    n = mm[1]
    vm.ex = variance.mean.rpkm.gene(log(T[n, ZT.ex]))
    vm.int = variance.mean.rpkm.gene(log(T[n, ZT.int]))
    vm.ex[, 2] = sqrt(sqrt(vm.ex[,2]))
    vm.int[, 2] = sqrt(sqrt(vm.int[,2]))
    
    plot(vm.ex, xlim=range(vm.ex[,1], vm.int[,1]), ylim=range(vm.ex[,2], vm.int[, 2]), log='', col='blue')
    points(vm.int, col='green')
    
  }
  
  xx = c();
  aa = cbind(alphas.genes[c(1:nrow(T)), ], alphas.genes[-c(1:nrow(T)), ])
  #for(n in 1:nrow(T))
  for(n in 1:500)
  {
    vm.ex = variance.mean.counts(T[n, ZT.ex], normalization = TRUE)
    vm.int = variance.mean.counts(T[n, ZT.int], normalization = TRUE)
    #vm.ex = cbind(vm.ex, (T$alpha.mRNA[n]*vm.ex[,2]^2 + vm.ex[,2]))
    #vm.int = cbind(vm.int, (T$alpha.premRNA[n]*vm.int[,2]^2 + vm.int[,2]))
    vm.ex = cbind(vm.ex, (aa[n, c(1:12)]*vm.ex[,2]^2 + vm.ex[,2]))
    vm.int = cbind(vm.int, (aa[n, c(13:24)]*vm.int[,2]^2 + vm.int[,2]))
    vm = rbind(vm.ex, vm.int);
    xx = rbind(xx, vm)
  }
  
  plot(xx[, c(3:4)], cex=0.25, log='xy', xlab='replicate variance', ylab='predicted variance NB')
  abline(0, 1, lwd=2.0, col='red')
  
  kk = which(T$qv.rpkm.mRNA<0.005)
  data = as.matrix(T[kk, grep('count.mRNA', colnames(T))])
  colnames(data) = rep(colnames(data)[c(1:12)], 4)
  
  pdfname = '/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myplots/Total/Total_counts/PCA_2vs3.pdf';
  pdf(pdfname, width=16, height=16)
  biplot(princomp(data[c(1:2000), ]), 2:3, cex = 0.8, arrow.len=0.05, xlim=c(-0.5, 0.5))
  dev.off()
  
  #### comparison between Deseq and my own method
  
}



estimate.dispersion.parameter.identify.outliers = function(T = T, gene.index = 1, zt = seq(0,94,by = 2), ZT.index = ZT.ex)
{
  #T = T; gene.index = mm[5];  zt = seq(0,94,by = 2); alpha.mRNA = TRUE; 
  #if(alpha.mRNA){
   # ZT.index = ZT.ex; 
  #}else{
  #  ZT.index = ZT.int; 
  #} 
  #xx = data.frame(alphas[c(1:nrow(T)),], alphas[-c(1:nrow(T)), ]);colnames(xx) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''));
  #yy = data.frame(alphas.genes[c(1:nrow(T)),], alphas.genes[-c(1:nrow(T)), ]);colnames(yy) = c(paste('alpha.mRNA.ZT', c(0:11)*2, sep=''),  paste('alpha.premRNA.ZT', c(0:11)*2, sep=''));
  
  reads = as.numeric(T[gene.index, ZT.index])
  outlier.index = c()
  nb.additonal = 1; 
  
  while(nb.additonal>0 & length(outlier.index)<=6)
  {
    res.fit = estimate.dispersion.parameter(reads = reads, zt = seq(0,94,by = 2), outlier.index = outlier.index)
    load(file='size_factors_48_samples.Rdata')
    #set.scaling.factors();
    #size.factors = scaling.factors;
    #reads.norm = reads/as.numeric(size.factors)
    
    mu = size.factors*compute.s.beta(t=zt, res.fit[1], res.fit[2], res.fit[3], res.fit[4]);
    alpha = 10^(res.fit[5]*(1+res.fit[6]*cos(w*(zt-res.fit[7]))));
    err1 = qnbinom(0.025, size = 1/alpha, mu=mu);
    err2 = qnbinom(0.975, size = 1/alpha, mu=mu);
    #additional.index = which(reads<err1 | reads>err2);
    loglike = -2*dnbinom(as.numeric(reads), size=1/alpha, mu=mu, log = TRUE)
    
    additional.index = index.outliers.loglike(loglike);
    cat('additional Outlier  ', additional.index,  '\n');
    
    if(length(outlier.index)==0)
    {
      outlier.index = additional.index;
    }else{
      additional.index = setdiff(additional.index, outlier.index);
      outlier.index = c(outlier.index, additional.index);
    }
    nb.additonal = length(additional.index);
    
    outlier.index = outlier.index[order(outlier.index)]
    cat('Outlier  ', outlier.index,  '\n');
    cat('parameters  ', res.fit[c(1:4)], '\n');
    cat('alpha  ',  res.fit[c(5:7)], '\n');
    
    ### NOT detect outliers
    nb.additonal = 0;
  }
  
  test.fit = FALSE
  if(test.fit)
  {
    estimate = rough.estimate.mean.alpha(reads);
    estimate[2, which(estimate[2,]<=0)] = 10^-6;
    
    xx = alphas;
    yy = alphas.genes
    aa = estimate[2,]
    par.init = res.fit;
    mu = size.factors*compute.s.beta(t=zt, par.init[1], par.init[2], par.init[3], par.init[4]);
    alpha = 10^(par.init[5]*(1+par.init[6]*cos(w*(zt-par.init[7]))));
    
    par(mfcol=c(1,2),cex=1.0)
    plot(zt, reads, log='y', cex=1.2, ylim=range(err1, err2), main=T[gene.index, 1])
    points(zt[outlier.index], reads[outlier.index], col='red', cex=1.2, pch=16)
    points(zt, mu, col='blue',type='l')
    points(zt, mu, col='blue',type='p', pch=20)
    arrows(zt, err1, zt, err2, length=0.05, angle=90, code=3, col='blue', lwd=1.2, cex=0.7)
    
    if(alpha.mRNA)
    {
      xx.des = xx[gene.index, grep('alpha.mRNA', colnames(xx))];
      yy.des = yy[gene.index, grep('alpha.mRNA', colnames(yy))];
    }else{
      xx.des = xx[gene.index, grep('alpha.premRNA', colnames(xx))];
      yy.des = yy[gene.index, grep('alpha.premRNA', colnames(yy))];
    }
    plot(c(0:11)*2, rep(estimate[2,], 1), ylim=range(aa, alpha, xx.des, yy.des), log='y')
    points(c(0:11)*2, alpha[c(1:12)], col='darkblue',type='l', lwd=2.0, cex=1.5)
    points(c(0:11)*2, xx.des, col='red', pch=0, cex=1.2)
    points(c(0:11)*2, yy.des, col='black', pch=16, cex=0.8)
  }
  
  return(list(param.fit = res.fit, outlier.index = outlier.index))
  
}

estimate.dispersion.parameter = function(reads, zt = seq(0,94,by = 2), outlier.index = c())
{
  #reads = as.numeric(T[gene.index, ZT.index])
  if(length(outlier.index)>0) reads[outlier.index] = NA;
  load(file='size_factors_48_samples.Rdata');
  reads.norm = reads/as.numeric(size.factors);
  estimate = rough.estimate.mean.alpha(reads);
  estimate[2, which(estimate[2,]<=0)] = 10^-6;
  stat.mm = f24_R2_alt2(reads.norm, t=c(0:47)*2)
  #stat.aa = f24_R2_alt2(estimate[2, ], t=c(0:11)*2)
  stat.aa = f24_R2_alt2(log10(estimate[2, ]), t=c(0:11)*2)
  if(stat.aa[4]<0) 
  {
    stat.aa[4] = - stat.aa[4]
    stat.aa[5] = (stat.aa[5]+12)%%24
  }
  #plot(c(0:11)*2, log10(estimate[2,]))
  #tt = c(0:11)*2;
  #points(c(0:11)*2, stat.aa[2]*(1+stat.aa[4]*cos(2*pi/24*(tt - stat.aa[5]))), type='l', col='blue')
  Nfit.alpha = 6
  for(model in c(2,4))
  {
    if(model==4)
    {
      PAR.INIT.alpha = cbind(rep(min(reads.norm[which(!is.na(reads.norm)==TRUE)]), Nfit.alpha), 
                            rep(stat.mm[3], Nfit.alpha), 
                            rep(stat.mm[5], Nfit.alpha), 
                            lseq(1, 5, length = Nfit.alpha), 
                            rep(mean(stat.aa[2]), Nfit.alpha),  
                            rep(stat.aa[4], Nfit.alpha), 
                            rep(stat.aa[5], Nfit.alpha))
      
      errors.fit = rep(NA, Nfit.alpha)
      for(fit.nb in 1:Nfit.alpha)
      {
        # fit.nb = 1;
        par.init = PAR.INIT.alpha[fit.nb,];
        opt = optim(par.init, f2min.alpha, reads=reads, zt = zt, model = 4, method = 'L-BFGS-B', 
                    lower = c(min(estimate[1,])/5, 0, 0, 1, -6, 0, 0), 
                    upper = c(max(estimate[1,])*5, max(c(1000, max(estimate[1,])-min(estimate[1,])))*2, 24, 5, 0, 1, 24))
        #print(opt.s$convergence)
        res.fit = opt$par
        errors.fit[fit.nb] = opt$value
        eval(parse(text = paste('res.fit.', fit.nb, ' = res.fit', sep = '')))
      }
      
      ## choose the best-fitting parameters
      imin.alpha = which.min(errors.fit); 
      eval(parse(text = paste('res.fit.m4 = res.fit.', imin.alpha, sep = '')))
      eval(parse(text = paste('error.m4 = errors.fit[', imin.alpha, ']', sep = ''))) 
    }
    if(model==2)
    {
      PAR.INIT.alpha = cbind(rep(min(reads.norm[which(!is.na(reads.norm)==TRUE)]), Nfit.alpha),
                             rep(stat.mm[3], Nfit.alpha), 
                             rep(stat.mm[5], Nfit.alpha), 
                             lseq(1, 5, length = Nfit.alpha), 
                             rep(stat.aa[2], Nfit.alpha))
      
      errors.fit = rep(NA, Nfit.alpha)
      for(fit.nb in 1:Nfit.alpha)
      {
        # fit.nb = 1;
        par.init = PAR.INIT.alpha[fit.nb,];
        opt = optim(par.init, f2min.alpha, reads=reads, zt = zt, model = 2, method = 'L-BFGS-B', 
                    lower = c(min(estimate[1,])/5, 0, 0, 1, -6), 
                    upper = c(max(estimate[1,])*5, max(c(1000, max(estimate[1,])-min(estimate[1,])))*2, 24, 5, 0))
        res.fit = opt$par
        errors.fit[fit.nb] = opt$value
        eval(parse(text = paste('res.fit.', fit.nb, ' = res.fit', sep = '')))
      }
      ## choose the best-fitting parameters
      imin.alpha = which.min(errors.fit); 
      eval(parse(text = paste('res.fit.m2 = res.fit.', imin.alpha, sep = '')))
      eval(parse(text = paste('error.m2 = errors.fit[', imin.alpha, ']', sep = ''))) 
    }
    if(model==3)
    {
      PAR.INIT.alpha = cbind(rep(mean(estimate[1,]), Nfit.alpha), sample(seq(0, 1, 0.05), Nfit.alpha), 
                             sample(seq(0, 24, 2), Nfit.alpha), lseq(1, 5, length = Nfit.alpha), 
                             rep(mean(estimate[2,]), Nfit.alpha),  rep(0, Nfit.alpha), sample(seq(0, 24, 2), Nfit.alpha));
    }
    if(model==1)
    {
      PAR.INIT.alpha = cbind(rep(mean(estimate[1,]), Nfit.alpha), rep(0, Nfit.alpha), 
                             sample(seq(0, 24, 2), Nfit.alpha), lseq(1, 5, length = Nfit.alpha), 
                             rep(mean(estimate[2,]), Nfit.alpha),  rep(0, Nfit.alpha), sample(seq(0, 24, 2), Nfit.alpha)); 
    }
  }
  
  ### model selection for flat or rhythmic alpha
  AIC.m2 = error.m2 + length(res.fit.m2)*2;
  AIC.m4 = error.m4 + length(res.fit.m4)*2;
  AIC = c(AIC.m2, AIC.m4)
  aic = AIC-min(AIC)
  prob.model = exp(-0.5*aic)
  prob.model = prob.model/sum(prob.model)
  if(prob.model[1]>=0.7)
  {
    #return(c(res.fit.m2, 0, 0));
    return(res.fit.m4);
  }else{
    return(res.fit.m4);
  }
}

f2min.alpha = function(par.init, reads, zt = seq(0,96,by = 2), model=4)
{
  w = 2*pi/24;
  #load(file='size_factors_48_samples.Rdata')
  set.scaling.factors();
  size.factors = scaling.factors;
  
  if(model==4)
  {
    mu = size.factors*compute.s.beta(t=zt, par.init[1], par.init[2], par.init[3], par.init[4]);
    #mu = size.factors*rep(estimate[1,], 4)
    #alpha = par.init[5]*(1+par.init[6]*cos(w*(zt-par.init[7])));
    alpha = 10^(par.init[5]*(1+par.init[6]*cos(w*(zt-par.init[7]))));
  }
  if(model==2)
  {
    mu = size.factors*compute.s.beta(t=zt, par.init[1], par.init[2], par.init[3], par.init[4]);
    #alpha = rep(par.init[5], length(zt));
    alpha = 10^rep(par.init[5], length(zt));
  }
  
  if(model==1)
  {
    mu = size.factors*rep(par.init[1], length(zt));
    alpha = rep(par.init[5], length(zt));
  }
  if(model==3)
  {
    mu = size.factors*rep(par.init[1], length(zt));
    alpha = par.init[5]*(1+par.init[6]*cos(w*(zt-par.init[7])));
  }
  
  kk = which(!is.na(reads)==TRUE)
  if(any(reads[kk]<0)|any(mu<=0))
  {
    error = 10^10;
  }else{
    error = -2*sum(dnbinom(reads[kk], size=1/alpha[kk], mu=mu[kk], log = TRUE))
  }
  #cat('param ...', par.init, '\n')
  #cat('error ...', error, '\n')
  return(error)
}

rough.estimate.mean.alpha = function(x, period=24, interval=2)
{
  # x = reads;period=24;interval=2;
  x = as.numeric(x)
  #load(file='size_factors_48_samples.Rdata')
  set.scaling.factors();
  size.factors = scaling.factors;
  x = x/as.numeric(size.factors)
  
  #if(length(outlier.index)>0) x[outlier.index] = NA;
  #kk = which(!is.na(x)==TRUE)
  index = period/interval;
  nn = length(x)/(period/interval)
  test = c()
  for(mm in 1:nn) test = cbind(test, x[c(1:index)+index*(mm-1)])
  test = data.frame(test)
  
  mm = apply(test, 1, mean, na.rm=TRUE)
  vv = apply(test, 1, var, na.rm = TRUE)
  #aa =  (vv - mm)/mm^2
  return(rbind(mean=mm, alpha = (vv - mm)/mm^2))
}

f24_R2_alt2=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
  
  kk = which(!is.na(x)==TRUE)
  x = x[kk]
  t = t[kk]
  n=length(x)
  #mu=mean(x)
  nb.timepoints=length(x)
  if(n<4)
  { 
    if(n==0) c(nb.timepoints=nb.timepoints, mean=NA, amp=NA, relamp=NA,phase=NA,pval=NA) 
    else 
    {
      c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)
    }
  }
  else
  {
    sig2=var(x)
    c=cos(2*pi*t/period)
    s=sin(2*pi*t/period)
    A = mean(x*c)-mean(x)*mean(c)
    B = mean(x*s)-mean(x)*mean(s)
    c1 = mean(c^2)-mean(c)^2
    c2 = mean(c*s)-mean(c)*mean(s)
    c3 = mean(s^2)-mean(s)^2
    b = (A*c2-B*c1)/(c2^2-c1*c3)
    a = (A-b*c2)/c1
    mu = mean(x)-a*mean(c)-b*mean(s)
    #	b=2*mean(x*s)
    x.hat=mu+a*c+b*s
    sig2.1=var(x-x.hat)
    if(is.na(a)||is.na(b)) {c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)}
    else
    {
      p=3
      R2=0
      if(sig2>0) R2=1-sig2.1/sig2
      # http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
      # I (Felix) checked that it works
      amp=max(x)-min(x)
      phase=period/(2*pi)*atan2(b, a)
      if(phase<0) phase=phase+period
      if(phase>period) phase=phase-period
      phase=(phase+offset)%%period
      pval = pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)
      
      c(nb.timepoints=nb.timepoints, mean=mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(mu),phase=phase, pval=pval)
    }
  }
}





##################################################################################################################
## ---------- finishing line of optimization --------------------------------
##################################################################################################################



################
######## Output 
###############
col.ramp.ex <- colorRampPalette(c('steelblue4', 'steelblue1' , 'white' , 'orangered' , 'orangered3'), space = "Lab")
col.ramp.int <- colorRampPalette(c('darkturquoise', 'turquoise' , 'white' , 'gold1' , 'darkgoldenrod2'), space = "Lab")
col.ramp.pval <- colorRampPalette(c('black', 'blue' , 'royalblue' , 'steelblue1' , 'lightsteelblue1', 'white'), space = "Lab")
col.ramp.median <- colorRampPalette(c('gray90', 'firebrick1' , 'black'), space = "Lab")

attribute.global.variable.colors = function(){
	col.ex <<- 'steelblue'; col.int <<- 'green3'; col.accepted.ex <<- 'steelblue2'; col.accepted.int <<- 'limegreen'; col.rejected <<- 'gray' ; zt <<- seq(0,46,by=2)
	col.phases <<- rainbow(n = 241,s = 0.9, v = 0.9); 
	col.deg <<- 'orangered'
	Ncolors <<- 21; col.rel.ex <<- col.ramp.ex(n = Ncolors); col.rel.int <<- col.ramp.int(n = Ncolors); 
	col.pval <<- col.ramp.pval(n = 11); col.median.level  <<- col.ramp.median(n = 15);col.kept <<- c('yellowgreen','tomato1')
	transparent.factor = 0.3
	transparent.gray = rgb(col2rgb('gray')[1]/255,col2rgb('gray')[2]/255,col2rgb('gray')[3]/255, transparent.factor)
	transparent.green3 = rgb(col2rgb('green3')[1]/255,col2rgb('green3')[2]/255,col2rgb('green3')[3]/255, transparent.factor)
	transparent.tomato = rgb(col2rgb('tomato')[1]/255,col2rgb('tomato')[2]/255,col2rgb('tomato')[3]/255, transparent.factor)
	transparent.black = rgb(col2rgb('black')[1]/255,col2rgb('black')[2]/255,col2rgb('black')[3]/255, transparent.factor)
	palette(c('gray','green3','tomato','black', transparent.gray ,transparent.green3, transparent.tomato, transparent.black))
}	

attribute.global.colors = function()
{
    col.ex <<- 'steelblue'; col.int <<- 'green3'; col.deg <<- 'orangered';
    col.m1 <<- 'gray'; col.m2 = 2; col.m3 = 3; col.m4 = 4;
}

contour.plot.objective.funciton.parameters = function(T = T, i=1, zt =  seq(0,46,by = 2), 
                                                      i.ex = ZT.ex, i.int = ZT.int, outliers = FALSE)
{
  # j = kk[jj3[1]]; zt =  seq(0,46,by = 2); i.ex = ZT.ex; i.int = ZT.int; outliers = FALSE;i=j;
  R.m = unlist(T[i, i.ex]) ## nb of reads for exon
  R.s = unlist(T[i, i.int]) ## nb of reads for intron
  L.m = T$length.mRNA[i];
  L.s = T$length.premRNA[i];
  alpha.m = rep(as.numeric(T[i, grep('alpha.mRNA.ZT', colnames(T))]), 4);
  alpha.s = rep(as.numeric(T[i, grep('alpha.premRNA.ZT', colnames(T))]), 4);
  
  MCMC = FALSE
  if(MCMC)
  {
    #sample loglike function using MCMC
    require(adaptMCMC)
    source('functions.R')
    
    #par= c(T$gamma.m3[i], T$eps.gamma.m3[i], T$phase.gamma.m3[i], T$splicing.k.m3[i], T$Min.int.m3[i])
    par = c(T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i])
    v0=c(T$gamma.stderr.m3[i], T$eps.gamma.stderr.m3[i], T$phase.gamma.stderr.m3[i], T$splicing.k.stderr.m3[i], T$Min.int.stderr.m3[i])^2*0.01
    v0[1] = 0.01
    v0[4] = 1
    #mc=MCMC(f2min.4mcmc, p.fixed, R.m, R.s, L.m, L.s, alpha.m, alpha.s, model=3, n=1000, 
    #        p0, v0, acc.rate=0.2, adapt=T,x=x,y=y[,channel],var=fit$value/length(y[,1]))
    mc=MCMC(f2min.4mcmc,n=10000, par, scale=v0, acc.rate=0.2, adapt=TRUE, gamma = 0.6,
            R.m=R.m, R.s=R.s, L.m=L.m, L.s=L.s, alpha.m=alpha.m, alpha.s=alpha.s, model=3)
    sampling = data.frame(mc$samples, mc$log.p, stringsAsFactors = FALSE)
    colnames(sampling) = c('gamma', 'eps.gamma', 'phase.gamma', 'splicing.k', 'Min.int', 'log.prob')
    
    par.true = c(T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i])
    par.es = c(T$gamma.m3[i], T$eps.gamma.m3[i], T$phase.gamma.m3[i], T$splicing.k.m3[i], T$Min.int.m3[i])
    #par.es = c(T$gamma.m3[i], T$eps.gamma.m3[i], T$phase.gamma.m3[i])
    
    #colfunc <- colorRampPalette(c("black", "white"))
    #colfunc(10)
    
    par(mfrow=c(1,3),cex=1.0)
    plot(sampling[, c(1,2)], cex=0.1, log='', xlim=c(log(2)/24, log(2)/(10/60)), ylim=c(0, 1), col='gray') 
    points(par.es[1], par.es[2], cex=1., col='red', pch=16)
    points(par.true[1], par.true[2], cex=1., col='blue', pch=16)
    plot(sampling[, c(1,3)], cex=0.1, log='', xlim=c(log(2)/24, log(2)/(10/60)), ylim=c(0, 24), col='gray')
    points(par.es[1], par.es[3], cex=1., col='red', pch=16)
    points(par.true[1], par.true[3], cex=1., col='blue', pch=16)
    plot(sampling[, c(2,3)], cex=0.1, xlim=c(0, 1), ylim=c(0, 24), col='gray')
    points(par.es[2], par.es[3], cex=1., col='red', pch=16)
    points(par.true[2], par.true[3], cex=1., col='blue', pch=16)
    #plot(sampling[, c(2,4)], cex=0.2)
    #points(par[1], par[4], cex=1.0, col='red')
    require('mcmcplots')
    denplot(mc)
    
    library(MASS)
    k <- kde2d(sampling[,1], sampling[,2], n=300, lims=c(c(log(2)/24, log(2)/(10/60)), c(0, 1)))
    image(k, col=r)
    
    ### contour plot
    #par(mfrow=c(2,2),cex=1.0)
    #test = sampling[, c(1:2, 6)]
    #xx = unique(test[,1]);
    #xx = xx[order(xx)]
    #yy = unique(test[,2])
    #yy = unique(yy)
    #vals = matrix(NA, nrow=length(xx), ncol=length(yy))
    #for(n1 in 1:length(xx))
    #{
    #  for(n2 in 1:length(yy))
    # {
    #    m1 = which(test[,1]==xx[n1] & test[,2] == yy[n2])
    #    vals[n1, n2] = test[m1, 3]
    #  }
    #}
    
    
  }
  #p.gamma = c(T$gamma[i], T$eps.gamma[i], T$phase.gamma[i])
  #p.gamma = c(T$gamma.m3[i], T$eps.gamma.m3[i], T$phase.gamma.m3[i])
  #p.fixed = c(T$splicing.k[i], T$Min.int[i])
  #par=p.gamma
  ### plot objective function
  Plot.objective.function = FALSE
  if(Plot.objective.function)
  {
    source('functions.R')
    nb.gamma = 100
    nb.eps = 20;
    nb.phase = 20;
    samples.gamma = lseq(log(2)/24, log(2)/(10/60), length=nb.gamma)
    samples.eps.gamma = seq(max(0, min(c(T$eps.gamma[i], T$eps.gamma.m3[i]))-0.1), min(max(c(T$eps.gamma[i], T$eps.gamma.m3[i]))+0.1, 1), length=nb.eps)
    samples.phase.gamma = seq((min(c(T$phase.gamma[i], T$phase.gamma.m3[i]))-1), (max(c(T$phase.gamma[i], T$phase.gamma.m3[i]))+1), length=nb.phase)
    par.true = c(T$gamma[i], T$eps.gamma[i], T$phase.gamma[i], T$splicing.k[i], T$Min.int[i])
    par.es = c(T$gamma.m3[i], T$eps.gamma.m3[i], T$phase.gamma.m3[i], T$splicing.k.m3[i], T$Min.int.m3[i])
    
    fvalues = as.matrix(NA, nrow=nb.gamma*nb.eps*nb.phase, ncol=4)
    for(nn1 in 1:nb.gamma)
    {
      for(nn2 in 1:nb.eps)
      {
        for(nn3 in 1:nb.phase)
        {
          fvalues[(nn1+nn2+nn3-2), ] = c(sample.gamma[nn1], samples.eps.gamma[nn2], samples.phase.gamma[nn3], 
                                         f2min.4mcmc(par=c(sample.gamma[nn1], samples.eps.gamma[nn2], samples.phase.gamma[nn3], par.true[4]/par.true[1]*samples.gamma[nn1]), par.true[5],   
                                                     R.m=R.m, R.s=R.s, L.m=L.m, L.s=L.s, alpha.m=alpha.m, alpha.s=alpha.s, model=3, zt = seq(0,94,by = 2)))
        }
      }
    }
    
    ptm <- proc.time()
    f2min.4mcmc(par=par, p.fixed=p.fixed, R.m=R.m, R.s=R.s, L.m=L.m, L.s=L.s, alpha.m=alpha.m, alpha.s=alpha.s, model=3, zt = seq(0,94,by = 2))
    proc.time() - ptm  
  }

}

f2min.4mcmc = function(par, R.m, R.s, L.m, L.s, alpha.m, alpha.s, outlier.m = c(), outlier.s = c(), model=4, zt = seq(0,94,by = 2), 
                 parametrization =c('cosine.beta'), debug = FALSE, absolute.signal = TRUE)
{
  #alpha.m = rep(alpha.m, 4);
  #alpha.s = rep(alpha.s, 4);
  #par = par.init; simulation.only=FALSE
  if(par[1]>log(2)/(10/60) | par[1]<log(2)/24 | par[2]<0 |par[2]>1 | par[3]<0 |par[3]>24)
  {
    err.fit = -Inf;
  }else{
    # initialization of the parameters
    #par = c(par, p.fixed)
    gamma = par[1];
    if(model==2)	
    {
      eps.gamma = 0.0;
      phase.gamma = 12.0; 
      splicing.k = par[2];
      param.synthesis.1 = par[3];
      param.synthesis.2 = par[4];  
      param.synthesis.3 = par[5]; 
      param.synthesis.4 = par[6]; 
    }
    if(model==3)
    {
      eps.gamma = par[2];
      phase.gamma = par[3];
      splicing.k = par[4];
      param.synthesis.1 = par[5];
      param.synthesis.2 = 0;  
      param.synthesis.3 = 0; 
      param.synthesis.4 = 1; 
    }
    if(model==4)
    {
      eps.gamma = par[2];
      phase.gamma = par[3];
      splicing.k = par[4];
      param.synthesis.1 = par[5];
      param.synthesis.2 = par[6];  
      param.synthesis.3 = par[7]; 
      param.synthesis.4 = par[8]; 
    }
    
    s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
    #print(par);
    #print(compute.m.beta(t = zt, gamma, eps.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4, simulation.only=TRUE))
    m = compute.m.beta(t = zt, gamma, eps.gamma, phase.gamma, splicing.k,
                       param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
    
    ### convert rpkm into mean of read numbers
    mu.m = convert.nb.reads(m, L.m);
    mu.s = convert.nb.reads(s, L.s);
    #print(m);
    err.fit = -0.5*NB.error(R.m = R.m, R.s = R.s, alpha.m = alpha.m, alpha.s = alpha.s, mu.m = mu.m, mu.s = mu.s, 
                       outlier.m = outlier.m, outlier.s = outlier.s, specie = 'both');
    #print(par);
    #print(mu.m)
    #print(mu.s)
    #print(err.fit)
  }
  #print('........')
  #print(err.fit)
  return(err.fit)
}

compare.parameters.RNA.seq.total.fake = function(T = T, select.method='BIC.best.model', parametrization = c('consine.beta'), quality = FALSE, plot.parameters = FALSE)
{
  T$best.model = eval(parse(text=paste('T$', select.method , sep='')))
  
  T = T[which(!is.na(T$best.model)==TRUE),]
  palette(c('gray','green3','tomato','black'))
  best.model = T$best.model;
  true.model=T$true.model;
  if(!plot.parameters){
    ## model selection
    Prediction.Quality = matrix(NA, nrow = 4, ncol =4)
    for(model in 1:4)
    {
      #h = hist(T$LRT.best.model[((model-1)*X+1):(model*X)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
      h = hist(T$best.model[which(T$true.model==model)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
      Prediction.Quality[model,] = h$counts/sum(h$counts)
    }
    barplot(t(Prediction.Quality), col = 1:4, names.arg = c('CS-CD','RS-CD','CS-RD','RS-RD'), border = NA, cex.main=0.8, main =paste(select.method))
    
    #### noise (sigma) estimation for premRNA and mRNA
    attribute.global.colors();
    lims = range(c(T$sigma.s, T$sigma.se))
    plot(T$sigma.s, T$sigma.se, ylim=lims, xlim=lims, col=col.int, type='p', cex=0.6, xlab='True values', ylab='Estimated values', cex.main=0.8, main='Sigma of premRNAs')
    abline(0, 1, col='gray', lwd=1.0)
    
    lims = range(c(T$sigma.m, T$sigma.me))
    plot(T$sigma.m, T$sigma.me, ylim=lims, xlim=lims, col=col.ex, type='p', cex=0.6, xlab='True values', ylab='Estimated values', 
         cex.main=0.8, main='Sigma of mRNAs');
    abline(0, 1, col='gray', lwd=1.0)
    #### nothing in this panel
    plot(T$sigma.m, T$sigma.me, type='n', main=NA, xlab=NA, ylab=NA, axes=FALSE)
  }else{
    mains = c('min.premRNA', 'amp.premRNA','phase.premRNA','beta.premRNA',
              'eps.degr.norm','phase.degr.norm','splicing.k/gamma', '', 
              'half-life (M2)', 'half-life (M3)', 'half-life (M4)', '')
    parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 
                       'eps.gamma','phase.gamma','splicing.k',rep('gamma', 5))
    models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), 
                     c(0,0,3,4),c(0,0,3,4), c(0,2,3,4), c(0,2,3,4),
                     c(0,2,0,0), c(0,0,3,0), c(0,0,0,4), c(0,2,3,4))
    
    bounds = set.bounds(model = 4, absolute.signal=TRUE);
    lower = bounds$lower[c(5:8, 2:4, rep(1, 5))]; 
    upper = bounds$upper[c(5:8, 2:4, rep(1, 5))];
    lower[1] = range(T$Min.int, na.rm = TRUE)[1];upper[1]=range(T$Min.int, na.rm = TRUE)[2];
    lower[7] = range(T$splicing.k, na.rm = TRUE)[1];upper[7]=range(T$splicing.k, na.rm = TRUE)[2];
    lower[8:12] = 10/60;upper[8:12]=24;
    #lower[4] = log(2)/(30/60);
    #upper[1] = log(2)/(1/60); lower[2] = 0.01;
    #lower[6] = 0.01;upper[6] = 20;
    #lower = lower[c(5:8, 1:4)]; upper = upper[c(5:8, 1:4)]; upper[1] = 10;
    
    ### 8 parpameters and 8 plots
    layout(matrix(c(1:12),nrow = 3, ncol = 4, byrow=TRUE))
    par(cex=0.7, mgp = c(1.0,0.3,0.), mar = c(1.2, 1.,1.,0.5), pty='s');
    for(p in 1:12)
    {
      #par(mfrow = c(3,4), cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      if(p==8|p==12){
        plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab=NA, ylab=NA, main=NA, axes=FALSE)
      }else{
        log = ''; if((p==1)|(p==2)|(p>=7)){log ='xy';}
        
        plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = NA, ylab = NA, cex.main=0.6, main = mains[p], log = log, axes = FALSE)
        
        cex = 0.6; pch.tt = 1; pch.ff = 2; cex.axis=0.7;
        for(m in c(2:4)) # for each model 2,3,4
        {
          eval(parse(text = paste('ttrue = T$', parameter.list[p],'[best.model == m & true.model==m]', sep ='')))
          eval(parse(text = paste('ftrue = T$', parameter.list[p],'[best.model == m & true.model!=m]', sep ='')))
          eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.m',m,'[best.model == m & true.model==m ]', sep ='')));
          eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.m',m,'[best.model == m & true.model!=m ]', sep ='')));
          
          if(models.p[m,p]>0){
            if(p<9){
              points(ttrue, estimated.t, pch = pch.tt, col = m, bg = m, cex = cex);
              points(ftrue, estimated.f, pch = pch.ff, col = m, bg = m, cex = cex);
            }else{
              points(log(2)/ttrue, log(2)/estimated.t, pch = pch.tt, col = m, bg = m, cex = cex);
              points(log(2)/ftrue, log(2)/estimated.f, pch = pch.ff, col = m, bg = m, cex = cex);
            }
           
          }
        }
        abline(a =0, b = 1, col = 'gray', lwd=1.5); box();
        if(p==3|p==6){
          axis(1,  at = seq(0, 24, by=6), las = 1,  cex.axis=cex.axis, tck=-0.02, labels = FALSE); 
          axis(1,at = seq(0, 24, by=6), las=1,cex.axis=cex.axis, lwd=0, line = -0.3)
          axis(2, at = seq(0, 24, by=6), las = 3, cex.axis=cex.axis, tck=-0.02, labels = FALSE); 
          axis(2, at = seq(0, 24, by=6),las=3, cex.axis=cex.axis, lwd=0, line = -0.2)
        }else{
          axis(1,  las = 1,  cex.axis=cex.axis, tck=-0.02, labels = FALSE); 
          axis(1, las=1,cex.axis=cex.axis, lwd=0, line = -0.3)
          axis(2,  las = 3, cex.axis=cex.axis, tck=-0.02, labels = FALSE); 
          axis(2, las=3, cex.axis=cex.axis, lwd=0, line = -0.2)
        }
        
      }
      
    }
  }  
	
}

compare.parameters.RNA.seq.total.fake.by.model = function(T = T, select.method='BIC.best.model', parametrization = c('consine.beta'), quality = FALSE)
{
  #### model selection
  palette(c('gray','green3','tomato','black'))
  Prediction.Quality = matrix(NA, nrow = 4, ncol =4)
  for(model in 1:4)
  {
    #h = hist(T$LRT.best.model[((model-1)*X+1):(model*X)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
    h = hist(T$best.model[which(T$true.model==model)], breaks = seq(0.5,4.5,by = 1), plot = FALSE)
    Prediction.Quality[model,] = h$counts/sum(h$counts)
  }
  
  for(model in 2:4)
  {
    barplot(t(Prediction.Quality), col = 1:4, names.arg = c('M1','M2','M3','M4'), border = NA, cex.main=0.8, main =paste(select.method))
    
    #### quality of model selection in function of amplitudes
    plot(c(0:10), c(0:10), type='n', main=paste('M', model, sep=''));
    #### noise (sigma) estimation for premRNA and mRNA
    attribute.global.colors();
    lims = range(c(T$alpha.mRNA.ZT0, T$a.m.ZT0.estimated), na.rm=TRUE)
    plot(T$alpha.mRNA.ZT0, T$a.m.ZT0.estimated, ylim=lims, xlim=lims, col=col.ex, type='p', cex=0.6, 
         xlab='True values', ylab='Estimated values', cex.main=0.8, main='alpha of mRNAs at ZT00', log='xy')
    abline(0, 1, col='gray', lwd=2.0)
    abline(log(2), 1, col='gray', lwd=2.0)
    abline(log(0.5), 1, col='gray', lwd=2.0)
    
    lims = range(c(T$alpha.premRNA.ZT0, T$a.s.ZT0estimated), na.rm=TRUE)
    plot(T$alpha.premRNA.ZT0, T$a.s.ZT0estimated, ylim=lims, xlim=lims, col=col.int, type='p', cex=0.6, 
         xlab='True values', ylab='Estimated values', cex.main=0.8, main='alpha of premRNAs at ZT00', log='xy')
    abline(0, 1, col='gray', lwd=2.0)
    abline(log(2), 1, col='gray', lwd=2.0)
    abline(log(0.5), 1, col='gray', lwd=2.0)
    #### nothing in this panel
    #plot(T$sigma.m, T$sigma.me, type='n', main=NA, xlab=NA, ylab=NA, axes=FALSE)
    
    best.model = T$best.model;
    true.model=T$true.model;
    
    mains = c('Min-premRNA', 'Amp.premRNA','Phase.premRNA','Beta.premRNA','Degr.Rate','Degr.eps.gamma','Degr.Phase','Splicing.K')
    parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 'gamma','eps.gamma','phase.gamma','splicing.k')
    models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), c(0,2,3,4),c(0,0,3,4),c(0,0,3,4),c(0,2,3,4))
    
    bounds = set.bounds(model = 4, absolute.signal=TRUE);
    lower = bounds$lower; 
    upper = bounds$upper;
    #lower[4] = log(2)/(30/60);
    #upper[1] = log(2)/(10/60);
    #lower[2] = 0.0;
    #lower[6] = 0.0;
    upper[6] = 20;
    lower = lower[c(5:8, 1:4)]
    upper = upper[c(5:8, 1:4)]
    upper[1] = 10^3;
    lower[1] = 0.01;
    
    ### 8 parpameters and 8 plots
    m = model;
    for(p in 1:8)
    {
      log = ''; 
      if((p==1)|(p==2)|(p==4)|(p==5)|(p==8)){log ='xy'}
      
      plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = 'True values', ylab = 'Estimated values', cex.main=0.8, main = mains[p], log = log)
      eval(parse(text = paste('ttrue = T$', parameter.list[p],'[best.model == m & true.model==m]', sep ='')))
      eval(parse(text = paste('ftrue = T$', parameter.list[p],'[best.model == m & true.model!=m]', sep ='')))
      
      eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.m',model,'[best.model == model & true.model==model ]', sep ='')));
      eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.m',model,'[best.model == model & true.model!=model ]', sep ='')));
      
      if(models.p[m,p]>0)
      {
        points(ttrue, estimated.t, pch = 21, col = model, bg = m, cex = 0.6);
        points(ftrue, estimated.f, pch = 2, col = m, bg = m, cex = 0.6);
        
      }
      abline(a = 0, b = 1, col = 'gray', lwd=2.0)
      
    } 
  }
  
}

#compute.non.identifiability.parameter.region.by.model = function(T = T, alphas = c(10^-3, 10^-2, 10^-1), nb.realization = 50, model=3)
#{
  # alphas = c(10^-3, 10^-2, 10^-1); nb.realization = 50; model=3
#}

################################################################################### finishing line for the model selection methods

## clean the results of model selection:
cleaning.model.selection.results = function(T=Tt, parametrization =c('sigmoid','cosine'),absolute.signal = FALSE, model=4)
{
## T=Tt; absolute.signal = FALSE; model=4; i = 1;Nfit = NA; debug = FALSE; parametrization =c('sigmoid');method =  c('integration','simulation'); ZT.int = grep('.rel.ampl.int', colnames(T));ZT.ex = grep('.rel.ampl.ex', colnames(T));	zt = seq(0,46,by = 2)
	
	parametrization = parametrization[1];
	
	bounds = set.bounds(model = model, parametrization = parametrization, absolute.signal = absolute.signal); 
	upper = bounds$upper; 
	lower = bounds$lower;
	
	model = c(2,3,4)
	
	LRT.best.model = T$LRT.best.model
	AIC.best.model = T$AIC.best.model
	AICc.best.model =T$AICc.best.model
	BIC.best.model = T$BIC.best.model
	
	cutoff = 0.000001
	
	for(n in 1:nrow(T))
	{
		cat(n, '.....\n')
		bool.2 =	abs(T$gamma.m2[n]-upper[1])<cutoff|abs(T$gamma.m2[n]-lower[1])<cutoff|
		abs(T$fold.change.int.m2[n]-upper[4])<cutoff|abs(T$fold.change.int.m2[n]-lower[4])<cutoff|
		abs(T$phase.int.m2[n]-upper[5])<cutoff|	abs(T$phase.int.m2[n]-lower[5])<cutoff|
		abs(T$up.time.int.m2[n]-upper[6])<cutoff|abs(T$up.time.int.m2[n]-lower[6])<cutoff|
		abs(T$down.time.int.m2[n]-lower[7])<cutoff|abs(T$down.time.int.m2[n]-upper[7])<cutoff
		
		bool.3 =	abs(T$gamma.m3[n]-upper[1])<cutoff|abs(T$gamma.m3[n]-lower[1])<cutoff|
		abs(T$eps.gamma.m3[n]-upper[2])<cutoff|abs(T$eps.gamma.m3[n]-lower[2])<cutoff|
		abs(T$phase.gamma.m3[n]-upper[3])<cutoff|abs(T$phase.gamma.m3[n]-lower[3])<cutoff	
		
		bool.4 =abs(T$gamma.m4[n]-upper[1])<cutoff|abs(T$gamma.m4[n]-lower[1])<cutoff|
		abs(T$eps.gamma.m4[n]-upper[2])<cutoff|abs(T$eps.gamma.m4[n]-lower[2])<cutoff|
		abs(T$phase.gamma.m4[n]-upper[3])<cutoff|abs(T$phase.gamma.m4[n]-lower[3])<cutoff|
		abs(T$fold.change.int.m4[n]-upper[4])<cutoff|abs(T$fold.change.int.m4[n]-lower[4])<cutoff|
		abs(T$phase.int.m4[n]-upper[5])<cutoff|	abs(T$phase.int.m4[n]-lower[5])<cutoff|
		abs(T$up.time.int.m4[n]-upper[6])<cutoff|abs(T$up.time.int.m4[n]-lower[6])<cutoff|abs(T$down.time.int.m4[n]-lower[7])<cutoff|abs(T$down.time.int.m4[n]-upper[7])<cutoff
		
		#parameters.model = c(c('gamma.m2', 'fold.change.int.m2', 'phase.int.m2', 'up.time.int.m2', 'down.time.int.m2'), c('gamma.m3', 'eps.gamma.m3', 'phase.gamma.m3'), 
		#					 c('gamma.m4', 'eps.gamma.m4', 'phase.gamma.m4', 'fold.change.int.m4', 'phase.int.m4', 'up.time.int.m4', 'down.time.int.m4'))
		#index = match(parameters.model, colnames(T))
		
		#bool.2;bool.3;bool.4
		#T[n,index];
		#upper;lower
		#LRT.best.model[n];AIC.best.model[n];AICc.best.model[n]
		
		if(bool.2)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==2) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==2) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==2) AICc.best.model[n]=NA
			}
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==2) BIC.best.model[n]=NA
			}
		}
		
		
		
		if(bool.3)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==3) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==3) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==3) AICc.best.model[n]=NA
			}
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==3) BIC.best.model[n]=NA
			}
			
		}
		
		
		if(bool.4)
		{
			if(!is.na(LRT.best.model[n])) 
			{
				if(LRT.best.model[n]==4) LRT.best.model[n]=NA
			}
			if(!is.na(AIC.best.model[n])) 
			{
				if(AIC.best.model[n]==4) AIC.best.model[n]=NA
			}
			if(!is.na(AICc.best.model[n]))
			{
				if(AICc.best.model[n]==4) AICc.best.model[n]=NA
			}
			if(!is.na(BIC.best.model[n])) 
			{
				if(BIC.best.model[n]==4) BIC.best.model[n]=NA
			}
		}
		
	}
	T$LRT.best.model.eff = LRT.best.model
	T$AIC.best.model.eff = AIC.best.model
	T$AICc.best.model.eff = AICc.best.model
	T$BIC.best.model.eff = BIC.best.model
	return(T)
	
	
}

### using Z score to make some selection when extracting estimated parameters
plot.figure.fit.results.zscore = function(T = T, best.model = best.model, gamma = gamma, gamma.zscore = gamma.zscore, gamma.fc = gamma.fc, gamma.fc.zscore = gamma.fc.zscore, phase.gamma = phase.gamma,phase.gamma.score = phase.gamma.score, phase.int = phase.int, phase.int.score = phase.int.score, 
up.time.int = up.time.int, up.time.int.score = up.time.int.score, P = P)
{
	#T =T; best.model = best.model; gamma = gamma; gamma.fc = gamma.fc; phase.gamma = phase.gamma; phase.int = phase.int; P = P;
	
	layout(mat = matrix(c(1,1:9),nrow = 2, ncol = 5, byrow = TRUE), widths =c(0.65,0.55,1,1,1.05))
	#layout(mat = matrix(c(1:8),nrow = 2, ncol = 4, byrow = TRUE), widths =c(1.2,1,1,1.05))
	par(cex = 0.6, mgp = c(0,0.5,0),las = 1, tcl = -0.3, cex.main = 1, cex.axis = 0.9)
	
	par(mar = c(2,0,2,2))
	
	ylim = c(-0.1,1.15)
	
	h = hist(best.model,breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	H = matrix(c(h$density, 0,h$counts[2:4]/sum(h$counts[2:4])), 4,2)
	space = 1
	b = barplot(H, beside = FALSE, border = FALSE, col = 1:4, space = space, axes = FALSE, xlim = c(-6.2,4), ylim = ylim)
	lines(x = b + c(space/2,-space/2), y = H[1,], col = 'gray', lty = '22')
	lines(x = b + c(space/2,-space/2), y = rep(1,2), col = 'gray', lty = '22')
	axis(4, las = 1)
	text(x = rep(0.7,4), y = c(0.25,0.7,0.85,0.96), c('Constitutive transcription\n& constant degradation', 'Rhythmic transcription', 'Rhythmic degradation','Rhythmic transcr. & degr.'), col = c('gray50',2:4), pos = 2, cex = 0.9, font = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.7, 0.9), col = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.85,0.98), col = 3)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.97,1), col = 4)
	
	mtext('A',side =3, at = -6, font = 2, cex = 0.8, line = 0.8)
	
	
	fractions.ex = matrix(NA, nrow = 14, ncol = 4); 
	numbers.ex = rep(0,14)
	fractions.int = matrix(NA, nrow = 14, ncol = 4); 
	numbers.int = rep(0,14)
	for(r in 1:14)
	{
		ok = (log2(T$exon.median)>(r-1))&(log2(T$exon.median)<=r)
		fractions.ex[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.ex[r] = sum(ok, na.rm = TRUE)
		ok = (log2(T$intron.median)>(r-1))&(log2(T$intron.median)<=r)
		fractions.int[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.int[r] = sum(ok, na.rm = TRUE)
	}
	
	par(mgp = c(0,0.5,0), mar = c(2,0,2,0))
	
	b = barplot(t(fractions.ex), col = c(1:4), border = c(1:4), ylim = ylim, space = 0.0, xlab = 'exonic probes expression level', axes = FALSE)
	#text(b,rep(1.005,14), numbers.ex, pos = 4, cex = 0.7, srt = 90, offset = 0.1)
	pos = c(b[1]-diff(b)[1],b)+diff(b)[1]/2
	text(pos[seq(1,15,by = 2)],rep(0,15)[seq(1,15,by = 2)], c(0:14)[seq(1,15,by = 2)], pos = 1)
	
	numbers.ex
	frac.ex = numbers.ex/max(numbers.ex)*0.13
	frac.ex = frac.ex+1.02
	polygon(rep(seq(b[1]-(b[2]-b[1])/2,b[length(b)]+(b[2]-b[1])/2,by = b[2]-b[1] ),each = 2), c(1.02,rep(frac.ex,each=2),1.02), border = NA, col = 'gray40')
	#text(0,1.08,'# of mRNA per\nexpr. level', col = 'gray40',pos = 4, cex = 0.8)	
	mtext('# of mRNA per\nexpr. level', line = -0.4 ,side = 3, at = 3,col = 'gray40', cex = 0.5 )
	mtext('B',side =3, at = -0.5, font = 2, cex = 0.8, line = 0.8)
	
		
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')
	
	# half-lives 
	j2 = best.model == 2
	j3 = best.model == 3
	j4 = best.model == 4

	#j2 = which(best.model == 2 & log(2)/T$gamma.m2^2*T$gamma.stderr.m2<2)
	#j3 = which(best.model == 3 & log(2)/T$gamma.m3^2*T$gamma.stderr.m3<2)
	#j4 = which(best.model == 4 & log(2)/T$gamma.m4^2*T$gamma.stderr.m4<2)
	
	breaks = lseq(log(2)/max(gamma,na.rm = TRUE),log(2)/min(gamma,na.rm = TRUE), len = 31)
		
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0))
	
	h2 = hist(log(2)/gamma[j2], breaks = breaks, plot = FALSE)
	h3 = hist(log(2)/gamma[j3], breaks = breaks, plot = FALSE)
	h4 = hist(log(2)/gamma[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h2$breaks[1], h2$breaks, h2$breaks[length(h2$breaks)]),c(0, h2$counts, h2$counts[length(h2$counts)],0), log = 'x', type = 'n', main = 'Maximal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h2$counts, h3$count, h4$counts))
	polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	
	mtext('C',side =3, at = 0.12, font = 2, cex = 0.8, line = 0.8)
	
	
	
	# MINIMAL half-lives
		
	minimal.h.l = pmax(log(2)/gamma/gamma.fc,5/60)
	
	breaks = lseq(min(minimal.h.l[j3|j4],na.rm = TRUE),max(minimal.h.l[j3|j4],na.rm = TRUE), len = 31)
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0.5))
	
	h3 = hist(minimal.h.l[j3], breaks = breaks, plot = FALSE)
	h4 = hist(minimal.h.l[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Minimal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h3$count, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	mtext('D',side =3, at = 0.05, font = 2, cex = 0.8, line = 0.8)
	
	
	# degradation fold changes
	par(mar = c(3,2,3,0.5))
	
	breaks = lseq(1,max(gamma.fc,na.rm = TRUE), len = 31)
	
	h3 = hist(gamma.fc[j3], breaks = breaks, plot = FALSE)
	h4 = hist(gamma.fc[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Degr. FC',  xlab = 'fold changes', ylab = '', ylim = range(0,h3$counts, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	
	mtext('E',side =3, at = 0.6, font = 2, cex = 0.8, line = 1.7)
	
	
	# fold.changes in exons  ## it does not work at all with Nacho results....!
	###
	ind = c(1:(ncol(P)-1))
	
	FC = apply(P[,ind],1,max)/apply(P[,ind],1,min); FC = pmax(FC,1)
	boxplot(FC~best.model,log = 'y', pch = 18, col = 5:8, border = 1:4, xlim = c(1.5,4.5), main = 'mRNA FC', axes = FALSE, names = c('CS-CD','RS-CD','CS-RD','RS-RD')) # , ylim = c(1,10)
	box(); axis(2)
	axis(1, at = 2, labels = 'RS-CD', col.axis =2, las = 2);
	axis(1, at = 3, labels = 'CS-RD', col.axis =3, las = 2);
	axis(1, at = 4, labels = 'RS-RD', col.axis =4, las = 2);
	mtext('F',side =3, at = 1, font = 2, cex = 0.8, line = 1.7)
	
	
	# degradation phases
	par(mgp = c(1.6,0.5,0), mar = c(1,1,3,0.5))
	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
	h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE);
	h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	h3 = hist(phase.gamma[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Phase of max. degr.', axes = FALSE, plot = TRUE, ylim = range(0,h3$counts, h4$counts))
	h4 = hist(phase.gamma[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE, plot = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
	# factor.rose.diag = 0.75
	# h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE)
	# h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	
	# pg3 = circular(phase.gamma[j3], type = 'angles',units = 'hours',template = 'clock24')
	# pg4 = circular(phase.gamma[j4], type = 'angles',units = 'hours',template = 'clock24')
	# rose.diag(pg3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main =  'Phase of max. degr.', asp = 1)
	# points.circular(pg3, bins = 24, col = 3+4)
	# rose.diag(pg4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
	# points.circular(pg4, bins = 24, col = 4+4)
	
	
	mtext('G',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
	
	
	par(mar = c(3,2,3,0.5))
	
	phase.ex = apply(P[,ind],1,which.max)/10
	diff.phase = (phase.ex-phase.gamma)%%24
	#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Delay between\nmax. mRNA & degr.', axes = FALSE)
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
	# h3 = hist(diff.phase[j3], breaks = breaks,  plot = FALSE)
	# h4 = hist(diff.phase[j4], breaks = breaks, plot = FALSE)
	# diff3 = circular(diff.phase[j3], type = 'angles',units = 'hours',template = 'clock24')
	# diff4 = circular(diff.phase[j4], type = 'angles',units = 'hours',template = 'clock24')
	# rose.diag(diff3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main = 'Delay between\nmax. mRNA & degr.', asp = 1)
	# points.circular(diff3, bins = 24, col = 3+4)
	# rose.diag(diff4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
	# points.circular(diff4, bins = 24, col = 4+4)
	
	mtext('H',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
	diff.phase = (phase.int-phase.gamma)%%24
	#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, xlab = '[h]', ylab = '', main = 'Delay between\nmax. pre-mRNA & degr.', axes = FALSE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	mtext('I',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
}


plot.figure.fit.results = function(T = T, best.model = T$model.MCMC, gamma = T$gamma.MCMC, gamma.fc = T$eps.gamma.MCMC, phase.gamma = T$phase.gamma.MCMC, phase.int = T$phase.int.MCMC, up.time.int = up.time.int, P = P)
{
	#T =T; best.model = best.model; gamma = gamma; gamma.fc = gamma.fc; phase.gamma = phase.gamma; phase.int = phase.int; P = P;
	
	layout(mat = matrix(c(1,1:9),nrow = 2, ncol = 5, byrow = TRUE), widths =c(0.65,0.55,1,1,1.05))
	#layout(mat = matrix(c(1:8),nrow = 2, ncol = 4, byrow = TRUE), widths =c(1.2,1,1,1.05))
	par(cex = 0.6, mgp = c(0,0.5,0),las = 1, tcl = -0.3, cex.main = 1, cex.axis = 0.9)
	
	par(mar = c(2,0,2,2))
	
	ylim = c(-0.1,1.15)
	
	h = hist(best.model,breaks = seq(0.5,4.5,by = 1), plot = FALSE)
	H = matrix(c(h$density, 0,h$counts[2:4]/sum(h$counts[2:4])), 4,2)
	space = 1
	b = barplot(H, beside = FALSE, border = FALSE, col = 1:4, space = space, axes = FALSE, xlim = c(-6.2,4), ylim = ylim)
	lines(x = b + c(space/2,-space/2), y = H[1,], col = 'gray', lty = '22')
	lines(x = b + c(space/2,-space/2), y = rep(1,2), col = 'gray', lty = '22')
	axis(4, las = 1)
	text(x = rep(0.7,4), y = c(0.25,0.7,0.85,0.96), c('Constitutive transcription\n& constant degradation', 'Rhythmic transcription', 'Rhythmic degradation','Rhythmic transcr. & degr.'), col = c('gray50',2:4), pos = 2, cex = 0.9, font = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.7, 0.9), col = 2)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.85,0.98), col = 3)
	lines(x = c(0.7 , b[1]-space/2), y =c(0.97,1), col = 4)
	
	mtext('A',side =3, at = -6, font = 2, cex = 0.8, line = 0.8)
	
	
	fractions.ex = matrix(NA, nrow = 14, ncol = 4); 
	numbers.ex = rep(0,14)
	fractions.int = matrix(NA, nrow = 14, ncol = 4); 
	numbers.int = rep(0,14)
	for(r in 1:14)
	{
		ok = (log2(T$exon.median)>(r-1))&(log2(T$exon.median)<=r)
		fractions.ex[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.ex[r] = sum(ok, na.rm = TRUE)
		ok = (log2(T$intron.median)>(r-1))&(log2(T$intron.median)<=r)
		fractions.int[r,] =  hist( best.model[ok], breaks = seq(0.5,4.5,by=1),plot = FALSE)$density
		numbers.int[r] = sum(ok, na.rm = TRUE)
	}
	
	par(mgp = c(0,0.5,0), mar = c(2,0,2,0))
	
	b = barplot(t(fractions.ex), col = c(1:4), border = c(1:4), ylim = ylim, space = 0.0, xlab = 'exonic probes expression level', axes = FALSE)
	#text(b,rep(1.005,14), numbers.ex, pos = 4, cex = 0.7, srt = 90, offset = 0.1)
	pos = c(b[1]-diff(b)[1],b)+diff(b)[1]/2
	text(pos[seq(1,15,by = 2)],rep(0,15)[seq(1,15,by = 2)], c(0:14)[seq(1,15,by = 2)], pos = 1)
	
	numbers.ex
	frac.ex = numbers.ex/max(numbers.ex)*0.13
	frac.ex = frac.ex+1.02
	polygon(rep(seq(b[1]-(b[2]-b[1])/2,b[length(b)]+(b[2]-b[1])/2,by = b[2]-b[1] ),each = 2), c(1.02,rep(frac.ex,each=2),1.02), border = NA, col = 'gray40')
	#text(0,1.08,'# of mRNA per\nexpr. level', col = 'gray40',pos = 4, cex = 0.8)	
	mtext('# of mRNA per\nexpr. level', line = -0.4 ,side = 3, at = 3,col = 'gray40', cex = 0.5 )
	mtext('B',side =3, at = -0.5, font = 2, cex = 0.8, line = 0.8)
	
	
	bounds = set.bounds(model = 4, parametrization = 'sigmoid')
	
	cutoff = 1.036
	
	# half-lives 
	j2 = best.model == 2
	j3 = best.model == 3
	j4 = best.model == 4
	
	#j2 = which(best.model == 2 & log(2)/T$gamma.m2^2*T$gamma.stderr.m2<2)
	#j3 = which(best.model == 3 & log(2)/T$gamma.m3^2*T$gamma.stderr.m3<2)
	#j4 = which(best.model == 4 & log(2)/T$gamma.m4^2*T$gamma.stderr.m4<2)
	
	
	breaks = lseq(log(2)/max(gamma,na.rm = TRUE),log(2)/min(gamma,na.rm = TRUE), len = 31)
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0))
	
	h2 = hist(log(2)/gamma[j2[which(gamma.zscore[j2]>cutoff)]], breaks = breaks, plot = FALSE)
	h3 = hist(log(2)/gamma[j3[which(gamma.zscore[j3]>cutoff)]], breaks = breaks, plot = FALSE)
	h4 = hist(log(2)/gamma[j4[which(gamma.zscore[j4]>cutoff)]], breaks = breaks, plot = FALSE)
	
	plot(c(h2$breaks[1], h2$breaks, h2$breaks[length(h2$breaks)]),c(0, h2$counts, h2$counts[length(h2$counts)],0), log = 'x', type = 'n', main = 'Maximal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h2$counts, h3$count, h4$counts))
	polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = 2+4, border = 2, lwd = 1)
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	
	mtext('C',side =3, at = 0.12, font = 2, cex = 0.8, line = 0.8)
	
	
	
	# MINIMAL half-lives
	
	minimal.h.l = pmax(log(2)/gamma/gamma.fc,5/60)
	
	breaks = lseq(min(minimal.h.l[j3|j4],na.rm = TRUE),max(minimal.h.l[j3|j4],na.rm = TRUE), len = 31)
	
	par(mgp = c(1.6,0.5,0), mar = c(3,2,2,0.5))
	
	h3 = hist(minimal.h.l[j3], breaks = breaks, plot = FALSE)
	h4 = hist(minimal.h.l[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Minimal half-lives',  xlab = 'half-lives [h]', ylab = '', ylim = range(h3$count, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	mtext('D',side =3, at = 0.05, font = 2, cex = 0.8, line = 0.8)
	
	
	
    # degradation fold changes
	par(mar = c(3,2,3,0.5))
	
	breaks = lseq(1,max(gamma.fc,na.rm = TRUE), len = 31)
	
	h3 = hist(gamma.fc[j3], breaks = breaks, plot = FALSE)
	h4 = hist(gamma.fc[j4], breaks = breaks, plot = FALSE)
	
	plot(c(h3$breaks[1],h3$breaks,h3$breaks[length(h3$breaks)]),c(0,h3$counts,h3$counts[length(h3$counts)],0), log = 'x', type = 'n', main = 'Degr. FC',  xlab = 'fold changes', ylab = '', ylim = range(0,h3$counts, h4$counts))
	polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$counts, each=2),0), col = 3+4, border = 3, lwd = 1)
	polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = 4+4, border = 4, lwd = 1)
	
	mtext('E',side =3, at = 0.6, font = 2, cex = 0.8, line = 1.7)
	
	
    # fold.changes in exons  ## it does not work at all with Nacho results....!
    ###
	ind = c(1:(ncol(P)-1))
	
	FC = apply(P[,ind],1,max)/apply(P[,ind],1,min); FC = pmax(FC,1)
	boxplot(FC~best.model,log = 'y', pch = 18, col = 5:8, border = 1:4, xlim = c(1.5,4.5), main = 'mRNA FC', axes = FALSE, names = c('CS-CD','RS-CD','CS-RD','RS-RD')) # , ylim = c(1,10)
	box(); axis(2)
	axis(1, at = 2, labels = 'RS-CD', col.axis =2, las = 2);
	axis(1, at = 3, labels = 'CS-RD', col.axis =3, las = 2);
	axis(1, at = 4, labels = 'RS-RD', col.axis =4, las = 2);
	mtext('F',side =3, at = 1, font = 2, cex = 0.8, line = 1.7)
	
	
    # degradation phases
    ##par(mgp = c(1.6,0.5,0), mar = c(1,1,3,0.5))
	breaks = seq(bounds$upper[3],bounds$lower[3], len = 25)
	h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE);
	h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	h3 = hist(phase.gamma[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Phase of max. degr.', axes = FALSE, plot = TRUE, ylim = range(0,h3$counts, h4$counts))
	h4 = hist(phase.gamma[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE, plot = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
    # factor.rose.diag = 0.75
    # h3 = hist(phase.gamma[j3], breaks = breaks, plot = FALSE)
    # h4 = hist(phase.gamma[j4], breaks = breaks, plot = FALSE)
	
# pg3 = circular(phase.gamma[j3], type = 'angles',units = 'hours',template = 'clock24')
# pg4 = circular(phase.gamma[j4], type = 'angles',units = 'hours',template = 'clock24')
# rose.diag(pg3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main =  'Phase of max. degr.', asp = 1)
# points.circular(pg3, bins = 24, col = 3+4)
# rose.diag(pg4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
# points.circular(pg4, bins = 24, col = 4+4)
	
	
	mtext('G',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
	
	
	par(mar = c(3,2,3,0.5))
	
	phase.ex = apply(P[,ind],1,which.max)/10
	diff.phase = (phase.ex-phase.gamma)%%24
#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j3], breaks = breaks, col = 3+4, border = 3, xlab = '[h]', ylab = '', main = 'Delay between\nmax. mRNA & degr.', axes = FALSE)
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, add = TRUE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	
# h3 = hist(diff.phase[j3], breaks = breaks,  plot = FALSE)
# h4 = hist(diff.phase[j4], breaks = breaks, plot = FALSE)
# diff3 = circular(diff.phase[j3], type = 'angles',units = 'hours',template = 'clock24')
# diff4 = circular(diff.phase[j4], type = 'angles',units = 'hours',template = 'clock24')
# rose.diag(diff3, bins = 24, col = 3+4, border = 3, prop = 1/sqrt(max(h3$density))* factor.rose.diag, main = 'Delay between\nmax. mRNA & degr.', asp = 1)
# points.circular(diff3, bins = 24, col = 3+4)
# rose.diag(diff4, bins = 24, col = 4+4, border = 4, prop = 1/sqrt(max(h3$density))* factor.rose.diag, add = TRUE)
# points.circular(diff4, bins = 24, col = 4+4)
	
	mtext('H',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
	diff.phase = (phase.int-phase.gamma)%%24
#diff.phase[which(diff.phase>12)] = diff.phase[which(diff.phase>12)]-24
	h = hist(diff.phase[j4], breaks = breaks, col = 4+4, border = 4, xlab = '[h]', ylab = '', main = 'Delay between\nmax. pre-mRNA & degr.', axes = FALSE)
	axis(1, at = seq(0,24,by = 6)); axis(2); box()
	mtext('I',side =3, at = -1, font = 2, cex = 0.8, line = 1.7)
	
}

plot.this.gene.full = function(pdf.folder = '',gene = 'Dbp', ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE, with.fit = FALSE, fitting.method = c('MCMC','Optim','both'), parametrization = c('sigmoid','cosine'))
{
    ylim.rel.type = ylim.rel.type[1]; parametrization = parametrization[1]; fitting.method = fitting.method[1]
    
    #global variable
    ylim.default <<- c(0,2.2); ylim.default.log2 <<- c(log2(0.1),log2(4))
    attribute.global.variable.colors()
    
    # create pdf and plots
    pdfname = paste(pdf.folder,'/',gene,'.pdf',sep = '')
    nrow = 4; ncol = 6; dim = 3
    pdf(pdfname, width = ncol/3*dim *1.7 , height = nrow* dim)
    if(fitting.method == 'both'){layout(matrix(c(rep(1,6),rep(2:5,each = 3), rep(6:8,each=2)),nrow = nrow,ncol = ncol, byrow = TRUE), widths = rep(1,ncol), heights = c(1.5,1,1,1))}
    else{layout(matrix(c(rep(1,6),rep(2:7,each = 3)),nrow = nrow,ncol = ncol, byrow = TRUE), widths = rep(1,ncol), heights = c(1.5,1,1,1))}
    par(mar = c(0.5, 4, 2, 0.5) + 0.1, mgp = c(2, 0, 0)); make.UCSC.plot( gene = gene) ; cat('\tUCSC plot done \n');
    par(mar = c(3, 3, 2, 0.5) + 0.1, mgp = c(2, 0.7, 0)); show.all.probesets(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel); cat('\tprobesets plots done \n')
    show.summarized.signal(gene = gene, ylim.rel.type = ylim.rel.type, ylim.rel = ylim.rel, log2.rel = log2.rel, with.fit = with.fit, fitting.method = fitting.method, parametrization = parametrization); cat('\tsummary plots done \n')
    
    dev.off()
}

###########
##### Show fitting results for individual examples
###########
show.summarized.signal = function(gene = gene, ylim.rel.type = c('fixed','adjusted','personalized'), ylim.rel = c(0,2), log2.rel = FALSE, with.cosine.fit = FALSE, with.fit = FALSE, fitting.method = 'MCMC', parametrization = c('sigmoid','cosine'), abline = FALSE, dt = 4, title = 'default'){
    ylim.rel.type = ylim.rel.type[1]; parametrization = parametrization[1]
    i = T$gene == gene;
    if(sum(i)>0){
        ex = unlist(T[i,grep('.rel.ampl.ex', colnames(T))]); int = unlist(T[i,grep('.rel.ampl.int', colnames(T))]);hline = 1; ylab = 'relative signal [lin scale]'
        if(log2.rel){ex = log2(ex);hline = 0; ylab = 'relative signal [log2 scale]'};
        
        if(ylim.rel.type == 'fixed'){ylim = ylim.default; if(log2.rel){ylim = ylim.default.log2}}else if(ylim.rel.type == 'adjusted'){ylim = range(ex,int)}else{ylim = ylim.rel}
        if(title == 'default'){title = 'summarized signal'; if(with.cosine.fit){title = paste(title,'+ cos. fits (2 comp.)')}; if(with.fit){}}
        plot(1,1, type = 'n', xlab = 'time (hours)', ylab = ylab, axes = FALSE, xlim = range(zt), ylim = ylim); box();
        time.axis(dt = dt); axis(2)
        if(abline){abline(h = 1, col = 'gray')}
        points(zt, ex, col = col.ex, type = 'p', lwd = 2, pch = 16); points(zt, int, col = col.int, type = 'p', lwd = 1, pch = 5)
        if(with.cosine.fit){
            zt.p = seq(0,48,by = 0.05)
            cos.ex = T$mean.ex[i] + min(1,T$rel.ampl.ex[i]) * cos(2*pi*(zt.p-T$phase.ex[i])/24) +  min(1,T$rel.ampl.12.ex[i]) * cos(2*pi*(zt.p-T$phase.12.ex[i])/12)
            cos.int = T$mean.int[i] + min(1,T$rel.ampl.int[i]) * cos(2*pi*(zt.p-T$phase.int[i])/24) +  min(1,T$rel.ampl.12.int[i]) * cos(2*pi*(zt.p-T$phase.12.int[i])/12)
            points(zt.p, cos.ex, type = 'l', lwd = 2, col = col.ex);points(zt.p, cos.int, type = 'l', lwd = 2, col = col.int)
            points(zt, ex, col = col.ex, type = 'l', lwd = 1);points(zt, int, col = col.int, type = 'l', lwd = 0.5)
        }else if(with.fit){
            zt.p = seq(0,48,by = 0.05)
            
            gamma = 0; eps.gamma = 0; phase.gamma = 0;
            rel.ampl.int = 0; phase.int = 0; rel.ampl.12.int = 0; phase.12.int = 0
            fold.change = 1; up.time = 12; down.time = 12;
            
            if((fitting.method == 'MCMC')&!is.null(T$model.MCMC[i]) && !is.na(T$model.MCMC[i])){
                best.model = T$model.MCMC[i]
                gamma = T$gamma.MCMC[i]; eps.gamma = T$eps.gamma.MCMC[i]; phase.gamma = T$phase.gamma.MCMC[i];
                rel.ampl.int = T$rel.ampl.int.MCMC[i]; phase.int = T$phase.int.MCMC[i]; rel.ampl.12.int = T$rel.ampl.12.int.MCMC[i]; phase.12.int = T$phase.12.int.MCMC[i]
                if(title == 'default'){title = paste(title,'-',fitting.method)}
                if(title == 'model'){title = c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]}
            }else if(!is.na(T$LRT.best.model[i])){
                best.model = T$LRT.best.model[i]; if(!is.null(T$LRT.eff.best.model[i])){best.model = T$LRT.eff.best.model[i]}
                eval(parse(text =paste('gamma = T$gamma.m',best.model,'[i]',sep = '')));eval(parse(text =paste('eps.gamma = T$eps.gamma.m',best.model,'[i]',sep = '')));eval(parse(text =paste('phase.gamma = T$phase.gamma.m',best.model,'[i]',sep = '')));
                if(parametrization == 'cosine'){
                    eval(parse(text = paste('rel.ampl.int = T$rel.ampl.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; rel.ampl.12.int = T$rel.ampl.12.int.m',best.model,'[i]; phase.12.int = T$phase.12.int.m',best.model,'[i]',sep = '')))}
                else{eval(parse(text = paste( 'fold.change = T$fold.change.int.m',best.model,'[i]; phase.int = T$phase.int.m',best.model,'[i]; up.time = T$up.time.int.m', best.model,'[i]; down.time = T$down.time.int.m', best.model,'[i]', sep = '')))}
                if(title == 'default'){title = paste(title,'- Optim')}
                if(title == 'model'){title = c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]}
            }
            
            cat('is.null(fold.change) :', is.null(fold.change),'\n')
            
            if(is.null(fold.change)){fold.change = 1; phase.int = 0; up.time = 12; down.time = 12}
            if(is.null(eps.gamma)){if(parametrization == 'cosine'){eps.gamma = 0; phase.gamma = 0}else{eps.gamma = 1; phase.gamma = 0}}
            
            
            
            cat('parametrization :', parametrization, '\n gamma = ',gamma,'\n eps.gamma = ',eps.gamma, '\n phase.gamma = ',phase.gamma,'\n fold.change = ',fold.change, '\n phase.int = ',phase.int, '\n up.time = ',up.time, '\n down.time = ', down.time, '\n' );
            real.gamma = gamma
            gamma = min(gamma, 5)
            
            profile.int = rep(1,length(zt.p));
            if(best.model > 1){
                if(parametrization == 'cosine'){profile.int = compute.s(t = zt.p, eps.24.S = rel.ampl.int, phase.24.S = phase.int, eps.12.S = rel.ampl.12.int, phase.12.S = phase.12.int  )}
                else{profile.int = compute.sigmoid(t = zt.p, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)}
            }
            
            profile.deg = rep(1,length(zt.p));
            if(best.model > 2){
                if(parametrization == 'cosine'){profile.deg = cos.deg + eps.gamma * cos(2*pi*(zt.p-phase.gamma)/24)}
                else{profile.deg = compute.sigmoid(t = zt.p, fold.change = eps.gamma, phase = phase.gamma, up.time = 12, down.time = 12)}
            }
            
            
            if(best.model == 1){profile.ex = rep(1,length(zt.p))}
            if(best.model > 1){
                if(parametrization == 'cosine'){profile.ex = compute.m.changing(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, eps.24.S = rel.ampl.int, phase.24.S = phase.int, eps.12.S = rel.ampl.12.int, phase.12.S = phase.12.int)}
                else{profile.ex = compute.m.sigmoid(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, fold.change = fold.change, phase = phase.int, up.time = up.time, down.time = down.time)}
            }
            
            #cat('length(zt.p) = ',length(zt.p), '\t length(profile.ex) = ',length(profile.ex),'\t length(profile.int) = ',length(profile.int),'\t length(profile.deg) = ',length(profile.deg),'\n' )
            
            points(zt.p, profile.ex, type = 'l', lwd = 2, col = col.ex);points(zt.p, profile.int, type = 'l', lwd = 2, col = col.int); points(zt.p, profile.deg, type = 'l', lwd = 2, col = col.deg);
            
            
            if(best.model>1){legend('topright',legend = paste('half-life =',round(log(2)/real.gamma, digits = 2),'h'), text.col = col.deg, bty = 'n', lty = 0)}
            if(title == 'default'){legend('topleft',legend = paste('Model =',c('CS-CD','RS-CD','CS-RD','RS-RD')[best.model]), text.col = c('gray','green4','tomato','black')[best.model], bty = 'n', lty = 0)}
            points(zt, ex, col = col.ex, type = 'l', lwd = 1);points(zt, int, col = col.int, type = 'l', lwd = 0.5)
            
        }else{points(zt, ex, col = col.ex, type = 'l', lwd = 2);points(zt, int, col = col.int, type = 'l', lwd = 2)}
    }else{plot(1,1,type = 'n', axes = FALSE, xlab = '', ylab = ''); text(1,1, 'no summarized signal for this gene')}
    col.title = 'black'; if(any(title == c('CS-CD','RS-CD','CS-RD','RS-RD'))){col.title = best.model}
    mtext(title, side = 3, line = 1, cex = 0.8, font = 2, col = col.title)
    
}
time.axis = function(dt = 4){axis(1, at = zt, labels = rep('',length(zt))); axis(1, at = seq(0,48,by = dt))}

###########################
###########################
####### PLots for Results
###########################
###########################
Transform.Scoring.parameters = function(T = T)
{
  kk = which(!is.na(T$BIC.best.model>1)==TRUE)
  keep = matrix(data=NA, nrow=length(kk), ncol=26)
  colnames(keep) = c('gene', 'index', 'best.model', 'prob.best.model', 'non.identifiability.gamma.L', 'non.identifiability.gamma.R',
                     'half.life', 'cv.half.life', 'splicing.time', 'cv.splicing.time',
                     'eps.gamma', 'cv.eps.gamma', 'phase.gamma', 'sem.phase.gamma', 
                     'Min.int', 'cv.Min.int','Amp.int', 'cv.Amp.int', 'phase.int', 'sem.phase.int', 'beta.int', 'cv.beta.int',
                     'Fried', 'Shar', 'Schwa', 'score.params')
  #keep[,4] = T$BIC.prob.[kk]
  keep = data.frame(T$gene[kk], keep[,-1], stringsAsFactors=FALSE);colnames(keep)[1] = 'gene'
  keep$index = kk; keep$best.model = T$BIC.best.model[kk];
  
  load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/mRNAs_half_lives_databases.Rdata')
  parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 'gamma','eps.gamma','phase.gamma','splicing.k')
  models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), c(0,2,3,4),c(0,0,3,4),c(0,0,3,4),c(0,2,3,4))
  
  ### extract parameter estimation from the fitting result
  for(n in 1:nrow(keep))
  #for(n in c(1:10))
  {
    #n = 1
    cat(n, '\n')
    model = keep$best.model[n]; index = keep$index[n];
    eval(parse(text = paste('keep$prob.best.model[n] = T$BIC.prob.m', model, '[index]', sep=''))); ## prob of best model
   
    if(model>1){
      ## analysis of non-identifiability
      eval(parse(text = paste('keep$non.identifiability.gamma.L[n] = T$non.identifiability.gamma.L.m', model, '[index]', sep=''))); 
      eval(parse(text = paste('keep$non.identifiability.gamma.R[n] = T$non.identifiability.gamma.R.m', model, '[index]', sep='')));
      ### 8 parameters
      w = 2*pi/24;
      for(n.p in 1:length(parameter.list))
      {
        if(models.p[model, n.p]>0)
        {
          # n.p = 8;
          par = parameter.list[n.p];
          
          if(n.p<=4){
            eval(parse(text = paste('keep$', par, '[n] = T$', par, '.m', model, '[index]', sep='')));
            if(n.p!=3){
              eval(parse(text = paste('keep$cv.', par, '[n] = as.numeric(T$', par, '.stderr.m', model, '[index])/T$', par, '.m', model, '[index]', sep='')));
            }else{
              ## standard error for estimated pre-mRNA phase
              eval(parse(text = paste('keep$sem.', par, '[n] = T$', par, '.stderr.m', model, '[index]', sep = '')));
            }
          }else{
            gamma = c(); delta.gamma=c(); eval(parse(text = paste('gamma = T$gamma.m', model, '[index]', sep='')));
            eval(parse(text = paste('delta.gamma = T$gamma.stderr.m', model, '[index]', sep='')));
            if(n.p==5){keep$half.life[n] = log(2)/gamma; keep$cv.half.life[n] = delta.gamma/gamma;}
            if(n.p==6){
              eval(parse(text = paste('epsilon = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.epsilon = T$', par, '.stderr.m', model, '[index]', sep='')));
              Y = epsilon*sqrt(w^2/gamma^2+1); 
              delta.Y = sqrt((sqrt(w^2/gamma^2+1)*delta.epsilon)^2 + (epsilon*w^2/gamma^3*delta.gamma/sqrt(w^2/gamma^2+1))^2);
              keep$eps.gamma[n] = Y; keep$cv.eps.gamma[n] = delta.Y/Y;
            }
            if(n.p==7){
              eval(parse(text = paste('phi = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.phi = T$', par, '.stderr.m', model, '[index]', sep='')));
              ## standard error for estimation of phase.gamma
              Y = (phi - atan2(w, gamma)/w)%%24; 
              delta.Y = sqrt((delta.phi)^2 + (delta.gamma/(gamma^2+w^2))^2);
              keep$phase.gamma[n] = Y; 
              keep$sem.phase.gamma[n] = delta.Y;
            }
            if(n.p==8){
              eval(parse(text = paste('a = T$', par, '.m', model, '[index]', sep='')));
              eval(parse(text = paste('delta.a = T$', par, '.stderr.m', model, '[index]', sep='')));
              Y = log(2)/a/gamma*60; delta.Y = log(2)*60*sqrt((delta.gamma/a/gamma^2)^2 + (delta.a/a^2/gamma)^2)
              keep$splicing.time[n] = Y; keep$cv.splicing.time[n] = delta.Y/Y;
            }
          }
        }
      } 
    }
    ## scoring the estimation of parameters
    #if(length(test)==length(which(!is.na(test)))) keep$score.params[n] = max(test)
  }
  ### add half-lived measured from other data sets
  mm = match(keep$gene, fried[,1]); keep$Fried = as.numeric(fried[mm, 2]);
  mm = match(keep$gene, shar[,1]); keep$Shar = as.numeric(shar[mm, 2]);
  mm = match(keep$gene, schwa[,1]); keep$Schwa = as.numeric(schwa[mm, 2]);
  
  return(keep);
}

scoring.estimation.params = function(test)
{
  # test = keep[1, ];
  test.pars1 = as.numeric(test[c(grep('cv.', names(test)))]);
  test.pars2 = as.numeric(test[c(grep('sem.', names(test)))])/2;
  test.pars = c(test.pars1, test.pars2);
  return(mean(as.numeric(test.pars), na.rm = TRUE));
}

plot.summary.splicing.half.life = function(folder='myplots/Total/', T=T)
{
    cols=c('bisque','olivedrab','palevioletred3','purple4')
    
    #### Extract results from Table T
    kk = which(T$BIC.best.model>1)
    keep = matrix(data=NA, nrow=length(kk), ncol=24)
    colnames(keep) = c('gene', 'index', 'best.model', 'prob.best.model', 'splicing.time', 'cv.splicing.time',
                        'half.life', 'cv.half.life', 'Fried', 'Shar', 'Schwa',
                        'gamma.phase', 'cv.gamma.phase', 'gamma.fc', 'cv.gamma.fc',
                        'phase.premRNA', 'amp.premRNA', 'pval.premRNA', 'qv.premRNA',
                        'phase.mRNA', 'amp.mRNA', 'pval.mRNA', 'qv.mRNA', 'score.params')
    keep[,1] = T$gene[kk]
    keep[,2] = kk
    keep[,3] = T$BIC.best.model[kk]
    keep[,4] = T$prob.best.model[kk]
    keep[,16] = T$phase.premRNA[kk]
    keep[,17] = T$amp.premRNA[kk]
    keep[,18] = T$pval.premRNA[kk]
    keep[,19] = T$qv.premRNA[kk]
    keep[,20] = T$phase.mRNA[kk]
    keep[,21] = T$amp.mRNA[kk]
    keep[,22] = T$pval.mRNA[kk]
    keep[,23] = T$qv.mRNA[kk]
    keep = data.frame(T$gene[kk], keep[,-1], stringsAsFactors=FALSE)
    colnames(keep)[1] = 'gene'
    
    load(file='myRdata/mRNAs_half_lives_databases.Rdata')
    
    for(n in 1:length(kk))
    {
        cat(n, '\n')
        test = c()
        best.model = T$BIC.best.model[kk[n]];
        
        ### splicing time
        xx = eval(parse(text = paste('log(2)/T$splicing.k.m', best.model, '[kk[n]]', sep='')))
        yy = eval(parse(text = paste('log(2)/(T$splicing.k.m', best.model, '[kk[n]])^2*T$splicing.k.stderr.m', best.model, '[kk[n]]', sep='')))
        yy = yy/xx
        keep$splicing.time[n] = xx*60 ## in mins
        keep$cv.splicing.time[n] = yy
        test = c(test, yy)
        
        ### half-life
        model = best.model
        if(model==2)
        {
            xx = eval(parse(text=paste('log(2)/T$gamma.m', model, '[kk[n]]', sep='')))
            yy = eval(parse(text=paste('log(2)*T$gamma.stderr.m', model, '[kk[n]]/(T$gamma.m', model, '[kk[n]])^2', sep='')))
        }else{
            xx = eval(parse(text=paste('log(2)/(T$gamma.m', model, '[kk[n]]+0.5*T$amp.gamma.m', model, '[kk[n]])', sep='')))
            yy = eval(parse(text=paste('log(2)*T$gamma.stderr.m', model, '[kk[n]]/(T$gamma.m', model, '[kk[n]]+0.5*T$amp.gamma.m', model, '[kk[n]])^2+log(2)*0.5*T$amp.gamma.stderr.m', model, '[kk[n]]/(T$gamma.m', model, '[kk[n]]+0.5*T$amp.gamma.m', model, '[kk[n]])^2', sep='')))
        }
        yy = yy/xx
        test = c(test, yy)
    
        keep$half.life[n] = xx;
        keep$cv.half.life[n] = yy;
        gg = T$gene[kk[n]];
        mm = match(gg, fried[,1])
        keep$Fried[n] = as.numeric(fried[mm, 2])
        mm = match(gg, shar[,1])
        keep$Shar[n] = as.numeric(shar[mm, 2])
        mm = match(gg, schwa[,1])
        keep$Schwa[n] = as.numeric(schwa[mm, 2])
        
        ### phases and fold-changes of rhythmic degradation
        if(model>2)
        {
            xx = eval(parse(text=paste('T$phase.gamma.m', model, '[kk[n]]', sep='')))
            yy = eval(parse(text=paste('T$phase.gamma.stderr.m', model, '[kk[n]]', sep='')))
            yy = yy/4;
            test = c(test, yy)
            keep$gamma.phase[n] = xx
            keep$cv.gamma.phase[n] = yy
            
            xx = eval(parse(text=paste('(T$gamma.m', model, '[kk[n]]+T$amp.gamma.m', model, '[kk[n]])/T$gamma.m', model, '[kk[n]]', sep='')))
            yy = eval(parse(text=paste('T$amp.gamma.m', model, '[kk[n]]*T$gamma.stderr.m', model, '[kk[n]]/(T$gamma.m', model, '[kk[n]])^2 + T$amp.gamma.stderr.m', model, '[kk[n]]/T$gamma.m', model, '[kk[n]]', sep='')))
            yy = yy/xx
            test = c(test, yy)
            keep$gamma.fc[n] = xx
            keep$cv.gamma.fc[n] = yy
        }
        if(model==2 | model==4)
        {
            test = c(test, eval(parse(text=paste('T$Min.int.stderr.m', model, '[kk[n]]', sep='')))/eval(parse(text=paste('T$Min.int.m', model, '[kk[n]]', sep=''))),
                           eval(parse(text=paste('T$Amp.int.stderr.m', model, '[kk[n]]', sep='')))/eval(parse(text=paste('T$Amp.int.m', model, '[kk[n]]', sep=''))),
                           eval(parse(text=paste('T$phase.int.stderr.m', model, '[kk[n]]', sep='')))/4,
                           eval(parse(text=paste('T$beta.int.stderr.m', model, '[kk[n]]', sep='')))/eval(parse(text=paste('T$beta.int.m', model, '[kk[n]]', sep=''))))
        }
        
        #test = test[which(!is.na(test))]
        if(length(test)==length(which(!is.na(test)))) keep$score.params[n] = max(test)
    }
    
    
    Save = FALSE
    if(Save)
    {
        save(keep, file ='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Total_models_parameters_results.Rdata')
        
        xx = data.frame(keep$gene, keep$best.model, keep$half.life, keep$gamma.fc, keep$gamma.phase, keep$score.params, stringsAsFactors=FALSE)
        kk = which(xx[,2]>2)
        xx = xx[kk,]
        colnames(xx) = c('gene', 'best.model', 'gamma', 'gamma.amp', 'gamma.phase', 'score.params')
        xx$gamma = log(2)/xx$gamma
        xx$gamma.amp = log2(xx$gamma.amp)
        
        write.table(xx, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Tables/Rhythmic_Decay_Amp_Phase.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        ss = kkeep
        kk = which(ss[,3]>2)
        ss = ss[kk, c(1, 3, 2, 7)]
        ss = ss[order(ss[,2]),]
        colnames(ss) = c('gene', 'best.model', 'half.life', 'hl.score')
        write.table(ss, file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/model_RNA_seq_total/Tables/Model_selection_results_genes_rhythmic_decay.txt', sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
    
    jj = which(keep$prob.best.model>0.6 & keep$score.params<0.5)
    length(which(keep$best.model[jj]==2))/length(jj)
    length(which(keep$best.model[jj]==3))/length(jj)
    length(which(keep$best.model[jj]==4))/length(jj)
    
    #### Plot results
    #plot(keep$splicing.time, keep$cv.splicing.time)
    jj = which(keep$prob.best.model>0.6 & keep$cv.splicing.time<0.5)
    #kk = which(keep$prob.best.model>0.6 & keep$score.params<0.5)
    
    pdfname = paste(folder, 'Splicing_time_distribution.pdf', sep='')
    pdf(pdfname, width=2.5, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    breaks = 40;
    lty = 2;
    hist(keep$splicing.time[jj], col='gold3', breaks=breaks, main=NA, xlab='splicing time [min]', ylab=NA)
    abline(v=median(keep$splicing.time[jj]), lwd=2.0, col='black', lty=lty)
    
    dev.off()
    
    #### Half-life
    cat('correlation between database')
    jj = which(!is.na(keep$Fried)==TRUE & !is.na(keep$Shar)==TRUE)
    cat(signif(cor(keep$Fried[jj], keep$Shar[jj]), d=2), '\n');
    
    jj = which(!is.na(keep$Fried)==TRUE & !is.na(keep$Schwa)==TRUE)
    cat(signif(cor(keep$Fried[jj], keep$Schwa[jj]), d=2), '\n');
    
    jj = which(!is.na(keep$Shar)==TRUE & !is.na(keep$Schwa)==TRUE)
    cat(signif(cor(keep$Shar[jj], keep$Schwa[jj]), d=2), '\n');
    
    cat('correlation between our estimated half-life and database')
    jj = which(keep$cv.half.life<0.5 & !is.na(keep$Fried)==TRUE)
    cat(signif(cor(keep$half.life[jj], keep$Fried[jj]), d=4), '\n');
    
    jj = which(keep$cv.half.life<0.5 & !is.na(keep$Shar)==TRUE)
    cat(signif(cor(keep$half.life[jj], keep$Shar[jj]), d=4), '\n');
    
    jj = which(keep$cv.half.life<0.5 & !is.na(keep$Schwa)==TRUE)
    cat(signif(cor(keep$half.life[jj], keep$Schwa[jj]), d=4), '\n');
    
    #par(mfrow=c(1,3))
    
    #lims = range(as.matrix(keep[, c(2, 4:6)]), na.rm=TRUE)
    #lims = as.numeric(lims)
    databases = c('Fried', 'Shar', 'Schwa')
    cex = 0.25;
    for(db in databases)
    {
      pdfname = paste(folder, 'half-life_comparison_', db, '.pdf', sep='')
      
      jj = eval(parse(text=paste('which(keep$cv.half.life<0.5 & keep$half.life>10/60 & !is.na(keep$', db, ')==TRUE)', sep='')))
      #jj = which(!is.na(keep[,kk])==TRUE)
      R = eval(parse(text=paste('signif(cor(keep$half.life[jj], keep$', db, '[jj]), d=2)', sep='')));
      
      pdf(pdfname, width=2.2, height=2.2)
      par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
      
      for(m in c(2:4))
      {
        kk = jj[which(keep$best.model[jj]==m)]
        yy1 = keep$half.life[kk]
        xx1 = eval(parse(text=paste('keep$', db, '[kk]', sep='')));
        
        if(m==2){
          plot(xx1, yy1, cex=cex, col=cols[m], cex.lab=0.7, cex.main=0.8, xlab=paste(db, ' [hr]', sep=''), ylab='estimated half-lives [hr]', main=paste('R = ', R, sep=''), log='xy', axes=FALSE)
        }else{
          points(xx1, yy1, cex=cex, col=cols[m])
        }
        #cor(xx1, yy1)
        #cat('model', m, '\n');
      }
      box()
      axis(2,at=c(1, 2, 5, 10, 24), cex.axis=0.8)
      axis(1, at=c(1, 2, 5, 10, 24, 50), las=1,cex.axis = 0.8)
      abline(0, 1, lwd=1., col='darkgray')
      dev.off()
      
    }
    
    #### Half-life distribution
    pdfname = paste(folder, 'Half_life_distribution.pdf', sep='')
    pdf(pdfname, width=2.5, height=2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    #breaks = 50
    #hist(keep[,2], breaks=breaks)
    
    h2 = hist(keep$half.life[which(keep$best.model==2 & keep$cv.half.life<0.5)], breaks = 80, plot = FALSE)
    #h3 = hist(keep[which(keep[,3]==3), 2], breaks = breaks, plot = FALSE)
    h4 = hist(keep$half.life[which(keep$best.model==4 & keep$cv.half.life<0.5)], breaks = 40, plot = FALSE)
    
    lims = range(h2$counts, h4$counts)
    plot(c(h2$breaks[1], h2$breaks, h2$breaks[length(h2$breaks)]),c(0, h2$counts, h2$counts[length(h2$counts)],0), type = 'n', main = NA, cex.lab=1.0, xlab = 'half-lives [hr]', ylab = '', xlim=c(0, 12), ylim = lims)
    polygon(c(rep(h2$breaks, each = 2)), c(0, rep(h2$counts, each=2),0), col = cols[2], border = cols[2], lwd = 1)
    
    #polygon(c(rep(h3$breaks, each = 2)), c(0, rep(h3$density, each=2),0), col = cols[3], border = 3, lwd = 1)
    polygon(c(rep(h4$breaks, each = 2)), c(0, rep(h4$counts, each=2),0), col = cols[4], border = cols[4], lwd = 1)
    
    abline(v=median(keep$half.life[which(keep$best.model==2 & keep$cv.half.life<0.5)]), lwd=1.5, col=cols[2], lty=2)
    abline(v=median(keep$half.life[which(keep$best.model==4 & keep$cv.half.life<0.5)]), lwd=1.5, col=cols[4], lty=2)
    
    dev.off()
    
    #########
    #### Rhythmic degradation (phases and fc for model 3 and 4)
    #########
    jj = which(keep$cv.gamma.phase<0.5 & keep$best.model==3)
    
    pdfname = paste(folder, 'Rhythmic_degradation_phases_Model_3.pdf', sep='')
    pdf(pdfname, width=2.2, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    breaks = seq(0, 24, by=1)
    cutoff = 0.5;
    #hist(keep[which(keep$phase.score<0.7)], 5)
    #h2 = hist(keep[which(keep[,3]==2), 2], breaks = 90, plot = FALSE)
    hist(keep$gamma.phase[jj], breaks = breaks, col=cols[3], axes=FALSE, main='Phase of Degr.', cex.main=1.0, xlab='[hr]', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 1.0)
    box()

    dev.off()
    
    jj = which(keep$cv.gamma.fc<0.5 & keep$best.model==3)
    pdfname = paste(folder, 'Rhythmic_degradation_fold_changes_Model_3.pdf', sep='')
    pdf(pdfname, width=2.2, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    breaks = seq(0, 5, by=0.2)
    fc = log2(keep$gamma.fc[jj])
    hist(fc, breaks = breaks, col=cols[3], axes=FALSE, main='FC of Degr.', cex.main=1.0, xlab='log2(peak/trough)', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 5, by=1), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    jj = which(keep$cv.gamma.phase<0.5 & keep$best.model==4)
    pdfname = paste(folder, 'Rhythmic_degradation_phases_Model_4.pdf', sep='')
    pdf(pdfname, width=2.2, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    breaks = seq(0, 24, by=1)
    cutoff = 0.5;
    #hist(keep[which(keep$phase.score<0.7)], 5)
    #h2 = hist(keep[which(keep[,3]==2), 2], breaks = 90, plot = FALSE)
    hist(keep$gamma.phase[jj], breaks = breaks, col=cols[4], axes=FALSE, main='Phase of Degr.', cex.main=1.0, xlab='[hr]', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=6), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    jj = which(keep$cv.gamma.fc<0.5 & keep$best.model==4)
    pdfname = paste(folder, 'Rhythmic_degradation_fold_changes_Model_4.pdf', sep='')
    pdf(pdfname, width=2.2, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    breaks = seq(0, 5, by=0.2)
    fc = log2(keep$gamma.fc[jj])
    hist(fc, breaks = breaks, col=cols[4], axes=FALSE, main='FC of Degr.', cex.main=1.0, xlab='log2(peak/trough)', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 5, by=1), las=1,cex.axis = 1.0)
    box()
    
    dev.off()
    
    #########
    #### Advantages of rhythmic degradation
    #########
    pdfname = paste(folder, 'Delays_RD_mRNA_Model_3.pdf', sep='')
    pdf(pdfname, width=2.1, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    jj = which(keep$cv.gamma.phase<0.5 & keep$best.model==3)
    breaks = seq(0, 24, by=1)
    dd = keep$gamma.phase[jj] - keep$phase.mRNA[jj]
    for(n in 1:length(dd))
    {
        if(dd[n]<0) dd[n] = dd[n]+24;
        if(dd[n]>24) dd[n] = dd[n]-24;
    }
    hist(dd, breaks = breaks, col=cols[3], axes=FALSE, main='Delay (RD-mRNA) ', xlim=c(0, 24), cex.main=1.0, xlab='[hr]', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    pdfname = paste(folder, 'Delays_RD_mRNA_Model_4.pdf', sep='')
    pdf(pdfname, width=2.1, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    jj = which(keep$cv.gamma.phase<0.5 & keep$best.model==4)
    breaks = seq(0, 24, by=1)
    dd = keep$gamma.phase[jj] - keep$phase.mRNA[jj]
    for(n in 1:length(dd))
    {
        if(dd[n]<0) dd[n] = dd[n]+24;
        if(dd[n]>24) dd[n] = dd[n]-24;
    }
    hist(dd, breaks = breaks, col=cols[4], axes=FALSE, main='Delay (RD-mRNA) ', xlim=c(0, 24), cex.main=1.0, xlab='[hr]', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    pdfname = paste(folder, 'Delays_RD_premRNA_Model_4.pdf', sep='')
    pdf(pdfname, width=2.1, height=2.)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    jj = which(keep$cv.gamma.phase<0.5 & keep$best.model==4)
    breaks = seq(0, 24, by=1)
    dd = keep$gamma.phase[jj] - keep$phase.premRNA[jj]
    for(n in 1:length(dd))
    {
        if(dd[n]<0) dd[n] = dd[n]+24;
        if(dd[n]>24) dd[n] = dd[n]-24;
    }
    hist(dd, breaks = breaks, col=cols[4], axes=FALSE, main='Delay (RD-premRNA) ', xlim=c(0, 24), cex.main=1.0, xlab='[hr]', ylab=NA, cex.lab=1.0)
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=4), las=1,cex.axis = 1.0)
    box()
    dev.off()
    
    ### test of advantages
    TEST = FALSE
    if(TEST)
    {
        ## relative amplitude of rhythmic degradation
        kk = which(keep$best.model==4)
        aa = keep$gamma.fc[kk]
        aa = (aa-1)/(aa+1)
        hist(aa, breaks=50)
        
        ## test of Michaelis-Menten
        tt = seq(0, 48, by=0.5)
        m = 1*(1+0.5*cos(2*pi/24*(tt-6)))
        m1 = 1*(1+0.5*cos(2*pi/24*(tt-6+12)))
        K = 0.1
        r = 1/(K+m)
        m = m/mean(m)
        m1 = m1/mean(m1)
        r = r/mean(r)
        plot(tt, r, type='l', col='red', ylim=range(r,m,m1))
        points(tt, m, type='l', col='green')
        points(tt, m1, type='l', col='blue')
        
        ####
        ## test other things
        m2 = which(keep$best.model==2)
        m3 = which(keep$best.model==3)
        m4 = which(keep$best.model==4)
        
        ii = grep('abs.mRNA', colnames(Tt))
        jj = grep('abs.premRNA', colnames(Tt))
        data1 = (as.matrix(Tt[keep$index,ii]))
        data2 = (as.matrix(Tt[keep$index,jj]))
        
        stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
        stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
        #stat1 = cbind(stat1, qvals(stat1[,6]))
        #stat2 = cbind(stat2, qvals(stat2[,6]))
        boxplot(stat1[,3] ~ keep$best.model)
        
        plot((keep$phase.mRNA[m4]-keep$phase.premRNA[m4]), (stat1[m4, 4]/stat2[m4, 4]), cex=0.5, col=cols[4], ylim=c(0, 5))
        points((keep$phase.mRNA[m2]-keep$phase.premRNA[m2]), (stat1[m2, 4]/stat2[m2, 4]), cex=0.5, col=cols[2], pch=2)
        abline(h=1, col='red', lwd=2.0)
        
        par(mfrow=c(1,2))
        for(m in c(2, 4))
        {
            #m = 4;
            #m = 2
            cutoff = 100
            m4 = which(keep$best.model==m & keep$cv.half.life< cutoff & keep$half.life>0.2)
            #if(m==2) m4 = m4[c(1:500)]
            
            #dd = keep$phase.mRNA[m4]-keep$phase.premRNA[m4]
            if(m==2) dd = T$phase.mRNA[keep$index[m4]] - T$phase.int.m2[keep$index[m4]]
            if(m==4) dd = T$phase.mRNA[keep$index[m4]] - T$phase.int.m4[keep$index[m4]]
            for(n in 1:length(dd))
            {
                if(dd[n]>24) dd[n] = dd[n]-24;
                if(dd[n]<0) dd[n] = dd[n]+24;
                
            }
            
            plot(keep$half.life[m4], dd, cex=0.4, log='x')
            xx = seq(min(keep$half.life[m4]), max(keep$half.life[m4]), length.out=200)
            yy = atan(2*pi/24*xx/log(2))/(2*pi)*24
            points(xx, yy, type='l', col='red', lwd=4.0)
            yy = atan(2*pi/24*xx/log(2))/(2*pi)*24
            points(xx, yy, type='l', col='red', lwd=4.0)
            abline(h=6, col='darkgray', lwd=4.0)
            #yy = atan(2*pi/24*xx/log(2))
            #points(xx, yy, type='l', col='red')
            #if(m==4)
            #{
            #   plot(keep$half.life[m4]*(1+keep$gamma.fc[m4])/2, dd, cex=0.4, log='x')
            #   xx = seq(min(keep$half.life[m4]), max(keep$half.life[m4]), length.out=200)
            #   yy = atan(2*pi/24*xx/log(2))/(2*pi)*24
                #points(xx, yy, type='l', col='red', lwd=4.0)
                #    abline(h=6, col='darkgray', lwd=4.0)
            
            #    plot(keep$half.life[m4]*(1+keep$gamma.fc[m4])/(2*keep$gamma.fc[m4]), dd, cex=0.4, log='x')
            #xx = seq(min(keep$half.life[m4]), max(keep$half.life[m4]), length.out=200)
            #   yy = atan(2*pi/24*xx/log(2))/(2*pi)*24
                #points(xx, yy, type='l', col='red', lwd=4.0)
                #   abline(h=6, col='darkgray', lwd=4.0)
                #}
        }
        
        ### Compare model 2, 3 and 4
        load(file ='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Total_models_parameters_results.Rdata')
        ii = grep('abs.mRNA', colnames(Tt))
        jj = grep('abs.premRNA', colnames(Tt))
        data1 = (as.matrix(Tt[keep$index,ii]))
        data2 = (as.matrix(Tt[keep$index,jj]))
        
        stat1 = t(apply(log2(data1),1, f24_R2_alt2, t=c(0:47)*2))
        stat2 = t(apply(log2(data2),1, f24_R2_alt2, t=c(0:47)*2))
        
        par(mfrow=c(4,2))
        boxplot(stat1[,2] ~ keep$best.model, main='Mean')
        boxplot(-log10(stat1[,6]) ~ keep$best.model, main='-log10(p-value)')
        boxplot(stat1[,3] ~ keep$best.model, main='log2(amplitude)', ylim=c(0, 2))
        boxplot(stat1[,4] ~ keep$best.model, main='relative amplitude', ylim=c(0, 1))
        hist(keep$phase.mRNA[which(keep$best.model==2)], breaks=seq(0, 24, by=1))
        hist(keep$phase.mRNA[which(keep$best.model==3)], breaks=seq(0, 24, by=1))
        hist(keep$phase.mRNA[which(keep$best.model==4)], breaks=seq(0, 24, by=1))
        
        ### Compare ratio of relative amplitudes in function of delay for model 2 and 4;
        ii = grep('abs.mRNA', colnames(Tt))
        jj = grep('abs.premRNA', colnames(Tt))
        data1 = (as.matrix(Tt[keep$index,ii]))
        data2 = (as.matrix(Tt[keep$index,jj]))
        
        stat1 = t(apply(log2(data1),1, f24_R2_alt2, t=c(0:47)*2))
        stat2 = t(apply(log2(data2),1, f24_R2_alt2, t=c(0:47)*2))
        stat1[, 4] = t(apply((data1),1, f24_R2_alt2, t=c(0:47)*2))[,4]
        stat2[, 4] = t(apply((data2),1, f24_R2_alt2, t=c(0:47)*2))[,4]
        
        cols=c('bisque','olivedrab','palevioletred3','purple4')
        
        xx = data.frame(keep$gene, keep$index, keep$best.model, keep$half.life, stat1[,2], stat2[, 2], stat1[,4], stat2[, 4], stat1[,5], stat2[, 5])
        colnames(xx) = c('gene', 'index', 'best.model', 'half.life', 'mean.mRNA', 'mean.premRNA', 'relamp.mRNA', 'relamp.premRNA', 'phase.mRNA', 'phase.premRNA')
        xx = xx[which(xx$best.model==2 |xx$best.model==4), ]
        
        kk = which(xx$half.life>0.2)
        yy = xx[kk,]
        dd = yy$phase.mRNA - yy$phase.premRNA
        for(n in 1:length(dd))
        {
            if(dd[n]>18) dd[n] = dd[n]-24;
            if(dd[n]<(-6)) dd[n] = dd[n]+24;
            
        }
        rr = yy$relamp.mRNA/yy$relamp.premRNA
        par(mfrow=c(2,2))
        hist(dd[which(yy$best.model==2)], breaks=24)
        hist(rr[which(yy$best.model==2)], breaks=20, xlim=c(0, 2))
        hist(dd[which(yy$best.model==4)], breaks=24)
        hist(rr[which(yy$best.model==4)], breaks=20, xlim=c(0, 2))
        
        ss = c()
        ss1 = c()
        kk = xx$index[which(xx$best.model==4 & xx$half.life>0.2)]
        zt = seq(0, 94, by=1)
        for(n in kk)
        {
            Min = T$Min.int.m4[n]
            Amp = T$Amp.int.m4[n]
            beta = T$beta.int.m4[n]
            phase = T$phase.int.m4[n];
            s2 = compute.s.beta(zt, Min, Amp, phase, beta);
            stat = f24_R2_alt2(s2, t=zt)
            b = stat[4]
            splicing.k = T$splicing.k.m4[n]
            gamma.mean = T$gamma.m4[n]+T$amp.gamma.m4[n]/2; ## mean half-life of 2h
            
            gamma = T$gamma.m4[n]
            for(m in 1:1)
            {
                phase.gamma = phase +12
                amp.gamma = T$amp.gamma.m4[n]
                
                m2 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
                stat = f24_R2_alt2(m2, t=zt)
                dd = stat[5] - phase
                if(dd>18) dd = dd-24;
                if(dd<(-6)) dd = dd+24;
                ss = rbind(ss, c(dd, stat[4]/b))
                
                amp.gamma = 0
                m2 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
                stat = f24_R2_alt2(m2, t=zt)
                dd = stat[5] - phase
                if(dd>18) dd = dd-24;
                if(dd<(-6)) dd = dd+24;
                ss1 = rbind(ss1, c(dd, stat[4]/b))
            }
        }
        
        xlim = c(-6, 18)
        ylim = range(xx$relamp.mRNA/xx$relamp.premRNA)
        ylim = c(0.1, 10)
        
        for(m in c(2, 4))
        {
            mm = which(xx$best.model==m & xx$half.life>0.2)
            #mm = which(keep$best.model==m & keep$cv.half.life< cutoff & keep$half.life>0.2)
            #if(m==2) m4 = m4[c(1:500)]
            
            #dd = keep$phase.mRNA[m4]-keep$phase.premRNA[m4]
            dd = xx$phase.mRNA[mm] - xx$phase.premRNA[mm]
            #if(m==4) dd = T$phase.mRNA[keep$index[m4]] - T$phase.int.m4[keep$index[m4]]
            for(n in 1:length(dd))
            {
                if(dd[n]>18) dd[n] = dd[n]-24;
                if(dd[n]<(-6)) dd[n] = dd[n]+24;
                
            }
            if(m==2)
            {
                plot(dd, xx$relamp.mRNA[mm]/xx$relamp.premRNA[mm], type='n', cex=0.5, log='y', xlim=xlim, ylim=ylim, col=cols[2], xlab='delay (mRNA-premRNA) [hr]', ylab='raito of relative amplitude')
                dd = seq(0,6, length.out=100)
                yy = cos(pi/12*dd)
                points(dd, yy, type='l', col='red', lwd=4.0)
                #yy = cos(pi/12*dd)+1
                #points(dd, yy, type='l', col='red', lwd=4.0)
                abline(h=1, col='darkgray', lwd=4.0)
                abline(v=6, col='darkgray', lwd=4.0)
                abline(v=0, col='darkgray', lwd=4.0)
            }
            #yy = atan(2*pi/24*xx/log(2))
            #points(xx, yy, type='l', col='red')
            if(m==4)
            {
                points(dd, xx$relamp.mRNA[mm]/xx$relamp.premRNA[mm], cex=0.5, col=cols[4], pch=16)
            }
           
        }
        points(ss[,1], ss[,2], cex=0.4, col='orange', pch=2)
        points(ss1[,1], ss1[,2], cex=0.4, col='green', pch=3)
        
        
        ss = c()
        ss1 = c()
        mm = which(xx$best.model==4 & xx$half.life>0.2)
        kk = xx$index[mm]
        zt = seq(0, 94, by=2)
        for(n in kk)
        {
            Min = T$Min.int.m4[n]
            Amp = T$Amp.int.m4[n]
            beta = T$beta.int.m4[n]
            phase = T$phase.int.m4[n];
            s2 = compute.s.beta(zt, Min, Amp, phase, beta);
            stat = f24_R2_alt2(s2, t=zt)
            b = stat[4]
            splicing.k = T$splicing.k.m4[n]
            gamma.mean = T$gamma.m4[n]+T$amp.gamma.m4[n]/2; ## mean half-life of 2h
            
            gamma = T$gamma.m4[n]
            for(m in 1:1)
            {
                phase.gamma = phase +12
                amp.gamma = T$amp.gamma.m4[n]
                
                m2 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
                stat = f24_R2_alt2(m2, t=zt)
                ss = rbind(ss, c(stat[5], stat[4]))
                
                amp.gamma = 0
                m2 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
                stat = f24_R2_alt2(m2, t=zt)
                ss1 = rbind(ss1, c(stat[5], stat[4]))
            }
        }
        
        par(mfrow=c(1,2))
        nb = 100
        dd = xx$phase.mRNA[mm] - xx$phase.premRNA[mm]
        for(n in 1:length(dd))
        {
            if(dd[n]>18) dd[n] = dd[n]-24;
            if(dd[n]<(-6)) dd[n] = dd[n]+24;
        }
        plot(xx$half.life[mm], dd, cex=0.6, ylim=c(-1, 6), col='blue', log='x')
        aa = xx$half.life[mm]
        bb = dd
        o1 = order(xx$half.life[mm])
        aa = runmean(aa[o1], nb)
        bb = runmean(bb[o1], nb)
        points(aa, bb, col='blue', type='l', lwd=4.0)
        #abline(0, 1, col='red', lwd=2.0)
        dd = ss[,1] - xx$phase.premRNA[mm]
        for(n in 1:length(dd))
        {
            if(dd[n]>18) dd[n] = dd[n]-24;
            if(dd[n]<(-6)) dd[n] = dd[n]+24;
        }
        points(xx$half.life[mm], dd, col='orange', cex=0.4)
        aa = xx$half.life[mm]
        bb = dd
        o1 = order(xx$half.life[mm])
        aa = runmean(aa[o1], nb)
        bb = runmean(bb[o1], nb)
        points(aa, bb, col='orange', type='l', lwd=4.0)
        
        dd = ss1[,1] - xx$phase.premRNA[mm]
        for(n in 1:length(dd))
        {
            if(dd[n]>18) dd[n] = dd[n]-24;
            if(dd[n]<(-6)) dd[n] = dd[n]+24;
        }
        points(xx$half.life[mm], dd, col='black', cex=0.4)
        aa = xx$half.life[mm]
        bb = dd
        o1 = order(xx$half.life[mm])
        aa = runmean(aa[o1], nb)
        bb = runmean(bb[o1], nb)
        points(aa, bb, col='black', type='l', lwd=4.0)
        
        plot(xx$half.life[mm], xx$relamp.mRNA[mm]/xx$relamp.premRNA[mm], cex=0.5, log='xy', col='blue', ylim=c(0.1, 10.0))
        #abline(0, 1, col='red', lwd=2.0)
        points(xx$half.life[mm], ss[,2]/xx$relamp.premRNA[mm], col='orange', cex=0.4)
        points(xx$half.life[mm], ss1[,2]/xx$relamp.premRNA[mm], col='black', cex=0.4)
        
        ### test another hypothesis
        zt = seq(0, 24, by=0.1)
        Min = 5
        Amp = 5.0
        beta = 1.0
        phase = 6;
        s2 = compute.s.beta(zt, Min, Amp, phase, beta);
        
        splicing.k = log(2)/(5/60);
        gamma = log(2)/2; ## mean half-life of 2h
        
        amp.gamma = 0.0;
        phase.gamma = 18;
        
        m20 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
        gamma = log(2)/4
        m21 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
        gamma = log(2)/20
        m22 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
        
        deg = gamma + amp.gamma*(1+cos(2*pi/24*(zt-phase.gamma)))/2
        
        s2 = s2/mean(s2)
        m20 = m20/mean(m20)
        m21 = m21/mean(m21)
        m22 = m22/mean(m22)
        deg = deg/mean(deg)
        
        lims = range(s2, m20, m21, m22, deg)
        
        pdf.name = paste("myplots/Total/Illustration_M2.pdf", sep='')
        pdf(pdf.name, width=2.4, height=2.0)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        plot(zt, s2, type='l', lwd=2.5, lty=1, col=col.int, ylim=lims, xlab=NA, ylab=NA, main=NA, axes=FALSE)
        points(zt, m20, type='l', lwd=2.5, col=col.ex, lty=1)
        points(zt, m21, type='l', lwd=2.5, col=col.ex, lty=2)
        points(zt, m22, type='l', lwd=2.5, col=col.ex, lty=3)
        #points(zt, deg, type='l', lwd=2.0, lty=1, col=col.deg)
        box()
        axis(2, cex.axis=1.0)
        axis(1, at=seq(0, 24, by=3), las=1,cex.axis = 1.0)
        abline(v=6, col='gray', lwd=2.0)
        abline(v=12, col='gray', lwd=2.0)
        
        dev.off()
        
        #######
        ### test phase shift, amplitude amplification and others
        #######
        ii = grep('abs.mRNA', colnames(Tt))
        jj = grep('abs.premRNA', colnames(Tt))
        data1 = (as.matrix(Tt[keep$index,ii]))
        data2 = (as.matrix(Tt[keep$index,jj]))
        
        stat1 = t(apply(log2(data1),1, f24_R2_alt2, t=c(0:47)*2))
        stat2 = t(apply(log2(data2),1, f24_R2_alt2, t=c(0:47)*2))
        stat1[, 4] = t(apply((data1),1, f24_R2_alt2, t=c(0:47)*2))[,4]
        stat2[, 4] = t(apply((data2),1, f24_R2_alt2, t=c(0:47)*2))[,4]
        #cols=c('bisque','olivedrab','palevioletred3','purple4')
        xx = data.frame(keep$gene, keep$index, keep$best.model, keep$half.life, keep$gamma.fc, stat1[,2], stat2[, 2], stat1[,4], stat2[, 4], stat1[,5], stat2[, 5])
        colnames(xx) = c('gene', 'index', 'best.model', 'half.life', 'gamma.fc', 'mean.mRNA', 'mean.premRNA', 'relamp.mRNA', 'relamp.premRNA', 'phase.mRNA', 'phase.premRNA')
        xx = xx[which(xx$half.life>0.2 & xx$best.model==4), ]
        xx$gamma.fc = (xx$gamma.fc-1)/(xx$gamma.fc+1)
        
        k = log(2)/(5/60)
        w = 2*pi/24;
        A1 = 0.2;
        phase1 = 3
        A2 = 0.5;
        delta = seq(0, 24, by=0.1)
        #phase2 = phase1
        
        ggamma = log(2)/seq(0.5, 5, by=0.5)
        for(n in 1:length(ggamma))
        {
            ss = c()
            gamma = ggamma[n]
            for(dd in delta)
            {
                C = A1*A2/2*gamma*(w*sin(w*dd)+gamma*cos(w*dd))
                shift = atan2((w*(A1-A2*cos(w*dd))-A2*gamma*sin(w*dd)), (gamma*(A1-A2*cos(w*dd))+A2*w*sin(w*dd)))
                ampl = gamma*sqrt(A1^2+A2^2-2*A1*A2*cos(w*dd))/(sqrt(gamma^2+w^2)-C)
                mean = k/gamma*(2*(gamma^2+w^2-C)/(2*(gamma^2+w^2)-gamma^2*A2^2))
                
                shift = shift/(2*pi)*24 - atan2(w, gamma)/(2*pi)*24
                #shift = shift/(2*pi)*24
                #if(shift>24) shift = shift -24
                #if(shift<0) shift = shift +24
                #shift = shift
                ampl = ampl/(gamma*A1/(sqrt(gamma^2+w^2)))
                #ampl = ampl/A1
                mean = mean/(k/gamma)
                
                ss = rbind(ss, c(shift, ampl, mean, mean*ampl))
            }
            colnames(ss) = c('phase.shift', 'ampl.amplification', 'mean', 'absolute.amp')
            ss[,1] = ss[,1]/6
            
            if(n==1)
            {
                plot(delta, ss[,1], col='darkred', type='l', lwd=2.0, ylim=c(-4, 4), xlab='Delay (Degr-Prod)')
            }else{
                points(delta, ss[,1], col='red', type='l', lwd=2.0)
            }
            points(delta, ss[,2], col='blue', type='l', lwd=2.0)
            #points(delta, ss[,3], col='orange', type='l', lwd=2.0)
            #points(delta, ss[,4], col='darkgray', type='l', lwd=2.0)
            #abline(h=(atan2(w, gamma)/(2*pi)*24)/6, lwd=2.0, col='black')
        }
        abline(v=c(6, 12, 18), lwd=2.0, lty=2, col='gray')
        abline(h=c(0, 1), lwd=2.0, lty=2, col='gray')



    }
    
}

compute.phase.amp.4premrna.mrna.with.fitting.res = function(index=1, T=T)
{
  #Tt = T; i = 4;
  zt.p = seq(0,96,by = 0.1)
  gene.index = index;
  i = index;
  cat(i, '\n');
  #### compute the fitting results for reads and also RPKM
  best.model = T$BIC.best.model[i]
  mm = grep(paste('.m', best.model, sep=''),colnames(T))
  mm = mm[-c(1, length(mm))];
  par = unlist(T[i, mm])
  gamma = par[1];
  model = best.model;
  if(model>1){
    if(model==2){
      eps.gamma = 0.0; phase.gamma = 0.0; splicing.k = par[2];
      param.synthesis.1 = par[3]; param.synthesis.2 = par[4]; param.synthesis.3 = par[5]; param.synthesis.4 = par[6];}
    if(model==3){
      eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4]; param.synthesis.1 = par[5];
      param.synthesis.2 = 0; param.synthesis.3 = 0; param.synthesis.4 = 1;}
    if(model==4){
      eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4];
      param.synthesis.1 = par[5]; param.synthesis.2 = par[6]; param.synthesis.3 = par[7]; param.synthesis.4 = par[8];}
    
    ### compute the m and s with fitted parameters
    w = 2*pi/24;
    #read.m = convert.nb.reads(compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
    #                                         param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4), L.m);
    #read.s = convert.nb.reads(compute.s.beta(zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4), L.s);
    m = compute.m.beta(t = zt.p, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                       param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
    s = compute.s.beta(t = zt.p, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
    
    m22 = compute.m.beta(t = zt.p, gamma, 0, (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                         param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
  }
  ### phases and relative amplitudes for m, s and m22
  if(model==1) {phase.amps = c(NA, 0, NA, 0, NA, 0)
  }else{
    m = m/mean(m); s = s/mean(s); m22 = m22/mean(m22);
    kk = which(zt.p<24); 
    zt.p = zt.p[kk];s = s[kk]; m = m[kk]; m22 = m22[kk];
    m.rel = (max(m)-min(m))/2; m.phi = zt.p[which(m==max(m))][1];
    if(model==2){
      s.rel = (max(s)-min(s))/2; s.phi = zt.p[which(s==max(s))][1];
      #m2.rel = (max(m22) - min(m22))/2; m2.phi = zt.p[which(m22==max(m22))];
      phase.amps = c(m.phi, m.rel, s.phi, s.rel, m.phi, m.rel);
    }
    if(model==3){
      phase.amps = c(m.phi, m.rel, NA, 0, NA, 0);
    }
    if(model==4){
      s.rel = (max(s)-min(s))/2; s.phi = zt.p[which(s==max(s))][1];
      m2.rel = (max(m22) - min(m22))/2; m2.phi = zt.p[which(m22==max(m22))];
      phase.amps = c(m.phi, m.rel, s.phi, s.rel, m2.phi, m2.rel);
    }
  }
  names(phase.amps) = c('phase.mRNA', 'relamp.mRNA', 'phase.premRNA', 'relamp.premRNA', 'phase.m22', 'relamp.m22');
  return(phase.amps)
}

plot.genes.examples = function(index=c(1), pdfname='myplots/gene_examples/Examples_Test.pdf', T=T, RF.KO = FALSE, Figure=FALSE, folder=folder, 
                               summary.plot=FALSE, RD.advantage = FALSE, shape.test.M3 = FALSE, outliers.removal = FALSE)
{
    #Tt = T; i = which(T$gene=='Per1');
    #### plot to check fitting with examples
    Tt = T;
    zt = seq(0,94,by = 2);
    zt.p = seq(0,96,by = 0.1)
    ZT.ex = grep('.count.mRNA', colnames(T))
    ZT.int = grep('.count.premRNA', colnames(T))
    #index = index[which(!is.na(Tt$BIC.best.model[index])==TRUE)]
    
    ### half-lives database
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/mRNAs_half_lives_databases.Rdata')
    load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata')
    bounds = set.bounds(model=4); upper = bounds$upper;lower = bounds$lower; gamma.max = upper[1]; gamma.min = lower[1];
    
    if(!Figure)
    {
        if(RF.KO){pdf(pdfname, width = 20, height = 8)
        }else{pdf(pdfname, width = 12, height = 8)}
    }
    for(i in index)
    {
      gene.index = i;
      R.m = as.numeric(unlist(T[gene.index, ZT.ex])); a.m = R.m;
      R.s = as.numeric(unlist(T[gene.index, ZT.int])); a.s = R.s;
      #a.m = as.numeric(rep(alphas[gene.index, c(1:12)], 4))
      #a.s = as.numeric(rep(alphas[gene.index, c(13:24)], 4))
      a.m = rep(as.numeric(T[gene.index, grep('alpha.mRNA.ZT', colnames(T))]), 4);
      a.s = rep(as.numeric(T[gene.index, grep('alpha.premRNA.ZT', colnames(T))]), 4);
      L.m = as.numeric(unlist(T$length.mRNA[gene.index]));
      L.s = as.numeric(unlist(T$length.premRNA[gene.index]));
      ### normalized RPKM
      M = norm.RPKM(R.m, L.m); S = norm.RPKM(R.s, L.s);
      outlier.M = as.numeric(unlist(strsplit(as.character(T$outlier.m[gene.index]), ';')));outlier.M = outlier.M[which(!is.na(outlier.M)==TRUE)];
      outlier.S = as.numeric(unlist(strsplit(as.character(T$outlier.s[gene.index]), ';')));outlier.S = outlier.S[which(!is.na(outlier.S)==TRUE)];
      #### condition of ad lib
      M0 = M; S0 = S; 
      #### compute the fitting results for reads and also RPKM
      best.model = T$BIC.best.model[i]
      mm = grep(paste('.m', best.model, sep=''),colnames(T))
      mm = mm[-c(1, length(mm))]
      #mm = mm[1:(length(mm)/2)]
      par = unlist(Tt[i, mm])
      gamma = par[1];
      model = best.model;
      if(model==1){
        read.m = convert.nb.reads(rep(mean(M), length(R.m)), L.m); read.s = convert.nb.reads(rep(mean(S), length(R.s)), L.s);
        m = rep(mean(M), length(zt.p)); s = rep(mean(M), length(zt.p));
      }else{
        if(model==2){
          eps.gamma = 0.0; phase.gamma = 0.0; splicing.k = par[2];
          param.synthesis.1 = par[3]; param.synthesis.2 = par[4]; param.synthesis.3 = par[5]; param.synthesis.4 = par[6];}
        if(model==3){
          eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4]; param.synthesis.1 = par[5];
          param.synthesis.2 = 0; param.synthesis.3 = 0; param.synthesis.4 = 1;}
        if(model==4){
          eps.gamma = par[2]; phase.gamma = par[3]; splicing.k = par[4];
          param.synthesis.1 = par[5]; param.synthesis.2 = par[6]; param.synthesis.3 = par[7]; param.synthesis.4 = par[8];}
        #print(c(param.synthesis.3, phase.gamma))
        w = 2*pi/24;
        mu.m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                               param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
        read.m = convert.nb.reads(mu.m, L.m);
        mu.s = compute.s.beta(zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
        read.s = convert.nb.reads(mu.s, L.s);
        m = compute.m.beta(t = zt.p, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                           param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
        s = compute.s.beta(t = zt.p, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
        m0 = m; s0 = s;
        #s = compute.s.beta(t = zt, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
        #print(par);
        #print(compute.m.beta(t = zt, gamma, eps.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4, simulation.only=TRUE))
        #w = 2*pi/24;
        #m = compute.m.beta(t = zt, gamma, eps.gamma*sqrt(1+w^2/gamma^2), (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
        #                   param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
        ### convert rpkm into mean of read numbers
        #mu.m = convert.nb.reads(m, L.m);
        #mu.s = convert.nb.reads(s, L.s);
        
        ### theoretical mRNA profiles without RD
        m22 = compute.m.beta(t = zt.p, gamma, 0, (phase.gamma-atan2(w, gamma)/w), splicing.k*gamma,
                             param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4);
      }
      
      ### check outliers for the fitted model
      if(outliers.removal)
      {
        #loglike.m = -2*dnbinom(as.numeric(R.m), size=1/a.m, mu=as.numeric(read.m), log = TRUE);outlier.M = index.outliers.loglike(loglike.m, c=1)
        #loglike.s = -2*dnbinom(as.numeric(R.s), size=1/a.s, mu=as.numeric(read.s), log = TRUE); outlier.S = index.outliers.loglike(loglike.s, c=1);
        rrs.m = (M-mu.m)^2;outlier.M = index.outliers.loglike(rrs.m);
        rrs.s = (S-mu.s)^2;outlier.S = index.outliers.loglike(rrs.s);
        cat('mRNA outliers : ', outlier.M, '\n')
        cat('premRNA outliers : ', outlier.S, '\n')
        #outlier.M = c(13:36) 
        #outlier.S = c(13:24, 37:48) 
        #outlier.m = intersect(index.outliers.loglike(loglike.m[3,]), intersect(index.outliers.loglike(loglike.m[1,]), index.outliers.loglike(loglike.m[2,])))
        #outlier.s = intersect(index.outliers.loglike(loglike.s[3,]), intersect(index.outliers.loglike(loglike.s[1,]), index.outliers.loglike(loglike.s[2,])))
      }
     
      ### transform back the normalized parameters for epsilon and phase of RD
      eps.gamma.r = eps.gamma*sqrt(1+w^2/gamma^2); phase.gamma.r = (phase.gamma-atan2(w, gamma)/w)%%24;
      eps.gamma = eps.gamma.r; phase.gamma = phase.gamma.r;
      
      ##### extract all information of fitting and data
      cat(i, '...'); cat(as.character(Tt$gene[i]));
      if(!is.na(Tt$BIC.best.model[i])){cat('\n');
      }else{cat('...NO model selection result');}
      gg = Tt$gene[i];
      mm = match(gg, shar[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl1 = signif(shar[mm,2], d=2)}else{hl1 = NA;}
      mm = match(gg, fried[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl2 = signif(fried[mm,2], d=2)}else{hl2 = NA;}
      mm = match(gg, schwa[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl3 = signif(schwa[mm,2], d=2)}else{hl3 = NA;}
      hl.ref = paste(hl1, 'h(shar); ', hl2, 'h(fried); ', hl3, 'h(schwar)',  sep='', collapse=';')
      
      kk = which(R.WT.RF$gene==gg)[1]
      if(length(kk)>0){
        M1 = unlist(R.WT.RF[kk, grep('.abs.mRNA', colnames(R.WT.RF))]);
        S1 = unlist(R.WT.RF[kk, grep('.abs.premRNA', colnames(R.WT.RF))]);
      }else{
        M1 = rep(NA, 24)
        S1 = rep(NA, 24)
      }
      kk = which(R.KO.RF$gene==gg)[1]
      if(length(kk)>0){
        M2 = unlist(R.KO.RF[kk, grep('.abs.mRNA', colnames(R.KO.RF))]);
        S2 = unlist(R.KO.RF[kk, grep('.abs.premRNA', colnames(R.KO.RF))]);
      }else{
        M2 = rep(NA, 12)
        S2 = rep(NA, 12)
      }
      #col.ex <<- 'steelblue'; col.int <<- 'green3'; col.accepted.ex <<- 'steelblue2'; col.accepted.int <<- 'limegreen'; col.rejected <<- 'gray' ; zt <<- seq(0,46,by=2)
      #col.phases <<- rainbow(n = 241,s = 0.9, v = 0.9);
      #col.deg <<- 'orangered'
      #Ncolors <<- 21; col.rel.ex <<- col.ramp.ex(n = Ncolors); col.rel.int <<- col.ramp.int(n = Ncolors);
      #col.pval <<- col.ramp.pval(n = 11); col.median.level  <<- col.ramp.median(n = 15);col.kept <<- c('yellowgreen','tomato1')
      if(Figure)
      {
        Old.version.plot4examples = FALSE ## NOT used
        if(Old.version.plot4examples)
        {
          pdfname = paste(folder, '/RPKM_fitting_results_', gg, '.pdf', sep='')
          pdf(pdfname, width = 1.6, height = 2.)
          par(cex = 0.7, las = 1, tcl = -0.3)
          #par(mfcol = c(2,1))
          layout(matrix(c(1,2),nrow = 2, ncol = 1), heights=c(1,1.3))
          #par(mgp = c(0.5,0.3,0.), mar = c(0.,1.2,0.7,0.5));
          Read.counts.fitting = FALSE
          if(Read.counts.fitting)
          {
            ### read counts plot with fitting
            ylim = range(c(R.m, R.s, read.s, read.m));xlim = c(0, 96)
            plot(zt, R.s, col = 'green3', type = 'l',lwd = 0.7, pch = 1, ylim=ylim, xlim=xlim, log='y', main=gg, cex.main=0.7, cex.lab=0.7, xlab=NA, 
                 ylab='read counts', axes=FALSE)
            points(zt, R.s, col='green3', type='p', cex=0.5)
            points(zt, R.m, type='l', lwd=0.7, col='steelblue')
            points(zt, R.m, type='p', cex=0.5, col='steelblue', pch=5)
            points(zt, read.s, type='l', lwd=1.5, col='limegreen')
            points(zt, read.m, type='l', lwd=1.5, col='steelblue2')
            axis(1,at=seq(0, 96, by=12),cex.axis =0.6)
            #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
            #lims = signif(lims, d=1)
            #by = signif((lims[2]-lims[1])/4,d=1)
            #print(gene)
            #print(lims)
            axis(2, las=1, cex.axis = 0.6)
            box();
          }
          #ylim = c(0, 2.5)
          ylim = range(c(S0, s0, M0, m0));xlim = c(0, 96)
          if(ylim[1]<10^-6){
            ylim = range(c(s0, m0));
          }
          #### absolute RPKM plot without fitting
          par(mgp = c(0.5,0.3,0.), mar = c(0.,1.2,0.7,0.5));
          
          plot(zt, S0, col = 'green3', type = 'l', lwd = 0.8, ylim=ylim, xlim=xlim, 
               main = gg, cex.main=0.7, cex.lab=0.7, xlab=NA, log='y', ylab='RPKM (log)', axes=FALSE);
          points(zt, S0, col='green3', type='p', cex=0.3, pch=5);
          points(zt, M0, type='l', lwd=0.8, col='steelblue');
          points(zt, M0, type='p', cex=0.4, col='steelblue', pch=16);
          points(zt.p, s0, type='l', lwd=1., col='limegreen')
          points(zt.p, m0, type='l', lwd=1., col='steelblue2')
          #axis(1,at=seq(0, 96, by=6), tck=-0.04, labels = NA)
          axis(2, las=3, cex.axis = 0.6, tck=-.03, labels = FALSE)
          axis(2, las=3, cex.axis =0.6, lwd=0, line = -0.3)
          box();
          #abline(h=mean(S), col='gray', lwd=1.0);abline(h=mean(M), col='gray', lwd=1.0);
          
          #### relative RPKM plot with fitting
          #S = S/mean(S);M = M/mean(M);s = s/mean(s);m = m/mean(m);
          #S = S/mean(s);M = M/mean(m);s = s/mean(s);m = m/mean(m);
          ylim = range(c(S, s, M, m)); xlim = c(0, 96);
          cols = c('gray','green3','tomato','black')
          mains = c('CS-CD M1', 'RS-CD M2', 'CS-RD M3', 'RS-RD M4')
          par(mgp = c(0.5,0.3,0), mar = c(1.0,1.2,1.0,0.5));
          plot(zt, S, col = 'green3', type = 'l', lwd = 0.5, pch = 1, ylim=ylim, xlim=xlim, 
               main = mains[model], col.main=cols[model], cex.main=0.5, cex.lab=0.7, xlab=NA, log='', ylab='rel.RPKM', axes=FALSE)
          points(zt, S, col='green3', type='p', cex=0.3, pch=5)
          points(zt, M, type='l', lwd=0.5, col='steelblue')
          points(zt, M, type='p', cex=0.4, col='steelblue', pch=16)
          abline(h=1, col='gray', lwd=1.0);
          points(zt.p, s, type='l', lwd=1., col='limegreen')
          points(zt.p, m, type='l', lwd=1., col='steelblue2')
          axis(1,at=seq(0, 96, by=6), tck=-0.03, labels = NA)
          axis(1,at=seq(0, 96, by=24),cex.axis =0.6, lwd=0, line = -0.4)
          lims = signif(lims, d=1)
          by = signif((lims[2]-lims[1])/4,d=1)
          axis(2,at=seq(0, 10, by=1), las=1,cex.axis = 0.6, tck=-.03, labels = FALSE)
          axis(2, las=3, cex.axis =0.6, lwd=0, line = -0.3)
          box();
          if(model>1){
            if(model==2){dd = rep(gamma, length(zt.p));}
            if(model==3|model==4){dd = gamma*(1+eps.gamma*(1+cos(2*pi/24*(zt.p-phase.gamma))))};
            dd = dd/mean(dd)
            points(zt.p, dd, type='l', lwd=1., col='red')
          }  
        }
        
        if(!RD.advantage){
          ### plot examples in Figure 2 to illustrate M2, M3 and M4
          pdfname = paste(folder, '/RPKM_fitting_results_', gg, '.pdf', sep='')
          pdf(pdfname, width = 2.2, height = 1.6)
          par(cex = 0.7, las = 1, tcl = -0.3)
          par(mgp = c(1.0,0.5,0), mar = c(1.8,1.4,1.5,1.4));
          
        }else{
          ### plot examples in Figure 5 to illustrate roles of RD: amplifiers and phase tuners
          pdfname = paste(folder, '/RD_advantages_RPKM_', gg, '.pdf', sep='')
          pdf(pdfname, width = 1.6, height = 1.2)
          par(cex = 0.7, las = 1, tcl = -0.3)
          par(mgp = c(1.0,0.5,0), mar = c(1.8,1.4,1.5,1.4));
          #par(mgp = c(0.5,0.3,0.), mar = c(0.,1.2,0.7,0.5));
          #layout(matrix(c(1,2),nrow = 2, ncol = 1), heights=c(1,1.3))
        }
        ### degradation rate
        if(model>1){
          if(model==2){dd = rep(gamma, length(zt.p));}
          if(model==3|model==4){dd = gamma*(1+eps.gamma.r*(1+cos(2*pi/24*(zt.p-phase.gamma.r))));}
          dd = dd/mean(dd)
        }
        if(outliers.removal){
          if(length(outlier.S)>0) S[outlier.S] = NA;
          if(length(outlier.M)>0) M[outlier.M] = NA;
        }
        
        ## choose to plot on top of each other or in the same scale
        S = S/mean(s, na.rm=TRUE)+2;M = M/mean(m, na.rm=TRUE)+4;s = s/mean(s)+2;m = m/mean(m)+4;m22 = m22/mean(m22)+4;
        #S = S/mean(s, na.rm=TRUE);M = M/mean(m, na.rm=TRUE);s = s/mean(s);m = m/mean(m);m22 = m22/mean(m22);
        
        plot.log.scale = FALSE
        if(plot.log.scale){
          S = log2(S0); M = log2(M0); s = log2(s0); m =log2(m0); m22 = log2(m22);
        }
        #test  = matrix(data=unlist(S), nrow=12, ncol=4, byrow=FALSE)
        kk = which(zt.p<24); zt.p = zt.p[kk];s = s[kk]; m = m[kk]; m22 = m22[kk];dd = dd[kk];
        zt = seq(0, 22, by=2);
        err.S = mean.err(S, interval = 2)[2,];mean.S = mean.err(S, interval = 2)[1,];
        err.M = mean.err(M, interval = 2)[2,];mean.M = mean.err(M, interval = 2)[1,];
        xlim = c(0, 24); 
        
        if(!RD.advantage){
          lwd=1.;cex=0.5;
          ylim = c(min(c(0)), max(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22, 6), na.rm = TRUE));
          #ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22))
          #if(ylim[1]<10^-6){ylim = range(c(s, m, dd));}
          #### relative RPKM with fitting for one day
          plot(zt, mean.S, col = 'green3', type = 'n', lwd = lwd, ylim=ylim, xlim=xlim, 
               main = gg, cex.main=1.0, xlab=NA, log='', ylab=NA, axes=FALSE);
          points(zt, mean.S, col='green3', type='p', cex=cex, pch=16);
          arrows(zt, mean.S-err.S, zt, mean.S+err.S, length=0.03, angle=90, code=3, col='green3', lwd=1.0, cex=1.0)
          #points(zt, M, type='l', lwd=0.8, col='steelblue');
          points(zt, mean.M, type='p', cex=cex, col='steelblue', pch=16);
          arrows(zt, mean.M-err.M, zt, mean.M+err.M, length=0.03, angle=90, code=3, col=col.ex, lwd=1.0, cex=1.0)
          points(zt.p, dd, type='l', lwd=1.5, col=col.deg)
          #axis(1,at=seq(0, 96, by=6), tck=-0.04, labels = NA)
          points(zt.p, s, type='l', lwd=lwd, col=col.accepted.int)
          points(zt.p, m, type='l', lwd=lwd, col=col.accepted.ex)
          #arrows(-1, max(m), zt.p[which(m==max(m))[1]], max(m), length = 0.0, angle = 90, code = 3, col=col.accepted.ex, cex=1.0, lwd=1.2, lty=2)
          #arrows(26, max(s), zt.p[which(s==max(s))[1]], max(s), length = 0.0, angle = 90, code = 3, col=col.accepted.int, cex=1.0, lwd=1.2, lty=2)
          
          if(shape.test.M3){
            m.rel = (max(m)-min(m))/2; m.phi = zt.p[which(m==max(m))][1];
            m3.cos = m.rel*cos(2*pi/24*(zt.p-m.phi))+mean(m);
            points(zt.p, m3.cos, type='l', lty=2, lwd=1.2, col='darkcyan')
          }
          abline(h=5, col=col.ex, lwd=1.0, lty=2);
          abline(h=c(2,4), col='black', lwd=1.0, lty=2)
          if(model==4 & RD.advantage) points(zt.p, m22, type='l', lwd=1., col='darkgray')
          if(model==2|model==4) abline(h=3, col=col.int, lwd=1.0, lty=2);
          if(model>2) {
            abline(h=1, col=col.deg, lwd=1.0, lty=2);
            #arrows(zt.p[which(dd==max(dd))[1]], mean(dd), zt.p[which(dd==max(dd))[1]], max(dd), length = 0.0, angle = 90, code = 3, col=col.deg, cex=1.0, lwd=1.2)
          }
          
          axis(1,at=seq(0, 24, by=4), tck=-0.02, labels = FALSE)
          axis(2, las=1, at = seq(0, 10, by=0.5), cex.axis = 1.0, tck=-.02, labels = FALSE)
          axis(4, las=1, at = seq(0, 10, by=0.5), cex.axis = 1.0, tck=-.02, labels = FALSE)
          axis(1,at=seq(0, 24, by=4),cex.axis = 1.0, lwd=0, line = -0.3)
          #axis(2, las=1,at = c, cex.axis =1.0, lwd=0, line = -0.2)
          mtext(c('0' ,'1', '2'), side = 2, line=0.3, at=c(0:2), col=col.deg, cex = 0.75);
          #mtext('0.5', side = 2, line=0.3, at=c(0.5), col=col.deg, cex = 0.5); mtext('1.5', side = 2, line=0.3, at=c(1.5), col=col.deg, cex = 0.5);
          mtext(c('0' ,'1', '2'), side = 2, line=0.3, at=c(4:6), col=col.ex, cex = 0.75);
          #mtext('0', side = 2, line=0.3, at=c(2), col=col.ex, cex = 0.6);mtext('2', side = 2, line=0.3, at=c(4), col=col.ex, cex = 0.6);
          mtext(c('0' ,'1', '2'), side = 4, line=0.3, at=c(2:4), col=col.int, cex = 0.75);
          #mtext('0', side = 4, line=0.3, at=c(1), col=col.int, cex = 0.7);mtext('2', side = 4, line=0.3, at=c(3), col=col.int, cex = 0.7);
          box();
          #abline(h=mean(S), col='gray', lwd=1.0);abline(h=mean(M), col='gray', lwd=1.0);
          dev.off()
          }else{
            lwd=1.5;cex=0.7;
            ylim = c(min(c(0)), max(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22, 6), na.rm = TRUE));
            #ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22))
            #if(ylim[1]<10^-6){ylim = range(c(s, m, dd));}
            #### relative RPKM with fitting for one day
            plot(zt, mean.S, col = 'green3', type = 'n', lwd = lwd, ylim=ylim, xlim=xlim, 
                 main = gg, cex.main=1.0, xlab=NA, log='', ylab=NA, axes=FALSE);
            points(zt, mean.S, col='green3', type='p', cex=cex, pch=16);
            arrows(zt, mean.S-err.S, zt, mean.S+err.S, length=0.03, angle=90, code=3, col='green3', lwd=1.0, cex=1.0)
            #points(zt, M, type='l', lwd=0.8, col='steelblue');
            points(zt, mean.M, type='p', cex=cex, col='steelblue', pch=16);
            arrows(zt, mean.M-err.M, zt, mean.M+err.M, length=0.03, angle=90, code=3, col=col.ex, lwd=1.0, cex=1.0)
            points(zt.p, dd, type='l', lwd=1.5, col=col.deg)
            #axis(1,at=seq(0, 96, by=6), tck=-0.04, labels = NA)
            points(zt.p, s, type='l', lwd=lwd, col=col.accepted.int)
            points(zt.p, m, type='l', lwd=lwd, col=col.accepted.ex)
            points(zt.p, m22, type='l', lwd=1., col='darkgray')
            #arrows(-1, max(m), zt.p[which(m==max(m))[1]], max(m), length = 0.0, angle = 90, code = 3, col=col.accepted.ex, cex=1.0, lwd=1.2, lty=2)
            #arrows(26, max(s), zt.p[which(s==max(s))[1]], max(s), length = 0.0, angle = 90, code = 3, col=col.accepted.int, cex=1.0, lwd=1.2, lty=2)
            
            if(shape.test.M3){
              m.rel = (max(m)-min(m))/2; m.phi = zt.p[which(m==max(m))][1];
              m3.cos = m.rel*cos(2*pi/24*(zt.p-m.phi))+mean(m);
              points(zt.p, m3.cos, type='l', lty=2, lwd=1.2, col='darkcyan')
            }
            abline(h=5, col=col.ex, lwd=1.0, lty=2);
            abline(h=c(2,4), col='black', lwd=1.0, lty=2)
            if(model==4 & RD.advantage) points(zt.p, m22, type='l', lwd=1., col='darkgray')
            if(model==2|model==4) abline(h=3, col=col.int, lwd=1.0, lty=2);
            if(model>2) {
              abline(h=1, col=col.deg, lwd=1.0, lty=2);
              #arrows(zt.p[which(dd==max(dd))[1]], mean(dd), zt.p[which(dd==max(dd))[1]], max(dd), length = 0.0, angle = 90, code = 3, col=col.deg, cex=1.0, lwd=1.2)
            }
            
            axis(1,at=seq(0, 24, by=4), tck=-0.02, labels = FALSE)
            axis(2, las=1, at = seq(0, 10, by=0.5), cex.axis = 1.0, tck=-.02, labels = FALSE)
            axis(4, las=1, at = seq(0, 10, by=0.5), cex.axis = 1.0, tck=-.02, labels = FALSE)
            axis(1,at=seq(0, 24, by=4),cex.axis = 1.0, lwd=0, line = -0.3)
            #axis(2, las=1,at = c, cex.axis =1.0, lwd=0, line = -0.2)
            mtext(c('0' ,'1', '2'), side = 2, line=0.3, at=c(0:2), col=col.deg, cex = 0.75);
            #mtext('0.5', side = 2, line=0.3, at=c(0.5), col=col.deg, cex = 0.5); mtext('1.5', side = 2, line=0.3, at=c(1.5), col=col.deg, cex = 0.5);
            mtext(c('0' ,'1', '2'), side = 2, line=0.3, at=c(4:6), col=col.ex, cex = 0.75);
            #mtext('0', side = 2, line=0.3, at=c(2), col=col.ex, cex = 0.6);mtext('2', side = 2, line=0.3, at=c(4), col=col.ex, cex = 0.6);
            mtext(c('0' ,'1', '2'), side = 4, line=0.3, at=c(2:4), col=col.int, cex = 0.75);
            #mtext('0', side = 4, line=0.3, at=c(1), col=col.int, cex = 0.7);mtext('2', side = 4, line=0.3, at=c(3), col=col.int, cex = 0.7);
            box();
            #abline(h=mean(S), col='gray', lwd=1.0);abline(h=mean(M), col='gray', lwd=1.0);
            dev.off()
            
            Old.plot = FALSE
            if(Old.plot)
            {
              ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M), (mean.M-err.M), m, m22), na.rm = TRUE);
              #ylim = range(c((mean.S+err.S), (mean.S-err.S), s, (mean.M+err.M),(mean.M-err.M), m, dd, m22))
              #if(ylim[1]<10^-6){ylim = range(c(s, m, dd));}
              #### relative RPKM with fitting for one day
              plot(zt, mean.S, col = 'green3', type = 'n', lwd = 1.0, ylim=ylim, xlim=xlim, 
                   main = gg, cex.main=0.8, xlab=NA, log='', ylab=NA, axes=FALSE);
              points(zt, mean.S, col='green3', type='p', cex=0.6, pch=16);
              arrows(zt, mean.S-err.S, zt, mean.S+err.S, length=0.03, angle=90, code=3, col='green3', lwd=1, cex=1.0)
              #points(zt, M, type='l', lwd=0.8, col='steelblue');
              points(zt, mean.M, type='p', cex=0.6, col='steelblue', pch=16);
              arrows(zt, mean.M-err.M, zt, mean.M+err.M, length=0.03, angle=90, code=3, col=col.ex, lwd=1.0, cex=1.0)
              #points(zt.p, dd, type='l', lwd=1.2, col=col.deg)
              #axis(1,at=seq(0, 96, by=6), tck=-0.04, labels = NA)
              points(zt.p, s, type='l', lwd=1.2, col=col.accepted.int)
              points(zt.p, m, type='l', lwd=1.2, col=col.accepted.ex)
              abline(h=3, col=col.ex, lwd=1.0, lty=3);
              if(model==4 & RD.advantage) points(zt.p, m22, type='l', lwd=1., col='darkgray')
              if(model==2|model==4) abline(h=2, col=col.int, lwd=1.0, lty=3);
              if(model>2) abline(h=1, col=col.deg, lwd=1.0, lty=3);
              
              axis(1,at=seq(0, 24, by=6), tck=-0.02, labels = FALSE)
              axis(1,at=seq(0, 24, by=6),cex.axis = 0.8, lwd=0, line = -0.3)
              axis(2, las=1, at = seq(0, 4, by=0.5), cex.axis = 0.7, tck=-.02, labels = FALSE)
              axis(4, las=1, at = seq(0, 4, by=0.5), cex.axis = 0.7, tck=-.02, labels = FALSE)
              #axis(2, las=1,at = c, cex.axis =1.0, lwd=0, line = -0.2)
              #mtext('1', side = 2, line=0.3, at=c(1), col=col.deg, cex = 0.7);
              #mtext('0.5', side = 2, line=0.3, at=c(0.5), col=col.deg, cex = 0.5); mtext('1.5', side = 2, line=0.3, at=c(1.5), col=col.deg, cex = 0.5);
              mtext('1', side = 2, line=0.3, at=c(3), col=col.ex, cex = 0.7);
              #mtext('2', side = 2, line=0.3, at=c(4), col=col.ex, cex = 0.7);
              mtext('0', side = 2, line=0.3, at=c(2), col=col.ex, cex = 0.6);
              #mtext('0', side = 4, line=0.3, at=c(1), col=col.int, cex = 0.7);
              mtext('1', side = 4, line=0.3, at=c(2), col=col.int, cex = 0.7);
              mtext('2', side = 4, line=0.3, at=c(3), col=col.int, cex = 0.7);
              box();
              #abline(h=mean(S), col='gray', lwd=1.0);abline(h=mean(M), col='gray', lwd=1.0);
              dev.off()
            }
            
        }
        
        ### summary plot
        if(summary.plot & model>1)
        {
          w = 2*pi/24; 
          #m.rel = T$rel.amp.rpkm.mRNA[i]; m.phi = T$phase.rpkm.mRNA[i];
          #m.rel = (max(m-2)/min(m-2)-1)/(max(m-2)/min(m-2)+1); m.phi = T$phase.rpkm.mRNA[i];
          m.rel = (max(m)-min(m))/2; m.phi = zt.p[which(m==max(m))][1];
          if(model==2){
            #s.rel = T$rel.amp.rpkm.premRNA[i]; s.phi = T$phase.rpkm.premRNA[i];
            s.rel = (max(s)-min(s))/2; s.phi = zt.p[which(s==max(s))][1];
            m2.rel = 0;m2.phi = 0;
            d.rel = 0; d.phi= 0;
          }
          if(model==3){
            s.rel = 0; s.phi = 0;
            m2.rel = s.rel/sqrt(w^2/gamma^2+1);m2.phi = (s.phi + atan2(w, gamma)/w)%%24;
            d.rel = eps.gamma; d.phi= phase.gamma;
          }
          if(model==4){
            #s.rel = T$rel.amp.rpkm.premRNA[i]; s.phi = T$phase.rpkm.premRNA[i];
            s.rel = (max(s)-min(s))/2; s.phi = zt.p[which(s==max(s))][1];
            m2.rel = s.rel/sqrt(w^2/gamma^2+1);m2.phi = (s.phi + atan2(w, gamma)/w)%%24;
            d.rel = eps.gamma; d.phi= phase.gamma;
          }
          
          cex.text = 0.6;
          if(RD.advantage){
            ### plot examples in Figure 2 to illustrate M2, M3 and M4
            pdfname = paste(folder, '/Summary_plot_', gg, '.pdf', sep='')
            pdf(pdfname, width = 1.2, height = 1.2)
            par(cex = 0.7, las = 1, tcl = -0.3)
            par(mgp = c(0,0,0), mar = c(0.,0.2,0.5,0.), pty='s');
            
            #ref = max(c(m.rel, m2.rel, s.rel));
            ref = max(c(m.rel, m2.rel, s.rel, d.rel));
            s.rel = s.rel/ref; m.rel = m.rel/ref; m2.rel = m2.rel/ref; d.rel = d.rel/ref;
            #ref2 = max(c(m.rel, m2.rel, s.rel));
            ref2 = max(c(m.rel, m2.rel, s.rel, d.rel));
            s.a = s.rel*cos(w*s.phi);s.b = s.rel*sin(w*s.phi);CC = (s.a -1i*s.b) * exp(1i*pi/2);s.a = Re(CC); s.b=Im(CC);
            m.a = m.rel*cos(w*m.phi);m.b = m.rel*sin(w*m.phi);CC = (m.a -1i*m.b) * exp(1i*pi/2);m.a = Re(CC); m.b=Im(CC);
            m2.a = m2.rel*cos(w*m2.phi);m2.b = m2.rel*sin(w*m2.phi);CC = (m2.a -1i*m2.b) * exp(1i*pi/2);m2.a = Re(CC); m2.b=Im(CC);
            
            d.a = d.rel*cos(w*d.phi);d.b = d.rel*sin(w*d.phi);CC = (d.a -1i*d.b) * exp(1i*pi/2);d.a = Re(CC); d.b=Im(CC);
            
            rmm=ref2+0.3; rr=c(-rmm,rmm); xlim = rr; ylim = rr;
            
            plot(s.a, s.b, main=NA, cex.main = 0.7, type='n', xlim=xlim, ylim=ylim, axes=FALSE, xlab='', ylab='')
            phi=seq(0,2*pi,len=1000); motif.cycle = ref2;lwd=1.5;
            lines(motif.cycle*cos(phi), motif.cycle*sin(phi), col='black', lwd=lwd);
            ticks.ll = 0.07;cex=1.0;
            
            text(0, 1, '0', col='black', cex = cex, offset=0.2, pos=3); arrows(0, 1, 0, (1-ticks.ll), col='black', lwd=1.2, length=0.0);
            text(1, 0, '6', col='black', cex = cex, offset=0.2, pos=4);arrows(1, 0, (1-ticks.ll), 0, col='black', lwd=1.2, length=0.0);
            text(-1, 0, '18', col='black', cex = cex, offset=0.2, pos=2);arrows(-1, 0, (-1+ticks.ll),0, col='black', lwd=1.2, length=0.0);
            text(0, -1, '12', col='black', cex = cex, offset=0.2, pos=1);arrows(0, -1, 0, (-1+ticks.ll), col='black', lwd=1.2, length=0.0);
            ll = 0.07; lwd=2.5;
            
            if(model==4 & RD.advantage) arrows(0, 0, m2.a, m2.b,  col='darkgray', lwd=lwd, lty=1, length=ll);
            if(model==2|model==4) arrows(0, 0, s.a, s.b,  col=col.accepted.int, lwd=lwd, length=ll);
            if(model>2) arrows(0, 0, d.a, d.b,  col=col.deg, lwd=lwd, length=ll);
            arrows(0, 0, m.a, m.b,  col=col.accepted.ex, lwd=lwd, length=ll);
            text(0.8, 1.15, paste('hl = ', signif(log(2)/gamma, d=2), 'h', sep=''), col='blue',  cex=0.8)
            #text(-0.8, 1.25, gg, cex=0.7)
            
            dev.off();
            
          }else{
            ### plot examples in Figure 5 to illustrate roles of RD: amplifiers and phase tuners
            pdfname = paste(folder, '/Summary_plot_', gg, '.pdf', sep='')
            pdf(pdfname, width = 1.2, height = 1.2)
            par(cex = 0.7, las = 1, tcl = -0.3)
            par(mgp = c(0,0,0), mar = c(0.,0.2,0.5,0.), pty='s');
            #ref = max(c(m.rel, m2.rel, d.rel, s.rel, 1))
            ref = max(c(m.rel, m2.rel, s.rel, d.rel));
            s.rel = s.rel/ref; m.rel = m.rel/ref; m2.rel = m2.rel/ref; d.rel = d.rel/ref;
            ref2 = max(c(m.rel, m2.rel, s.rel, d.rel));
            s.a = s.rel*cos(w*s.phi);s.b = s.rel*sin(w*s.phi);CC = (s.a -1i*s.b) * exp(1i*pi/2);s.a = Re(CC); s.b=Im(CC);
            m.a = m.rel*cos(w*m.phi);m.b = m.rel*sin(w*m.phi);CC = (m.a -1i*m.b) * exp(1i*pi/2);m.a = Re(CC); m.b=Im(CC);
            m2.a = m2.rel*cos(w*m2.phi);m2.b = m2.rel*sin(w*m2.phi);CC = (m2.a -1i*m2.b) * exp(1i*pi/2);m2.a = Re(CC); m2.b=Im(CC);
            d.a = d.rel*cos(w*d.phi);d.b = d.rel*sin(w*d.phi);CC = (d.a -1i*d.b) * exp(1i*pi/2);d.a = Re(CC); d.b=Im(CC);
            
            #rmm=1.0+0.3; rr=c(-rmm,rmm); xlim = rr; ylim = rr;
            rmm=ref2+0.3; rr=c(-rmm,rmm); xlim = rr; ylim = rr;
            plot(s.a, s.b, main=NA, cex.main = 0.7, type='n', xlim=xlim, ylim=ylim, axes=FALSE, xlab='', ylab='')
            phi=seq(0,2*pi,len=1000); motif.cycle = 1;lwd=1.5;
            lines(motif.cycle*cos(phi), motif.cycle*sin(phi), col='black', lwd=lwd);
            ticks.ll = 0.07;cex=1.0;
            text(0, 1, '0', col='black', cex = cex, offset=0.2, pos=3); arrows(0, 1, 0, (1-ticks.ll), col='black', lwd=1.2, length=0.0);
            text(1, 0, '6', col='black', cex = cex, offset=0.2, pos=4);arrows(1, 0, (1-ticks.ll), 0, col='black', lwd=1.2, length=0.0);
            text(-1, 0, '18', col='black', cex = cex, offset=0.2, pos=2);arrows(-1, 0, (-1+ticks.ll),0, col='black', lwd=1.2, length=0.0);
            text(0, -1, '12', col='black', cex = cex, offset=0.2, pos=1);arrows(0, -1, 0, (-1+ticks.ll), col='black', lwd=1.2, length=0.0);
            ll = 0.07; lwd=2.5;
            if(model==4 & RD.advantage) arrows(0, 0, m2.a, m2.b,  col='darkgray', lwd=lwd, lty=1, length=ll);
            if(model==2|model==4) arrows(0, 0, s.a, s.b,  col=col.accepted.int, lwd=lwd, length=ll);
            if(model>2) arrows(0, 0, d.a, d.b,  col=col.deg, lwd=lwd, length=ll);
            arrows(0, 0, m.a, m.b,  col=col.accepted.ex, lwd=lwd, length=ll);
            if(gg!='Fus' & gg!='Arntl' & gg!='Cry1')
            text(0.8, 1.15, paste('hl = ', signif(log(2)/gamma, d=2), 'h', sep=''), col='blue',  cex=cex)
            #text(-0.8, 1.25, gg, cex=0.7)
            
            dev.off();
          }
          
        }
      }else{
        ### Figure = FALSE
        if(!RF.KO){
          par(mfcol = c(2,2))
        }else{
          par(mfcol = c(2,3))
        }
        ### absolute plots
        xlim = c(0, 96)
        ylim = range(c(S, M))
        plot(zt, S, col = 'green3', type = 'l', lwd = 1.5, pch = 1, ylim=ylim, xlim=xlim, log='y', main=Tt$gene[i], ylab='absolute signals(log)')
        points(zt, S, col='green3', type='p', cex=0.8)
        points(zt, M, type='l', lwd=1.5, col='steelblue')
        points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
        
        ylim = range(c(S, M, s, m))
        plot(zt, S, col = 'green3', type = 'l', lwd = 0.7, pch = 1, ylim=ylim, xlim=xlim, log='y', main=paste('Model ', best.model, ' ( ', signif(Tt$prob.best.model[i]*100, 2), '%)', sep=''), ylab='absolute signals(log)')
        points(zt, S, col='green3', type='p', cex=0.8)
        points(zt.p, s, type='l', lwd=2.0, col='limegreen')
        points(zt, M, type='l', lwd=0.7, col='steelblue')
        points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
        points(zt.p, m, type='l', lwd=2.0, col='steelblue2')
        
        ### relative plots
        S = S/mean(S)
        M = M/mean(M)
        s = s/mean(s)
        m = m/mean(m)
        
        xlim = c(0, 96)
        ylim = range(c(S, M))
        plot(zt, S, col = 'green3', type = 'l', lwd = 1.5, pch = 1, ylim=ylim, xlim=xlim, log='', ylab='relative signals(linear)')
        points(zt, S, col='green3', type='p', cex=0.8)
        points(zt, M, type='l', lwd=1.5, col='steelblue')
        points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
        
        ylim = c(0, 2.5)
        #ylim = range(c(S, M, s, m))
        #if(ylim[1]>0) ylim[1]=0;
        #if(ylim[2]<2.5) ylim[2]=2.5;
        plot(zt, S, col = 'green3', type = 'l', lwd = 0.7, pch = 1, ylim=ylim, xlim=xlim, log='', ylab='relative signals(linear)')
        points(zt, S, col='green3', type='p', cex=0.8)
        points(zt.p, s, type='l', lwd=2.0, col='limegreen')
        points(zt, M, type='l', lwd=0.7, col='steelblue')
        points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
        points(zt.p, m, type='l', lwd=2.0, col='steelblue2')
        if(model==2)
        {
          text(18, ylim[2]-0.1, paste('splicing = ', signif(log(2)/splicing.k*60, d=2), 'min'), col='black')
          if(abs(gamma-gamma.max)<0.00001 | abs(gamma-gamma.min)<0.00001){
            text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/gamma, d=2), ' h'), col='red')
          }else{
            text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/gamma, d=2), ' h'), col='black')
          }
          text(32, ylim[2]-0.5, paste(hl.ref), col='blue')
        }
        if(model==3|model==4)
        {
          dd = gamma(1+ eps.gamma*(1+cos(2*pi/24*(zt.p-phase.gamma))))
          dd = dd/mean(dd)
          points(zt.p, dd, type='l', lwd=2.0, col='red')
          
          text(18, ylim[2]-0.1, paste('splicing = ', signif(log(2)/splicing.k*60, d=2), ' min'), col='black')
          
          if(abs(gamma-gamma.max)<0.00001 | abs(gamma-gamma.min)<0.00001){
            text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/(gamma+amp.gamma), d=2),  '-', signif(log(2)/gamma, d=2), 'h'), col='red')
          }else{
            text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/(gamma+amp.gamma), d=2),  '-', signif(log(2)/gamma, d=2), 'h'), col='black')
          }
          text(32, ylim[2]-0.5, paste(hl.ref), col='blue')
        }
        
        if(RF.KO)
        {
          xlim = c(0, 96)
          ylim = range(c(S0, S1, S2))
          plot(zt, S0, col = 'green3', type = 'l', lwd = 1.5, pch = 1, ylim=ylim, xlim=xlim, log='y', main='premRNA', ylab='absolute signals(log)')
          points(zt, S0, col='green3', type='p', cex=0.8)
          points((c(0:23)*4+2), S1, type='l', lwd=1.2, col='green3', lty=2)
          points((c(0:23)*4+2), S1, type='p', cex=0.7, col='green3', lty=2)
          points((c(12:23)*4+2), S2, type='l', lwd=1.2, col='black', lty=2)
          points((c(12:23)*4+2), S2, type='p', cex=0.7, col='black', lty=2)
          cols = c('green3', 'green3', 'black')
          legend('topright', legend = c('WT.Ad','WT.RF', 'KO.RF') , cex=1.2, pch=c(1,2,2), col=cols, pt.cex=1.0, pt.lwd=1.2, pt.bg= cols, border = NA, bty = 'n')
          
          xlim = c(0, 96)
          ylim = range(c(M0, M1, M2))
          plot(zt, M0, col = 'steelblue', type = 'l', lwd = 1.5, pch = 1, ylim=ylim, xlim=xlim, log='y', main='mRNA', ylab='absolute signals(log)')
          points(zt, M0, col='steelblue', type='p', cex=0.8)
          points((c(0:23)*4+2), M1, type='l', lwd=1.2, col='steelblue', lty=2)
          points((c(0:23)*4+2), M1, type='p', cex=0.7, col='steelblue', lty=2)
          points((c(12:23)*4+2), M2, type='l', lwd=1.2, col='darkgray', lty=2)
          points((c(12:23)*4+2), M2, type='p', cex=0.7, col='darkgray', lty=2)
          cols = c('steelblue', 'steelblue', 'darkgray')
          legend('topright', legend = c('WT.Ad','WT.RF', 'KO.RF') , cex=1.2, pch=c(1,2,2), col=cols, pt.cex=1.0, pt.lwd=1.2, pt.bg= cols, border = NA, bty = 'n')
          
        }
      }

    }
    if(!Figure)
    {
      dev.off()
    }
}

mean.nona = function(x)
{	
  kk = which(!is.na(x))
  if(length(kk)>0) 
  {
    xx = mean(x[which(!is.na(x)==TRUE)]);
  }else{
    xx = NA;
  }
  return(xx)
}
sme.nona = function(x)
{
  kk = which(!is.na(x))
  if(length(kk)>1)
  {
    return(sd(x[kk])/sqrt(length(kk)))
  }else{
    return(NA)
  }
}
mean.err = function(x, period=24, interval=3)
{
  x = as.numeric(x)
  kk = which(!is.na(x)==TRUE)
  index = period/interval;
  nn = length(x)/(period/interval)
  test = c()
  for(mm in 1:nn) test = cbind(test, x[c(1:index)+index*(mm-1)])
  test = data.frame(test)
  
  mean = apply(test, 1, mean.nona)
  err = apply(test, 1, sme.nona)
  
  return(rbind(mean=mean, err=err))
}

mean.replicates = function(x, period=24, interval=2)
{
  x = as.numeric(x)
  kk = which(!is.na(x)==TRUE)
  index = period/interval;
  nn = length(x)/(period/interval)
  test = c()
  for(mm in 1:nn) test = cbind(test, x[c(1:index)+index*(mm-1)])
  test = data.frame(test)
  
  xx = apply(test, 1, mean.nona)
  xx = (xx-mean(xx))/sd(xx)
  #err = apply(test, 1, sme.nona)
  return(xx)
}

plot.half.life.comparison = function(index, pdfname, T)
{
    index = index[which(!is.na(T$BIC.best.model[index])==TRUE)]
    ### half-lives database
    load(file='myRdata/mRNAs_half_lives_databases.Rdata')
    #load(file='myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata')
    
    pdf(pdfname, width=12.0, height=12.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    par(mfrow=c(3,3))
    
    for(model in c(2:4))
    {
        kk = index[which(T$BIC.best.model==model)]
        gg = T$gene[kk]
        if(model==2)
        {
            xx = eval(parse(text=paste('log(2)/T$gamma.m', model, '[kk]', sep='')))
        }else{
            xx = eval(parse(text=paste('log(2)/(T$gamma.m', model, '[kk]+0.5*T$amp.gamma.m', model, '[kk])', sep='')))
        }
        
        mm = match(gg, fried[,1])
        jj = which(!is.na(mm)==TRUE)
        mm = mm[which(!is.na(mm)==TRUE)]
        xx1 = xx[jj]
        yy1 = as.numeric(fried[mm, 2])
        
        cex = 0.7
        plot(xx1, yy1, cex=cex, xlab='half-lives estimated', ylab='fried', main=paste('model', model, ', R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        #cor(xx1, yy1)
        abline(0, 1, lwd=2.0, col='red')
        
        mm = match(gg, shar[,1])
        jj = which(!is.na(mm)==TRUE)
        mm = mm[which(!is.na(mm)==TRUE)]
        xx1 = xx[jj]
        yy1 = as.numeric(shar[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='half-lives estimated', ylab='shar', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        #cor(xx1, yy1)
        abline(0, 1, lwd=2.0, col='red')
        
        mm = match(gg, schwa[,1])
        jj = which(!is.na(mm)==TRUE)
        mm = mm[which(!is.na(mm)==TRUE)]
        xx1 = xx[jj]
        yy1 = as.numeric(schwa[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='half-lives estimated', ylab='schwa', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        abline(0, 1, lwd=2.0, col='red')
        #cor(xx1, yy1)
    }
    
    dev.off()


}

plot.genes.examples.conditions = function(index=c(1),pdfname='Examples.pdf', Tt=T, RF.KO = TRUE, Figure=FALSE, folder=folder)
{
    #pdfname='Examples_Test.pdf'; Tt=T; RF.KO = TRUE; Figure=FALSE; folder=folder
    #i = 11; Tt=T;
    #### plot to check fitting with examples
    zt = seq(0,94,by = 2);
    zt.p = seq(0,96,by = 0.5)
    index = index[which(!is.na(Tt$BIC.best.model[index])==TRUE)]
    
    ### half-lives database
    load(file='myRdata/mRNAs_half_lives_databases.Rdata')
    load(file='myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata')
    bounds = set.bounds(model=4)
    upper = bounds$upper
    lower = bounds$lower
    gamma.max = upper[1]
    gamma.min = lower[1];
    
    if(!Figure)
    {
        pdfname = paste(folder, '/', 'Examples.pdf', sep='')
        
        if(RF.KO)
        {
            pdf(pdfname, width = 12, height = 8)
        }else{
            pdf(pdfname, width = 12, height = 8)
        }
    }
    
    for(i in index)
    {
        ##### extract all information of fitting and data
        cat(i, '...');
        cat(as.character(Tt$gene[i]), '\n');
        gg = Tt$gene[i];
        mm = match(gg, shar[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl1 = signif(shar[mm,2], d=2)}else{hl1 = NA;}
        mm = match(gg, fried[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl2 = signif(fried[mm,2], d=2)}else{hl2 = NA;}
        mm = match(gg, schwa[,1]);mm = mm[which(!is.na(mm))];if(length(mm)>0){ hl3 = signif(schwa[mm,2], d=2)}else{hl3 = NA;}
        
        hl.ref = paste(hl1, 'h(shar); ', hl2, 'h(fried); ', hl3, 'h(schwar)',  sep='', collapse=';')
        
        best.model = Tt$BIC.best.model[i]
        mm = grep(paste('.m', best.model, sep=''),colnames(Tt))
        mm = mm[-c(1, length(mm))]
        #mm = mm[1:(length(mm)/2)]
        par = unlist(Tt[i, mm])
        M = unlist(Tt[i, grep('.abs.mRNA', colnames(Tt))]);
        S = unlist(Tt[i, grep('.abs.premRNA', colnames(Tt))]);
        M0 = M;
        S0 = S;
        
        kk = which(R.WT.RF$gene==gg)[1]
        if(length(kk)>0) {
            M1 = unlist(R.WT.RF[kk, grep('.abs.mRNA', colnames(R.WT.RF))]);
            S1 = unlist(R.WT.RF[kk, grep('.abs.premRNA', colnames(R.WT.RF))]);
        }else{
            M1 = rep(NA, 24)
            S1 = rep(NA, 24)
        }
        kk = which(R.KO.RF$gene==gg)[1]
        if(length(kk)>0) {
            M2 = unlist(R.KO.RF[kk, grep('.abs.mRNA', colnames(R.KO.RF))]);
            S2 = unlist(R.KO.RF[kk, grep('.abs.premRNA', colnames(R.KO.RF))]);
        }else{
            M2 = rep(NA, 12)
            S2 = rep(NA, 12)
        }
        
        gamma = par[1];
        model = best.model;
        if(model==1)
        {
            s = rep(exp(mean(log(S))), length(zt.p));
            m = rep(exp(mean(log(M))), length(zt.p));
        }else{
            if(model==2)
            {
                amp.gamma = 0.0;
                phase.gamma = 0.0;
                splicing.k = par[2];
                param.synthesis.1 = par[3];
                param.synthesis.2 = par[4];
                param.synthesis.3 = par[5];
                param.synthesis.4 = par[6];
            }
            if(model==3)
            {
                amp.gamma = par[2];
                phase.gamma = par[3];
                splicing.k = par[4];
                param.synthesis.1 = par[5];
                param.synthesis.2 = 0;
                param.synthesis.3 = 0;
                param.synthesis.4 = 1;
            }
            if(model==4)
            {
                amp.gamma = par[2];
                phase.gamma = par[3];
                splicing.k = par[4];
                param.synthesis.1 = par[5];
                param.synthesis.2 = par[6];
                param.synthesis.3 = par[7];
                param.synthesis.4 = par[8];
            }
            
            #print(c(param.synthesis.3, phase.gamma))
            s = compute.s.beta(t = zt.p, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
            m = compute.m.beta(t = zt.p, gamma, amp.gamma, phase.gamma, splicing.k, param.synthesis.1, param.synthesis.2, param.synthesis.3, param.synthesis.4)
            
        }
        
        if(!Figure)
        {
            if(!RF.KO){
                par(mfcol = c(2,1))
            }else{
                par(mfcol = c(2,2))
            }
            ### absolute plots
            xlim = c(0, 96)
            ylim = range(c(S, M))
            plot(zt, S, col = 'green3', type = 'l', lwd = 1.5, pch = 1, ylim=ylim, xlim=xlim, log='y', main=Tt$gene[i], ylab='absolute signals(log)')
            points(zt, S, col='green3', type='p', cex=0.8)
            points(zt, M, type='l', lwd=1.5, col='steelblue')
            points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
            
            ### relative plots
            S = S/mean(S)
            M = M/mean(M)
            s = s/mean(s)
            m = m/mean(m)
            
            ylim = c(0, 2.5)
            #ylim = range(c(S, M, s, m))
            #if(ylim[1]>0) ylim[1]=0;
            #if(ylim[2]<2.5) ylim[2]=2.5;
            plot(zt, S, col = 'green3', type = 'l', lwd = 0.7, pch = 1, ylim=ylim, xlim=xlim, log='', ylab='relative signals(linear)')
            points(zt, S, col='green3', type='p', cex=0.8)
            points(zt.p, s, type='l', lwd=2.0, col='limegreen')
            points(zt, M, type='l', lwd=0.7, col='steelblue')
            points(zt, M, type='p', cex=0.8, col='steelblue', pch=5)
            points(zt.p, m, type='l', lwd=2.0, col='steelblue2')
            if(model==2)
            {
                text(18, ylim[2]-0.1, paste('splicing = ', signif(log(2)/splicing.k*60, d=2), 'min'), col='black')
                if(abs(gamma-gamma.max)<0.00001 | abs(gamma-gamma.min)<0.00001){
                    text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/gamma, d=2), ' h'), col='red')
                }else{
                    text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/gamma, d=2), ' h'), col='black')
                }
                text(32, ylim[2]-0.5, paste(hl.ref), col='blue')
            }
            if(model==3|model==4)
            {
                dd = gamma + amp.gamma*(1+cos(2*pi/24*(zt.p-phase.gamma))/2)
                dd = dd/mean(dd)
                points(zt.p, dd, type='l', lwd=2.0, col='red')
                
                text(18, ylim[2]-0.1, paste('splicing = ', signif(log(2)/splicing.k*60, d=2), ' min'), col='black')
                
                if(abs(gamma-gamma.max)<0.00001 | abs(gamma-gamma.min)<0.00001){
                    text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/(gamma+amp.gamma), d=2),  '-', signif(log(2)/gamma, d=2), 'h'), col='red')
                }else{
                    text(18, ylim[2]-0.3, paste('half-life = ', signif(log(2)/(gamma+amp.gamma), d=2),  '-', signif(log(2)/gamma, d=2), 'h'), col='black')
                }
                text(32, ylim[2]-0.5, paste(hl.ref), col='blue')
            }
            
            if(RF.KO)
            {
                xlim = c(0, 96)
                ylim = range(c(S0, S1, S2))
                plot(zt, S0, col = 'green3', type = 'l', lwd = 1.5, lty=2, ylim=ylim, xlim=xlim, log='y', main='premRNA', ylab='(log2)')
                points(zt, S0, col='green3', type='p', cex=0.8)
                points((c(0:23)*4+2), S1, type='l', lwd=1.5, col='red')
                points((c(0:23)*4+2), S1, type='p', cex=0.8, col='red')
                points((c(12:23)*4+2), S2, type='l', lwd=1.5, col='black')
                points((c(12:23)*4+2), S2, type='p', cex=0.8, col='black')
                cols = c('green3', 'green3', 'black')
                legend('topright', legend = c('WT.Ad','WT.RF', 'KO.RF') , cex=1.2, pch=c(1,1,1), col=c('green3', 'red', 'black'), pt.cex=1.0, pt.lwd=1.2, pt.bg= cols, border = NA, bty = 'n')
                
                xlim = c(0, 96)
                ylim = range(c(M0, M1, M2))
                plot(zt, M0, col = 'steelblue', type = 'l', lwd = 1.5, pch = 1, lty=2, ylim=ylim, xlim=xlim, log='y', main='mRNA', ylab='(log2)')
                points(zt, M0, col='steelblue', type='p', cex=0.8)
                points((c(0:23)*4+2), M1, type='l', lwd=1.5, col='red')
                points((c(0:23)*4+2), M1, type='p', cex=0.8, col='red')
                points((c(12:23)*4+2), M2, type='l', lwd=1.5, col='black')
                points((c(12:23)*4+2), M2, type='p', cex=0.8, col='black')
                cols = c('steelblue', 'steelblue', 'darkgray')
                legend('topright', legend = c('WT.Ad','WT.RF', 'KO.RF') , cex=1.2, pch=c(1,1,1), col=c('steelblue', 'red', 'black'), pt.cex=1.0, pt.lwd=1.2, pt.bg= cols, border = NA, bty = 'n')
                
            }
            
        }else{
            
            pdfname = paste(folder, '/fitting_results_', gg, '.pdf', sep='')
            pdf(pdfname, width = 2.5, height = 3.2)
            par(cex = 0.7, las = 1, mgp = c(1.2,0.4,0), mar = c(2.2,2.5,1.,0.8)+0.1, tcl = -0.3)
            par(mfcol = c(2,1))
            
            ### absolute plot
            xlim = c(0, 96)
            ylim = range(c(S, M))
            plot(zt, S, col = 'green3', type = 'l', lwd = 1., pch = 1, ylim=ylim, xlim=xlim, log='y', main=gg, cex.main=0.7, cex.lab=0.7, xlab=NA, ylab='absolute signals(log)', axes=FALSE)
            points(zt, S, col='green3', type='p', cex=0.5)
            points(zt, M, type='l', lwd=1., col='steelblue')
            points(zt, M, type='p', cex=0.5, col='steelblue', pch=5)
            
            axis(1,at=seq(0, 96, by=12),cex.axis =0.6)
            #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
            lims = signif(lims, d=1)
            by = signif((lims[2]-lims[1])/4,d=1)
            #print(gene)
            #print(lims)
            axis(2,las=1,cex.axis = 0.6)
            box();
            
            
            #### relative plot with fitting
            S = S/mean(S)
            M = M/mean(M)
            s = s/mean(s)
            m = m/mean(m)
            
            ylim = c(0, 2.5)
            #ylim = range(c(S, M, s, m))
            #if(ylim[1]>0) ylim[1]=0;
            #if(ylim[2]<2.5) ylim[2]=2.5;
            cols = c('gray','green3','tomato','black')
            mains = c('CS-CD', 'RS-CD', 'CS-RD', 'RS-RD')
            
            plot(zt, S, col = 'green3', type = 'l', lwd = 0.7, pch = 1, ylim=ylim, xlim=xlim, main = mains[model], col.main=cols[model], cex.main=0.7, cex.lab=0.7, xlab='ZT [hr]', log='', ylab='relative signals(linear)', axes=FALSE)
            points(zt, S, col='green3', type='p', cex=0.5)
            points(zt.p, s, type='l', lwd=1.5, col='limegreen')
            points(zt, M, type='l', lwd=0.7, col='steelblue')
            points(zt, M, type='p', cex=0.5, col='steelblue', pch=5)
            points(zt.p, m, type='l', lwd=1.5, col='steelblue2')
            cex.text = 0.5;
            start.text = 72
            if(model==2)
            {
                text(start.text, ylim[2]-0.1, paste('splicing = ', signif(log(2)/splicing.k*60, d=2), 'min'), col='orange', cex=cex.text)
                text(start.text, ylim[2]-0.35, paste('half-life = ', signif(log(2)/gamma, d=2), ' h'), col='orange', cex=cex.text)
                
            }
            if(model==3|model==4)
            {
                dd = gamma + amp.gamma*(1+cos(2*pi/24*(zt.p-phase.gamma))/2)
                dd = dd/mean(dd)
                points(zt.p, dd, type='l', lwd=2.0, col='red')
                
                text(start.text, ylim[2]-0.1, paste('splicing.time = ', signif(log(2)/splicing.k*60, d=2), ' min'), col='orange',  cex=cex.text)
                text(start.text, ylim[2]-0.35, paste('half.life = ', signif(log(2)/(gamma+amp.gamma), d=2),  '-', signif(log(2)/gamma, d=2), 'h'), col='red',  cex=cex.text)
            }
            
            axis(1,at=seq(0, 96, by=12),cex.axis =0.6)
            #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
            lims = signif(lims, d=1)
            by = signif((lims[2]-lims[1])/4,d=1)
            #print(gene)
            #print(lims)
            axis(2,las=1,cex.axis = 0.6)
            box();
            
            
            dev.off()
        }
        
    }
    
    if(!Figure)
    {
        dev.off()
    }
    
}


plot.summary.comparison.AL.RF.KO = function(T, model=4, qv = 1, folder=folder)
{
    #pdfname='Examples_Test.pdf'; Tt=T; RF.KO = TRUE; Figure=FALSE; folder=folder; qv = 0.05;
    #i = 11; Tt=T; model = 4;
    #zt = seq(2,94,by = 4);
    #load(file='/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata')
    load(file=paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/Cry_Bmal_WT_AD_NRF_Bmal_KO_RF.Rdata', sep=''))
    
    index = which(T$BIC.best.model==model & T$qv.rpkm.mRNA<qv);
    ggs = T$gene[index];
    Recompute.statistics = FALSE
    if(Recompute.statistics)
    {
      ii = grep('abs.mRNA', colnames(T))
      ii = ii[which(c(1:length(ii))%%2==0)]
      jj = grep('abs.premRNA', colnames(T))
      jj = jj[which(c(1:length(jj))%%2==0)]
      data1 = (as.matrix(T[index,ii]))
      data2 = (as.matrix(T[index,jj]))
      
      stat1 = t(apply(log2(data1),1, f24_R2_alt2, t=zt))
      stat2 = t(apply(log2(data2),1, f24_R2_alt2, t=zt))
      stat1[, 4] = t(apply((data1),1, f24_R2_alt2, t=zt))[,4]
      stat2[, 4] = t(apply((data2),1, f24_R2_alt2, t=zt))[,4]
      
      kk = match(T$gene[index], R.WT.RF$gene)
      ii = grep('abs.mRNA', colnames(R.WT.RF))
      jj = grep('abs.premRNA', colnames(R.WT.RF))
      data3 = (as.matrix(R.WT.RF[kk,ii]))
      data4 = (as.matrix(R.WT.RF[kk,jj]))
      
      stat3 = t(apply(log2(data3),1, f24_R2_alt2, t=zt))
      stat4 = t(apply(log2(data4),1, f24_R2_alt2, t=zt))
      stat3[, 4] = t(apply((data3),1, f24_R2_alt2, t=zt))[,4]
      stat4[, 4] = t(apply((data4),1, f24_R2_alt2, t=zt))[,4]
    }
    
    pdfname = paste(folder, '/', 'Comparison_AL_RF_KO_summary.pdf', sep='')
    pdf(pdfname, width = 14, height = 6)
    
    par(mfcol = c(2,4), pty='s')
    ## mRNA
    #plot(stat1[,2], stat3[,2], cex=0.4, main='Mean Expression', ylab='mRNA (NRF)', xlab='ALF')
    #abline(0, 1, lwd=2.0, col='red')
    index = match(ggs, R.WT.AD$gene);
    index = index[which(!is.na(index)==TRUE)];
    plot((R.WT.AD$phase.mRNA[index]), (R.WT.RF$phase.mRNA[index]), cex=0.4, 
         main='mRNA phases', ylab='NRF', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    #abline(v=2, lwd=2.0, col='red');abline(h=2, lwd=2.0, col='red')
    plot((R.WT.AD$phase.premRNA[index]), (R.WT.RF$phase.premRNA[index]), cex=0.4, 
         main='premRNA phase', ylab='NRF', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    plot((R.WT.AD$rel.amp.mRNA[index]), (R.WT.RF$rel.amp.mRNA[index]), cex=0.4, log='xy',
         main='mRNA rel.amp', ylab='NRF', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    plot((R.WT.AD$rel.amp.premRNA[index]), (R.WT.RF$rel.amp.premRNA[index]), cex=0.4, log='xy',
         main='premRNA rel.amp', ylab='NRF', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    ### WT vs KO in NRF
    plot((R.WT.RF$phase.mRNA[index]), (R.KO.RF$phase.mRNA[index]), cex=0.4, 
         main='mRNA phases',ylab='KO.NRF', xlab='WT.NRF')
    abline(0, 1, lwd=2.0, col='red')
    #abline(v=2, lwd=2.0, col='red');abline(h=2, lwd=2.0, col='red')
    plot((R.WT.RF$phase.premRNA[index]), (R.KO.RF$phase.premRNA[index]), cex=0.4, 
         main='premRNA phase', ylab='KO.NRF', xlab='WT.NRF')
    abline(0, 1, lwd=2.0, col='red')
    plot((R.WT.RF$rel.amp.mRNA[index]), (R.KO.RF$rel.amp.mRNA[index]), cex=0.4, log='xy',
         main='mRNA rel.amp', ylab='KO.NRF', xlab='WT.NRF')
    abline(0, 1, lwd=2.0, col='red')
    plot((R.WT.RF$rel.amp.premRNA[index]), (R.KO.RF$rel.amp.premRNA[index]), cex=0.4, log='xy',
         main='premRNA rel.amp', ylab='KO.NRF', xlab='WT.NRF')
    abline(0, 1, lwd=2.0, col='red') 
   
    dev.off()
}

model.sel.wt.ko.allModel = function(data.wt.ko, time.wt=c(0:11)*4, time.ko=c(0:11)*4, period=24)
{
  #index =which(keep$gene=='Per2'); data.wt.ko =  D[index, ]; time.wt=c(0:11)*4; time.ko=c(0:11)*4;period=24;
  wt = as.numeric(data.wt.ko[c(1:length(time.wt))])
  ko = as.numeric(data.wt.ko[-c(1:length(time.wt))])
  if(length(which(is.na(data.wt.ko)==TRUE))>0)
  {
    kk = which(!is.na(wt)==TRUE); wt = wt[kk]; time.wt = time.wt[kk];
    jj = which(!is.na(ko)==TRUE); ko = ko[jj]; time.ko = time.ko[jj];
  }
  
  if(length(wt)<=3|length(ko)<=3)
  {
    model.prob = c(NA, NA);
  }else{
    ### centralize the data and do not take into account of intercepts
    wt = wt - mean(wt);
    ko = ko - mean(ko);
    c.wt=cos(2*pi*time.wt/period)
    s.wt=sin(2*pi*time.wt/period)
    c.ko=cos(2*pi*time.ko/period)
    s.ko=sin(2*pi*time.ko/period)
    c.all=cos(2*pi*c(time.wt, time.ko)/period)
    s.all=sin(2*pi*c(time.wt, time.ko)/period)
    ## model 1 both flat
    rss1 = sum(wt^2 + ko^2);
    ## model 2 wt rhythmic and ko flat
    rss2 =  sum(lm(wt ~ 0+c.wt+s.wt)$residuals^2) + sum(ko^2)
    ## model 3 wt flat and ko rhytmic
    rss3 =  sum(wt^2) + sum(lm(ko ~ 0+c.ko+s.ko)$residuals^2)  
    ## model 4 wt and ko both rhytmic with the same parameters
    rss4 =  sum(lm(c(wt, ko) ~ 0+c.all+s.all)$residuals^2)
    ## model 5 wt and ko both rhytmic with different parameters
    rss5 =  sum(lm(wt ~ 0+c.wt+s.wt)$residuals^2) + sum(lm(ko ~ 0+c.ko+s.ko)$residuals^2) 
    
    ### BIC
    n = length(wt) + length(ko);
    BIC1 = n*log(rss1/n) + 0*log(n);
    BIC2 = n*log(rss2/n) + 2*log(n);
    BIC3 = n*log(rss3/n) + 2*log(n);
    BIC4 = n*log(rss4/n) + 2*log(n);
    BIC5 = n*log(rss5/n) + 4*log(n);
    
    BIC = c(BIC1, BIC2, BIC3, BIC4, BIC5)
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    prob.model = c(which(prob.model==max(prob.model)), prob.model[which(prob.model==max(prob.model))])
    #prob.wt.ko = prob.model
  }
  #print(pval.wt.ko)
  names(prob.model) = c('best.model', 'prob.best.model')
  
  return(prob.model)
}

library(plotrix)
library("circular")
make_circ_coord = function(t,x,ttot=24) {
  dt=(t[2]-t[1])*.45
  a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
  h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
  list(angles=a,heights=h)
}
circular_phase24H_histogram<-function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5, breaks=seq(0,24, by=1)){
  #color_DHS = rgb(0.6,0,0.2)
  par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
  #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
  br=breaks;
  h=hist(x, br=br,plot=FALSE)
  co=make_circ_coord(br[-1],h$counts)
  radial.plot(co$heights,co$angle,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
}

