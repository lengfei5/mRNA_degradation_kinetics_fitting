####
#### We work only on the absolute value of RNA-seq and cosine.beta parametrization
####
library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)
#library(plotrix)
library(numDeriv)

###### splicing demonstration and intron selection
selection = function(vect, cutoff=-4)
{
    length(which(vect>cutoff))>=6
}

Introns.signal.correction = function(vect, ss)
{
    #if((n%%(1000))==0){cat(round(100*n/nrow(Introns), digits = 1),'%\n')} ### indicating the progress
    test = vect
    test = unlist(strsplit(as.character(unlist(test)),","))
    nb.introns = length(test)/24
    rpkm.introns = matrix(0, nb.introns, 24)
    test = unlist(strsplit(as.character((test)),"/"))
    ii = c(1:nb.introns)
    jj = c(1:24)
    for(m in jj) rpkm.introns[,m] = as.numeric(test[c(ii*2+(m-1)*2*nb.introns)])
    
    test = test[ii*2-1]
    test = unlist(strsplit(as.character(unlist(test)),":"))
    test = test[ii*2]
    test = unlist(strsplit(as.character(unlist(test)),"-"))
    length.introns = as.numeric(test[ii*2]) - as.numeric(test[ii*2-1])
    kk = which(length.introns>200)
    if(length(kk)>0)
    {
        length.introns = length.introns[kk] ## select the introns with length>200bp
        rpkm.introns = rpkm.introns[kk,] ## select the introns with length>200bp
        
        norm.introns = rpkm.introns
        if(length(kk)>1)
        {
            for(m in 1:ncol(rpkm.introns))
            {
                norm.introns[,m] = rpkm.introns[,m]*10^9/ss[m]/length.introns
            }
            norm.introns = log2(norm.introns+10^-10)
            sel = apply(norm.introns, 1, selection) ## select introns with log2(rpkm)>-6
            length.sel = length.introns[sel]
            rpkm.sel = rpkm.introns[sel,]
            
        }else{
            norm.introns = log2(rpkm.introns*10^9/ss[m]/length.introns + 10^-10)
            #norm.introns = log2(norm.introns+10^-10)
            sel = selection(norm.introns)
            length.sel = length.introns[sel]
            rpkm.sel = rpkm.introns[sel]
        }
        
        if(sum(sel)>1) {
            nb.reads = apply(rpkm.sel, 2, sum)
            return(log2(nb.reads*10^9/sum(length.sel)/ss+10^-10));
        }else if(sum(sel)==1){
            return(log2(rpkm.sel*10^9/sum(length.sel)/ss+10^-10));
        }else{
            return(rep(NA, 24));
        }
        
    }else{
        return(rep(NA, 24));
    }
    #print(c(mean(keep[n,]), mean(unlist(compare[n,]))))
}

Introns.signal.correction.v2 = function(vect, ss)
{
    #if((n%%(1000))==0){cat(round(100*n/nrow(Introns), digits = 1),'%\n')} ### indicating the progress
    test = vect[c(1:24)]
    strand = as.numeric(vect[25])
    test = unlist(strsplit(as.character(unlist(test)),","))
    nb.introns = length(test)/24
    
    nb.reads = matrix(0, nb.introns, 24)
    test = unlist(strsplit(as.character((test)),"/"))
    ii = c(1:nb.introns)
    jj = c(1:24)
    for(m in jj) nb.reads[,m] = as.numeric(test[c(ii*2+(m-1)*2*nb.introns)])
    
    test = test[ii*2-1]
    test = unlist(strsplit(as.character(unlist(test)),":"))
    test = test[ii*2]
    test = unlist(strsplit(as.character(unlist(test)),"-"))
    length.introns = as.numeric(test[ii*2]) - as.numeric(test[ii*2-1])
    introns.start = as.numeric(test[ii*2-1])
    
    kk = which(length.introns>0)
    
    if(length(kk)>=3)
    {
        length.introns = length.introns[kk]
        introns.start = introns.start[kk]
        nb.reads = nb.reads[kk,]
        
        if(strand>0)
        {
            o1 = order(introns.start)
        }else{
            o1 = order(-introns.start)
        }
        nb.reads = nb.reads[o1,]
        length.introns = length.introns[o1]
    
    
        jj = c((length(kk)-2), (length(kk)-1), length(kk))
        length.introns = length.introns[jj] ## select the last three introns with length>200bp
        nb.reads= nb.reads[jj,]
        
        norm.introns = log2((apply(nb.reads, 2, sum))/sum(length.introns)/ss*10^9)
        
        if(length(which(norm.introns>(-4.8)))<18)
        {
            return(rep(NA, 24));
        }else{
            return(norm.introns)
        }
        
    }else{
        return(rep(NA, 24));
    }
    
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

Introns.splicing = function(vect, ss)
{
    #if((n%%(1000))==0){cat(round(100*n/nrow(Introns), digits = 1),'%\n')} ### indicating the progress
    test = vect[c(1:24)]
    strand = as.numeric(vect[25])
    test = unlist(strsplit(as.character(unlist(test)),","))
    nb.introns = length(test)/24
    rpkm.introns = matrix(0, nb.introns, 24)
    test = unlist(strsplit(as.character((test)),"/"))
    ii = c(1:nb.introns)
    jj = c(1:24)
    for(m in jj) rpkm.introns[,m] = as.numeric(test[c(ii*2+(m-1)*2*nb.introns)])
    
    test = test[ii*2-1]
    test = unlist(strsplit(as.character(unlist(test)),":"))
    test = test[ii*2]
    test = unlist(strsplit(as.character(unlist(test)),"-"))
    length.introns = as.numeric(test[ii*2]) - as.numeric(test[ii*2-1])
    introns.start = as.numeric(test[ii*2-1])
    if(nb.introns>1)
    {
        if(strand>0)
        {
            o1 = order(introns.start)
        }
        else{
            o1 = order(-introns.start)
        }
        rpkm.introns = rpkm.introns[o1,]
        length.introns = length.introns[o1]
    }
    
    kk = which(length.introns>200)
    
    if(length(kk)>0)
    {
        length.introns = length.introns[kk] ## select the introns with length>200bp
        rpkm.introns = rpkm.introns[kk,] ## select the introns with length>200bp
        
        norm.introns = rpkm.introns
        if(length(kk)>1)
        {
            for(m in 1:ncol(rpkm.introns))
            {
                norm.introns[,m] = rpkm.introns[,m]*10^9/ss[m]/length.introns
            }
            norm.introns = log2(norm.introns+10^-10)
            return(c(length(kk), packing(length.introns), t(apply(norm.introns, 2, packing))))
            
        }else{
            norm.introns = log2(rpkm.introns*10^9/ss[m]/length.introns + 10^-10)
            #norm.introns = log2(norm.introns+10^-10)
            return(c(length(kk), length.introns, norm.introns))
        }
        
    }else{
        return(c(0, rep(NA, 25)));
    }
    
}
######## prepare tables of WT-NF and KO-NF;
prepare.NF.KO.tables = function()
{
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    
    R = read.table('/Users/jiwang/RNA_seq_Data/Total_RNA_Seq_AL_RF_KO_Tables/WT_RF_Intron_Exon_RFP.txt',sep='\t', header=TRUE)
    R = R[, -1]
    kk = grep('RFP', colnames(R))
    R = R[, -kk]
    
    ii = grep('WT_RF_Exon', colnames(R))
    jj = grep('WT_RF_Intron', colnames(R))
    ## change mRNA and premRNA data into linear scale and statistics in linear scale
    data1 = 2^R[,ii]
    data2 = 2^R[,jj]
    
    tt = c(0:23)*4+2
    stat1 = t(apply(data1,1, f24_R2_alt2, t=tt))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=tt))
    
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    R = data.frame(R[,1], data1, data2, stat1, stat2, stringsAsFactors=FALSE)
    colnames(R)[1] = 'gene'
    colnames(R)[2:25] = paste('ZT', tt, '.abs.mRNA', sep='')
    colnames(R)[26:49] = paste('ZT', tt, '.abs.premRNA', sep='')
    #stat = colnames(R)[50:56]
    #stat = unlist(strsplit(as.character(stat), '.'))
    colnames(R)[50:56] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.mRNA', sep='')
    colnames(R)[57:63] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.premRNA', sep='')
    #kk = which(R$mean.mRNA>R$mean.premRNA)
    #R = R[kk,]
    R.WT.RF = R
    
    R = read.table('/Users/jiwang/RNA_seq_Data/Total_RNA_Seq_AL_RF_KO_Tables/KO_RF_Intron_Exon_RFP.txt',sep='\t', header=TRUE)
    R = R[, -1]
    kk = grep('RFP', colnames(R))
    R = R[, -kk]
    
    ii = grep('KO_RF_Exon', colnames(R))
    jj = grep('KO_RF_Intron', colnames(R))
    ## change mRNA and premRNA data into linear scale and statistics in linear scale
    data1 = 2^R[,ii]
    data2 = 2^R[,jj]
    
    tt = c(0:11)*4+2
    stat1 = t(apply(data1,1, f24_R2_alt2, t=tt))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=tt))
    
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    R = data.frame(R[,1], data1, data2, stat1, stat2, stringsAsFactors=FALSE)
    colnames(R)[1] = 'gene'
    colnames(R)[2:13] = paste('ZT', tt, '.abs.mRNA', sep='')
    colnames(R)[14:25] = paste('ZT', tt, '.abs.premRNA', sep='')
    #stat = colnames(R)[50:56]
    #stat = unlist(strsplit(as.character(stat), '.'))
    colnames(R)[26:32] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.mRNA', sep='')
    colnames(R)[33:39] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.premRNA', sep='')
    #kk = which(R$mean.mRNA>R$mean.premRNA)
    #R = R[kk,]
    R.KO.RF = R
    
    save(R.WT.RF, R.KO.RF, file=paste('myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata', sep=''))
}

######## FUNCTIONS FOR SIMULATED DATA and for the optimization and model selection
set.global.sigma = function() ## this function is to used to generate fake data with different standard deviation
{
	sigma.s <<- 0.2065864
	sigma.m <<- 0.1013937
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

#############################################################################################################################
######################################################  Generate simulated data
#############################################################################################################################
generate.fake.data = function(T = T[1:100,], X = 100, model = 1, zt = seq(0,46,by = 2), parametrization = 'cosine.beta', absolute.signal=TRUE)
{
	#model = 4; T = F[((model-1)*X+1):(model*X),]; absolute.signal=TRUE; parametrization = 'cosine.beta'; zt = seq(0, 94, by=2);
	set.seed(42)
	ZT.int = grep('.abs.premRNA', colnames(T)); 
	ZT.ex = grep('.abs.mRNA', colnames(T));
    
	if(absolute.signal)
	{
		###### This is part to simulate RNA-seq absolute levels 
		T$gene = paste('fake_m',model,'_',1:X,sep = '')
		
		### Define the limits of parameters: beta, eps.gamma, gamma
		bounds = set.bounds(model = 4, parametrization = parametrization, absolute.signal = TRUE); 
		lower = bounds$lower; 
		upper = bounds$upper;
		
        ### parameter of premRNAs
		Min.int.min = lower[5];
		Min.int.max = upper[5];
		
		Amp.int.min = 0.1;
		Amp.int.max = upper[6];
		
		phase.int.min = lower[7];
		phase.int.max = upper[7];
		
		beta.int.min = lower[8];
		beta.int.max = upper[8];
        
        ### parameters of mrna degradation
        gamma.min = lower[1]
        gamma.max = log(2)/(3/60);
        
        amp.gamma.min = 0.1;
        amp.gamma.max = upper[2];
        
        phase.gamma.min = lower[3];
        phase.gamma.max = upper[3];
        
        splicing.k.min = lower[4];
        splicing.k.max = upper[4];
		
		### Shared parameters for all models
        T$Min.int = pmin(pmax(2^rnorm(X, mean = -2, sd = 1.5), Min.int.min), Min.int.max)
		T$splicing.k = pmin(pmax(log(2)/rnorm(X, mean = 5/60, sd = 2/60), splicing.k.min), splicing.k.max)
		T$gamma = pmin(pmax(log(2)/exp(rnorm(X, mean = log(5), sd = 1.0)), gamma.min), gamma.max)
		
		### Parameters specific for certain models
		if((model == 2)|(model == 4)) ## randomly generate parameters for rhythmic pre-mrna
		{
            rel.amp.int= pmax(pmin(rnorm(n = X, mean = 0.4, sd = 0.2), 1.0), 0.05)
            fc.int = (1+rel.amp.int)/(1-rel.amp.int);
            T$Amp.int = T$Min.int * fc.int ## assign premRNA amplitudes according to relative amplitudes
			T$phase.int =  sample(seq(phase.int.min ,phase.int.max,by = 0.1),X, replace = TRUE)
			T$beta.int = pmax(pmin(rnorm(n = X, mean = 1, sd = 1), beta.int.max), beta.int.min)
		
		}else{
			T$Amp.int = rep(0,X); 
			T$phase.int = 0; 
			T$beta.int = rep(1,X); 
		}
		if(model>2) ## parameters for rhythmic degradation
		{
            rel.amp.gamma= pmax(pmin(rnorm(n = X, mean = 0.4, sd = 0.2), 1.0), 0.05)
            fc.gamma = (1+rel.amp.int)/(1-rel.amp.int);
            T$amp.gamma = T$gamma*fc.gamma;
            
            if(model==3)
            {
                T$phase.gamma = sample(seq(phase.gamma.min, phase.gamma.max,by = 0.1),X, replace = TRUE)
            }
            if(model==4)
            {
                phase.diff = sample(seq(10, 20, by=0.1), X, replace=TRUE)
                T$phase.gamma = (T$phase.int + phase.diff)%%24;
            }
		}else{
			T$amp.gamma = rep(0,X)
			T$phase.gamma = rep(0,X)
		}
		
        abs.int = T[, ZT.int]*0
		abs.ex = T[,ZT.ex]*0
		for(i in 1:X)
		{
			## absolute values of premRNAs
			abs.int[i,] = compute.s.beta(t = zt, Min = T$Min.int[i], Amp = T$Amp.int[i], phase = T$phase.int[i], beta = T$beta.int[i]);
						
			## absolute values of mRNAs
			if(model == 1)
			{
				abs.ex[i,] = T$splicing.k[i] * T$Min.int[i]/T$gamma[i];
			} 
			if(model>1)
			{
				
				abs.ex[i,] = compute.m.beta(t = zt, 
											gamma = T$gamma[i], amp.gamma = T$amp.gamma[i], phase.gamma = T$phase.gamma[i], splicing.k = T$splicing.k[i],
											Min = T$Min.int[i], Amp = T$Amp.int[i], phase = T$phase.int[i], beta = T$beta.int[i])
			}
		}
        
		log.abs.int = log(abs.int);
		log.abs.ex = log(abs.ex);
		
        ### set the noise
        SS = as.matrix(log(T[, ZT.int]))
        MM = as.matrix(log(T[, ZT.ex]))
        sigma.S = apply(SS, 1, variance.estimate.replicates, nb.replicates=4)
        sigma.M = apply(MM, 1, variance.estimate.replicates, nb.replicates=4)
        T$sigma.s = sigma.S
        T$sigma.m = sigma.M
		noise.int = matrix(rnorm(X*length(ZT.int), mean = 0, sd = sigma.S),nrow = X, ncol = length(ZT.int))
		noise.ex = matrix(rnorm(X*length(ZT.ex), mean = 0, sd = sigma.M),nrow = X, ncol = length(ZT.ex))
		
		log.noise.abs.int = log.abs.int + noise.int
		log.noise.abs.ex = log.abs.ex + noise.ex
				
		noise.abs.int = exp(log.noise.abs.int)
		noise.abs.ex = exp(log.noise.abs.ex)
		
		T[,ZT.int] = noise.abs.int;
		T[,ZT.ex] =  noise.abs.ex
        #T$mean.ex = apply(noise.abs.ex,1,mean)
        
		return(T)
	}
}

##################################################################################################################
##################################################################################################################
################################################################################################################## Functions for the optimization
##################################################################################################################
##################################################################################################################

################
### Main Function for model selection and fitting
################
make.fits.with.all.models.for.one.gene = function(T = T, gene.index = 1, debug = FALSE, parametrization = c('cosine.beta'), zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int, absolute.signal = TRUE)
{
	#T = T; gene.index = j; debug = TRUE; parametrization = 'cosine.beta';  zt = zt; i.ex = ZT.ex; i.int = ZT.int; absolute.signal = TRUE
    
    #####  Prediction of variance (sd here) by LOESS
    Variance.estimate.loess = TRUE
    if(Variance.estimate.loess)
    {
        S = unlist(T[gene.index, i.int])
        M = unlist(T[gene.index, i.ex])
        
        load(file=paste('my_genes_RNA_seq_analysis_variance_mean_logRPKM_loess_v1.Rdata', sep=''))
        
        #if(nb.replicates==2) xx = rbind(x[c(1:12)], x[c(13:24)])
        x.s = as.numeric(log(S))
        xx.s = cbind(x.s[c(1:12)], x.s[c(13:24)], x.s[25:36], x.s[37:48])
        sigma.ss = rep((predict(var.mean.intron.loe, data.frame(mean = as.numeric(apply(xx.s, 1, mean)))))^2, 4)
        
        x.m = as.numeric(log(M))
        xx.m = cbind(x.m[c(1:12)], x.m[c(13:24)], x.m[25:36], x.m[37:48])
        sigma.mm = rep((predict(var.mean.exon.loe, data.frame(mean = as.numeric(apply(xx.m, 1, mean)))))^2, 4)
        
        sigma.ss = as.numeric(sigma.ss)
        sigma.mm = as.numeric(sigma.mm)
        #sigma.ss = variance.estimate.loess(log(S), nb.replicates=4)
        #sigma.mm = variance.estimate.loess(log(M), nb.replicates=4)
    }

    
    if(debug){cat('\t\t sigma of premRNA \n ', sigma.ss[c(1:12)], '\n')}
    if(debug){cat('\t\t sigma of mRNA \n ', sigma.mm[c(1:12)], '\n')}
    
	for(model in 1:4)
	{
		if(debug){cat('\t starting model ',model,'\n');}
		param.fit = make.fit.spec.model(T = T, gene.index = gene.index, model = model, debug = debug, zt = zt, i.ex = i.ex, i.int = i.int, sigma.ss=sigma.ss, sigma.mm=sigma.mm);
		
		if(model == 1)
		{
			Param.fit.for.gene = param.fit
		}else{
			Param.fit.for.gene = c(Param.fit.for.gene, param.fit)
		}
		if(debug){cat('\t model ',model,' finished \n')};
	}
	## model selction for one gene
	#BIC = my.BIC.this.gene(Param.fit.for.gene)
	#Param.fit.for.gene = c(Param.fit.for.gene, BIC)
	return(Param.fit.for.gene)
}

## BIC model-selection for one gene
my.BIC.this.gene.loglike = function(param.fits.results, nb.data=96, absolute.signal=TRUE)
{
    set.nb.param(absolute.signal=TRUE);
    index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), names(param.fits.results))
    error.m1 = param.fits.results[index[1]]
    error.m2 = param.fits.results[index[2]]
    error.m3 = param.fits.results[index[3]]
    error.m4 = param.fits.results[index[4]]
    
    ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
    BIC.m1 = log(nb.data)*n.param[1] + (error.m1)
    BIC.m2 = log(nb.data)*n.param[2] + (error.m2)
    BIC.m3 = log(nb.data)*n.param[3] + (error.m3)
    BIC.m4 = log(nb.data)*n.param[4] + (error.m4)
    
    BIC = c(BIC.m1, BIC.m2, BIC.m3, BIC.m4)
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    
    BIC = c(BIC, which.min(BIC), prob.model[which.min(BIC)])
    names(BIC) = c('BIC.m1', 'BIC.m2', 'BIC.m3', 'BIC.m4', 'BIC.best.model', 'prob.best.model')
    
    return(BIC)								
}	


make.fit.spec.model = function(T = T, gene.index = 1, model = 1, debug = FALSE, parametrization = c('cosine.beta'),zt = seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int,
sigma.ss=rep(0.2, 48), sigma.mm=rep(0.1, 48), absolute.signal = TRUE)
{
	param.fit = NA
	S = unlist(T[gene.index, i.int])
	M = unlist(T[gene.index, i.ex])
    
	## Model 1 premrna and mran are both flat
	if(model == 1)
	{
		s = rep(exp(mean(log(S))),length(S))
		m = rep(exp(mean(log(M))), length(M))
        
		err = error(S,s,M,m, sigma.ss=sigma.ss, sigma.mm=sigma.mm)
		param.fit = err;
		names(param.fit) = 'error.m1';
	}
	## parameter estimations for model 2,3,4
	if(model > 1)
	{
		param.fit = make.optimization(T = T, i = gene.index, model = model, Nfit = NA, debug = debug, zt = zt, i.ex = i.ex, i.int = i.int, sigma.ss=sigma.ss, sigma.mm=sigma.mm)
	}
	return(param.fit)
}

variance.estimate.replicates = function(x, nb.replicates=4, log=TRUE)
{
    if(!log)
    {
        x = log(x)
    }
    x = as.numeric(x)
    #xx = matrix(x, ncol=4)
    if(nb.replicates==4) xx = rbind(x[c(1:12)], x[c(13:24)], x[25:36], x[37:48])
    if(nb.replicates==2) xx = rbind(x[c(1:12)], x[c(13:24)])
    
    sd = mean(apply(xx, 2, sd))
    
    return(sd)
}

variance.mean.counts = function(x)
{
    # x = exon
    #var.mean = matrix(NA, ncol=3, nrow=(nrow(x)*12))
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

variance.mean.log = function(x)
{
    # x = log(exon[1, ])
    #var.mean = matrix(NA, ncol=3, nrow=(nrow(x)*12))
    var.mean = c()
    #mm = c()
    #vv = c()
    for(n in 1:12)
    {
        index = c(n, (n+12), (n+24), (n+36))
        #mm = c(mm, mean(x[index]))
        mean = apply(x[, index], 1, mean)
        var = apply(x[, index], 1, var)
        #vv = c(vv, var(x[index]))
        var.mean = rbind(var.mean, cbind(rep(n, nrow(x)), mean, var))
    }
    #sd = mean(apply(xx, 2, sd))
    #mm.vv = c(mm, vv)
    #return(mm.vv)
    return(var.mean)
}

variance.mean.log.replicates = function(x)
{
    mm = c()
    vv = c()
    for(n in 1:12)
    {
        index = c(n, (n+12), (n+24), (n+36))
        mm = c(mm, mean(x[index]))
        #mean = apply(x[, index], 1, mean)
        #var = apply(x[, index], 1, var)
        vv = c(vv, var(x[index]))
        #var.mean = rbind(var.mean, cbind(rep(n, nrow(x)), mean, var))
    }
    #sd = mean(apply(xx, 2, sd))
    mm.vv = c(mm, vv)
    return(mm.vv)

}

####################################################
#### Main Optimization Function to estimat parameters for model 2, 3, 4
####################################################
make.optimization = function(T = T, i = 1, model = 2, Nfit = NA, debug = FALSE, parametrization =c('cosine.beta'), zt =  seq(0,46,by = 2), i.ex = ZT.ex, i.int = ZT.int,
sigma.ss=rep(0.2, 48), sigma.mm=rep(0.1, 48), absolute.signal = TRUE)
{
	#i = j; zt =  seq(0,94,by = 2); i.ex = ZT.ex; i.int = ZT.int;absolute.signal = TRUE; Nfit=NA; debug = TRUE;
	param.fit = NA
	S = unlist(T[i, i.int])
	M = unlist(T[i, i.ex])
	
	fitting.factor = 1; ### define the nb of initial values in optimization
	
	#######
	####### FITTING the absolute signal
	#######
	if(absolute.signal) 
	{
		bounds = set.bounds(model = model); 
		upper = bounds$upper; 
		lower = bounds$lower;
        
        ### Here we want first to fit the pre-mRNA profile to identify parameters mean.int, fold.change.int, phase.int, beta.int
		## Prefit premRNA for model2 and model4 in order to have the good initial values for premRNA parameters
		prefit = TRUE 
		if((model==2|model==4) & prefit) 
		{
			Nfit.S = 4*fitting.factor;
			index = i
			
			### Choose different initial values of parameters for S fitting
			Min.init = rep(min(S), Nfit.S)
			Amp.init = rep((max(S)-min(S)), Nfit.S)
 			phase.init = zt[which.max(S)]
			phase.init = (rep(phase.init, Nfit.S)+rnorm(Nfit.S,sd = 3))%%24
			beta.min = 1;
			beta.max = 5;
			beta.init = lseq(beta.min, beta.max, length = Nfit.S)
			
			PAR.INIT.S = cbind(Min.init,Amp.init, phase.init, beta.init)
	
			errors.fit.s = rep(NA, Nfit.S)
			
			limit.factor = 5
			for(fit.nb.s in 1:Nfit.S)
			{
				par.init.s = PAR.INIT.S[fit.nb.s,]
				opt.s = optim(par.init.s, f2min.int, S=S, sigma.ss=sigma.ss, zt = zt, method = 'L-BFGS-B', lower = c(min(S)/limit.factor, (max(S)-min(S))/limit.factor, 0, 1), upper = c(max(S), (max(S)-min(S))*limit.factor, 24, 5))
				#print(opt.s$convergence)
				res.fit.s = opt.s$par
				errors.fit.s[fit.nb.s] = opt.s$value
				eval(parse(text = paste('res.fit.s.', fit.nb.s, ' = res.fit.s', sep = '')))
			}
			## choose the best-fitting parameters for S
			imin.s = which.min(errors.fit.s); 
			eval(parse(text = paste('res.fit.s = res.fit.s.', imin.s, sep = '')))			
		}
        
		###
		### Now to fit mRNA and premRNA all together for model2, model3 and model4
		###
		if(is.na(Nfit))
		{
			if(model==2) Nfit = fitting.factor*12;
			if(model==3) Nfit = fitting.factor*10;
			if(model==4) Nfit = fitting.factor*12;
		}
		
		###
		### Define the initial values for degradation and splicing parameters
		max.half.life = 24;
		min.half.life = 10/60;
		gamma.init = c(rep(log(2)/lseq(max.half.life, min.half.life, length = Nfit), Nfit%/%Nfit), rep(log(2)/5,Nfit%%Nfit))
		splicing.k.init = c(rep(log(2)/lseq(0.5/60, 0.5, length = Nfit), Nfit%/%Nfit), rep(log(2)/(5/60),Nfit%%Nfit))
		eps.gamma.init = sample(seq(0.01, 1, length = 20), Nfit, replace = TRUE)
		phase.gamma.init = sample(rep(c(6,12, 18),100), size=Nfit, replace=TRUE)
		
		if(model==2)
		{
			Min.init = rep(res.fit.s[1], Nfit)
			Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
			phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
			beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
			
			PAR.INIT = cbind(gamma.init, splicing.k.init, Min.init, Amp.init, phase.init, beta.init); 
			colnames(PAR.INIT) = c('gamma','splicing.k','Min.int','Amp.int','phase.int','beta.int');
            upper[2] = log(2)/(1/60/60) ## maximum splicing rate of 1s
			lower[3] = min(S)/5;
			upper[3] = max(S);
			lower[4] = (max(S)-min(S))/5;
			upper[4] = (max(S)-min(S))*5;
		}
		if(model==3)
		{	
			Min.init = rep(min(S), Nfit)
			PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, splicing.k.init, Min.init)
			colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int');
			upper[4] = log(2)/(1/60/60) ## maximum splicing rate of 1s
			lower[5] = min(S)/5;
			upper[5] = max(S);
		}
		if(model==4)
		{
			Min.init = rep(res.fit.s[1], Nfit)
			Amp.init = seq(res.fit.s[2]/2, res.fit.s[2]*2, length=Nfit)
			phase.init = (rep(res.fit.s[3],Nfit)+rnorm(Nfit,sd = 2))%%24
			beta.init = lseq(max(1, (res.fit.s[4]-2)), min((res.fit.s[4]+2), 10), length=Nfit)
			
			PAR.INIT = cbind(gamma.init, eps.gamma.init, phase.gamma.init, splicing.k.init, Min.init, Amp.init, phase.init, beta.init);
			colnames(PAR.INIT) = c('gamma','eps.gamma','phase.gamma','splicing.k', 'Min.int','Amp.int','phase.int','beta.int')
			upper[4] = log(2)/(1/60/60) ## maximum splicing rate of 1s
            lower[5] = min(S)/5;
			upper[5] = max(S);
			lower[6] = (max(S)-min(S))/5;
			upper[6] = (max(S)-min(S))*5;
		}
        
		#################################
		#### Fit the S and M all together 
		#################################
		errors.fit = rep(NA, Nfit)
        
        #fit.number = 1;
        #while(fit.number<(Nfit+1))
        for(fit.number in 1:Nfit)
		{
			##cat(fit.number, '\n');
			par.init = PAR.INIT[fit.number,]
            #if(debug){cat('\t\t starting optimization # ',fit.number,' for model ', model,':',par.init,'\n')}
			
            maxit = 100;
            limit = 1000;
            converg = 1;
			
            ptm <- proc.time()
			while(converg!=0 && maxit<limit)
			{
				#compute.m.beta(zt, par.init[1], par.init[2], par.init[3], par.init[4], par.init[5], par.init[6], par.init[7], par.init[8], simulation.only=FALSE)
				opt1 = optim(par.init, f2min, M = M, S = S, sigma.ss=sigma.ss, sigma.mm=sigma.mm, model = model, debug = debug, zt = zt, method = 'L-BFGS-B', lower = lower, upper = upper, control = list(maxit=maxit,trace=0))
				converg = opt1$convergence;
				maxit = maxit + 200;
                #cat(opt$par, '\n')
			}
            proc.time() - ptm
            
            opt = opt1;
            res.fit = opt$par
            errors.fit[fit.number] = opt$value
			
            eval(parse(text = paste('res.fit.', fit.number, ' = res.fit', sep = '')))
            #eval(parse(text = paste('res.fit.stderr.', fit.number, ' = res.fit.stderr', sep = '')))
            if(debug){cat('\t\t optimization # ',fit.number,' finished : \t\t\t',res.fit,'\t',opt$value,'\n')}
            #fit.number = fit.number+1;
            
        }
        
		## choose the best estimated parameters
		imin = which.min(errors.fit); 
		eval(parse(text = paste('res.fit = res.fit.', imin, sep = '')))
        #eval(parse(text = paste('res.fit.stderr = res.fit.stderr.', imin, sep = '')))
        
        #### calculate standard errors for parameters
        ttry = try(sqrt(diag(solve(0.5*hessian(func=f2min, M = M, S = S, sigma.ss=sigma.ss, sigma.mm=sigma.mm, model = model, debug = debug, zt = zt, x=res.fit)))), silent = TRUE)
        if(!inherits(ttry, "try-error")){
            #opt$stderr = sqrt(diag(2*opt$value/(48-n.param[model])*solve(opt$hessian)))
            res.fit = c(res.fit, sqrt(diag(solve(0.5*hessian(func=f2min, M = M, S = S, sigma.ss=sigma.ss, sigma.mm=sigma.mm, model = model, debug = debug, zt = zt, x=res.fit)))))
        }else{
            #opt$stderr = rep(NA, length(opt$par))
            res.fit = c(res.fit, rep(NA, length(res.fit)))
        }
        
        param.fit = c(errors.fit[imin], res.fit);
        names(param.fit) = paste(c('error', colnames(PAR.INIT), paste(colnames(PAR.INIT), '.stderr', sep='')),'.m',model,sep = '')
        #names(param.fit) = paste(c('error', colnames(PAR.INIT)),'.m', model,sep = '')
	}
	
	return(param.fit);
	if(debug){cat('\t\t\t optimization finished\n')}
}

#################################
## Set the boundary of parameters for each model
#################################
set.bounds = function(model = 4, parametrization =c('cosine.beta'), absolute.signal = TRUE)
{
	### minimum and maximum of gamma
	gamma.min = log(2)/24
	gamma.max = log(2)/(10/60);

	### relative amplitudes of gamma
	eps.gamma.min = 0;
    eps.gamma.max = 1;
    
	phase.gamma.min = 0;
	phase.gamma.max = 24; 
	
	splicing.k.max = log(2)/(1/60/60) ## 1s of splicing time;
	splicing.k.min = log(2)/(0.5) ## 30 min of splicing time;

	Min.int.max = 2^13; ## intron signal mean 2^9;
	Min.int.min = 2^(-22);

	Amp.int.min = 0;
	Amp.int.max = 200;
	
	phase.int.min = 0;
	phase.int.max = 24; 
	
	beta.int.min = 1;
    beta.int.max = 5;
	
	param.synthesis.upper = c(Min.int.max, Amp.int.max, phase.int.max, beta.int.max)
	param.synthesis.lower = c(Min.int.min, Amp.int.min, phase.int.min, beta.int.min)
		
	if(absolute.signal)
	{
		if(model==1)
		{
			upper = c(gamma.max, splicing.k.max, mean.int.max)
			lower = c(gamma.min, splicing.k.min, mean.int.min)
		}
		if(model==2)
		{
			upper = c(gamma.max, splicing.k.max, param.synthesis.upper)
			lower = c(gamma.min, splicing.k.min, param.synthesis.lower)
		}
		if(model==3)
		{
			upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max, Min.int.max)
			lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min, Min.int.min)
		}
		if(model==4)
		{
			upper = c(gamma.max, eps.gamma.max, phase.gamma.max, splicing.k.max,  param.synthesis.upper)
			lower = c(gamma.min, eps.gamma.min, phase.gamma.min, splicing.k.min,  param.synthesis.lower)
		}
	}
	return(list(lower = lower, upper = upper))
}

###
### ERROR function for premRNA profile
###
f2min.int = function(par.int, S, sigma.ss=rep(0.2, length(S)), parametrization ='cosine.beta', debug = FALSE, zt = seq(0,46,by = 2), absolute.signal = TRUE)
{
	S = as.numeric(unlist(S)); 
		
	if(absolute.signal)
	{
		s = compute.s.beta(t = zt, Min = par.int[1], Amp = par.int[2], phase = par.int[3], beta= par.int[4])
	}	
		
	err = sum((log(S)-log(s))^2/sigma.ss^2);
	return(err)
}

####
#### ERROR function for premRNA and mRNA
####
f2min = function(par, M, S, sigma.ss=rep(0.2, length(S)), sigma.mm=rep(0.1, length(M)), model, parametrization =c('cosine.beta'), debug = FALSE, zt = seq(0,46,by = 2), absolute.signal = TRUE)
{
	#par = par.init;
    #simulation.only=FALSE
	# initialization of the parameters
	if(absolute.signal)
	{
		gamma = par[1];
		
		if(model==2)	
		{
			eps.gamma = 0.0;
			phase.gamma = 0.0; 
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
		
		#print(m);
		err.fit = error(S=S, s=s, M=M, m=m, sigma.ss=sigma.ss, sigma.mm=sigma.mm)
		#print(err.fit)
	}
	
	#print('........')
	return(err.fit)
}

error = function(S,s,M,m, sigma.ss=rep(0.2, length(S)), sigma.mm=rep(0.1, length(M)), log = TRUE, intense.debug=FALSE)
{
    S = as.numeric(unlist(S));
    s = as.numeric(unlist(s));
    M =  as.numeric(unlist(M));
    m =  as.numeric(unlist(m));
    
    if(log) ## calculate separately the errors of premRNAs and mRNAs
    {
        if(any(S<=0)|any(s<=0)){error.S = 10^10;
        }else{
            error.S = sum((log(S)-log(s))^2/sigma.ss^2);
        }
        if(any(M<=0)|any(m<=0)){
            error.M = 10^10
        }else{
            error.M = sum((log(M)-log(m))^2/sigma.mm^2);
        }
    }
    
    ### Here we estimate gene-specific variances with replicates
    error = error.M + error.S ## -2loglike
    
    if(intense.debug){cat('______________________________-------- error.S = ',error.S,'\n')}
    if(intense.debug){cat('______________________________-------- error.M = ',error.M,'\n')}
    return(error)  
}


############
############ Absolute premRNA profile
############
compute.s.beta = function(t = zt, Min = 1, Amp = 2, phase = 12, beta = 1)
{
	w <<- 2*pi/24;
	s = Min+Amp*((1+cos(w*(t-phase)))/2)^beta
	return(s)
}

compute.s.beta.v1 = function(t = zt, mean = 1, fold.change = 2, phase = 12, beta = 1)
{
	w <<- 2*pi/24;
	#s = 1/(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)*(1+(fold.change-1)*((1+cos(w*(t-phase)))/2)^beta)
	s = mean/(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)*(1+(fold.change-1)*((1+cos(w*(t-phase)))/2)^beta)
	return(s)
}

############
############ Absolute mRNA profile with two methods: Integration and Simulation
############
#### Integration 
integrate.m = function(t=seq(0, 46, by=2), gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, splicing.k=log(2)/(5/60), Min = 0.5, Amp=2.0, phase=12, beta=1)
{
	w <<- 2*pi/24;
	m = rep(0, length(t))
	
	Tstable = 24*(ceiling(log(2000)/gamma/24) +1)
	for(i in 1:length(t))
	{
		time = t[i]; 
		m[i] = splicing.k*exp(-Gamma(t = Tstable+time, gamma = gamma, eps.gamma=eps.gamma, phase.gamma=phase.gamma)) * integrate(f2integrate, lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, Min, Amp, phase, beta))$value
	}
	return(m)
}
f2integrate = function(t, par)
{
	gamma = par[1]; 
	eps.gamma= par[2];
	phase.gamma= par[3]; 
	Min = par[4];
	Amp = par[5];
	phase = par[6];
	beta = par[7];
	return(exp(Gamma(t, gamma, eps.gamma, phase.gamma)) * compute.s.beta(t, Min, Amp, phase, beta))
}
#### This is the integral of degradation function
Gamma = function(t = 0, gamma = log(2)/3, eps.gamma = 0.2, phase.gamma = 0)
{
	w <<- 2*pi/24;
	Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
    #Gamma = gamma*t + amp.gamma/2*(t+ sin(w*(t-phase.gamma))/w)
	return(Gamma)
}	


#### Simulation
dmdt = function(t, y, par)
{
	m = y	
	gamma = par[1];
	eps.gamma = par[2];
	phase.gamma = par[3];
	splicing.k = par[4];
	param.synthesis.1 = par[5]; 
	param.synthesis.2 = par[6]; 
	param.synthesis.3 = par[7]; 
	param.synthesis.4 = par[8];
	
	s.t = compute.s.beta(t=t, Min = param.synthesis.1, Amp = param.synthesis.2, phase = param.synthesis.3, beta =  param.synthesis.4)
		
	gamma.t = gamma *(1 + eps.gamma*(cos(2*pi/24*(t-phase.gamma))))
	
	dmdt = splicing.k * s.t - gamma.t * m
	
	list(dmdt,NULL)
}

simulate.m = function(t, par)
{
	gamma = par[1];
	
	Tstable =  24*(ceiling(log(2000)/gamma/24) + 5) ## burning time
	t.res = 2; 
	if(length(t)!=1){t.res = (t[2]-t[1])}
	t.sup = seq(0, Tstable+max(t),by= t.res)
	
	soln = lsoda(dmdt,
				 times= t.sup, ## times
				 y = 0, #init.conditions
				 par=par, rtol = 1e-6, atol = 1e-6) ## parameter values
	
	#soln[match(t+48, soln[,1]),2]
	i.last = nrow(soln); 
	i.keep = seq(i.last - length(t)+1,i.last,by = 1)
	#plot(soln[,1],soln[,2], type = 'b',col='red')
	m = soln[i.keep,2]; 
	#mean = mean(m);m = m/mean
	#cat(soln[,2],'\n')
	#cat(parametrization,'\n')
	return(m)
}

########### Main function to calculate mRNA profile
compute.m.beta = function(t=seq(0, 94, by=2), gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, splicing.k=log(2)/(5/60), Min = 0.5, Amp=5, phase=12, beta=1, simulation.only=FALSE)
{
	#gamma=log(2)/5; eps.gamma=0.25; phase.gamma=12; splicing.k=log(2)/(5/60); mean = 10; fold.change=10; phase=6; beta=1;t = seq(0,46, by=2);simulation.only=FALSE
    #gamma=4.159; eps.gamma=0.3; phase.gamma=24; splicing.k=43.135; Min = 0.01418; Amp=0.215; phase=6.3; beta=5;t = seq(0,94, by=2);simulation.only=FALSE
	w <<- 2*pi/24;
	zt = t;
	
	#### because m is periodic, so just compute the first period and then repeat it.
	if(((max(t)-min(t))>24) & (max(t)%%24 == 24-t[2]+t[1]))
	{
        nb.period = (max(t)-min(t)+t[2]-t[1])/24
		zt = t[1:(length(t)/nb.period)]
	}
	
	if(!simulation.only)
	{
		m.try = try(integrate.m(t = zt, gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta), silent = TRUE)
		if(!inherits(m.try, "try-error")){
			m = m.try
		}else{
			
            #cat('error in integration and start simulation !\n')
            #cat(zt, '\n')
            #cat(c(gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta), '\n')
			#print(c(gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta));
			m = simulate.m(t = zt, par = c(gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta));
		}
	}else{		
		#cat('HERE')
		m = simulate.m(t = zt, par = c(gamma, eps.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta));
	}
	
	if(length(zt)!=length(t))
    {
        nb.period = (max(t)-min(t)+t[2]-t[1])/24
        m = rep(m, nb.period)
    }
	
	return(m)
}

##################################################################################################################
##################################################################################################################
################################################################################################################## finishing line of optimization
##################################################################################################################
##################################################################################################################

######################
# Functions for model selecion
######################
set.nb.param = function(absolute.signal=TRUE)
{
	n.param <<- c(0,6,5,8);
}

fdr2pval = function(P,fdr = 0.1){
	N = length(P)
	sorted_P = sort(P)
	in_fdr = which((sorted_P - c(1:N)*fdr/N)<=0)
	x = length(in_fdr)
	if(x>0){
		in_fdr = in_fdr[length(in_fdr)]; 
		p_th = sorted_P[in_fdr]
	}else{p_th = 1}
	return(p_th)
}
passes.fdr = function(t = t, fdr = 0.1){
	pval.col = grep('.pval', colnames(t)); X = length(pval.col);
	pval.lim = apply(t[,pval.col],2,fdr2pval,0.1)
	pass = t( t(t[,pval.col])<= pval.lim )
	return(apply(pass,1,sum)>(X/2))
}



##### LAURA OLD codes for the model selection parts
Likelihood.Ratio.Test.LAURA = function(T = T, pval = 0.05, pval.stringent = 10^(-15), V1 = FALSE, combined = FALSE, FDR = FALSE, fdr.value = 0.2, direct = FALSE, zt = seq(0,46,by = 2), absolute.signal=TRUE)
{
	model.orders = list(prod = c(1,2,4), deg = c(1,3,4), all = c(1,4));
	set.nb.param(absolute.signal=absolute.signal);
	LRP = matrix(NA, nrow = nrow(T), ncol = length(unlist(model.orders)) - length(model.orders))
	i = 1
	
	for(P in 1:length(model.orders)){
		order = unlist(model.orders[P])
		for(m in 1:(length(order)-1)){
			m1 = order[m]; m2 = order[m+1]
			eval(parse(text = paste('err1 = T$error.m',m1,sep = '')))
			eval(parse(text = paste('err2 = T$error.m',m2,sep = '')))
			LR = length(zt)*(log(err1)-log(err2));#
			df = n.param[m2]-n.param[m1]
			pvals.LR = pchisq(LR,df = df ,lower.tail = FALSE);
			LRP[,i] = pvals.LR; i = i+1 	
			
		}
		
	}
	LRP = data.frame(LRP, stringsAsFactors = FALSE); colnames(LRP) = paste('LRT.pval',c('p','P','d','D','A'),sep = '.')
	LRP$LRT.pval.C = pchisq(-2*log(LRP$LRT.pval.P * LRP$LRT.pval.D), df = 4, lower.tail = FALSE)
	LOOSE = LRP<pval; colnames(LOOSE) =  c('p','P','d','D','A','C'); LOOSE = data.frame(LOOSE, stringsAsFactors = FALSE)
	STRINGENT = LRP<pval.stringent;  colnames(STRINGENT) = c('p','P','d','D','A','C'); STRINGENT = data.frame(STRINGENT, stringsAsFactors = FALSE)
	
	
	if(V1)
	{
		LRP$LRT.best.model = NA
		LRP$LRT.best.model[!LOOSE$P &!LOOSE$D] = 1
		LRP$LRT.best.model[!LOOSE$P & LOOSE$D] = 2
		LRP$LRT.best.model[LOOSE$P & !LOOSE$D] = 3
		LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
	}else{
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)] =apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[LOOSE$P & LOOSE$D] = 4
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(direct){
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)],LRP$LRT.pval.A[(LOOSE$p|LOOSE$d|LOOSE$A)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(combined){
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)] = apply(cbind(LRP$LRT.pval.p[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)], LRP$LRT.pval.d[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)],LRP$LRT.pval.C[(LOOSE$p|LOOSE$d|LOOSE$C)&!is.na(LOOSE$p)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}
	if(FDR){
		pval.m2 = fdr2pval(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], fdr = fdr.value); MODEL2 = LRP$LRT.pval.p<=pval.m2
		pval.m3 = fdr2pval(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], fdr = fdr.value); MODEL3 = LRP$LRT.pval.d<=pval.m3
		pval.m4 = fdr2pval(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], fdr = fdr.value); MODEL4 = LRP$LRT.pval.A<=pval.m4
		
		fdr.m2 = fdrtool(LRP$LRT.pval.p[!is.na(LRP$LRT.pval.p)], statistic = 'pvalue')
		fdr.m3 = fdrtool(LRP$LRT.pval.d[!is.na(LRP$LRT.pval.d)], statistic = 'pvalue')
		fdr.m4 = fdrtool(LRP$LRT.pval.A[!is.na(LRP$LRT.pval.A)], statistic = 'pvalue')
		FDR.M2 = LRP$LRT.pval.p; FDR.M2[!is.na(LRP$LRT.pval.p)] = fdr.m2$qval
		FDR.M3 = LRP$LRT.pval.p; FDR.M3[!is.na(LRP$LRT.pval.p)] = fdr.m3$qval
		FDR.M4 = LRP$LRT.pval.p; FDR.M4[!is.na(LRP$LRT.pval.p)] = fdr.m4$qval
		
		MODEL2 = FDR.M2<= fdr.value
		MODEL3 = FDR.M3<= fdr.value
		MODEL4 = FDR.M4<= fdr.value
		
		LRP$LRT.best.model = 1
		LRP$LRT.best.model[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)] = apply(cbind(LRP$LRT.pval.p[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)], LRP$LRT.pval.d[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)],LRP$LRT.pval.A[(MODEL2| MODEL3 | MODEL4)&!is.na(MODEL2)]),1, which.min)+1
		LRP$LRT.best.model[is.na(LRP$LRT.pval.p)] = NA
	}	
	
	ok = !is.na(LRP$LRT.best.model)
	
	LRP$LRT.quality = NA;
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 1)] = -log10(1-apply(cbind(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 1)], LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 1)]),1,min))
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 2)] = -log10(LRP$LRT.pval.p[ok & (LRP$LRT.best.model == 2)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 3)] = -log10(LRP$LRT.pval.d[ok & (LRP$LRT.best.model == 3)])
	LRP$LRT.quality[ok & (LRP$LRT.best.model == 4)] = -log10(LRP$LRT.pval.A[ok & (LRP$LRT.best.model == 4)])
	
#	LRP$LRT.quality = NA;
#	LRP$LRT.quality[LRP$LRT.best.model == 1] = apply(cbind(LRP$LRT.pval.p[LRP$LRT.best.model == 1], LRP$LRT.pval.d[LRP$LRT.best.model == 1]),1,min)
#	LRP$LRT.quality[LRP$LRT.best.model == 2] = 1-LRP$LRT.pval.p[LRP$LRT.best.model == 2]
#	LRP$LRT.quality[LRP$LRT.best.model == 3] = 1-LRP$LRT.pval.d[LRP$LRT.best.model == 3]
#	LRP$LRT.quality[LRP$LRT.best.model == 4] = 1-apply(cbind(LRP$LRT.pval.P[LRP$LRT.best.model == 4], LRP$LRT.pval.D[LRP$LRT.best.model == 4]),1,max)
	
	
	#colnames(LRP) = paste(colnames(LRP), '.Laura', sep='')
	return(LRP)
}

AIC.LAURA = function(T = T, correction = FALSE, number.of.data.point = 48, fact = 2, absolute.signal=TRUE)
{
	set.nb.param(absolute.signal=absolute.signal);
	N = number.of.data.point;
	if(correction){corr = 1}else{corr = 0}
	AIC.m1 = 2*n.param[1] + N*log(T$error.m1) + corr*2*n.param[1]*(n.param[1]+1)/(N-n.param[1]-1) # + N*log(2*pi*exp(1)*T$error.m1/N)
	AIC.m2 = 2*n.param[2] + N*log(T$error.m2) + corr*2*n.param[2]*(n.param[2]+1)/(N-n.param[2]-1) # + N*log(2*pi*exp(1)*T$error.m2/N)
	AIC.m3 = 2*n.param[3] + N*log(T$error.m3) + corr*2*n.param[3]*(n.param[3]+1)/(N-n.param[3]-1) # + N*log(2*pi*exp(1)*T$error.m3/N)
	AIC.m4 = 2*n.param[4] + N*log(T$error.m4) + corr*2*n.param[4]*(n.param[4]+1)/(N-n.param[4]-1) # + N*log(2*pi*exp(1)*T$error.m4/N) 
	
	AIC = data.frame(AIC.m1, AIC.m2, AIC.m3, AIC.m4, stringsAsFactors = FALSE)
	AIC$best.model = as.numeric(apply(AIC,1,which.min))
	
	if(correction){
		colnames(AIC) = c('AICc.m1', 'AICc.m2','AICc.m3','AICc.m4', 'AICc.best.model');
		#colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
		
	}else{
		colnames(AIC)[5] = 'AIC.best.model';
		#colnames(AIC) = paste(colnames(AIC), '.Laura', sep='');
	}	
	
	return(AIC)								
}

BIC.LAURA = function(T = T, fact = 2, absolute.signal=TRUE)
{
	set.nb.param(absolute.signal=absolute.signal);
	nb.data = 48
	BIC.m1 = log(nb.data)*n.param[1] + nb.data*log(T$error.m1)
	BIC.m2 = log(nb.data)*n.param[2] + nb.data*log(T$error.m2)
	BIC.m3 = log(nb.data)*n.param[3] + nb.data*log(T$error.m3)
	BIC.m4 = log(nb.data)*n.param[4] + nb.data*log(T$error.m4)
	
	BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
	BIC$best.model = as.numeric(apply(BIC,1,which.min))
	
	colnames(BIC)[5] = 'BIC.best.model'
	#colnames(BIC) = paste(colnames(BIC), '.Laura', sep='')
	return(BIC)								
}	

prob.best.model = function(BIC)
{
    bic = BIC-min(BIC)
    prob.model = exp(-0.5*bic)
    prob.model = prob.model/sum(prob.model)
    
    #prob.model = exp(-0.5*BIC)
    #prob.model = prob.model/sum(prob.model)
    return(prob.model[which.min(BIC)])
}

my.BIC.loglike = function(T, nb.data=96, absolute.signal=TRUE)
{
    set.nb.param(absolute.signal)
    
    index = match(c('error.m1', 'error.m2', 'error.m3', 'error.m4'), colnames(T))
    error.m1 = T[,index[1]]
    error.m2 = T[,index[2]]
    error.m3 = T[,index[3]]
    error.m4 = T[,index[4]]
    
    ## the formula of BIC used here is chi-square+k*ln(n)==error/sigma^2+k(ln(n)) in which sigma of noise is supported to be known.
    BIC.m1 = (error.m1) + log(nb.data)*n.param[1]
    BIC.m2 = (error.m2) + log(nb.data)*n.param[2]
    BIC.m3 = (error.m3) + log(nb.data)*n.param[3]
    BIC.m4 = (error.m4) + log(nb.data)*n.param[4]
    
    BIC = data.frame(BIC.m1, BIC.m2, BIC.m3, BIC.m4, stringsAsFactors = FALSE)
    BIC$best.model = as.numeric(apply(BIC,1, which.min))
    BIC$prob.best.model = as.numeric(apply(BIC[, c(1:4)], 1, prob.best.model))
    
    colnames(BIC) = c('BIC.m1', 'BIC.m2', 'BIC.m3', 'BIC.m4', 'BIC.best.model', 'prob.best.model')

    return(BIC)
}

########
######## Output 
########
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
}	

attribute.global.colors = function()
{
    col.ex <<- 'steelblue';
    col.int <<- 'green3';
    col.deg <<- 'orangered';
}

compare.parameters.RNA.seq.total.fake = function(T = T, select.method='BIC.best.model', parametrization = c('consine.beta'), quality = FALSE)
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
    barplot(t(Prediction.Quality), col = 1:4, names.arg = c('CS-CD','RS-CD','CS-RD','RS-RD'), border = NA, cex.main=0.8, main =paste(select.method))
    
    #### noise (sigma) estimation for premRNA and mRNA
    attribute.global.colors();
    lims = range(c(T$sigma.s, T$sigma.se))
    plot(T$sigma.s, T$sigma.se, ylim=lims, xlim=lims, col=col.int, type='p', cex=0.6, xlab='True values', ylab='Estimated values', cex.main=0.8, main='Sigma of premRNAs')
    abline(0, 1, col='gray', lwd=1.0)
    
    lims = range(c(T$sigma.m, T$sigma.me))
    plot(T$sigma.m, T$sigma.me, ylim=lims, xlim=lims, col=col.ex, type='p', cex=0.6, xlab='True values', ylab='Estimated values', cex.main=0.8,, main='Sigma of mRNAs')
    abline(0, 1, col='gray', lwd=1.0)
    
    #### nothing in this panel
    plot(T$sigma.m, T$sigma.me, type='n', main=NA, xlab=NA, ylab=NA, axes=FALSE)
    
    best.model = T$best.model;
    true.model=T$true.model;
    
	mains = c('Min-premRNA', 'Amp.premRNA','Phase.premRNA','Beta.premRNA','Min.Degr.Rate','Degr.Amp','Degr.Phase','Splicing.K')
	parameter.list = c('Min.int','Amp.int','phase.int','beta.int', 'gamma','eps.gamma','phase.gamma','splicing.k')
	models.p = cbind(c(0,2,3,4), c(0,2,0,4), c(0,2,0,4),c(0,2,0,4), c(0,2,3,4),c(0,0,3,4),c(0,0,3,4),c(0,2,3,4))
	
    bounds = set.bounds(model = 4, absolute.signal=TRUE);
	lower = bounds$lower; 
	upper = bounds$upper;
	#lower[4] = log(2)/(30/60);
    upper[1] = log(2)/(1/60);
    lower[2] = 0.01;
    lower[6] = 0.01;
    upper[6] = 20;
    lower = lower[c(5:8, 1:4)]
    upper = upper[c(5:8, 1:4)]
    upper[1] = 10;
    
	### 8 parpameters and 8 plots
	for(p in 1:8)
	{
		log = ''; 
		if((p==1)|(p==2)|(p==5)|(p==6)|(p==8)){log ='xy'}
        
        plot(c(lower[p], upper[p]),c(lower[p],upper[p]),type = 'n', xlab = 'True values', ylab = 'Estimated values', cex.main=0.8, main = mains[p], log = log)
        
        for(m in c(2:4)) # for each model 2,3,4
        {
            eval(parse(text = paste('ttrue = T$', parameter.list[p],'[best.model == m & true.model==m]', sep ='')))
            eval(parse(text = paste('ftrue = T$', parameter.list[p],'[best.model == m & true.model!=m]', sep ='')))
            
            eval(parse(text = paste('estimated.t = T$', parameter.list[p],'.m',m,'[best.model == m & true.model==m ]', sep ='')));
            eval(parse(text = paste('estimated.f = T$', parameter.list[p],'.m',m,'[best.model == m & true.model!=m ]', sep ='')));
            
            if(models.p[m,p]>0)
            {
                points(ttrue, estimated.t, pch = 21, col = m, bg = m, cex = 0.6);
                points(ftrue, estimated.f, pch = 2, col = m, bg = m, cex = 0.6);
                
            }
        }
        abline(a =0, b = 1, col = 'gray', lwd=1.5)
        
	}
	
}




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
            cat('model', m, '\n');
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


plot.genes.examples = function(index=c(1),pdfname='myplots/gene_examples/Examples_Test.pdf', Tt=T, RF.KO = FALSE, Figure=FALSE, folder=folder)
{
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
        if(RF.KO)
        {
            pdf(pdfname, width = 20, height = 8)
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
        
        
        #col.ex <<- 'steelblue'; col.int <<- 'green3'; col.accepted.ex <<- 'steelblue2'; col.accepted.int <<- 'limegreen'; col.rejected <<- 'gray' ; zt <<- seq(0,46,by=2)
        #col.phases <<- rainbow(n = 241,s = 0.9, v = 0.9);
        #col.deg <<- 'orangered'
        #Ncolors <<- 21; col.rel.ex <<- col.ramp.ex(n = Ncolors); col.rel.int <<- col.ramp.int(n = Ncolors);
        #col.pval <<- col.ramp.pval(n = 11); col.median.level  <<- col.ramp.median(n = 15);col.kept <<- c('yellowgreen','tomato1')
        if(!Figure)
        {
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


plot.summary.comparison.AL.RF = function(T, model=4, folder=folder)
{
    #pdfname='Examples_Test.pdf'; Tt=T; RF.KO = TRUE; Figure=FALSE; folder=folder
    #i = 11; Tt=T; model = 4;
    zt = seq(2,94,by = 4);
    load(file='myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata')
    
    index = which(T$BIC.best.model==model)
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

    pdfname = paste(folder, '/', 'Comparison_AL_RF_summary.pdf', sep='')
    pdf(pdfname, width = 12, height = 6)
    
    par(mfrow = c(2,4))
    
    ## mRNA
    plot(stat1[,2], stat3[,2], cex=0.4, main='Mean Expression', ylab='mRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    
    plot(-log10(stat1[,6]), -log10(stat3[,6]), cex=0.4, main='Rhythmicity (-log10 p-value)', ylab='mRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    abline(v=2, lwd=2.0, col='red')
    abline(h=2, lwd=2.0, col='red')
    
    plot(stat1[,5], stat3[,5], cex=0.4, main='Phase', ylab='mRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    
    plot(stat1[,4], stat3[,4], cex=0.4, main='Relative Amplitude', ylab='mRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    
    ## premRNA
    plot(stat2[,2], stat4[,2], cex=0.4, main='Mean Expression', ylab='premRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    
    plot(-log10(stat2[,6]), -log10(stat4[,6]), cex=0.4, main='Rhythmicity (-log10 p-value)', ylab='premRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    abline(v=2, lwd=2.0, col='red')
    abline(h=2, lwd=2.0, col='red')
    
    plot(stat2[,5], stat4[,5], cex=0.4, main='Phase', ylab='premRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')
    
    plot(stat2[,4], stat4[,4], cex=0.4, main='Relative Amplitude', ylab='premRNA (NRF)', xlab='ALF')
    abline(0, 1, lwd=2.0, col='red')

    dev.off()
}


