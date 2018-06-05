####################### Front:: HERE we started again this project and we really hope that we can finish it this time !!!
#############
####### Import total RNA-seq data
#############
RNA.Seq = FALSE
if(RNA.Seq)
{
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    #### previously used polyA-RNA-seq
    data.version = '_v3'
    load(file=paste('my_genes_RNA_seq_analysis',data.version, '.Rdata', sep=''))
    #source("f24_modified_1.0.r")
    polya = R
    
    R = read.table('/Users/jiwang/RNA_seq_Data/Total_RNA_Seq_AL_RF_KO_Tables/WT_AL_Intron_Exon_RFP.txt',sep='\t', header=TRUE)
    R = R[, -1]
    kk = grep('RFP', colnames(R))
    R = R[, -kk]
    
    ii = grep('WT_AL_Exon', colnames(R))
    jj = grep('WT_AL_Intron', colnames(R))
    ## change mRNA and premRNA data into linear scale and statistics in linear scale
    data1 = 2^R[,ii]
    data2 = 2^R[,jj]

    stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
    
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    
    R = data.frame(R[,1], data1, data2, stat1, stat2, stringsAsFactors=FALSE)
    colnames(R)[1] = 'gene'
    colnames(R)[2:49] = paste('ZT', c(0:47)*2, '.abs.mRNA', sep='')
    colnames(R)[50:97] = paste('ZT', c(0:47)*2, '.abs.premRNA', sep='')
    #stat = colnames(R)[50:56]
    #stat = unlist(strsplit(as.character(stat), '.'))
    colnames(R)[98:104] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.mRNA', sep='')
    colnames(R)[105:111] = paste(c('nb.timepoints', 'mean', 'amp', 'rel.amp', 'phase','pval', 'qv'), '.premRNA', sep='')
    
    kk = which(R$mean.mRNA>R$mean.premRNA)
    R = R[kk,]
    
    data.version = '_total_v1'
    save(R, file=paste('myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))

    examples = c('Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
    mm = match(examples, R$gene)
    
    #load(file=paste('myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))
    prepare.NF.KO.tables();
    
}

#### a favor for friend
CLDN.check = FALSE
if(CLDN.check)
{
    data.version = '_total_v1'
    load(file=paste('myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))
    load(file=paste('myRdata/Cry_Bmal_WT_Bmal_KO_RF.Rdata', sep=''))
    
    mm = grep('Cldn', R$gene)
    
    pdf('myplots/CLDNs_mRNAs_expression_log2.pdf', width=6, height=4)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    for(n in mm)
    {
        gg = R$gene[n]
        y0 = R[n, grep('abs.mRNA', colnames(R))]
        y1 = R.WT.RF[which(R.WT.RF$gene==gg), grep('abs.mRNA', colnames(R.WT.RF))]
        y2 = R.KO.RF[which(R.KO.RF$gene==gg), grep('abs.mRNA', colnames(R.KO.RF))]
        lims = range(y0, y1, y2)
        plot(0, 1, xlim=c(0, 96), ylim=lims, type='n', xlab='time (ZT)', ylab='log2(RPKM).mRNA', main=gg)
        points(c(0:47)*2, y0, type='b', lwd=1.0, col='blue')
        points(c(0:23)*4, y1, type='b', lwd=1.0, col='orange')
        points(c(0:11)*4, y2, type='b', lwd=1.0, col='red')
        legend('topright', lty=c(1), legend = c('WT.AD','WT.RF', 'KO.RF') , cex=1.0, pt.cex=1.5, pt.lwd=1.0, col= c('blue', 'orange', 'red'), border = NA, bty = 'n')

    }
    
    dev.off()

}

########################################
########### some validation for the total RNA-seq data by comparing with qPCR, ExonArray, nascent-RNA-seq
########################################
Data.Control = FALSE
if(Data.Control)
{
    data.version = '_total_v1'
    load(file=paste('myRdata/my_genes_RNA_seq_analysis', data.version, '.Rdata', sep=''))
    
    source("f24_modified_1.0.r")
    source('model_RNA_seq_total/functions.R')
    
    ###### Summary of data: expression of premRNAs and mRNAs, phases of rhythmic premRNAs and mRNAs
    Tt = R;
   
    ii = grep('abs.mRNA', colnames(Tt))
    jj = grep('abs.premRNA', colnames(Tt))
    data1 = log2(as.matrix(Tt[,ii]))
    data2 = log2(as.matrix(Tt[,jj]))
    
    stat1 = t(apply(data1,1, f24_R2_alt2, t=c(0:47)*2))
    stat2 = t(apply(data2,1, f24_R2_alt2, t=c(0:47)*2))
    stat1 = cbind(stat1, qvals(stat1[,6]))
    stat2 = cbind(stat2, qvals(stat2[,6]))
    
    
    sel = which(stat1[,2]>-1)
    Tt = Tt[sel, ]
    
    data.version = '_total_v1'
    save(R, Tt, file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    
    
    pdf('myplots/premRNAs_expression_log2.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    x = stat2[sel,2]
    hist(x,breaks=25,xlab=NA,ylab=NA, main=NA,col='gray',cex.lab=1.0,xlim=c(-5, 15), axes=FALSE)
    abline(v=median(x), lwd=1.5, col='navajowhite4')
    axis(1,at=seq(-5,15,by=3),cex.axis =0.7)
    axis(2,at = seq(0,1000, by=200),las=1,cex.axis = 0.7)
    
    dev.off()
    
    pdf('myplots/mRNAs_expression_log2.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    x = stat1[sel,2]
    hist(x,breaks=50,xlab=NA,ylab=NA, main=NA,col='gray',cex.lab=1.0,xlim=c(-5, 15), axes=FALSE)
    abline(v=median(x), lwd=1.5, col='navajowhite4')
    axis(1,at=seq(-5,15,by=3),cex.axis =0.7)
    axis(2,at = seq(0,800, by=200),las=1,cex.axis = 0.7)
    
    dev.off()
    
    
    library(plotrix)
    library("circular")
    make_circ_coord = function(t,x,ttot=24)
    {
        dt=(t[2]-t[1])*.45
        a=(rep(t,rep(4,length(t)))+rep(c(-dt,-dt,dt,dt),length(t)))*2*pi/ttot
        h=rep(x,rep(4,length(x)))*rep(c(0,1,1,0),length(t))
        list(angles=a,heights=h)
    }
    circular_phase24H_histogram = function(x,color_hist = rgb(0.6,0,0.2), cex.axis=0.5, cex.lab=0.5, lwd=0.5)
    {
        #color_DHS = rgb(0.6,0,0.2)
        par(lwd=lwd,cex.axis=cex.axis, cex.main=0.1,cex.lab=cex.lab)
        #par(mfrow=c(1,1),mar=c(4.5,4.5,1,.5)+.1,las=1)
        br=0:24
        h=hist(x, br=br,plot=FALSE)
        co=make_circ_coord(br[-1],h$counts)
        radial.plot(co$heights,co$angles,br[-1]-br[2], clockwise=TRUE,start=pi/2,main=NA, rp.type='p',poly.col=color_hist)
    }
    
    kk = which(stat2[,7]<0.05)
    phases = stat2[kk, 5]
    
    pdf('myplots/Phases_rhythmic_premRNAs.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.0,0.5,0.2), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    kk = which(stat1[,7]<0.05)
    phases = stat1[kk, 5]
    
    pdf('myplots/Phases_rhythmic_RNAs.pdf', width=2.2, height=2.2)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    circular_phase24H_histogram(phases, col=rgb(0.7,0.,0.2), cex.axis=0.7, cex.lab=0.01, lwd=0.5)
    
    dev.off()
    
    ### compare with half-lives
    half.life.comparison = FALSE
    if(half.life.comparison)
    {
        shar = read.table('/Users/jiwang/Degradation_Liver/_x_MicroArray/complementary_data/half-lives/Sharova.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
        fried = read.table('/Users/jiwang/Degradation_Liver/_x_MicroArray/complementary_data/half-lives/Friedel.txt', stringsAsFactors = FALSE, sep = '\t', header = TRUE)
        schwa = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Annotations/Schwarnhausser_2011.csv', header=TRUE)
        ii = match(c("Gene.Names", "mRNA.half.life.replicate..h."), colnames(schwa))
        schwa = schwa[,ii]
        schwa = schwa[which(!is.na(schwa[,2])==TRUE), ]
        gene = c()
        
        for(n in 1:nrow(schwa))
        {
            test = schwa[n, 1]
            test = unlist(strsplit(as.character(test), ';'))
            mm = match(test, Tt$gene)
            if(length(which(!is.na(mm)==TRUE))>0) {gene = c(gene, as.character(test[which(!is.na(mm)==TRUE)][1]))
            }else{
                gene = c(gene, as.character(schwa[n, 1]))
            }
        }
        schwa[, 1] = gene
        
        colnames(shar) = c('gene', 'half.life')
        colnames(fried) = c('gene', 'probeset.ID','hlmin','half.life')
        shar = shar[which(shar[,2]>0),]
        fried = fried[which(fried[,4]>0),]
        fried = fried[, c(1,4)]
        colnames(schwa) = c('gene', 'half.life')
        #save(shar, fried, schwa, file='myRdata/mRNAs_half_lives_databases.Rdata')
        
        overlap2 = intersect(Tt$gene, shar[,1])
        overlap2 = intersect(overlap2, fried[,1])
        overlap2 = intersect(overlap2, schwa[,1] )
        
        ii = match(overlap2, Tt$gene)
        
        pdf('myplots/Compare_half-lives_polyA_Total_Exon_Array_Nascent.pdf', width=16.0, height=6.0)
        par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        par(mfrow=c(1,3))
        
        xx1 = as.numeric(Tt$mean.mRNA[ii]/Tt$mean.premRNA[ii])
        mm = match(overlap2, fried[,1])
        yy1 = as.numeric(fried[mm, 2])
        
        cex = 0.7
        plot(xx1, yy1, cex=cex, xlab='total', ylab='fried', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        mm = match(overlap2, shar[,1])
        yy1 = as.numeric(shar[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='total', ylab='shar', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        mm = match(overlap2, schwa[,1])
        yy1 = as.numeric(schwa[mm, 2])
        plot(xx1, yy1, cex=cex, xlab='total', ylab='schwa', main=paste('R = ', signif(cor(xx1, yy1), d=2), sep=''), log='xy')
        cor(xx1, yy1)
        
        dev.off()
        
    }
    
}

####################
########## Determine the noise in the premRNA and mRNA data
####################
Noise.Determination = FALSE
if(Noise.Determination)
{
    
    cat('ESTIMATION OF NOISE for RNA-seq data\n')
    
    data.version = '_total_v1'
    load(file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    source("f24_modified_1.0.r")
    source('model_RNA_seq_total/functions.R')
    T = Tt;
    
    exon = as.matrix(T[,c(2:49)])
    intron = as.matrix(T[,c(50:97)])
    res.ex = as.matrix(T[,c(98:104)])
    res.int = as.matrix(T[,c(105:111)])
    
    ######
    ###### variances vs mean
    ######
    mean.s = apply(log(intron), 1, mean)
    var.s = apply(log(intron), 1, var)
    mean.m = apply(log(exon), 1, mean)
    var.m = apply(log(exon), 1, var)
    genes = T$gene
    
    non.rhythmic.sel = TRUE
    if(non.rhythmic.sel)
    {
        cutoff = 0.1
        kk = which(res.ex[,6]>cutoff & res.int[,6]>cutoff)
        
        pdf('myplots/Noise_estimation_Total.pdf', width=6.0, height=3.0)
        par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
        
        par(mfrow=c(1,2))
        
        cex = 0.1
        
        xlim = range(c(mean.s, mean.m))
        ylim = range(c(var.s, var.m))
        #ylim = c(0, 0.5)
        plot(mean.s, var.s, xlim=xlim, ylim=ylim, cex=cex, col='blue', xlab='mean', ylab='var')
        points(mean.m, var.m, cex=cex, col='black')
        
        plot(var.m, var.s, cex=cex, log='xy')
        ss = seq(min(var.m), max(var.m), length.out=100)
        points(ss, ss, cex=cex, col='red', type='l', lwd=2.0)
        #points(ss, 2*ss, cex=cex, col='green', type='l')
        #points(ss, 4*ss, cex=cex, col='orange', type='l')
        #abline(0, 0, col='red', lwd=1.5)
        #abline(0, log(2), col='red', lwd=1.5)
        #abline(0, log(4), col='red', lwd=1.5)
        dev.off()
        
    }
    
    ######
    ###### estimate the gene-specific variances
    ######
    source('model_RNA_seq_total/functions.R')
    zt = seq(0, 94, by=2)
    Min = 5
    Amp = 1.0
    beta = 5.0
    phase = 6;
    s1 = compute.s.beta(zt, Min, Amp, phase, beta);
    
    fold.change = Amp/Min+1;
    mean = Min*(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)
    s2 = compute.s.beta.v1(zt, mean, fold.change, phase, beta)
    
    gamma = 4.1588831;
    amp.gamma = 0.2;
    phase.gamma = 15;
    splicing.k = log(2)/(5/60);
    m1 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta, simulation.only=FALSE)
    mm1 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta, simulation.only=TRUE)
    
    S0 = s1
    M0 = m1
    
    source('model_RNA_seq_total/functions_28_Jan_2015.R')
    eps.gamma = amp.gamma/(2*gamma+amp.gamma)
    m2 = compute.m.beta(zt, gamma/(1-eps.gamma), eps.gamma, phase.gamma, splicing.k, mean, fold.change, phase, beta)
    
    PLOT.Comparison = FALSE
    if(PLOT.Comparison)
    {
        s1 = s1/mean(s1)
        s2 = s2/mean(s2)
        m1 = m1/mean(m1)
        m2 = m2/mean(m2)
        
        plot(zt, s1, type='l', col='red', lwd=2.0, ylim=range(s1, s2, m1, m2)) ################ GOOD, compute.s.beta works perfectly
        points(zt, s2, type='l', col='darkred', lwd=1.5)
        points(zt, m1, type='l', col='blue', lwd=1.5)
        points(zt, m2, type='l', col='darkblue', lwd=1.5)
        abline(h=1, col='gray', lwd=3.0)
        abline(h=0, col='gray', lwd=3.0)
    }
    
    #### Add noise
    test = lseq(0.001, 2, length.out=1000)
    keep = c()
    S0 = S0[1:24]
    M0 = M0[1:24]
    for(n in 1:length(test))
    {
        log.abs.int = log(S0);
        log.abs.ex = log(M0);
        
        #set.global.sigma();
        sigma.s = test[n]
        sigma.m = test[n]
        
        noise.int = rnorm(length(S0), mean = 0, sd = sigma.s)
        noise.ex =  rnorm(length(M0), mean = 0, sd = sigma.m)
        
        log.noise.abs.int = log.abs.int + noise.int
        log.noise.abs.ex = log.abs.ex + noise.ex
        
        S = exp(log.noise.abs.int)
        M = exp(log.noise.abs.ex)
        source('model_RNA_seq_total/functions.R')
        test.s = variance.estimate.replicates(log(S), nb.replicates=2)
        test.m = variance.estimate.replicates(log(M), nb.replicates=2)
        #print(c(sigma.s, test.s))
        #print(c(sigma.m, test.m))
        keep = rbind(keep, c(test.s, test.m))
    }
    
    #### check genome-wide noise
    exon = log(as.matrix(T[,c(2:49)]))
    intron = log(as.matrix(T[,c(50:97)]))
    
    test = cbind(apply(exon, 1, variance.estimate.replicates, nb.replicates=4), apply(intron, 1, variance.estimate.replicates, nb.replicates=4))
    
    par(mfrow=c(1,3))
    hist(test[,1], breaks=100, main='Exon estimated Sigma')
    abline(v=median(test[,1]), col='red')
    hist(test[,2], breaks=100, main='Exon estimated Sigma', col='gray')
    abline(v=median(test[,2]), col='red')
    plot(test[,1], test[,2], xlab='Exon', ylab='Intron', log='', cex=0.5)
    abline(0, 1, col='red', lwd=2.0)
    
    pdf('myplots/Noise_estimation_by_replicates.pdf', width=5.0, height=5.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    
    cex = 0.7
    
    xlim = range(keep, test)
    ylim = xlim
    plot(test, keep[,1], xlim=xlim, ylim=ylim, cex=cex, col='blue', xlab='true sigma', ylab='estimated sigma', log='xy')
    points(test, keep[,2], cex=cex, col='black')
    abline(0, 1, col='red', lwd=2.0)
    
    dev.off()

    plot(zt, S, ylim=range(S, M), type='p', cex=1.0, col='red')
    points(zt, M, type='p', cex=1.0, col='blue')
    points(zt, S0, type='l', lwd=2.0, col='red')
    points(zt, M0, type='l', lwd=2.0, col='blue')
    
}

##################
######## Generate simulated FITS with real data, fit the model and check the parameter estimatioin
##################
my.do.step = TRUE
data.version = '_total_v1'
fake.version = '_total_s1'

if(my.do.step)
{
    cat('GENERATE THE FAKE DATA\n')
    
    load(file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    source("f24_modified_1.0.r")
    T = Tt;
    
    #### 100 fake genes for each model
    X = 100
    zt = seq(0, 94, by=2)
    ZT.int = grep('.abs.premRNA', colnames(T));
    ZT.ex = grep('.abs.mRNA', colnames(T));
    
    F = T[1:(4*X),];
    source('model_RNA_seq_total/functions.R')
    
    for(model in c(1:4))
    {
        cat('model : ',model, '\n')
        Fm = generate.fake.data(T = F[((model-1)*X+1):(model*X),],X = X, model = model, zt=zt)
        if(model == 1){FF = Fm}
        else{FF = rbind(FF,Fm)}
    }
    
    F = FF;
    
    ### statistics
    exon = as.matrix(log2(F[,ZT.ex]))
    intron = as.matrix(log2(F[,ZT.int]))
    
    res.ex = t(apply(exon,1, f24_R2_alt2, t=c(0:47)*2))
    res.int = t(apply(intron,1, f24_R2_alt2, t=c(0:47)*2))
    res.ex = cbind(res.ex, qv=qvals(res.ex[,6]))
    res.int = cbind(res.int, qv=qvals(res.int[,6]))
    
    F[,c(98:104)] = res.ex
    F[,c(105:111)] = res.int
    
    #### check fake data (rhythmicity, ampitudes for model 2, 3 and 4, especially for model 4)
    sigma.se = apply(log(F[, ZT.int]), 1, variance.estimate.replicates, nb.replicates=4)
    sigma.me = apply(log(F[, ZT.ex]), 1, variance.estimate.replicates, nb.replicates=4)
    
    F$sigma.se = sigma.se
    F$sigma.me = sigma.me
    
    pdf('myplots/Total/simulated_data_control.pdf', width=6, height=6)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,3,0.8)+0.1, tcl = -0.3)

    par(mfrow=c(2,2))
    lims = range(c(T$sigma.s, T$sigma.se))
    plot(F$sigma.s, sigma.se, ylim=lims, xlim=lims, type='p', cex=0.7, xlab='true value', ylab='estimated value', main='sigma of premRNAs')
    abline(0, 1, col='red', lwd=2.0)
    
    lims = range(c(T$sigma.m, T$sigma.me))
    plot(F$sigma.m, sigma.me, ylim=lims, xlim=lims, type='p', cex=0.7, main='sigma of mRNAs')
    abline(0, 1, col='red', lwd=2.0)
    
    kk = c(101:400)
    model = c(rep(2, 100), rep(3, 100), rep(4, 100))
    boxplot(-log10(F$pval.mRNA[kk]) ~ model, col=c('green', 'red', 'gray'), xlab='Model', ylab='-log10(pval)')
    boxplot((F$amp.mRNA[kk]) ~ model, col=c('green', 'red', 'gray'), xlab='Model', ylab='log2 fold change')
    
    dev.off()
    
    write.table(F, file = paste('simulated_data/my_simulated_data_RNA_seq',fake.version,'.txt',sep = ''), quote = FALSE, sep = '\t', row.names = FALSE)
    #eval(parse(text = paste('T',100*sd,'= F',sep ='')))
    
    save(F, file=paste('myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
    
    #############
    ### test the model selection and parameter estimation
    #############
    Test.Fitting = FALSE
    if(Test.Fitting)
    {
        cat('TEST THE FIT on THE FAKE DATA\n')
        fake.version = '_total_s1'
        load(file=paste('myRdata/my_simulated_data', fake.version, '.Rdata', sep=''))
        source("f24_modified_1.0.r")
        T = F;
        
        source('model_RNA_seq_total/functions.R')
        #best.model = model.check
        ZT.int = grep('.abs.premRNA', colnames(T))
        ZT.ex = grep('.abs.mRNA', colnames(T))
        zt = seq(0,94,by = 2)
        
        j = 276
        
        ptm <- proc.time()
        param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int);
        proc.time() - ptm
        
        #c(param.fits.results, my.BIC.this.gene(param.fits.results))
        c(param.fits.results, my.BIC.this.gene.loglike(param.fits.results))
    }
    
    Fitting.results  = FALSE
    if(Fitting.results)
    {
        #### Model selection for simulated data
        res.version = '_total_fake_s1';
        load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
        source('functions_RPKM_gamma_amp_v2.R')
        
        BIC = my.BIC.loglike(T = T)
        T = data.frame(T, BIC, stringsAsFactors = FALSE);
        T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
        Tt = T
        
        length(which(Tt$BIC.best.model==1))
        length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4))
        
        #### check parameter estimation and find good combination of parameters to fit (for NB model)
        Test.parameter.combination = FALSE
        if(Test.parameter.combination)
        {
          kk = which(T$BIC.best.model==3)
          gamma.true = T$gamma[kk]+0.5*T$amp.gamma[kk];
          gamma.estimate = (T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk])
          
          eps.true = T$amp.gamma[kk]/(T$gamma[kk]+0.5*T$amp.gamma[kk])/2;
          eps.estimate = T$amp.gamma.m3[kk]/(T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk])/2;
          #eps.true = T$amp.gamma[kk]/2;
          #eps.estimate = T$amp.gamma.m3[kk]/2
          ratio1 = gamma.estimate/gamma.true;
          ratio2 = eps.estimate/eps.true;
          w = 2*pi/24
          par(mfrow=c(1,2),cex=1.0)
          plot(gamma.true, eps.true, log='x')
          jj1 = which(ratio1>10 | ratio1<0.1)
          points(gamma.true[jj1], eps.true[jj1], col='red')
          ii1 = which(ratio1>10 | ratio1<0.1)
          abline(v=w)
          
          plot(gamma.true, eps.true, log='x')
          jj2 = which(ratio2>2 | ratio2<0.5)
          points(gamma.true[jj2], eps.true[jj2], col='red')
          abline(v=w)
          
          par(mfrow=c(2,2),cex=1.0)
          
          plot((T$gamma[kk]), (T$gamma.m3[kk]), cex=0.5, col='red', log='xy')
          abline(0, 1, lwd=2.0, col='gray')
          plot(T$amp.gamma[kk], T$amp.gamma.m3[kk], cex=0.5, col='red', log='xy')
          abline(0, 1, lwd=2.0, col='gray')
          
          #par(mfrow=c(1,2),cex=1.0)
          plot((T$gamma[kk]+0.5*T$amp.gamma[kk]), (T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk]), cex=0.5, col='red', log='xy')
          abline(0, 1, lwd=2.0, col='gray')
          
          #par(mfrow=c(1,2),cex=1.0)
          plot(T$amp.gamma[kk]/(T$gamma[kk]+0.5*T$amp.gamma[kk])/2, T$amp.gamma.m3[kk]/(T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk])/2, cex=0.5, 
               col='red', log='', xlim=c(0, 1), ylim=c(0, 1))
          abline(0, 1, lwd=2.0, col='gray')
          
          #plot(T$phase.gamma[kk], T$phase.gamma.m3[kk], cex=0.5, col='red')
          #abline(0, 1, lwd=2.0, col='gray')
          w = 2*pi/24
          plot(T$amp.gamma[kk]/(T$gamma[kk]+0.5*T$amp.gamma[kk])/2/sqrt(1+(w/(T$gamma[kk]+0.5*T$amp.gamma[kk]))^2), 
               T$amp.gamma.m3[kk]/(T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk])/2/sqrt(1+(w/(T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk]))^2), cex=0.5, col='red', log='')
          abline(0, 1, lwd=2.0, col='gray')
          
          kk = which(T$BIC.best.model==3 & T$rel.amp.mRNA>0.4)
          gamma.true = T$gamma[kk]+0.5*T$amp.gamma[kk];
          gamma.estimate = (T$gamma.m3[kk]+0.5*T$amp.gamma.m3[kk])
          
          par(mfrow=c(1,3),cex=1.0)
          plot(T$phase.gamma[kk], T$phase.gamma.m3[kk], cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          plot((T$phase.gamma[kk]+ atan2(w, gamma.true)/w), (T$phase.gamma.m3[kk] + atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          plot(atan2(w, gamma.true)/w, atan2(w, gamma.estimate)/w, col='red')
          #plot((T$phase.gamma[kk]+ atan2(w, gamma.true)/w), (T$phase.gamma.m3[kk] + 0*atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          #abline(0, 1, lwd=2.0, col='gray')
          #plot(T$phase.gamma[kk], (T$phase.gamma.m3[kk]-atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          
          
          
          kk = which(T$BIC.best.model==4)
          gamma.true = T$gamma[kk]+0.5*T$amp.gamma[kk];
          gamma.estimate = (T$gamma.m4[kk]+0.5*T$amp.gamma.m4[kk])
          par(mfrow=c(1,3),cex=1.0)
          plot(T$phase.gamma[kk], T$phase.gamma.m4[kk], cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          plot((T$phase.gamma[kk]+ atan2(w, gamma.true)/w), (T$phase.gamma.m4[kk] + atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          #plot(atan2(w, gamma.true)/w, atan2(w, gamma.estimate)/w, col='red')
          #plot((T$phase.gamma[kk]+ atan2(w, gamma.true)/w), (T$phase.gamma.m3[kk] + 0*atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          #abline(0, 1, lwd=2.0, col='gray')
          plot(T$phase.gamma[kk], (T$phase.gamma.m4[kk]-atan2(w, gamma.estimate)/w), cex=0.7, col='red')
          abline(0, 1, lwd=2.0, col='gray')
          
        }
        
        #### Plot selected mdoels and estimated parameters
        Plot = FALSE
        if(Plot)
        {
            source('model_RNA_seq_total/functions.R')
            library(stringr)
            #selecting.approach = c('LRT.best.model','LRT.best.model.FDR','BIC.best.model','AICc.best.model','AIC.best.model')
            selecting.approach = c('BIC.best.model')
            #sd.noise = c(0.1, 0.25, 0.5)
            #set.global.sigma();
            #filter = TRUE
            
            X = 100
            pdf.name = paste(file = 'myplots//Total/my_simulated_data_with_fits_results_BIC_AIC_LRT', res.version, '.pdf', sep='')
            pdf(pdf.name, width = 10, height = 7.5)
            par(mfrow = c(3,4), cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
            
            for(select.method in selecting.approach)
            {
                select.method = 'BIC.best.model'
                T$true.model = ceiling(c(1:400)/100)
                eval(parse(text = paste('T$best.model=T$', select.method, sep='')))
                kk = which(!is.na(T$best.model)==TRUE)
                T = T[kk,]
                
                compare.parameters.RNA.seq.total.fake(T, select.method)
            }
            
            dev.off()
            
        }


    }
    
    ###### Check individual example
    check.individual.example = FALSE
    if(check.individual.example)
    {
        absolute = TRUE
        sd = 0.5
        cat(sd, '\n')
        eval(parse(text = paste('T = T',100*sd,sep ='')))
        
        ZT.int = grep('.abs.premRNA', colnames(T))
        ZT.ex = grep('.abs.mRNA', colnames(T))
        ### make plots for each fake gene
        pdf.name = paste('myplots/Plots_fake_data_check_log', fake.version,'_', sd*100, '.pdf', sep='')
        pdf(pdf.name, width = 12, height = 7)
        
        zt.p = seq(0, 46, by=2)
        par(mfcol=c(1,2),cex=1.0)
        for(i in 1:nrow(T))
        {
            #i = 400
            eval(parse(text = paste('gamma = T$gamma[i];amp.gamma = T$amp.gamma[i];phase.gamma = T$phase.gamma[i];splicing.k = T$splicing.k[i]; Min=T$Min.int[i];Amp = T$Amp.int[i]; phase = T$phase.int[i]; beta = T$beta.int[i]', sep = '')))
            
            mt = compute.m.beta(zt.p, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
            #mmt = compute.m.beta(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, splicing.k = splicing.k, mean = mean, fold.change = fold.change, phase = phase.int, beta = beta, simulation.only=TRUE)
            #mmt = mmt/mean(mmt)
            st = compute.s.beta(zt.p, Min, Amp, phase, beta)
            st = (st)
            mt = (mt)
            Mt = (as.numeric(T[i,ZT.ex]))
            St = (as.numeric(T[i,ZT.int]))
            plot(zt, St, col='darkblue',type='b',lty=1, main=paste(T$gene[i], ', premRNA, phase=', signif(T$phase.premRNA[i],d=3), ', rel.ampl= ', signif(T$rel.amp.premRNA[i],d=3), sep=''))
            points(zt.p, st, col='darkred', type='l',lty=1, lwd=2.0)
            abline(h=mean(St),col='gray',lwd=2.0)
            abline(h=log(Min),col='gray50',lwd=2.0)
            print(mean(st)-mean(St));
            #print(mean(noise.int[i%%100,]))
            
            plot(zt, Mt, col='darkgreen', type='b',lty=1, main=paste(T$gene[i], ', mRNA, beta=', signif(T$beta.int[i],d=2), ', rel.ampl= ', signif(T$rel.amp.mRNA[i],d=3), sep=''))
            points(zt.p, mt, col='darkred', type='l',lty=1,lwd=2.0)
            #points(zt.p, mmt, col='darkorange', type='l',lty=1,lwd=2.0)
            abline(h=mean(Mt),col='gray',lwd=2.0)
            abline(h=(splicing.k*Min/gamma),col='gray50',lwd=2.0)
            
        }
        dev.off()
        
        ####
        #### check the weired behavior profiles because of beta parameter
        ####
        i = 302
        par(mfcol=c(1,2),cex=1.0)
        eval(parse(text = paste('gamma = T$gamma[i];eps.gamma = T$eps.gamma[i];phase.gamma = T$phase.gamma[i];splicing.k = T$splicing.k[i]; mean=T$mean.int[i];fold.change = T$fold.change.int[i]; phase.int = T$phase.int[i]; beta = T$beta.int[i]', sep = '')))
        beta = 1
        eps.gamma = 1.0
        cat(c(log(2)/gamma, eps.gamma, phase.gamma, splicing.k, mean, fold.change,phase.int, beta))
        mt = compute.m.beta(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, splicing.k = splicing.k, mean = mean, fold.change = fold.change, phase = phase.int, beta = beta)
        mt = mt/mean(mt)
        mmt = compute.m.beta(t = zt.p, gamma = gamma, eps.gamma = eps.gamma, phase.gamma = phase.gamma, splicing.k = splicing.k, mean = mean, fold.change = fold.change, phase = phase.int, beta = beta, simulation.only=TRUE)
        mmt = mmt/mean(mmt)
        st = compute.s.beta(t = zt.p, mean = mean, fold.change = fold.change, phase = phase.int, beta = beta)
        st = st/mean(st)
        
        plot(zt, as.numeric(T[i,ZT.int])/mean(as.numeric(T[i,ZT.int])), col='darkblue',type='b',lty=1, ylim=c(0.1,3.0), main=paste(T$gene[i], ', premRNA', sep=''))
        points(zt.p, st, col='darkred', type='l',lty=1, lwd=2.0)
        abline(h=1.0,col='gray',lwd=2.0)
        
        plot(zt, as.numeric(T[i,ZT.ex])/mean(as.numeric(T[i,ZT.ex])), col='darkgreen', type='b',lty=1, ylim=range(0.1, 3.0), main=paste(T$gene[i], ', mRNA', sep=''))
        points(zt.p, mt, col='darkred', type='l',lty=1,lwd=2.0)
        points(zt.p, mmt, col='darkorange', type='l',lty=1,lwd=2.0)
        abline(h=1.0,col='gray',lwd=2.0)
    }
    
}

########
######## FITS WITH REAL DATA (test part)
########
Real.Data.Fitting = FALSE
if(Real.Data.Fitting)
{
    cat('TEST THE FIT on THE REAL DATA\n')
    data.version = '_total_v1'
    load(file=paste('myRdata/my_genes_RNA_seq_analysis_sel', data.version, '.Rdata', sep=''))
    source("f24_modified_1.0.r")
    T = Tt;
    
    source('model_RNA_seq_total/functions.R')
    #best.model = model.check
    ZT.int = grep('.abs.premRNA', colnames(T))
    ZT.ex = grep('.abs.mRNA', colnames(T))
    zt = seq(0,94,by = 2)
    
    j = which(T$gene=='Dyrk1b')
    
    ptm <- proc.time()
    param.fits.results = make.fits.with.all.models.for.one.gene(T = T, gene.index = j, debug = TRUE, zt = zt, i.ex = ZT.ex, i.int = ZT.int);
    proc.time() - ptm
    
    #c(param.fits.results, my.BIC.this.gene(param.fits.results))
    c(param.fits.results, my.BIC.this.gene.loglike(param.fits.results))
    #T[j, 64:ncol(T)]
    
}

Fitting.Results = FALSE
if(Fitting.Results)
{
    Borrow.names = FALSE
    if(Borrow.names)
    {
        #res.version = '_v6'
        #load(file = paste('myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
        #names = colnames(T)[64:105]
        #mm = grep('stderr', names)
        names = names[-mm]
    }
    
    #res.version = '_total_test_v4'
    res.version = '_total_all_v1';
    load(file = paste('/Users/jiwang/Degradation_Liver/Main_Code_Degradation/myRdata/my_genes_fits_optim_RNA_seq_beta_cos', res.version, '.Rdata',sep = ''))
    source('model_RNA_seq_total/functions.R')
    
    BIC = my.BIC.loglike(T = T)
    T = data.frame(T, BIC, stringsAsFactors = FALSE);
    T = T[which(!is.na(T$BIC.best.model)==TRUE), ]
    
    Tt = T
    
    #### Result-check part
    fitting.check = FALSE
    if(fitting.check)
    {
        length(which(Tt$BIC.best.model==1))
        length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4))
        
        length(which(Tt$BIC.best.model==2))/length(which(Tt$BIC.best.model!=1))
        length(which(Tt$BIC.best.model==3 | Tt$BIC.best.model==4))/length(which(Tt$BIC.best.model!=1))
        
        length(which(Tt$BIC.best.model==3))/length(which(Tt$BIC.best.model!=1))
        length(which(Tt$BIC.best.model==4))/length(which(Tt$BIC.best.model!=1))
        
        #### Check Half-lives, half-life amplitudes, splicing.k
        source('model_RNA_seq_total/functions.R')
        bounds = set.bounds(model=4)
        upper = bounds$upper
        lower = bounds$lower
        
        gamma.max = log(2)/(1/60)
        gamma.min = log(2)/24;
        splicing.k.max = upper[4]
        splicing.k.min = lower[4]
        amp.gamma.max = upper[2]
        
        ### 10 min half-life boundary
        length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
        #### 24 h half-life boundary
        length(which(Tt$BIC.best.model==2 & (abs(Tt$gamma.m2-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3 & (abs(Tt$gamma.m3-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4 & (abs(Tt$gamma.m4-gamma.min)<0.00001)))/length(which(Tt$BIC.best.model==4))
        #### amplitude of gamma boundary
        length(which(Tt$BIC.best.model==3 & (abs(Tt$amp.gamma.m3-amp.gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4 & (abs(Tt$amp.gamma.m4-amp.gamma.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
        #### 1s splicing.k boundary
        length(which(Tt$BIC.best.model==2 & (abs(Tt$splicing.k.m2-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3 & (abs(Tt$splicing.k.m3-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4 & (abs(Tt$splicing.k.m4-splicing.k.max)<0.00001)))/length(which(Tt$BIC.best.model==4))
        #### 30min splicing.k boundary
        length(which(Tt$BIC.best.model==2 & (abs(Tt$splicing.k.m2-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==2))
        length(which(Tt$BIC.best.model==3 & (abs(Tt$splicing.k.m3-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==3))
        length(which(Tt$BIC.best.model==4 & (abs(Tt$splicing.k.m4-splicing.k.min)<0.00001)))/length(which(Tt$BIC.best.model==4))
        
        T[which(T$gene=='Arhgef18'), ]
        
        ##### Check individual examples
        if(!file.exists('myplots/Total/gene_examples/')){dir.create('myplots/Total/gene_examples/')}
        
        source('model_RNA_seq_total/functions.R')
        
        for(m in c(2:4))
        {
            index = which(T$BIC.best.model==m)
            index = index[order(T$pval.mRNA[index])]
            #index = kk[order(T$BIC.best.model[kk])]
            #index = c(1:10)
            pdfname = paste('myplots/Total/gene_examples/fitting_results_all_model_', m, '.pdf', sep='')
            plot.genes.examples(index, pdfname, T)
        
        }
        
        #### Check half-lives of mRNAs
        source('model_RNA_seq_total/functions.R')
        pdfname = 'myplots/Total/gene_examples/Comparison_half_lives_v2.pdf'
        plot.half.life.comparison(index, pdfname, T)
        
        ############
        ### check the noise for the genes with bad fitting
        ############
        noise.check.4fitting = FALSE
        if(noise.check.4fitting)
        {
            exon = log(as.matrix(T[,c(2:49)]))
            intron = log(as.matrix(T[,c(50:97)]))
            test = cbind(apply(exon, 1, variance.estimate.replicates, nb.replicates=4), apply(intron, 1, variance.estimate.replicates, nb.replicates=4))
            
            plot(test[,1], test[,2], xlab='Exon sigma', ylab='Intron sigma', log='', cex=0.5)
            mm = match(c('Hlf', 'Ppara'), T$gene)
            points(test[mm, 1], test[mm, 2], cex=0.7, col='red')
            abline(0, 1, col='red', lwd=2.0)
        }

    }
    
    #######################################
    #######################################
    ####### Result-plot part
    #######################################
    #######################################
    
    ##### Different senarios
    source('model_RNA_seq_total/functions.R')
    
    zt = seq(0, 24, by=0.1)
    Min = 5
    Amp = 5.0
    beta = 1.0
    phase = 6;
    s2 = compute.s.beta(zt, Min, Amp, phase, beta);
    
    amp.gamma = 0.0;
    phase.gamma = 18;
    splicing.k = log(2)/(5/60);
    gamma = log(2)/0.5;
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
    
    Amp = 0.0;
    gamma = log(2)/4;
    amp.gamma = gamma;
    s3 = compute.s.beta(zt, Min, Amp, phase, beta);
    m3 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    deg = gamma + amp.gamma*(1+cos(2*pi/24*(zt-phase.gamma)))/2
    s3 = s3/mean(s3)
    m3 = m3/mean(m3)
    deg = deg/mean(deg)
    
    lims = range(s3, m3, deg)
    
    pdf.name = paste("myplots/Total/Illustration_M3.pdf", sep='')
    pdf(pdf.name, width=2.4, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(zt, s3, type='l', lwd=2.5, lty=1, col=col.int, ylim=lims, xlab=NA, ylab=NA, main=NA, axes=FALSE)
    points(zt, m3, type='l', lwd=2.5, col=col.ex, lty=1)
    #points(zt, deg, type='l', lwd=2.5, lty=1, col=col.deg)
    box()
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=3), las=1,cex.axis = 1.0)
    #abline(v=6, col='gray', lwd=2.0)
    #abline(v=12, col='gray', lwd=2.0)
    
    dev.off()
    
    Amp = Min*0.5
    #gamma = log(2)/4;
    gamma = log(2)/0.5;
    amp.gamma = gamma*5;
    phase.gamma = 2;
    deg = gamma + amp.gamma*(1+cos(2*pi/24*(zt-phase.gamma)))/2
    #print(deg/mean(deg))
    m20 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    gamma = log(2)/4
    
    
    amp.gamma = gamma;
    deg = gamma + amp.gamma*(1+cos(2*pi/24*(zt-phase.gamma)))/2
    #print(deg/mean(deg))
    m21 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    gamma = log(2)/20
    amp.gamma = gamma;
    m22 = compute.m.beta(zt, gamma, amp.gamma, phase.gamma, splicing.k, Min, Amp, phase, beta)
    
    deg = gamma + amp.gamma*(1+cos(2*pi/24*(zt-phase.gamma)))/2
    #print(deg/mean(deg))
    
    s2 = s2/mean(s2)
    m20 = m20/mean(m20)
    m21 = m21/mean(m21)
    m22 = m22/mean(m22)
    deg = deg/mean(deg)
    
    lims = range(s2, m21, deg)
    
    pdf.name = paste("myplots/Total/Illustration_M4_2.pdf", sep='')
    pdf(pdf.name, width=2.4, height=2.0)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,2,0.8)+0.1, tcl = -0.3)
    plot(zt, s2, type='l', lwd=2.5, lty=1, col=col.int, ylim=lims, xlab=NA, ylab=NA, main=NA, axes=FALSE)
    #points(zt, m20, type='l', lwd=2.0, col=col.ex, lty=1)
    points(zt, m21, type='l', lwd=2.5, col=col.ex, lty=1)
    #points(zt, m22, type='l', lwd=2.0, col=col.ex, lty=3)
    #points(zt, deg, type='l', lwd=2.0, lty=1, col=col.deg)
    box()
    axis(2, cex.axis=1.0)
    axis(1, at=seq(0, 24, by=3), las=1,cex.axis = 1.0)
    abline(v=6, col='gray', lwd=2.0)
    abline(v=12, col='gray', lwd=2.0)
    
    dev.off()
    
    
    #### gene examples with summary of fitting
    source('model_RNA_seq_total/functions.R')
    for(m in c(2:4))
    {
        folder = paste('myplots/Total/gene_examples/Model_', m, sep='')
        if(!file.exists(folder)){dir.create(folder)}
        
        index = which(T$BIC.best.model==m)
        index = index[which(T$qv.mRNA[index]<0.005)]
        #if(m==2){index = index[which(T$qv.premRNA<0.005)]}
        #else{}
        #index = index[order(T$pval.mRNA[index])]
        #index = kk[order(T$BIC.best.model[kk])]
        #index = c(1:10)
        #pdfname = paste('myplots/gene_examples/fitting_results_all_model_', m, '.pdf', sep='')
        plot.genes.examples(index, T, Figure=TRUE, folder=folder)
        
    }
    
    ##########
    ## partition of models and percentages of mRNA post-transcriptionally regulated
    ##########
    source('model_RNA_seq_total/functions.R')
    length(which(Tt$BIC.best.model==1))/nrow(T);
    length(which(Tt$BIC.best.model==2))/nrow(T);
    length(which(Tt$BIC.best.model==3))/nrow(T);
    length(which(Tt$BIC.best.model==4))/nrow(T);
    
    cols=c('bisque','olivedrab','palevioletred3','purple4')
    counts = matrix(NA, nrow=3, ncol=3)
    counts[1,1] = length(which(Tt$BIC.best.model==2))
    counts[2,1] = length(which(Tt$BIC.best.model==3))
    counts[3,1] = length(which(Tt$BIC.best.model==4))
    
    counts[1,2] = length(which(Tt$BIC.best.model==2 & Tt$qv.mRNA<0.05))
    counts[2,2] = length(which(Tt$BIC.best.model==3 & Tt$qv.mRNA<0.05))
    counts[3,2] = length(which(Tt$BIC.best.model==4 & Tt$qv.mRNA<0.05))
    
    counts[1,3] = length(which(Tt$BIC.best.model==2 & Tt$qv.mRNA<0.05 & Tt$rel.amp.mRNA>0.1))
    counts[2,3] = length(which(Tt$BIC.best.model==3 & Tt$qv.mRNA<0.05 & Tt$rel.amp.mRNA>0.1))
    counts[3,3] = length(which(Tt$BIC.best.model==4 & Tt$qv.mRNA<0.05 & Tt$rel.amp.mRNA>0.1))
    
    #counts <- table(mtcars$vs, mtcars$gear)
    
    pdf.name = paste("myplots/Total/Partition_models.pdf", sep='')
    pdf(pdf.name, width=3, height=3)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,1.5,0.8)+0.1, tcl = -0.3)
    
    barplot(counts, main=NA, ylab= 'nb of rhythmic mRNAs', col=cols[c(2:4)], legend = c('RS-CD', 'CS-RD', 'RS-RD'))
    
    dev.off()
    
    pdf.name = paste("myplots/Total/Partition_models.pdf", sep='')
    pdf(pdf.name, width=2.5, height=2.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    pie.sales= c(length(which(Tt$BIC.best.model==1))/nrow(T),
    length(which(Tt$BIC.best.model==2))/nrow(T),
    length(which(Tt$BIC.best.model==3))/nrow(T),
    length(which(Tt$BIC.best.model==4))/nrow(T))
    
    names(pie.sales) <- c("49%", "33%","9%", "10%")
    cols=c('gray','green3','tomato','black')
    pie(pie.sales,col=cols)
    #legend(0.65, 1.25, legend = c(' ',' ', ' ', ' ') , cex=0.9, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()
    
    pdf.name = paste("myplots/Total/Pecentages_rhythmic_mRNAs.pdf", sep='')
    pdf(pdf.name, width=2.5, height=2.5)
    par(cex = 0.7, las = 1, mgp = c(1.6,0.5,0), mar = c(1.0,1.0,1.0,3)+0.1, tcl = -0.3)
    pie.sales= c(length(which(Tt$BIC.best.model==2))/length(which(T$BIC.best.model>1)),
    length(which(Tt$BIC.best.model==3))/length(which(T$BIC.best.model>1)),
    length(which(Tt$BIC.best.model==4))/length(which(T$BIC.best.model>1)))
    
    names(pie.sales) <- c("64%","17%", "19%")
    cols=c('green3','tomato','black')
    pie(pie.sales,col=cols)
    #legend(0.65, 1.25, legend = c(' ',' ', ' ', ' ') , cex=0.9, pt.cex=1.5, pt.lwd=1.0, fill= cols, border = NA, bty = 'n')
    #title(main="Distribution of proteins (4016 proteins of 5827) for each models", cex.main=1.8, font.main=2)
    #text(0, -0.98, paste("M.1: Constant mRNA and Constant Protein \n M.2: Rhythmic mRNA and Constant Protein \n M.3: Constant mRNA and Rhythmic Protein \n M.4: Rhythmic mRNA and Rhythmic Protein"))
    dev.off()

    #### Splicing time, Half-life and Rhythmic Degradation (phase, amplitude and delays compared with transcription and mRNAs)
    source('model_RNA_seq_total/functions.R')
    plot.summary.splicing.half.life(folder='myplots/Total/')
    
    #### Compare effects of different conditions
    source('model_RNA_seq_total/functions.R')
    for(m in c(2:4))
    {
        folder = paste('myplots/Total/gene_examples/Conditions_Model_', m, sep='')
        if(!file.exists(folder)){dir.create(folder)}
        
        index = which(T$BIC.best.model==m)
        index = index[which(T$qv.mRNA[index]<0.001)]
        index = index[order(T$qv.mRNA[index])]
        #index = index[c(1:10)]
        #if(m==2){index = index[which(T$qv.premRNA<0.005)]}
        #else{}
        #index = index[order(T$pval.mRNA[index])]
        #index = kk[order(T$BIC.best.model[kk])]
        #index = c(1:10)
        #pdfname = paste('myplots/gene_examples/fitting_results_all_model_', m, '.pdf', sep='')
        plot.genes.examples.conditions(index, T, Figure=FALSE, folder=folder)
    }
    
    source('model_RNA_seq_total/functions.R')
    for(m in c(1:4))
    {
        folder = paste('myplots/Total/gene_examples/Conditions_Model_', m, sep='')
        if(!file.exists(folder)){dir.create(folder)}
        plot.summary.comparison.AL.RF(T, model=m, folder=folder)
    }

}







