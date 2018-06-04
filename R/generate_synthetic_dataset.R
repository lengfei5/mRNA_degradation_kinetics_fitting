#############################################################################################################################
######################################################  Generate simulated read counts
#############################################################################################################################
generate.fake.read.counts = function(T = T[1:200,], Tt=Tt, X = 200, model = 1, zt = seq(0,94,by = 2), test.noise.identifiability=FALSE,
                                     parametrization = 'cosine.beta', absolute.signal=TRUE)
{
  # model = 4; T = F[((model-1)*X+1):(model*X),]; absolute.signal=TRUE; parametrization = 'cosine.beta'; zt = seq(0, 94, by=2);test.noise.identifiability=TRUE;
  set.seed(42)
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  
  if(absolute.signal)
  {
    ###### This is part to simulate RNA-seq absolute levels 
    T$gene = paste('fake_m',model,'_',1:X,sep = '')
    
    ### Define the limits of parameters: beta, eps.gamma, gamma
    bounds = set.bounds(model = 4, parametrization = parametrization, absolute.signal = TRUE); 
    lower = bounds$lower;	upper = bounds$upper;
    
    ### parameter of premRNAs
    Min.int.min = lower[5]; Min.int.max = upper[5];
    Amp.int.min = 0.01; Amp.int.max = upper[6];
    phase.int.min = lower[7]; phase.int.max = upper[7];
    beta.int.min = lower[8]; beta.int.max = upper[8];
    
    ### parameters of mrna degradation
    gamma.min = lower[1];
    #gamma.max = log(2)/(10/60);
    gamma.max = upper[1];   
    eps.gamma.min = 0.05; eps.gamma.max = upper[2];
    
    phase.gamma.min = lower[3]; phase.gamma.max = upper[3];
    splicing.k.min = lower[4]; splicing.k.max = upper[4];
    
    ### Shared parameters for all models
    T$Min.int = pmin(pmax(2^rnorm(X, mean = -2.5, sd = 2.0), Min.int.min), Min.int.max)
    T$gamma = pmin(pmax(log(2)/exp(rnorm(X, mean = log(3), sd = 1.)), gamma.min), gamma.max)
    T$splicing.k = pmin(pmax(T$gamma*T$mean.rpkm.mRNA/T$mean.rpkm.premRNA, splicing.k.min), splicing.k.max)
    #T$splicing.k = pmin(pmax(log(2)/rnorm(X, mean = 5/60, sd = 2/60), splicing.k.min), splicing.k.max)
    
    #T$Min.int = sample(lseq(Amp.int.min, Amp.int.max, length = 100), X, replace = TRUE);
    #T$splicing.k = sample(lseq(splicing.k.min, splicing.k.max, length=100), X, replace = TRUE)
    #T$gamma = sample(lseq(gamma.min, gamma.max, length = 100), X, replace=TRUE);
    ### Parameters specific for certain models
    if((model == 2)|(model == 4)) ## randomly generate parameters for rhythmic pre-mrna
    {
      #kk = which(Tt$qv.rpkm.premRNA<0.05); hist(Tt$rel.amp.rpkm.premRNA[kk])
      rel.amp.int= pmax(pmin(sample(Tt$rel.amp.rpkm.premRNA[which(Tt$qv.rpkm.premRNA<0.05)], X, replace = TRUE), 0.95), 0.05)
      #rel.amp.int = sample(seq(0.01, 0.9, length=100), X, replace=TRUE);
      fc.int = (1+rel.amp.int)/(1-rel.amp.int);
      T$Amp.int = T$Min.int * (fc.int-1) ## assign premRNA amplitudes according to relative amplitudes
      T$phase.int =  sample(seq(phase.int.min, phase.int.max,by = 0.1), X, replace = TRUE)
      #T$beta.int = pmax(pmin(rnorm(n = X, mean = 1, sd = 1), beta.int.max), beta.int.min)
      T$beta.int = pmax(pmin(rnorm(n = X, mean = 1, sd = 2), beta.int.max), beta.int.min)
    }else{
      ### parameters of pre-mRNA for model 1 and model 3
      T$Amp.int = rep(0,X); 
      T$phase.int = rep(0,X); 
      T$beta.int = rep(1,X); 
    }
    if(model>2) ## parameters for rhythmic degradation
    {
      rel.amp.gamma= pmax(pmin(rnorm(n = X, mean = 0.2, sd = 0.2), 0.8), 0.05)
      #rel.amp.gamma = sample(seq(0.05, 0.8, length=100), X, replace=TRUE);
      #fc.gamma = (1+rel.amp.int)/(1-rel.amp.int);
      #T$eps.gamma = T$gamma*fc.gamma;
      T$eps.gamma = rel.amp.gamma;
      if(model==3){T$phase.gamma = sample(seq(phase.gamma.min, phase.gamma.max,by = 0.1),X, replace = TRUE);}
      if(model==4){
        #phase.diff = sample(seq(0, 24, by=0.1), X, replace=TRUE)
        phase.diff = rnorm(n=X, mean=12, sd=6)%%24;
        T$phase.gamma = (T$phase.int + phase.diff)%%24;
      }
    }else{
      T$eps.gamma = rep(0,X)
      T$phase.gamma = rep(0,X)
    }
    abs.int = T[, ZT.int]*0
    abs.ex = T[, ZT.ex]*0
    for(i in 1:X)
    {
      ### convert rpkm into mean of read numbers
      L.m = T$length.mRNA[i]; L.s = T$length.premRNA[i];
      #alpha.m = T$alpha.mRNA[i];
      #alpha.s = T$alpha.premRNA[i];
      alpha.m = rep(as.numeric(T[i, grep('alpha.mRNA.ZT', colnames(T))]), 4); alpha.s = rep(as.numeric(T[i, grep('alpha.premRNA.ZT', colnames(T))]), 4);
      
      ## absolute values of premRNAs and mRNAs
      s = compute.s.beta(t = zt, Min = T$Min.int[i], Amp = T$Amp.int[i], phase = T$phase.int[i], beta = T$beta.int[i]);
      if(model == 1){
        m = T$splicing.k[i] * T$Min.int[i]/T$gamma[i];
      }else{
        m = compute.m.beta(t = zt, 
                           gamma = T$gamma[i], T$eps.gamma[i], phase.gamma = T$phase.gamma[i], splicing.k = T$splicing.k[i], 
                           T$Min.int[i], T$Amp.int[i], T$phase.int[i], T$beta.int[i])
      }
      
      mu.m = convert.nb.reads(m, L.m);
      mu.s = convert.nb.reads(s, L.s);
      
      R.m = c()
      R.s = c()
      if(!test.noise.identifiability)
      {
        for(jj in 1:48)
        {
          R.m = c(R.m, rnbinom(n=1, 1/alpha.m[jj], mu=mu.m[jj]))
          R.s = c(R.s, rnbinom(n=1, 1/alpha.s[jj], mu=mu.s[jj]))
        }
      }else{
        alpha.m = 10^-6;
        alpha.s = 10^-6;
        T[, grep('alpha.mRNA.ZT', colnames(T))] = alpha.m;
        T[, grep('alpha.premRNA.ZT', colnames(T))] = alpha.s;
        for(jj in 1:48)
        {
          R.m = c(R.m, rnbinom(n=1, 1/0.0001, mu=mu.m[jj]))
          R.s = c(R.s, rnbinom(n=1, 1/0.0001, mu=mu.s[jj]))
        }
      }
      abs.ex[i,] = R.m;
      abs.int[i,] = R.s;
    }
    T[,ZT.ex] =  abs.ex
    T[,ZT.int] = abs.int
    #T$mean.ex = apply(noise.abs.ex,1,mean)
  }
  return(T)
}

generate.fake.read.counts.4identifiability.test = function(T = T, alpha = 10^-5, nb.realization = 100, model = 3, zt=seq(0,94,by = 2))
{
  # model = 3; T = T; alpha=10^-5; gene.template = 'Per1'; nb.realization = 100; zt = seq(0, 94, by=2); Ts = T;
  set.seed(42)
  ZT.int = grep('.count.premRNA', colnames(T))
  ZT.ex = grep('.count.mRNA', colnames(T))
  index = which(T$gene=='Per1')
  
  ### define grids of parameter space
  gamma.grids = log(2)/lseq(10/60, 24, length.out = 20)
  eps.gamma.grids = seq(0.05, 0.6, by=0.05)
  Ts  = matrix(NA, nrow = length(gamma.grids)*length(eps.gamma.grids)*nb.realization, ncol=123)
  colnames(Ts) = colnames(T)[1:123]
  Ts = data.frame(Ts)
  
  ### define gene.name, length, noise of premRNA, mRNA and all parameters but gamma and eps.gamma with gene template
  Ts$length.mRNA = T$length.mRNA[index];
  Ts$length.premRNA = T$length.premRNA[index];
  L.m = T$length.mRNA[index];
  L.s = T$length.premRNA[index];
  
  params = c(0.2, 1.0, 13, 1, 0.4, 0.5, 20, 10)
  Ts$Min.int = params[1]; Ts$Amp.int = params[2]; Ts$phase.int = params[3]; Ts$beta.int = params[4];
  Ts$gamma = params[5]; Ts$eps.gamma = params[6]; Ts$phase.gamma = params[7]; Ts$k.gamma.ratio = params[8];
  
  if(model==2){params[6] = 0; params[7]=0; Ts$eps.gamma = 0; Ts$phase.gamma = 0; Ts$gamma = NA;}
  if(model==3){params[2] = 0; params[3] = 0; Ts$Amp.int = 0; Ts$phase.int = 0;Ts$gamma = NA;Ts$eps.gamma=NA;}
  
  Ts[, grep('alpha', colnames(Ts))] = alpha;
  
  ### gene name and parameters and also simulate read counts for premRNA and mRNAs
  nb = 1
  s = compute.s.beta(t = zt, Min = params[1], Amp = params[2], phase = params[3], beta = params[4]); ## always the same
  mu.s = convert.nb.reads(s, L.s);
  for(m1 in 1:length(gamma.grids))
  {
    for(m2 in 1:length(eps.gamma.grids))
    {
      ## m1 = 1; m2 = 1;
      kk = ((nb-1)*nb.realization + c(1:nb.realization))
      
      Ts$gene[kk] = paste('fake_m',model,'_',nb, '_r', c(1:nb.realization),  sep = '')
      Ts$gamma[kk] = gamma.grids[m1];
      Ts$eps.gamma[kk] = eps.gamma.grids[m2];
      
      ### compute the mean of read counts with given gamma and eps.gamma
      if(model == 1){
        m =  params[8]*params[1];
      }else{
        m = compute.m.beta(t=zt, gamma.grids[m1], eps.gamma.grids[m2], params[7], splicing.k=params[8]*gamma.grids[m1], 
                           params[1], params[2], params[3], params[4])
      }
      mu.m = convert.nb.reads(m, L.m);
      
      for(jj in 1:48)
      {
        R.m = rnbinom(n=nb.realization, 1/alpha, mu=mu.m[jj])
        R.s = rnbinom(n=nb.realization, 1/alpha, mu=mu.s[jj])
        Ts[kk, ZT.ex[jj]] = R.m;
        Ts[kk, ZT.int[jj]] = R.s;
      }
      
      nb = nb + 1;
    }
  }
  
  #test = c(1:300) + 1800
  #cbind(Ts$gene[test], Ts$gamma[test], Ts$eps.gamma[test])
  
  return(Ts)
}
