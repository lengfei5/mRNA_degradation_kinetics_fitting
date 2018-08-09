##########################################################################
##########################################################################
## Project: Rpackage for kinetic model fitting for mRNA degradation
## Script purpose: the parametrization of kinetic model and calculation of pre-mRNA and mRNA for given parameters
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Jun  4 12:15:17 2018
##########################################################################

######################################
######################################
## Section: main functions to calculating the temproal profiles for pre-mRAN and mRNA 
######################################
######################################
compute.s.beta = function(t = zt, Min = 1, Amp = 2, phase = 12, beta = 1)
{
  w = 2*pi/24;
  s = Min+Amp*((1+cos(w*(t-phase)))/2)^beta
  return(s)
}

compute.s.beta.v1 = function(t = zt, mean = 1, fold.change = 2, phase = 12, beta = 1)
{
  w = 2*pi/24;
  #s = 1/(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)*(1+(fold.change-1)*((1+cos(w*(t-phase)))/2)^beta)
  s = mean/(1+(fold.change-1)/4^beta*gamma(1+2*beta)/gamma(1+beta)^2)*(1+(fold.change-1)*((1+cos(w*(t-phase)))/2)^beta)
  return(s)
}

compute.m.beta = function(t=seq(0, 94, by=2), gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, 
                          splicing.k=log(2)/(5/60), Min = 0.5, Amp=5, phase=12, beta=1, simulation.only=FALSE)
{
  # gamma=log(2)/5; eps.gamma=0.25; phase.gamma=12; splicing.k=log(2)/(5/60); mean = 10; 
  # fold.change=10; phase=6; beta=1;t = seq(0,46, by=2);simulation.only=FALSE
  # gamma=4.159; eps.gamma=0.3; phase.gamma=24; splicing.k=43.135; Min = 0.01418; Amp=0.215; phase=6.3; beta=5;
  # t = seq(0,94, by=2);simulation.only=FALSE
  w = 2*pi/24;
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

######################################
######################################
## Section: low-level functions for the main functions to calculating temporal profiles of pre-mRNA and mRNAs
######################################
######################################
###############################
# compute m by integreation 
###############################
integrate.m = function(t=seq(0, 46, by=2), 
                       gamma=log(2)/5, eps.gamma=0.2, phase.gamma=12, splicing.k=log(2)/(5/60), 
                       Min = 0.5, Amp=2.0, phase=12, beta=1)
{
  w = 2*pi/24;
  m = rep(0, length(t))
  
  Tstable = 24*(ceiling(log(2000)/gamma/24) +1)
  for(i in 1:length(t))
  {
    time = t[i]; 
    m[i] = splicing.k*exp(-Gamma(t = Tstable+time, gamma = gamma, eps.gamma=eps.gamma, phase.gamma=phase.gamma)) *integrate(f2integrate, lower = 0, upper = Tstable+time, par = c(gamma, eps.gamma, phase.gamma, Min, Amp, phase, beta))$value
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
  w = 2*pi/24;
  Gamma = gamma*(t + eps.gamma/w * sin(w*(t-phase.gamma)))
  #Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
  #Gamma = gamma*(t+ eps.gamma/w * sin(w*(t-phase.gamma)))
  #Gamma = gamma*t + amp.gamma/2*(t+ sin(w*(t-phase.gamma))/w)
  return(Gamma)
}

###############################
# compute m by simulation
###############################
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
  
  w = 2*pi/24;
  
  s.t = compute.s.beta(t=t, Min = param.synthesis.1, Amp = param.synthesis.2, phase = param.synthesis.3, beta =  param.synthesis.4)
  
  gamma.t = gamma * ( 1 + eps.gamma*(cos(w*(t-phase.gamma))))
  
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
