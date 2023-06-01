 function [nN,pdf] = ThermaldecayGauss(t,T,Et,sigma,s)

 kB = 8.617343e-5; % eV/K
 dEa=0.001; Ea=[0:dEa:10]; 
 distr=exp(-0.5*((Ea-Et)./sigma).^2)/(sigma*sqrt(2*pi));
 isoT=T(1,:);
 
 for j=1:length(isoT);
 time=t(:,j);
 for ix=1:length(time);
 integrand= distr.*exp(-time(ix)*s*exp(-Ea./kB/(isoT(j)+273.15)));
 decay(ix) = sum(integrand)*dEa;
 end
 decay_sum(:,j)=decay;
 end
 nN = decay_sum;