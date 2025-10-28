%% Returns nN for thermal decays using a GAUSS distribution %%
%%% Inputs Time (t) and Temp (T) have identical size
%%% (typically 2D grids that match the experiment)

% Renske Lambert [renske.lambert@unil.ch], August 2015


function nN = ThermaldecayGAUSS(t,T,muEt,sigmaEt,s)

%%% Define constant
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]
isoT=T(1,:);

%%% Probability distribution of Ea
dEa=0.001;
Ea=0:dEa:9;
distr = exp(-0.5*((Ea-muEt)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));             % Eq. 14, p. 2-24 in Lambert, 2018

%%% Calculate the detrapping rate (s_tun) = inverse of the lifetime
K =@(x,Ea) s*exp(-Ea./kb./(x+273.15));

%%% Compute the isothermal holding data
for j=1:length(isoT)
    time=t(:,j);
    for ix=1:length(time)
        integrand = distr.*exp(-time(ix)*K(isoT(j),Ea));
        decay(ix) = sum(integrand)*dEa;
    end
    decay_sum(:,j) = decay;
end

nN = decay_sum;