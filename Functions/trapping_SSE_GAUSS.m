% This model uses a single saturating exponential function for trapping.
% It assumes a single trap population.
% Thermal decay assumes a gaussian distribution of activation energy,
% i.e., of thermal lifetimes
% Athermal fading is not accoutned for.
%%maximebernard%%
% last update: 27/01/2025

function [nNf] = trapping_SSE_GAUSS(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
Et = kparams.Et(1);
sigmaEt = kparams.sigmaEt(1);

% Define constances
kb = 8.617343e-5; 
magic_ratio = ddot/D0;
nrp = 100;
nstep = length(time);
% dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Create variables for GAUSS model Lambert et al (Submitted)
Ea=[(5/nrp):(5/nrp):5]; nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;

% Define rprime range and tau athermic
rprime = linspace(0.0,2.50,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
pr = ones(size(pr)); %Maxime
npr = sum(pr);

% computes nN for the random Tt path
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = s*exp(-(Ea')./(kb.*((T(i-1)+T(i))/2)));
    dt = (time(i-1)-time(i))*Ma;
    alpha = magic_ratio+inv_tauth;
	nN(:,i) = (nN(:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha));
	nNf(i) = pEa*nN(:,i);
end
nNf = nNf./npEa;
