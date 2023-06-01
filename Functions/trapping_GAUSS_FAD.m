function [nNf,nNtest] = trapping_GAUSS_FAD(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
rhop = 10.^kparams.rhop10(1);
rhop(rhop<0)=0; %added for negative rhop GK 23.11.2016
Et = kparams.Et(1);
sigmaEt = kparams.sigmaEt(1);

% Define constances
kb = 8.617343e-5; Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
dEa=0.01; 
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Create variables for GAUSS model Lambert et al (Submitted)
Ea=[(5/nrp):(5/nrp):5]; nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;

% Define rprime range and tau athermic
rprime = linspace(0.0,2.50,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
npr = sum(pr);

inv_tauath = ones(nEa,1)*(Hs*exp(-(rhop.^-(1./3)).*rprime)); %combine eq 1 and 3 from Kars et al 2008, convert to Ma

% computes nN for the random Tt path
nN = zeros(nEa,nrp,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = s*exp(-(Ea')./(kb.*T(i-1)))*ones(1,nEa);

    alpha = magic_ratio+inv_tauth+inv_tauath;
	nN(:,:,i) = (nN(:,:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha));
	nNf(i) = pEa*nN(:,:,i)*pr;
end
nNf = nNf./npEa./npr;
