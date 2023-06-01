function [nNf,nNtest] = trapping_BTM_FAD(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
rhop = 10.^kparams.rhop10(1); rhop(rhop<0)=0; %added for negative rhop GK 23.11.2016
Et = kparams.Et(1);
Eu = kparams.Eu(1);

% Define constances
kb = 8.617343e-5; Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
nEb = 100;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Define rprime range and tau athermic
rprime = linspace(0.0,2.50,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
npr = sum(pr);

inv_tauath = ones(nEb,1)*(Hs*exp(-(rhop.^-(1./3)).*rprime)); %combine eq 1 and 3 from Kars et al 2008, convert to Ma

% Create variables for lili2013 BTM
%Eb = linspace(0.01,Et+0.01,nEb)';
Eb = linspace(0,Et,nEb)';
pEb = exp(-Eb'/Eu);
npEb = sum(pEb);
T = temp+273.15;

% computes nN for the random Tt path
nN = zeros(nEb,nrp,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = s*exp(-(Et-Eb)./(kb.*T(i-1)))*ones(1,nrp);
    alpha = magic_ratio+inv_tauth+inv_tauath;
	nN(:,:,i) = (nN(:,:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha));
	nNf(i) = pEb*nN(:,:,i)*pr;
end
nNf = nNf./npEb./npr;
