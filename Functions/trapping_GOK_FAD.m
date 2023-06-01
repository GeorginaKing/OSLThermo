 function [nNf] = trapping_GOK_FAD(time,temp,kparams);
 
Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
rhop = 10.^kparams.rhop10(1);
E = kparams.Et(1);
b = kparams.GOK_b(1);
a = kparams.GOK_a(1);

% Define constants
kb = 8.617343e-5 ;
Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Define rprime range and tau athermic
rprime = linspace(0.0,2.50,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
npr = sum(pr);

inv_tauath = 1*(Hs*exp(-(rhop.^-(1./3)).*rprime)); %combine eq 1 and 3 from Kars et al 2008, convert to Ma
T = temp+273.15;

% computes nN for the random Tt path (Euler integration method)
nN=zeros(nstep,nrp);
    for j=2:nstep
            inv_tauth=s*exp(-E/(kb*(temp(j-1)+273.15))); %convert to Ma
            xkd=-a*magic_ratio*(1-nN(j-1,:)).^(a-1)-b*inv_tauth.*(nN(j-1,:)).^(b-1)-inv_tauath(:)';          
            xk=magic_ratio*(1-nN(j-1,:)).^a-inv_tauth.*(nN(j-1,:)).^b-inv_tauath(:)'.*nN(j-1,:);
            nN(j,:) = nN(j-1,:)+dt*xk./(1-dt*xkd);
            nNf(j) = sum(pr'.*nN(j,:)); 
    end

 
   nNf = nNf./npr ;



