%%Fits SAR data to a GOK model%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%
%%maxime Bernard%%
% last update: 22/01/2025
% Here, we assume athermal detrapping during lab irradiation
% Kinetic parameters extracted from fit here are the parameters for the
% un-faded curve

function out = FsarfitGK_GOK(beta0, x);
global rhop10;
global labDdot
global index_aliquot
global Hs % huntley's escape frequency factor

rhop=10.^rhop10;
D0=beta0(1); alpha=beta0(2);
c = beta0(3:end); % pre-exponential factors

out=[];
t = x;
if length(x) <100 % assumes x is time read from L measurements else time is modvec
    idx = index_aliquot;
    A_fac = c(idx);
    Ddot = labDdot(idx);
else 
    idx = 1;
    A_fac = mean(c(:));
    Ddot = mean(labDdot);
end
Fading=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3); % include fading during lab irradiation [Huntley, 2006]
Fgrowth = A_fac.*Fading.*(1-(1+(Ddot./D0).*alpha.*t).^(-1/alpha));
Fgrowth(isnan(Fgrowth))=0; 
out=[out; Fgrowth];



