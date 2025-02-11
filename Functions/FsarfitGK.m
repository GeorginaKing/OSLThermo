%%Fits SAR data to a single saturating exponential%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%
%%maxime bernard%%
% Last update: 24/01/2025

function out = FsarfitGK(beta0, x)
global rhop10 index_aliquot Hs labDdot

rhop=10.^rhop10;
D0=beta0(1); a=beta0(2:end);

out=[];
t = x;
if length(x) <100 % assumes x is time read from L measurements else time is modvec
    idx = index_aliquot;
    A_fac = a(idx);
else 
    idx = 1;
    A_fac = mean(a(:));
end

kars=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3);
Fgrowth = A_fac.*kars.*(1-exp(-labDdot(idx)./D0.*t));
Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[out;Fgrowth];

