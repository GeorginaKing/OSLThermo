%%Fits fading data to derive rho' after Huntley (2006)%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%
%%maxime bernard%%
% last update: 22/01/2025
% Here we assume no trapping and no thermal detrapping during room
% temperature isothermal holding

function out1 = FfitGK(beta0,t)
global index_aliquot Hs

rhop=10.^beta0(1); c=beta0(2:end)';

if length(t) <100 % assumes x is time read from L measurements else time is Fvec
    idx = index_aliquot;
else 
    idx = 1;
end

fade = c(idx).*exp(-rhop.*(log(1.8.*Hs.*t)).^3);
out1=fade;

