%%Fits SAR data to a single saturating exponential%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@uni-kolen.de%%

function out = SarfitGK(beta0, x)
a=beta0(1); D0=beta0(2);
Growth = a.*(1-exp(-x./D0));
Growth(isnan(Growth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Growth];

