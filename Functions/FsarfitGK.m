%%Fits SAR data to a single saturating exponential%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%

function out = FsarfitGK(beta0, x)
global rhop10
rhop=10.^rhop10;
a=beta0(1); D0=beta0(2);
kars=exp(-rhop*log(1.8*3e15.*(0.5*x)).^3);
Fgrowth = a.*kars.*(1-exp(-x./D0));
Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Fgrowth];

