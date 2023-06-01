%%Fits SAR data to a GOK model%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%

function out = FsarfitGK_GOK(beta0, x);
global rhop10;
rhop=10.^rhop10;
a=beta0(1); D0=beta0(2); alpha=beta0(3);
alpha(alpha>=2)=2;
kars=exp(-rhop*log(1.8*3e15.*(0.5*x)).^3);
Fgrowth = a.*kars.*(1-(1+(x./D0)*alpha).^(-1/alpha));
Fgrowth(isnan(Fgrowth))=0; 
out=[Fgrowth];



