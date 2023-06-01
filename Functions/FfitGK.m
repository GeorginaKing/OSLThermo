%%Fits fading data to derive rho' after Huntley (2006)%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%

function out1 = FfitGK(beta0,t)
c=beta0(1); rhop=10.^beta0(2);
fade = c.*exp(-rhop.*(log(1.8.*3e15.*t)).^3);
out1=fade;

