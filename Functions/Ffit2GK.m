%%Fits fading data to derive g2days after Huntley and Lamothe (2001)%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%

function out1 = Ffit2GK(beta0,t)
I=beta0(1); m=beta0(2);
fit = m*log10(t)+I;
out1=fit;


