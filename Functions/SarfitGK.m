%%Fits SAR data to a single saturating exponential%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@uni-kolen.de%%
%%maxime bernard%%
% Last update: 24/01/2025

function out = SarfitGK(beta0, x)
global rhop10 Hs natDdot
D0=beta0(1); a=beta0(2);

% x input is natural doses (Gy)
t = x./natDdot(1).*(3600*24*365.25*1000); %time in s

Fading=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3); % include fading during lab irradiation [Huntley, 2006]

Growth = a.*Fading.*(1-exp(-x./D0));
Growth(isnan(Growth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Growth];

