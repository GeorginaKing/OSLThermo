%%Fits SAR data to a GOK fit%%
%%Solves for all aliquots at once, Georgina King, August 2015%%
%%georgina.king@unil.ch%%
%%maxime bernard%%
% Last update: 24/01/2025

function out = SarfitGK_GOK(beta0, x);
global rhop10 Hs natDdot

% x input is natural doses (Gy)
t = x./natDdot(1).*(3600*24*365.25*1000); %time in s

rhop=10.^rhop10;
Ma = 1e6*365*24*3600;                       % 1 Ma in seconds
ka = Ma./1000; 
a=beta0(3); D0=beta0(1); alpha=beta0(2);
Fading=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3); % include fading during lab irradiation [Huntley, 2006]
Fgrowth = a.*Fading.*(1-(1+(x./D0)*alpha).^(-1/alpha));
Fgrowth(isnan(Fgrowth))=0; 
out=[Fgrowth];
