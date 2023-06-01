function [nN] = ThermaldecayGOK_GK(t,T,E,b,s)
%% Returns nN for thermal decays starting at nN=1
kB=8.617343e-5; %eV/K
K=@(x) s.*exp(-E./kB./(x+273.15));
nN=(1-(1-b).*s.*K(T).*t).^(1./(1-b)); %plots the isothermal holding data

