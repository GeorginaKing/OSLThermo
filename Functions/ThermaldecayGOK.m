%% Returns nN for thermal decays using a GOK distribution %%
%%% Inputs Time (t) and Temp (T) have identical size 
%%% (typically 2D grids that match the experiment)

% Georgina King [georgina.king@unil.ch], 2015
% Chloé Bouscary [chloebouscary@gmail.com], 2025


function nN = ThermaldecayGOK(t,T,E,b,s)

%%% Define constant
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]

%%% Calculate the detrapping rate (s_tun) = inverse of the lifetime
K =@(x) s.*exp(-E./kb./(x+273.15));

%%% Compute the isothermal holding data
nN = (1-(1-b).*K(T).*t).^(1./(1-b));       