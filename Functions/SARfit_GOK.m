%% Fits SAR data to a GOK model %%
%%% Solves for all aliquots at once

% Georgina King [georgina.king@unil.ch], August 2015


function out = SARfit_GOK(beta0,t)

global rhop10

rhop=10.^rhop10;
D0=beta0(1); 
alpha=beta0(2);
a=beta0(3);                                                                 % pre-exponential factors

Growth = a.*(1-(1+(t./D0)*alpha).^(-1/alpha));

Growth(isnan(Growth))=0;                                                    % removes NaN data

out = Growth;