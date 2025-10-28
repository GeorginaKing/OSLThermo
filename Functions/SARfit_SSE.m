%% Fits SAR data to a single saturating exponential, SSE model %%
%%% Solves for all aliquots at once

% Georgina King [georgina.king@unil.ch], August 2015


function out = SARfit_SSE(beta0,t)

D0=beta0(1); 
a=beta0(2);

Growth = a.*(1-exp(-t./D0));

Growth(isnan(Growth))=0;                                                    % removes NaN data

out = Growth;