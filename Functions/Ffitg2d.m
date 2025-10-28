%% Fits fading data to derive g2days after Huntley and Lamothe (2001) %%
%%% Solves for all aliquots at once

% Georgina King [georgina.king@unil.ch], August 2015

function out = Ffitg2d(beta0,t)

I=beta0(1); m=beta0(2);

fit = m*log10(t)+I;

out = fit;