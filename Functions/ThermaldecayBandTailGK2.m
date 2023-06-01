function [nN,pdf] = ThermaldecayBandTailGK2(t,T,Et,Eu,s)
%% Returns nN for thermal decays starting at nN=1
%
% Inputs Time and Temp have identical size 
% (typically 2D grids that match the experiment)
%
% Input pdf_type selects a probabiliy density distribution type:
% 'lin' 'sqr' 'exp' or 'gauss'
% Add more types and corresponding parameters beyond Eu as needed.
%
% renske.lambert@unil.ch
% August 2015

distr = @(Eb) exp(-Eb./Eu); % Eb/Eu is scale factor

normdistr = quad(distr, 0, Et); % Normalize distribution
pdf = @(Eb) distr(Eb)./normdistr;

kB = 8.617343e-5; % eV/K
K = @(x,Eb) s*exp(-(Et-Eb)./kB./(x + 273.15));

decay = @(Eb) pdf(Eb).*exp(-t.*K(T,Eb));
decay0 = @(Eb) pdf(Eb).*exp(-0.*K(T,Eb));
nN = quadv(decay, 0, Et)./quadv(decay0, 0, Et);     

