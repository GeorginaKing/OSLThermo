% This model assumes a log-normal distribution of D0 and a gaussian
% distribution of activation energy, i.e., of thermal lifetimes
% Athermal fading is accounted for and s_ath = 3e15 (Huntley, 2006)
%%maximebernard%%
% last update: 27/01/2025

function [nNf,nNtest] = trapping_GAUSS_GAUSS_FAD(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %
D0 = log10(kparams.D0(1)); % Input D0 is in lin scale
sigmaD0 = kparams.sigmaD0(1)/(log(10)*10^D0); % Input sigmaD0 is in lin scale
s = 10.^kparams.s10(1);
rhop = 10.^kparams.rhop10(1);
rhop(rhop<0)=0; %added for negative rhop GK 23.11.2016
Et = kparams.Et(1);
sigmaEt = kparams.sigmaEt(1);

% Define constances
nrp = 100; nmr = 30;
kb = 8.617343e-5; Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
nstep = length(time);

% Create variables for detrapping GAUSS model Lambert et al (Submitted)
Ea= linspace(5 / nrp, 5, nrp); nEa=length(Ea);
pEa = exp(-0.5 * ((Ea-Et) ./ sigmaEt).^2) / (sigmaEt * sqrt(2 * pi));
npEa = sum(pEa);
T = temp+273.15;

% Define rprime range and tau athermic
rprime = linspace(0.0,2.50,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
npr = sum(pr);

% Precompute D0 probability distribution
xx = linspace(D0 - 3 * sigmaD0, D0 + 3 * sigmaD0, 30)';
D0_p = normpdf(xx, D0, sigmaD0);
D0_p = D0_p ./ sum(D0_p);
nD0 = length(xx);

nN = zeros(nEa, nrp, nD0, nstep);
nNf = zeros(1, nstep);

% Precompute magic ratio and athermic tau
magic_ratio = reshape(ddot ./ 10.^xx, [1, 1, nD0]); % Reshape to (1 x 1 x nD0) for broadcasting
inv_tauath = Hs * exp(-(rhop.^-(1/3)) .* rprime);    % (nrp x 1)

% Loop over time steps
for i = 2:nstep
    dt = (time(i-1) - time(i)) * Ma;

    % Compute inv_tauth for the current time step (nEa x 1)
    inv_tauth = s * exp(-Ea' ./ (kb * T(i-1))); % Vectorized computation for all Ea

    % Compute alpha for all combinations (nEa x nrp x nD0)
    alpha = magic_ratio + ...                              % (1 x 1 x nD0)
            reshape(inv_tauth, [nEa, 1, 1]) + ...          % (nEa x 1 x 1)
            reshape(inv_tauath, [1, nrp, 1]);             % (1 x nrp x 1)

    % Update nN for this time step using vectorized operations
    nN(:,:,:,i) = (nN(:,:,:,i-1) + magic_ratio .* dt) ./ (1 + dt .* alpha);

    % Compute nNf for this time step using matrix operations
    % Combine dimensions: (nEa x nrp x nD0) to scalar
    nNf(i) = sum(sum(sum((reshape(pEa, [nEa, 1, 1]) .* nN(:,:,:,i)) .* ...
                        reshape(D0_p, [1, 1, nD0]) .* ...
                        reshape(pr, [1, nrp, 1]))));
end

% Normalize nNf
nNf = nNf ./ npEa ./ npr;