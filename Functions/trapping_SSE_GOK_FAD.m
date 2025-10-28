%% Model luminescence trapping with SSE and GOK %%
%%% Assumes a single trap population (Single saturating exponential) for trapping.
%%% Thermal detrapping assumes a general order kinetic decay. A slower than
%%% exponential thermal detrapping is explained by electron retrapping or
%%% distant-dependent probabilities of detrapping.
%%% Athermal fading is accounted for, with s_ath = 3e15 (Huntley, 2006).

% Georgina King [georgina.king@unil.ch], August 2015


function [nNf] = trapping_SSE_GOK_FAD(time,temp,kparams)

%%% Extract parameters
ddot = kparams.natDdot(1);                                                  % environmental radiation dose rate [Gy/s]
D0 = kparams.D0(1);                                                         % fading corrected characteristic dose of saturation [Gy]
s = 10.^kparams.s10(1);                                                     % thermal frequency factor [s-1]
rhop = 10.^kparams.rhop10(1); rhop(rhop<0)=0;                               % density of recombination centers [-], added for negative rhop
E = kparams.Et(1);                                                          % activation energy = trap depth below the conduction band [eV]
b = kparams.GOK_b(1);
a = 1;                                                                      % SSE is like a first order kinetic model

%%% Define constants
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]
Hs = 3e15;                                                                  % Hs=s_ath value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
Ma = 1e6*365.25*24*3600;                                                    % 1 [Myr] in [s]
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;                                  % time step in [Myr]

%%% Define rprime range to calculate tau athermal
rprime = linspace(0.0,2.50,nrp);                                            % creates vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.);                                      % calculates p(r'), eq.3 in Kars et al. (2008), Rad. Meas.
npr = sum(pr);

%%% Compute inv_tauath (1 x npr)
inv_tauath = (Hs*exp(-(rhop.^-(1./3)).*rprime));                            % inverse of athermal lifetime [s-1], combines eq.1, eq.2, and eq.4 from Kars et al. (2008), Rad. Meas.

%%% Precompute nN for the random Tt path (Euler integration method)
nN = zeros(nstep,nrp);
nNf = zeros(1,nstep);
T = temp+273.15;                                                            % transforms temperatures from [°C] to [K]

%%% Loop through the time steps
for i=2:nstep

    %%% Compute inv_tauth for the current time step (1 x 1)
    inv_tauth0 = s*exp(-E/(kb*T(i-1)));                                     % inverse of thermal lifetime [s-1], GOK model

    %%% Compute xkd and xk for all combinations of rprime for the current time step (1 x nrp)
    xkd = -a*magic_ratio*(1-nN(i-1,:)).^(a-1)-b*inv_tauth0.*(nN(i-1,:)).^(b-1)-inv_tauath(:)';
    xk = magic_ratio*(1-nN(i-1,:)).^a-inv_tauth0.*(nN(i-1,:)).^b-inv_tauath(:)'.*nN(i-1,:);

    %%% Update nN for this time step using vectorized operations (1 x nrp)
    nN(i,:) = nN(i-1,:)+dt*xk./(1-dt*xkd);

    %%% Compute nNf accounting for the probabilities of rprime (1 x nrp)
    nNf(i) = sum(nN(i,:).*pr');
end

%%% Normalize nNf by the sum of the probabilities of rprime
nNf = nNf./npr ;