%% Model luminescence trapping with GOK and GDE %%
%%% Assumes a General order kinetic (GOK) for trapping. GOK assumes
%%% trapping slower than SSE, due to e.g., Coulomb repulsive forces. Also,
%%% the effective trap lifetime varies with time (Guralnik et al., 2015).
%%% Thermal detrapping assumes a gaussian distribution of activation energies,
%%% i.e., of thermal lifetimes (Lambert, 2018).
%%% Athermal fading is accounted for, with s_ath = 3e15 (Huntley, 2006).

% Georgina King [georgina.king@unil.ch], August 2015


function [nNf] = trapping_GOK_GDE_FAD(time,temp,kparams)

%%% Extract parameters
ddot = kparams.natDdot(1);                                                  % environmental radiation dose rate [Gy/s]
D0 = kparams.D0(1);                                                         % fading corrected characteristic dose of saturation [Gy]
s = 10.^kparams.s10(1);                                                     % thermal frequency factor [s-1]
rhop = 10.^kparams.rhop10(1); rhop(rhop<0)=0;                               % density of recombination centers [-], added for negative rhop
Et = kparams.Et(1);                                                         % activation energy = trap depth below the conduction band [eV]
sigmaEt = kparams.sigmaEt(1);
a = kparams.GOK_a(1);

%%% Define constants
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]
Hs = 3e15;                                                                  % Hs=s_ath value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
Ma = 1e6*365.25*24*3600;                                                    % 1 [Myr] in [s]
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;                                  % time step in [Myr]

%%% Create variables for detrapping GDE model Lambert et al. (Submitted)
Ea = (5/nrp):(5/nrp):5;
nEa = length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;                                                            % transforms temperatures from [°C] to [K]

%%% Define rprime range to calculate tau athermal
rprime = linspace(0.0,2.50,nrp);                                            % creates vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.);                                      % calculates p(r'), eq.3 in Kars et al. (2008), Rad. Meas.
npr = sum(pr);

%%% Compute inv_tauath (nEa x nrp)
inv_tauath = ones(nEa,1)*(Hs*exp(-(rhop.^-(1./3)).*rprime));                % inverse of athermal lifetime [s-1], combines eq.1, eq.2, and eq.4 from Kars et al. (2008), Rad. Meas.

%%% Precompute nN for the random Tt path
nN = zeros(nEa,nrp,nstep);
nNf = zeros(1,nstep);

%%% Loop through the time steps
for i=2:nstep

    %%% Compute inv_tauth for the current time step (nEa x nrp)
    inv_tauth = s*exp(-(Ea')./(kb.*T(i-1)))*ones(1,nrp);                    % inverse of thermal lifetime [s-1], GAUSS model

    %%% Compute xkd and xk for all combinations of Ea and rprime for the current time step (nEa x nrp)
    xkd = -a*magic_ratio*(1-nN(:,:,i-1)).^(a-1)-inv_tauth-inv_tauath;
    xk = magic_ratio*(1-nN(:,:,i-1)).^a-inv_tauth.*(nN(:,:,i-1))-inv_tauath.*(nN(:,:,i-1));

    %%% Update nN for this time step using vectorized operations (nEa x nrp)
    nN(:,:,i) = nN(:,:,i-1)+dt*xk./(1-dt*xkd);

    %%% Compute nNf for this time step using matrix operations ([1 x nEa] * [nEa x nrp] * [npr x 1] > 1 value of nNf)
    nNf(i) = pEa*nN(:,:,i)*pr;
end

%%% Normalize nNf by the sum of the probabilities of Ea and rprime
nNf = nNf./npEa./npr;