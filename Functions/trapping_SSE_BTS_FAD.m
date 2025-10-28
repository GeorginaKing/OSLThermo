%% Model luminescence trapping with SSE and BTS %%
%%% Assumes a single trap population (Single saturating exponential) for trapping.
%%% Thermal detrapping assumes a band-tail states model (Li and Li, 2013).
%%% Athermal fading is accounted for, with s_ath = 3e15 (Huntley, 2006).

% Georgina King [georgina.king@unil.ch], August 2015


function [nNf] = trapping_SSE_BTS_FAD(time,temp,kparams)

%%% Extract parameters
ddot = kparams.natDdot(1);                                                  % environmental radiation dose rate [Gy/s]
D0 = kparams.D0(1);                                                         % fading corrected characteristic dose of saturation [Gy]
s = 10.^kparams.s10(1);                                                     % thermal frequency factor [s-1]
rhop = 10.^kparams.rhop10(1); rhop(rhop<0)=0;                               % density of recombination centers [-], added for negative rhop
Et = kparams.Et(1);                                                         % activation energy = trap depth below the conduction band [eV]
Eu = kparams.Eu(1);                                                         % width of the Urbach tail [eV]

%%% Define constants
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]
Hs = 3e15;                                                                  % Hs=s_ath value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
nEb = 100;
Ma = 1e6*365.25*24*3600;                                                    % 1 [Myr] in [s]
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;                                  % time step in [Myr]

%%% Define rprime range to calculate tau athermal
rprime = linspace(0.0,2.50,nrp);                                            % creates vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3);                                       % calculates p(r'), eq.3 in Kars et al. (2008), Rad. Meas.
npr = sum(pr);

%%% Compute inv_tauath (nEb x nrp)
inv_tauath = ones(nEb,1)*(Hs*exp(-(rhop.^-(1./3)).*rprime));                % inverse of athermal lifetime [s-1], combines eq.1, eq.2, and eq.4 from Kars et al. (2008), Rad. Meas.

%%% Create variables for detrapping Li and Li (2013), BTS
Eb = linspace(0,Et,nEb)';                                                   % distribution of band-tail states energy levels [eV]
aLi = 1;                                                                    % pre-exponential multiplier, = rho(0), the density of the sub-conduction band states (Eb=Ec=0), Li and Li (2013), J. of Lumi
pEb = aLi*exp(-Eb'/Eu);                                                     % Li and Li (2013), J. of Lumi.
npEb = sum(pEb);
T = temp+273.15;                                                            % transforms temperatures from [°C] to [K]

%%% Precompute nN for the random Tt path
nN = zeros(nEb,nrp,nstep);
nNf = zeros(1,nstep);

%%% Loop through the time steps
for i = 2:nstep

    %%% Compute inv_tauth for the current time step (nEb x nrp)
    inv_tauth = (s*exp(-(Et-Eb)./(kb.*T(i-1))))*ones(1,nrp);                % inverse of thermal lifetime [s-1], BTS model, eq.2 of Li and Li (2013), J. of Lumi.

    %%% Compute alpha for all combinations of rprime and Eb for the current time step (nEb x nrp)
    alpha = magic_ratio+inv_tauth+inv_tauath;

    %%% Update nN for this time step using vectorized operations (nEb x nrp)
    nN(:,:,i) = (nN(:,:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha));

    %%% Compute nNf for this time step using matrix operations ([1 x nEb] * [nEb x nrp] * [npr x 1] > 1 value of nNf)
    nNf(i) = pEb*nN(:,:,i)*pr;
end

%%% Normalize nNf by the sum of the probabilities of Eb and rprime
nNf = nNf./npEb./npr;