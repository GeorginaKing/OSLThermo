%% Returns nN for thermal decays using a BTS distribution %%
%%% Inputs Time (t) and Temp (T) have identical size
%%% (typically 2D grids that match the experiment)
%%% Input pdf_type selects a probabiliy density distribution type:
%%% 'lin' 'sqr' 'exp' or 'gauss'
%%% Add more types and corresponding parameters beyond Eu as needed

% Renske Lambert [renske.lambert@unil.ch], August 2015


function nN = ThermaldecayBTS(t,T,Et,Eu,s)

%%% Define constant
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]

%%% Probability distribution of Eb
distr =@(Eb) exp(-Eb./Eu);                                                  % Eb/Eu is scale factor
normdistr = quad(distr,0,Et);                                               % normalize distribution
pdf =@(Eb) distr(Eb)./normdistr;

%%% Calculate the detrapping rate (s_tun) = inverse of the lifetime
K =@(x,Eb) s*exp(-(Et-Eb)./kb./(x+273.15));

%%% Compute the isothermal holding data
decay =@(Eb) pdf(Eb).*exp(-t.*K(T,Eb));
decay0 =@(Eb) pdf(Eb).*exp(-0.*K(T,Eb));

nN = quadv(decay,0,Et)./quadv(decay0,0,Et);