%% Fits SAR data to a GOK model %%
%%% Solves for all aliquots at once
%%% Here, athermal detrapping during lab irradiation is assumed
%%% Kinetic parameters are extracted from the fit
%%% These are the parameters for the un-faded curve

% Georgina King [georgina.king@unil.ch], August 2015
% Maxime Bernard [maxime.bernard@unil.ch], 22/01/2025


function out = FSARfit_GOK(beta0,t)

global rhop10 Hs 
global Ddot index_aliquotSAR

rhop=10.^rhop10;
D0=beta0(1); 
alpha=beta0(2);
a=beta0(3:end);                                                             % pre-exponential factors

if length(t) < 100                                                          % assumes t is time read from L measurements, otherwise t is modvec or natdose
    idx = index_aliquotSAR;
    if length(a)>1
        A_fac = a(idx);
        Dr = Ddot(idx)';                                                    % dose rate in [Gy/s]
    else
        idx(1:end) = a;
        A_fac = idx';
        Dr(1:length(idx)) = Ddot';                                          % dose rate in [Gy/s]
        Dr = Dr';
    end
else 
    A_fac = mean(a);
    Dr = mean(Ddot);                                                        % dose rate in [Gy/s]
end

%%% Anomalous fading
kars = exp(-rhop.*log(1.8*Hs.*(0.5.*t)).^3);                                % fading correction factor (Kars et al., 2008), includes fading during lab irradiation (Huntley, 2006). 0.5*t to have irradiation times strating from the mid-irradiation point

%%% Compute the growth matrix
Fgrowth = A_fac.*kars.*(1-(1+(Dr./D0.*t).*alpha).^(-1/alpha));              % accumulation of luminescence signal, corrected for fading

Fgrowth(isnan(Fgrowth))=0;                                                  % removes NaN data

out = Fgrowth;