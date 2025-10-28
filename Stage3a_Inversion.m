%%%% STAGE 3a, Inversion %%%%
%%%% Inverts data to determine cooling histories
%%%% Requires output of Stage 2a
%%%% Follows method of King et al., 2016, QG
%%%% Script generates a random time-temperature histories
%%%% n/N values are predicted using the thermal histories

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arnaud Duverger, 2015         arnaud.duverger@ens.fr %
% Frédéric Herman, 2016         frederic.herman@unil.ch %
% Georgina King, 2016           georgina.king@unil.ch %
% Chloé Bouscary, 2022          chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.mat']);
nt = length(records);                                                       % number of signals

%%% Import kinetic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [records.params];                                                      % extract kinetic parameters

for k=1:nt
    NNnat(k,:) = [records(k).params.nNnat];
    nNnat = NNnat(:,1); snNnat = NNnat(:,2);
end

%%% Extract signal parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nt
    TypeMeasurement(k) = records(k).typeMeasurement;
    TypeSignal(k) = records(k).typeSignal;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model parameters (TO CHANGE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Model parameters
Tmax = 150;                             % maximum temperature in [°C]
if isnan(records(1).params.natT(1))
    Tsurf = 10; TsurfErr = 5;           % average modern surface or sample temperature in [°C], and its uncertainty
else
    Tsurf = records(1).params.natT(1);  % average modern surface or sample temperature in [°C], and its uncertainty, as given in the input excel file
    TsurfErr = records(1).params.natT(2);
end
tmax = 1;                               % maximum time in [Ma]
tmin = 0;                               % minimum time in [Ma]
nstep = 501;                            % discretisation of the model in time
time = linspace(tmin,tmax,nstep);

%%% Number of iterations
niter = 5000;                           % number of random realisations / Monte Carlo iterations

%%% Redefine uncertainties on nNnat (OPTIONAL)
snNnat = 0.05*nNnat;                                                        % set uncertainties to 5% of nNnat
% snNnat(snNnat./nNnat.*100<5) = 0.05*nNnat;                                  % if snNnat < 5% of nNnat, redefined snNnat to 5% of nNnat
% snNnat(snNnat<0.05) = 0.05;                                                 % if snNnat < 0.05, redefine snNnat to 0.05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Compute vectors for fitting the data with the different models
misOUT = zeros(niter,1);
Time = zeros(niter,nstep);
Temp = zeros(niter,nstep);
nNmod = zeros(niter,nstep,nt);

%% Run the INVERSION in parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through the different Monte Carlo iterations (different Tt paths)
parfor i = 1:niter
    Tmin = Tsurf-(TsurfErr*(-1+2*rand));                                    % Tmin plus or minus uncertainty

    %%% Create a random path
    tT = randpathAD([tmin Tmax],[tmax Tmin]);

    %%% Patch to keep model strictly monotonic, added GK 27.07.2016
    for ii=1:length(tT)-1; if (tT(ii+1)==tT(ii) || tT(ii+1)<tT(ii));  tT(ii+1)=tT(ii)+1e-7; end; end

    %%% Interpolate the random path on a regular grid
    temp = interp1(tT(:,1),tT(:,2),time,'linear');

    %%% Compute the output vectors for the loop
    nNf = zeros(nt,1); v = zeros(nstep,nt); Residuals=zeros(1,nt);

    %%% Loop through the different traps
    for k = 1:nt

        %%% Call forward model for one realisation of each trap
        if SAR_fittype==1 % SSE fit
            if ITH_fittype==1     % BTS
                v(:,k) = trapping_SSE_BTS_FAD(time,temp,kp(k));
            elseif ITH_fittype==2 % GOK
                v(:,k) = trapping_SSE_GOK_FAD(time,temp,kp(k));
            elseif ITH_fittype==3 % GAUSS
                v(:,k) = trapping_SSE_GAUSS_FAD(time,temp,kp(k));
            end

        elseif SAR_fittype==2 % GOK fit
            if ITH_fittype==1     % BTS
                v(:,k) = trapping_GOK_BTS_FAD(time,temp,kp(k));
            elseif ITH_fittype==2 % GOK
                v(:,k) = trapping_GOK_GOK_FAD(time,temp,kp(k));
            elseif ITH_fittype==3 % GAUSS
                v(:,k) = trapping_GOK_GAUSS_FAD(time,temp,kp(k));
            end

        end

        %%% Save the final nNmod value
        nNf(k,1) = v(end,k);

        %%% Misfit between the final nNmod value (here nNf) and nNnat
        Residuals(k) = (1/2*(nNnat(k)/snNnat(k))*(log(nNnat(k)./nNf(k)))).^2; % logarithmic misfit for one signal

    end

    %%% Save model outputs
    nNmod(i,:,:) = v;                                                       % saves modelled nNmod for each thermal history
    misOUT(i) = sum(Residuals)/nt;                                          % saves mistfit
    Time(i,:) = time;                                                       % saves time
    Temp(i,:) = temp;                                                       % saves temperature
end

%%% Sort data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedmisOUT,IX] = sort(misOUT);
sortedTime = Time(IX,:);
sortedTemp = Temp(IX,:);
sortednNmod = nNmod(IX,:,:);


%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tt.misfit           = sortedmisOUT;
Tt.time             = sortedTime;
Tt.temp             = sortedTemp;
Tt.nNmod            = sortednNmod;
Tt.nNnat            = nNnat;
Tt.snNnat           = snNnat;

for k=1:nt
    AgeOSL(k,:) = [records(k).params.AgeOSL];
    maxAgeOSL(k,:) = [records(k).params.maxAgeOSL];
end
Tt.AgeOSL           = AgeOSL;
Tt.maxAgeOSL        = maxAgeOSL;

Tt.TypeMeasurement  = TypeMeasurement;
Tt.TypeSignal       = TypeSignal;

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_Tt.mat'], 'Tt');


%%% Running time
tEnd = toc(tStart);
fprintf('Stage3a_Inversion took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));