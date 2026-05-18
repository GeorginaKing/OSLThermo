%%%% STAGE 4a, InversionExh %%%%
%%%% Inverts data to determine exhumation histories
%%%% Requires output of Stage 2a
%%%% Follows method of Biswas et al., 2018, EPSL
%%%% Script generates a random time-depth history then calculates a thermal
%%%% history based on the exhumation rate and geothermal gradient constraints.
%%%% n/N values are predicted using the thermal history
%%%% computed from the exhumation history.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rabiul Biswas, 2018           biswasrabiul@gmail.com %
% Nadja Stalder, 2020           nadjafranziska.stalder@unil.ch %
% Frédéric Herman, 2020         frederic.herman@unil.ch %
% Georgina King, 2022           georgina.king@unil.ch %
% Chloé Bouscary, 2025          chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.mat']);

nt = length(records);                                                       % number of traps
ntmin=1; ntmax=max(nt);                                                     % change to restrict nt

%%% Import kinetic data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [records.params];                                                      % extract kinetic parameters
kp = kp(ntmin:ntmax);                                                       % restrict kp to nt
ntrap = length(kp);                                                         % define number of traps for analysis

for k = 1:nt
    NNnat(k,:) = [records(k).params.nNnat];
end
NNnat = NNnat(ntmin:ntmax,:);
nNnat = NNnat(:,1);
snNnat = NNnat(:,2);

%%% Extract signal parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:nt
    TypeMeasurement(k) = records(k).typeMeasurement;
    TypeSignal(k) = records(k).typeSignal;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model parameters (TO CHANGE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Physics
kappa0 = 30;                            % thermal diffusivity [km2/Myr]
Lz     = 15;                            % depth [km], use 25 km for highT thermochronometers
G      = 30;                            % initial temperature gradient [°C/km]
if isnan(records(1).params.natT(1))
    Tsurf=10; TsurfErr=5;               % surface temperature [°C] and its uncertainty
else
    Tsurf=records(1).params.natT(1);    % surface temperature [°C] and its uncertainty, as given in the input excel file
    TsurfErr=records(1).params.natT(2);
end
%%% Numerics: discretization, time steps and plotting frequency
nz    = 31;                             % number of nodes
dz    = Lz/(nz-1);                      % node size (km)

%%% Choose the depth time history (i.e. exhumation rate history)
Zmax        = Lz;                       % maximum depth in [km]
Zmin        = 0;                        % surface, i.e. Z=0 km
tmax        = 5;                        % time in [Ma] (today)
tmin        = 0;
nstep       = 6000;                     % timestep, to ensure high resolution over Quaternary timescales
dt          = (tmax-tmin)./(nstep-1);   % calculate timestep, dt (Myr)
tvec        = (tmin:dt:tmax);           % create time vector

%%% Number of iterations
niter       = 10000;                    % number of iterations

%%% Redefine uncertainties on nNnat (OPTIONAL)
snNnat = 0.05*nNnat;                            % set uncertainties to 5% of nNnat
% snNnat(snNnat./nNnat.*100<5) = 0.05*nNnat;     % if snNnat < 5% of nNnat, redefined snNnat to 5% of nNnat
% snNnat(snNnat<0.05) = 0.05;                    % if snNnat < 0.05, redefine snNnat to 0.05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Luminescence model time vector
timemaxTT   = tmax;                                                         % fix luminescence data to end of model, i.e. tmax = present day
nstepTT     = timemaxTT*1000;                                               % ensure high resolution nstep
TT_length   = 0.5;                                                          % time in [Ma] over which luminescence is modelled
DtimeTT     = timemaxTT/(nstepTT-1);                                        % discretize
timeM       = timemaxTT-TT_length:DtimeTT:timemaxTT;                        % create high resolution time vector [Ma]

%%% Compute vectors for fitting the data in the loop
G_end  = zeros(niter,1);                                                    % matrix for saving the final geothermal gradient
Misfitage = 0;                                                              % define for parfor loop

%% Run the INVERSION in parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through the different iterations
parfor i=1:niter
    Zvec=[];
    Tb_top = Tsurf-(TsurfErr*(-1+2*rand));                                  % surface temperature plus or minus uncertainty
    Tb_bot = Tb_top+G*Lz;                                                   % maximum temperature at maximum depth (Lz), dependent on Geothermal gradient (G)

    tZ = randpathAD([1 0],[0 1]);                                           % generate a random path of depth and time (values 0-1)
    Z=tZ(:,1)*(Zmax-Zmin)+Zmin; Z=[Zmax; Z; Zmin];                          % reinterpolate onto actual depths (km)
    time=tZ(:,2)*tmax; time=[0; time; tmax];                                % reinterpolate onto actual times (Ma)

    for ii=1:length(time)-1; if (time(ii+1)==time(ii) || time(ii+1)<time(ii));  time(ii+1)=time(ii)+1e-7; end; end % patch to make sure time is strictly monotonic

    Zvec=interp1(time,Z,tvec,'linear');                                     % reinterpolate depth vector to length = tvec
    vz=-diff(Zvec)./diff(tvec); vz=max(vz,1e-8);                            % calculate exhumation velocity, truncate to remove values <1e-8 km/Myr

    %%% Visualise modelling
    %figure(1);
    %subplot(1,2,1); plot(time,Z,'-'); ylabel('Depth (km)'); xlabel('Time (Ma)'); axis square;
    %subplot(1,2,2); plot(tvec(1:end-1),vz); ylabel('Exhumation velocity (km/Ma)'); xlabel('Time (Ma)'); axis square;

    %%% Solve the 1D heat equation to create a thermal history
    Te    = zeros(1,nz);                                                    % temperatures
    zc    = zeros(1,nz);                                                    % depths
    qz    = zeros(1,nz-1);                                                  % flux with depth
    dTedt = zeros(1,nz);                                                    % discretize Temperature with time

    for iz=1:nz                                                             % initial conditions for temperature field
        zc(iz) = (iz-1)*dz;                                                 % depths
        Te(iz) = Tb_bot+(Tb_top-Tb_bot)*zc(iz)/Lz;                          % initial temperatures
    end

    depth = Lz-Zvec(1);                                                     % start at surface
    tmodel=0; it=0; Tpath=[]; tpath=[]; zpath=[]; ZPath=[]; TPath=[]; t_path=[];

    while (tmodel<tmax)                                                     % while time is less than the maximum time
        it=it+1;

        %%% Choose the velocity
        index = floor(tmodel/dt)+1;                                         % tmodel increases with ndt, but our exhumation rate paths are defined in 1 ka windows
        Vzm   = vz(index);                                                  % vz = exhumation rate of corresponding time tmodel

        %%% Compute diffusion
        qz   =  diff(Te(1:end),1,2)/dz;                                     % compute the flux
        dTedt(2:end-1) = kappa0*(diff(qz)/dz);                              % update result

        %%% Compute advection
        dTedt(2:end) = dTedt(2:end) - Vzm*diff(Te)/dz;

        %%% Check time step
        dtn    = min(1e-4,min(dz^2/kappa0/2.1,dz/100/Vzm));                 % timestep depends on exhumation rate
        Te = Te + dtn*dTedt;                                                % update temperature
        tmodel = tmodel+dtn;                                                % update time

        %%% Set boundary conditions
        Te(nz)=Tb_top;                                                      % surface boundary condition
        Te(1)=Te(2)+G*dz;                                                   % flux lower boundary condition
        %Te(1)=Tb_bot;                                                      % fixed lower boundary condition

        %%% Track rock and save thermal history
        %%% Tracking using second order runge-kutta to solve diff. equation
        indexV = floor((tmodel+dtn/2)/dt)+1;                                % time of model/timestep of path (= 1 ka)
        if index < nstep-1
            Vzeff = vz(indexV);                                             % choose velocity corresponding to timestep of model
            depth = depth+Vzeff*dtn;                                        % depth starts at 0

            if (Lz-depth) < Zvec(indexV+1)                                  % condition if model overshoots because timestep of Zvec>for depth
                depth = Lz-Zvec(indexV+1);
            end
        else
            depth = Lz;                                                     % rock needs to be at surface
        end

        %%% Interpolate temperature
        index2=floor(depth/dz)+1;                                           % depth of model/change in depth(= 1 km)
        if (index2 > 0)
            if (index2 < nz)
                u=(depth-zc(index2))/(zc(index2+1)-zc(index2));             % calculate change in depth
                TPath(it)=(1-u)*Te(index2)+u*Te(index2+1);                  % calculate change in temperature
            else
                TPath(it)=Te(end);
            end
        else
            TPath(it)=Te(1);
        end
        t_path(it)=tmodel;
        ZPath(it)=Lz-depth;
    end

    zpath=ZPath; Tpath=TPath; tpath=t_path; Tpath(end)=Tb_top;

    %%% Allocate matrix for parfor loop
    G_end_sum = zeros(1,nz);                                                % matrix for saving the geothermal gradient

    G_end(i) = (Te(nz-1)-Te(nz))./dz;                                       % calculate final Geothermal gradient
    for kk=nz:-1:2; G_end_sum(kk) = (Te(kk-1)-Te(kk))./dz; end              % output geothermal gradient

    %%% Visualise modelling
    %figure(2);
    %subplot(1,3,1); plot(t_path,TPath,'-'); ylabel('Temperature (°C)'); xlabel('Time (Ma)'); axis square;
    %subplot(1,3,2); plot(ZPath,TPath); ylabel('Temperature (°C)'); xlabel('Depth (km)'); axis square;
    %subplot(1,3,3); plot([2:nz],G_end_sum(i,2:end),'o'); ylabel('Geothermal gradient (°C/km)'); xlabel('Time (Ma)'); axis square;

    %%% Compute nN data
    tempM=interp1(tpath,Tpath,timeM,'linear');                              % interpolate Tt path to a higher resolution to ensure that it is stable in trapping model
    Tpath_out = interp1(tpath,Tpath,tvec,'linear');

    %%% Confirming interpolation of data
    %figure(3);
    %subplot(1,2,1); plot(timeM,tempM);
    %subplot(1,2,2); plot(tpath,Tpath); plot(tvec,Tpath_out,'ro');

    %%% Allocate matrixes for parfor loop
    Residuals   = zeros(1,length(nNnat));
    v           = zeros(ntrap,TT_length*1000);
    v2          = zeros(ntrap,1);
    nNf         = [];

    for k = 1:ntrap

        if SAR_fittype==1 % SSE fit
            if ITH_fittype==1     % BTS
                nNf = trapping_SSE_BTS_FAD(timeM,tempM,kp(k));
            elseif ITH_fittype==2 % GOK
                nNf = trapping_SSE_GOK_FAD(timeM,tempM,kp(k));
            elseif ITH_fittype==3 % GDE
                nNf = trapping_SSE_GDE_FAD(timeM,tempM,kp(k));
            end


        elseif SAR_fittype==2 % GOK fit
            if ITH_fittype==1     % BTS
                nNf = trapping_GOK_BTS_FAD(timeM,tempM,kp(k));
            elseif ITH_fittype==2 % GOK
                nNf = trapping_GOK_GOK_FAD(timeM,tempM,kp(k));
            elseif ITH_fittype==3 % GDE
                nNf = trapping_GOK_GDE_FAD(timeM,tempM,kp(k));
            end

        end

        Residuals(k) = (1/2*(nNnat(k)/snNnat(k))*(log(nNnat(k)./nNf(end)))).^2;
        v(k,:) = nNf; v2(k) = nNf(end);

    end

    residuals = Residuals;
    Misfit = sum(Residuals)/ntrap;

    Out_nNmax(i,:) = v2;                                                    % newline for parfor loop
    Out_misfit(i,1) = Misfit;                                               % output misfit values
    Out_ZPath(i,:) = Zvec;                                                  % output depth vector
    Out_tempM(i,:) = tempM;                                                 % output temperature path
    Out_G(i,:) = G_end_sum;                                                 % output fial geothermal gradient

    fprintf('Path%i   Misfit=%f \n ',i,Misfit);

    Misfit=0.;Misfitage=0.;residuals=0.;

end

%%% Sort data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('sorting table')
[sortedmisOUT,IX]    = sort(Out_misfit(:,1));
sortedZt             = Out_ZPath(IX,:);
sortedMaxnN          = Out_nNmax(IX,:,:);
sortedTLumi          = Out_tempM(IX,:);
sortedGend           = G_end(IX,:);
sortedGendsum        = Out_G(IX,:);

%%% Apply rejection algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob = exp(-sortedmisOUT);                                                  % prob = likelihood
scale = max(prob); prob = prob/scale;                                       % normalisation by the maximum of the likelihood, max(prob)

R = rand(length(sortedmisOUT),1); test = prob>R;                            % rejection criteria = for every position of prob: 1 if true (prob>R), 0 if wrong (prob<=R)

idefix = find(test);                                                        % gives the name of the position for which test is true (=1)
idefix = idefix(end:-1:1);                                                  % invert the data for the line before, to order them in decreasing order (now worse fit to best fit)
movea = length(idefix);                                                     % length of the selected data (nb of data that passed the rejection criteria)

%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter data for rejection algorithm
zt.misfit           = sortedmisOUT(idefix);
zt.temp             = sortedZt(idefix,:);
zt.time             = tvec;
zt.timeLumi         = timeM;
zt.TempLumi         = sortedTLumi(idefix,:);
zt.nNnat            = nNnat;
zt.snNnat           = snNnat;
zt.maxnNsort        = sortedMaxnN(idefix,:);
zt.G_end            = sortedGend(idefix);
zt.G_end_sum        = sortedGendsum(idefix,:);
zt.Lz               = Lz;
zt.dz               = dz;
zt.TypeMeasurement  = TypeMeasurement;
zt.TypeSignal       = TypeSignal;

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_zt.mat'], 'zt');


%%% Running time
tEnd = toc(tStart);
fprintf('Stage4a_InversionExh took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));