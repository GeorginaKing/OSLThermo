%%%% STAGE 4a, InversionExh %%%%
%%%% Inverts data to determine exhumation histories
%%%% Requires output of Stage 2a
%%%% Follows method of Biswas et al., 2018, EPSL
%%%% Script generates a random time-depth history then calculates a thermal
%%%% history based on the exhumation rate and geothermal gradient constraints.
%%%% n/N values are predicted using the thermal history computed from the exhumation history.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rabiul Biswas, 2018, biswasrabiul@gmail.com %
% Nadja Stalder, 2020, nadjafranziska.stalder@unil.ch %
% Frederic Herman, 2020, frederic.herman@unil.ch %
% Georgina King, 2022, georgina.king@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR; 
close all; clc;

load(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.mat']);

nt = length(records);           % number of traps
ntmin=1; ntmax=max(nt);         % change to restrict nt

%% Import Kinetic data
kp = [records.params];          % extract kinetic parameters
kp=kp(ntmin:ntmax);             % restrict kp to nt
ntrap=length(kp);               % define nt for analysis
NNnat = reshape([records.nNnat],2,nt)'; 
NNnat=NNnat(ntmin:ntmax,:); 
nN = NNnat(:,1); 
sigmanNnat = NNnat(:,1).*0.05; % define unceratinty on n/Nnat values

%% physics
kappa0 = 30.;                   % thermal diffusivity (W/m.K)
Lz     = 15.;                   % depth (km)
G      = 30.;                   % Initial temperature gradient

%% numerics: discretization, time steps and plotting frequency
nz    = 31;                     % number of nodes
dz    = Lz/(nz-1);              % node size (km)

%% choose the depth time history (i.e. exhumation rate history)
Zmax        =Lz;                        % Z in km
Zmin        =0;                         % Surface, i.e. Z=0 km
tmax        =5;                         % time im Ma
tmin        =0;             
nstep       =6000;                      % timestep, to ensure high resolution over Quaternary timescales      
dt          =(tmax-tmin)./(nstep-1);    % calculate timestep, dt
tvec        =(tmin:dt:tmax);            % create time vector
niter       =50000;                     % number of iterations

%% Trapped-charge model time vector 
timemaxTT   =tmax;                                      % fix trapped-charge data to end of model, i.e. 5 Ma = present day
nstepTT     =timemaxTT*1000;                            % ensure high resolution nstep
TT_length   =0.5;                                       % time in Ma over which TL modelled
DtimeTT      =timemaxTT/(nstepTT-1);                    % discretize
timeM        =timemaxTT-TT_length:DtimeTT:timemaxTT;    % create high resolution time vector (Ma) 
    
G_end  = zeros(niter,1);                % Matrix for saving the final geothermal gradient
Misfitage=0;                            % Define for parfor loop


%% Main Program
parfor i=1:niter
    Zvec=[]; 
    Tb_top = 15-(5*rand);                           % Surface temperature
    Tb_bot = Tb_top+G*Lz;                           % Maximum temperature at maximum depth (Lz), dependent on Geothermal gradient (G)
    
    tZ = randpath([1 0],[0 1]);                     % Generate a random path of depth and time (values 0-1)
    Z=tZ(:,1)*(Zmax-Zmin)+Zmin; Z=[Zmax; Z; Zmin];  % Reinterpolate onto actual depths (km)
    time=tZ(:,2)*tmax; time=[0; time; tmax];        % Reinterpolate onto actual times (Ma)  
    
    for ii=1:length(time)-1; if (time(ii+1)==time(ii) || time(ii+1)<time(ii));  time(ii+1)=time(ii)+1e-7; end;end; % patch to make sure time is strictly monotonic 
    
    Zvec=interp1(time,Z,tvec,'linear');             % Reinterpolate depth vector to length = tvec 
    vz=-diff(Zvec)./diff(tvec); vz=max(vz,1e-8);    % Calculate exhumation velocity, truncate to remove values <1e-8 km/Ma
    
    % visualise modelling
%     figure(1); subplot(1,2,1); plot(time,Z,'-'); ylabel('Depth [km]'); xlabel('Time [Ma]'); axis square; subplot(1,2,2); plot(tvec(1:end-1),vz); ylabel('Exhumation velocity [km/Ma]'); xlabel('Time [Ma]'); axis square; 

    %% solve the 1D heat equation to create a thermal history
    Te    = zeros(1,nz);                            % Temperatures
    zc    = zeros(1,nz);                            % depths
    qz    = zeros(1,nz-1);                          % flux with depth
    dTedt = zeros(1,nz);                            % discretize Temperature with time
    kappa = kappa0;                                 % Thermal diffusivity (W/Ma)
    
    for iz=1:nz                                     % initial conditions for temperature field
        zc(iz) = (iz-1)*dz;                         % depths
        Te(iz) = Tb_bot+(Tb_top-Tb_bot)*zc(iz)/Lz;  % initial temperatures
    end
   
    depth = Lz-Zvec(1);                             % start at surface
    tmodel=0; it=0; Tpath=[]; tpath=[]; zpath=[]; ZPath=[]; TPath=[]; t_path=[];                  
    
    while (tmodel<tmax)                             % while time is less than the maximum time
    it=it+1;
    
    % choose the velocity
    index = floor(tmodel/dt)+1;                     % tmodel increases with ndt, but our exhumation rate paths are defined in 1 ka windows
    Vzm   = vz(index);                              % vz = exhumation rate of corresponding time tmodel
   
    % compute diffusion
    qz   =  diff(Te(1:end),1,2)/dz;                 % compute the flux
    dTedt(2:end-1) = kappa*(diff(qz)/dz);           % update result
      
    % compute advection
    dTedt(2:end) = dTedt(2:end) - Vzm*diff(Te)/dz;
      
    % check time step
    dtn    = min(1e-4,min(dz^2/kappa0/2.1,dz/100/Vzm)); % timestep depends on exhumation rate
    Te = Te + dtn*dTedt;                            % update temperature
    tmodel = tmodel+dtn;                            % update time 
    
    % set boundary conditions
    Te(nz)=Tb_top;                                  % surface boundary condition
    Te(1)=Te(2)+G*dz;                               % flux lower boundary condition
    % Te(1)=Tb_bot;                                 % fixed lower boundary condition
    
    % track rock and save thermal history
    % tracking using second order runge-kutta to solve diff. equation
    
    indexV = floor((tmodel+dtn/2)/dt)+1;            % time of model/timestep of path(= 1 ka)
    if index < nstep-1                              
        Vzeff   = vz(indexV);                       % choose velocity corresponding to timestep of model
        depth = depth+Vzeff*dtn;                    % depth starts at 0
        
        if (Lz-depth)<Zvec(indexV+1);               % condition if model overshoots because timestep of Zvec>for depth
             depth=Lz-Zvec(indexV+1);
        end
    else
        depth = Lz;                                 % rock needs to be at surface
    end
    
    % interpolate temperature
    index2=floor(depth/dz)+1;                       % depth of model/change in depth(= 1 km)
    if (index2 > 0)
        if (index2 < nz)
          u    = (depth-zc(index2))/(zc(index2+1)-zc(index2));  % calculate change in depth
          TPath(it)=(1-u)*Te(index2)+u*Te(index2+1);            % calculate change in temperature
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
    
    % Allocate matrix for parfor loop
    G_end_sum = zeros(1,nz);                                    % Matrix for saving the geothermal gradient

    G_end(i) = (Te(nz-1) - Te(nz))./dz;                         % calculate final Geothermal gradient
    for kk=nz:-1:2; G_end_sum(kk) = (Te(kk-1)-Te(kk))./dz; end  % output geothermal gradient
    
% visualise modelling
%     figure(2); subplot(1,3,1); plot(t_path,TPath,'-'); ylabel('Temperature [C]'); xlabel('Time [Ma]'); axis square; 
%     subplot(1,3,2); plot(ZPath,TPath); ylabel('Temperature [C]'); xlabel('Depth [km]'); axis square;
%     subplot(1,3,3); plot([2:nz],G_end_sum(i,2:end),'o'); ylabel('Geothermal gradient'); xlabel('Time [Ma]'); axis square;
    
    
%% compute nN data
    tempM=interp1(tpath,Tpath,timeM,'linear');                  % interpolate Tt path to a higher resoultion to ensure that it is stable in trapping_BTM_FAD
    Tpath_out = interp1(tpath,Tpath,tvec,'linear');
    
%confirm interpolation of data%
%     figure(10); subplot(1,2,1);
%     plot(timeM,tempM);
%     subplot(1,2,2); plot(tpath,Tpath); plot(tvec,Tpath_out,'ro')

% Allocate matrixes for parfor loop
Residuals   = zeros(1,length(nN)); 
v           = zeros(ntrap,TT_length*1000); 
v2          = zeros(ntrap,1); 
nNf         = [];

for k=1:ntrap;
    
    nNf = trapping_BTM_FAD(timeM,tempM,kp(k));
    Residuals(k)=(1/2*(nN(k)/sigmanNnat(k))*(log(nN(k)./nNf(end)))).^2;
    v(k,:)=nNf; v2(k)=nNf(end); 
    
end
  
    residuals=Residuals; 
    Misfit=sum(residuals)/ntrap; 
    
    Out_nNmax(i,:)=v2;      % Output nN values
    Out_misfit(i,1)=Misfit; % Output misfit values
    Out_ZPath(i,:)=Zvec;    % Output depth vector
    Out_tempM(i,:)=tempM;   % Output Temperature Path (Lumi)
    Out_G(i,:) = G_end_sum; % Output final geothermal gradient

    fprintf('Path%i   Misfit=%f \n ',i,Misfit);

    Misfit=0.;Misfitage=0.;resdiduals=0.;
    
end

% Sort data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('sorting table')
[sortedmisOUT,IX]    = sort(Out_misfit(:,1));
sortedZt             = Out_ZPath(IX,:);
sortedMaxnN          = Out_nNmax(IX,:,:);
sortedTLumi          = Out_tempM(IX,:);
sortedGend           = G_end(IX,:);
sortedGendsum        = Out_G(IX,:);

% Apply rejection algorithm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = rand(length(sortedmisOUT),1); prob=exp(-sortedmisOUT); 
scale=max(prob); 
prob=prob/scale; test=prob>R;
idefix=find(test); idefix=idefix(end:-1:1); movea=length(idefix);

% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter data for rejection algorithm
zt.misfit   = sortedmisOUT(idefix);
zt.temp     = sortedZt(idefix,:);
zt.time     = tvec;
zt.timeLumi = timeM;
zt.TempLumi = sortedTLumi(idefix,:);
zt.nNnat    = nN;
zt.snNnat   = sigmanNnat;
zt.maxnNs   = sortedMaxnN(idefix,:);
zt.G_end    = sortedGend(idefix);
zt.G_end_sum= sortedGendsum(idefix,:);
zt.Lz       = Lz;
zt.dz       = dz;

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_zt.mat'], 'zt'); 

toc