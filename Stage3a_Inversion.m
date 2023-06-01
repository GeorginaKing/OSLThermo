%%%% STAGE 3a, Inversion %%%%
%%%% Inverts data to determine cooling histories
%%%% Requires output of Stage 2a
%%%% Follows method of King et al., 2016, QG 
%%%% Script generates a random time, temperature history
%%%% n/N values are predicted using the thermal history

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arnaud Duverger, 2015 arnaud.duverger@ens.fr %
% Frederic Herman, 2016 frederic.herman@unil.ch %
% Georgina King, 2016 georgina.king@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % records time of execution
clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR;  
close all; clc;

load(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.mat']);

nt = length(records);   % number of traps
niter = 2000;           % number of random realisations 

%% Import Kinetic data
kp = [records.params];  % extract kinetic parameters
NNnat = reshape([records.nNnat],2,nt)';
nNnat = NNnat(:,1); snNnat = NNnat(:,2);

%% Model parameters (TO CHANGE)
Tmax = 200;                 % Maximum temp in C
Tsurf = 10;  TsurfErr = 5;  % Average modern surface or sample temp in C, and error
tmax = 1;                   % Max time in Ma
tmin = 0;                   % Min time in Ma
nstep = 500;                % discretisation of the model in time
time = linspace(tmin,tmax,nstep);

%% Define uncertainties (OPTIONAL)
% pererr = [snNnat-nNnat<=0.05];
% snNnat(pererr) = 0.05;
snNnat = 0.05*nNnat;

misOUT = zeros(niter,1);
Time = zeros(niter,nstep);
Temp = zeros(niter,nstep);
nNmod = zeros(niter,nstep,nt);

parfor i = 1:niter 
    Tmin = Tsurf-(TsurfErr*rand);
	
    % create a random path
	tT = randpathAD([tmin Tmax],[tmax Tmin]); 
    
    % patch to keep model strictly monotonic, added GK 27.07.2016
    for ii=1:length(tT)-1; if (tT(ii+1)==tT(ii) || tT(ii+1)<tT(ii));  tT(ii+1)=tT(ii)+1e-7; end; end;
    
    % interpolates the random path on a regular grid
	temp = interp1(tT(:,1),tT(:,2),time,'linear');  
      
	nNf = zeros(nt,1); v = zeros(nstep,nt); Residuals=zeros(1,nt)
	
    %loop through different traps
    for k = 1:nt 
		
        if SAR_model==1
            if ITL_model==1
                v(:,k) = trapping_EXP_GOK_FAD(time,temp,kp(k));
            elseif ITL_model==2
                v(:,k) = trapping_GAUSS_FAD(time,temp,kp(k));        
            elseif ITL_model==3
                v(:,k) = trapping_BTM_FAD(time,temp,kp(k));
            end
        end
    
        if SAR_model==2
            if ITL_model==1
                v(:,k) = trapping_GOK_FAD(time,temp,kp(k));
            elseif ITL_model==2
                v(:,k) = trapping_GOK_GAUSS_FAD(time,temp,kp(k));
            elseif ITL_model==3
                v(:,k) = trapping_GOK_BTS_FAD(time,temp,kp(k));
            end
        end
        nNf(k,1) = v(end,k);
        Residuals(k) = (1/2*(nNnat(k)/snNnat(k))*(log(nNnat(k)./nNf(k)))).^2;
                
    end
    
    nNmod(i,:,:)=v;
    
	misOUT(i) = sum(Residuals)/nt; % save mistfit	
	Time(i,:) = time; % saves time 
	Temp(i,:) = temp; % save temp
end

% Sort data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sortedmisOUT,IX] = sort(misOUT);
sortedTime = Time(IX,:);
sortedTemp = Temp(IX,:);
sortednNmod = nNmod(IX,:,:);

% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tt.misfit = sortedmisOUT;
Tt.time = sortedTime;
Tt.temp = sortedTemp;
Tt.nNmod = sortednNmod;
Tt.nNnat = nNnat;
Tt.snNnat = snNnat;

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_Tt.mat'], 'Tt'); 
toc
