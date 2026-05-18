%%%% OSLThermo %%%%
%%%% Loader of scripts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arnaud Duverger, 2016     arnaud.duverger@ens.fr %
% Georgina King, 2018       georgina.king@unil.ch %
% Chloé Bouscary, 2024      chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

tStartall=tic;                                                              % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


addpath('./Functions/');
addpath('./Data/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SAMPLE PARAMETERS ==> TO CHANGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% List of file names %%%
filenamevec = {'KRG104'};

%%% Number of measurements %%%
nSARvec     =   [3;3];          % number of aliquots measured for SAR
nITHvec     =   [7;7];          % number of isothermal decay temperatures
nFADvec     =   [3;3];          % number of aliquots measured for athermal fading

%%% Select models %%%
SAR_fittype     =   2;          % 1=SSE; 2=GOK
ITH_fittype     =   3;          % 1=BTS; 2=GOK; 3=GDE

if SAR_fittype==1; SAR_MODEL='SSE'; elseif SAR_fittype==2; SAR_MODEL='GOK'; end
if ITH_fittype==1; ITH_MODEL='BTS'; elseif ITH_fittype==2; ITH_MODEL='GOK'; elseif ITH_fittype==3; ITH_MODEL='GDE'; end

if ~ismember(SAR_fittype, [1 2])
    disp('ERROR: Invalid SAR_fittype');
    return
end
if ~ismember(ITH_fittype, [1 2 3])
    disp('ERROR: Invalid ITH_fittype');
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNNING THE MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f=1:size(filenamevec,1)                                                 % loops through the different filenames

    filename = cell2mat(filenamevec(f,:));
    nSAR = nSARvec(f);
    nITH = nITHvec(f);
    nFAD = nFADvec(f);
    
    %%% Convert the Excel file to a .mat format %%%
    run Stage1_ExcelToStruct 

    %%% Fit the data using the selected model %%%
    run Stage2a_Fitparameters

    %%% Plot the results of the fits %%%
    run Stage2b_PlotFit

    %%% Invert the data following King et al. (2016, QG) %%%
    run Stage3a_Inversion

    %%% Plot the results of the inversion %%%
    run Stage3b_PlotTt

    %%% Invert the data following Biswas et al. (2018, EPSL) %%%
    % run Stage4a_InversionExh

    %%% Plot the results of the inversion for exhumation %%%
    % run Stage4b_PlotExh

end


%%% Running time
tEnd = toc(tStartall);
fprintf('OSLThermo ran for a total duration of %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));