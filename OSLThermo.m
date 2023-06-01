%%%% OSLThermo 
%%%% Loader of scripts 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arnaud Duverger, 2016, arnaud.duverger@ens.fr %
% Georgina King, 2018, georgina.king@unil.ch %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic                                             
clear all; close all; clc;
addpath('./Functions/'); 
addpath('./Data/');

%define global parameters for the other scripts
global SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SAMPLE PARAMETERS ==> TO CHANGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% List of file names %%%
filenamevec={'KRG104'};
 
%%% Number of measurements %%%
NITLvec = [7];  %number of isothermal decay temperatures
nSAR = [3];     %number of aliquots measured for SAR

%%%Select models%%%
SAR_model=[1]; %1 = SSE; %2 = GOK
ITL_model=[3]; %1 = GOK; 2 = GAUSS; 3 = BTS

if SAR_model==1; SAR_MODEL='SSE'; elseif SAR_model==2; SAR_MODEL='GOK'; end;
if ITL_model==1; ITL_MODEL='GOK'; elseif ITL_model==2; ITL_MODEL='Gauss'; elseif ITL_model==3; ITL_MODEL='BTS'; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RUNNING THE MODEL %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:size(filenamevec,1);    % Loops through the different filenames

    filename=cell2mat(filenamevec(j,:));
    NITL=NITLvec(j);
    
    %%% Convert the excel file to a .mat format %%
     run Stage1_ExcelToStruct 
 
     %%% Fit the data using the selected model %%%
%     run Stage2a_Fitparameters
 
%     %%% Plot the results of the fits %%%
%          run Stage2b_PlotFit
 
    %%% Invert the data following King et al. (2016, QG) %%%
%      run Stage3a_Inversion

    %% Plot the results of the inversion %%%
%     run Stage3b_PlotTt
    
    %%% Invert the data following Biswas et al. (2018, EPSL) %%%
%     run Stage4a_InversionExh

    %%% Plot the results of the inversion %%%
%      run Stage4b_PlotExh
    
    
end
toc