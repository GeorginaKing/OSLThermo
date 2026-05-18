%%%% STAGE 1, ExcelToStructure %%%%
%%%% Converts .xlsx files to .MAT data for input to subsequent scripts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benny Guralnik, 2014      benny.guralnik@gmail.com %
% Arnaud Duverger, 2016     arnaud.duverger@ens.fr %
% Chloé Bouscary, 2024      chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Define constant
ka = 1e3*365.25*24*3600;                                                    % [kyr] in [s]

%%% Extract data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[exist,sheets] = xlsfinfo([filename '.xlsx']);                              % sheets = Excel sheet names (e.g. Al-center, Ti-center), exist = file type (Excel spreadsheet)

for i=1:length(sheets)                                                      % loops through the different Excel sheets

    %%% Select the data on the Excel sheet
    [~,txt,rawExcel] = xlsread([filename '.xlsx'],i);                       % read data
    rawdata=rawExcel;
    rawdata(~cellfun(@isfloat,rawdata)) = {NaN};                            % creates cells for the raw data to be stored within
    rawArray=cell2mat(rawdata);                                             % transforms rawdata(cells) in rawArray(array double) to remove the NaN values at the end of the lines and columns

    while all(isnan(rawArray(:,end))) == true                               % selects the columns that have only NaN at the end of the file
        rawArray(:,end) = [];                                               % removes all these columns
    end

    while all(isnan(rawArray(end,:))) == true                               % selects the rows that have only NaN at the end of the file
        rawArray(end,:)=[];                                                 % removes all these rows
    end
    rawdata=num2cell(rawArray);                                             % transforms back rawArrayt(array double) to rawdata(cells)


    %%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    records(i).id = sheets{i};                                              % name of the measured signal, as indicated by the name of the Excel sheets (e.g. IRSL50, IRSL100 etc.)

    if isnan(rawExcel{1,7})
        records(i).typeMeasurement = {'LUMI'};                              % default value for type of measured signal, if data is missing in Excel
    else
        records(i).typeMeasurement = {txt{1,7}};                            % type of measured signal (IRSL, OSL etc.), as inputted on the Excel file
    end

    if isnan(rawdata{1,8})
        records(i).typeSignal = 9999;                                       % default value for center type, if data is missing in Excel
    else
        records(i).typeSignal = [rawdata{1,8}];                             % temperature of the signal measured in [°C] (50, 100 etc.), as inputted on the Excel file
    end

    if isnan(rawdata{2,4})
        records(i).params.natT = [NaN,NaN];                                 % default value (NaN), if data is missing in Excel
    else
        records(i).params.natT = [rawdata{2,4:5}];                          % sample’s natural temperature [°C], as inputted on the Excel file
    end

    records(i).params.natDdot = [rawdata{3,4:5}]./ka;                       % sample’s environmental dose rate, as inputted on the Excel file, and converted from [Gy/ka] to [Gy/s]

    firstrow = 5;                                                           % fixed start limit of the raw data
    lastrow = size(rawdata,1);                                              % maximum limit of the raw data
    for k=firstrow:2:lastrow                                                % allows to extract the data of time (or dose) and associated response, as they are in pairs in the Excel file
        j = (k-firstrow)/2+1;
        records(i).rawdata(j).T = rawdata{k,2};                             % temperature of measurement [°C]
        records(i).rawdata(j).labDdot = rawdata{k+1,2};                     % instrument dose rate [Gy/s]
        records(i).rawdata(j).t = [rawdata{k,4:end}].*1e3;                  % measurement time (irradiation or delay time), converted from [ks] to [s]
        records(i).rawdata(j).L = [rawdata{k+1,4:end}];                     % measured data (not normalised)
    end
end

save(['./ComputeData/' filename '.mat'],'records');                         % exports the data output to a .MAT file in the "Computed Data" folder, for reading into other stages


%%% Running time
tEnd = toc(tStart);
fprintf('Stage1_ExcelToStruct took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));