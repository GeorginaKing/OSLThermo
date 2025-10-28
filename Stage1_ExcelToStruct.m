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
    [~,txt,raw] = xlsread([filename '.xlsx'],i);                            % read data
    raw(~cellfun(@isfloat,raw)) = {NaN};                                    % creates cells for the raw data to be stored within
    raw1=cell2mat(raw);                                                     % transforms raw(cells) in raw1(array double) to remove the NaN values at the end of the lines and columns
    s=size(raw1,2);
    rawcol=raw1(:,s);
    while all(isnan(rawcol)) == true                                        % selects the columns that have only NaN at the end of the file, and remove them
        rawcol(:,all(isnan(rawcol),1)) = [];                                % removes all rows with all NaN
        raw1(:,s)=[];
        s=s-1;
        rawcol=raw1(:,s);
    end
    l=length(raw1);
    rawlin=raw1(l,:);
    while all(isnan(rawlin)) == true                                        % selects the rows that have only NaN at the end of the file, and remove them
        rawlin(all(isnan(rawlin),2),:) = [];                                % removes all rows with all NaN
        raw1(l,:)=[];
        l=l-1;
        rawlin=raw1(l,:);
    end
    raw=num2cell(raw1);                                                     % transforms back raw1(array double) to raw(cells)


    %%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    records(i).id = sheets{i};                                              % name of the measured signal, as indicated by the name of the Excel sheets (e.g. IRSL50, IRSL100 etc.)
    records(i).typeMeasurement = {txt{1,7}};                                % type of measured signal (IRSL, OSL etc.), as inputted on the Excel file
    records(i).typeSignal = [raw{1,8}];                                     % temperature of the signal measured in [°C] (50, 100 etc.), as inputted on the Excel file
    records(i).params.natT = [raw{2,4:5}];                                  % sample’s natural temperature [°C], as inputted on the Excel file
    records(i).params.natDdot = [raw{3,4:5}]./ka;                           % sample’s environmental dose rate, as inputted on the Excel file, and converted from [Gy/ka] to [Gy/s]

    firstrow = 5;                                                           % fixed start limit of the raw data
    lastrow = size(raw,1);                                                  % maximum limit of the raw data
    for k=firstrow:2:lastrow                                                % allows to extract the data of time (or dose) and associated response, as they are in pairs in the Excel file
        j = (k-firstrow)/2+1;
        records(i).rawdata(j).T = raw{k,2};                                 % temperature of measurement [°C]
        records(i).rawdata(j).labDdot = raw{k+1,2};                         % instrument dose rate [Gy/s]
        records(i).rawdata(j).t = [raw{k,4:end}].*1e3;                      % measurement time (irradiation or delay time), converted from [ks] to [s]
        records(i).rawdata(j).L = [raw{k+1,4:end}];                         % measured data (not normalised)
    end
end

save(['./ComputeData/' filename '.mat'],'records');                         % exports the data output to a .MAT file in the "Computed Data" folder, for reading into other stages


%%% Running time
tEnd = toc(tStart);
fprintf('Stage1_ExcelToStruct took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));