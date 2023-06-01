%%%% STAGE1, ExcelToStructure %%%%
%%%% Converts .xlsx files to .MAT data for input to subsequent scripts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benny Guralnik, 2014 benny.guralnik@gmail.com %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic %records time of execution

clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR

ka = 1e3.*365.*24.*3600; % ka in seconds
[exist,sheets] = xlsfinfo([filename '.xlsx']); % sheets = excel sheet names (e.g. IRSL50, IRSL100), exist = file type (Excel spreadsheet)

for i=1:length(sheets) 
	[~,~,raw] = xlsread([filename '.xlsx'],i,'A1:Q30'); % (A1:Q30) = extent of the excel data
    raw(~cellfun(@isfloat,raw)) = {NaN};                % applies the function @isfloat to the contents of each cell, one cell at a time

	records(i).id = sheets{i};
	records(i).params.natT = [raw{2,4:5}];              % Natural temperature
	records(i).params.natDdot = [raw{3,4:5}]./ka;       % Natural dose rate, converted to Gys
	lastrow = length(raw);
	firstrow = 5;   

	for k=firstrow:2:lastrow
		j = (k-firstrow)/2+1;
		records(i).rawdata(j).T = raw{k,2};             % Temperature
		records(i).rawdata(j).Ddot = raw{k+1,2};        % Instrument dose rate
		records(i).rawdata(j).t = [raw{k,4:end}].*1e3;  % Measurement time (irradiation or delay time)
		records(i).rawdata(j).L = [raw{k+1,4:end}]./max([raw{k+1,4:end}]); %normalises the luminescence data to the maximum signal intensity
	end
end

save(['./ComputeData/' filename '.mat'],'records') %exports MAT file to ComputData folder
toc
