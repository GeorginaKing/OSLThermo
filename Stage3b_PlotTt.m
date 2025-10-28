%%%% STAGE 3b, PlotTt %%%%
%%%% Plot cooling histories - requires output of Stage 3a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Frédéric Herman, 2015         frederic.herman@unil.ch %
% Arnaud Duverger, 2015         arnaud.duverger@ens.fr %
% Chloé Bouscary, 2024          chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Define figures parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultFigureColormap',jet());
set(0,'defaultfigurecolor',[1 1 1]);
set(gca,'FontSize',12);
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

%%% Load inversion output (Stage3a) and original n/N data %%%%%%%%%%%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_Tt.mat']);
misfit = Tt.misfit;
Time = Tt.time;
Temp = Tt.temp;
nNs = Tt.nNmod; maxnNs = max(nNs,[],2);
nNnat = Tt.nNnat; snNnat = Tt.snNnat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting parameters (TO CHANGE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define parameters for plotting
nAvM = 100;                                                                 % resolution of pdf

TypeMeasurement = Tt.TypeMeasurement;
TypeSignal = Tt.TypeSignal;

%%% Window of interest
ToI = 0.5;          % time of interest, last x [Myr]
tempoI = 120;       % temperatures of interest, below x [°C]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Restrict size
[m,n,nt] = size(maxnNs);
maxnNs = reshape(maxnNs,m,nt);

%%% To compute the PDF for the likelihood, the Tt paths are reinterpolated onto a grid (Av_matrix)
time = Time(1,:);
temp = Temp(1,:);
time_max = max(time);
time_min = min(time);
dt = (time_max-time_min)/(nAvM-1);
Temp_max = max(temp);
Temp_min = 0;
dT = (Temp_max-Temp_min)/(nAvM-1);
vec_time = time_min:dt:time_max;
vec_Temp = Temp_min:dT:Temp_max;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data selection %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Rejection algorithm
prob = exp(-misfit);                                                        % prob = likelihood
scale = max(prob); prob = prob/scale;                                       % normalisation by the maximum of the likelihood, max(prob)

R = rand(m,1); test = prob>R;                                               % rejection criteria = for every position of prob: 1 if true (prob>R), 0 if wrong (prob<=R)

idefix = find(test);                                                        % gives the name of the position for which test is true (=1)
idefix = idefix(end:-1:1);                                                  % invert the data for the line before, to order them in decreasing order (now worse fit to best fit)
movea = length(idefix);                                                     % length of the selected data (nb of data that passed the rejection criteria)

%%% Define the color scheme
max_likelihood = max(prob); min_likelihood = min(prob);                     % color scheme bounds
map = colormap(jet);
index_color = max(1,floor((length(map)-1)*(prob-min_likelihood)/(max_likelihood-min_likelihood))+1); % scaling the value of misfit to the colorbar scale

MTEMP = ones(movea,1)*TypeSignal(1:nt);
INDEX = index_color(idefix)*ones(1,nt);
MaxnNs = maxnNs(idefix,:);

Av_matrix = zeros(nAvM);

%%% Add accepted model to a matrix to compute the PDF, this is based on the rejection algorithm
for k = idefix'
    vec_temp = interp1(time,Temp(k,:),vec_time,'linear');
    Tpath = (0:nAvM-1)*nAvM+round((vec_temp-Temp_min)/dT)+1;
    Av_matrix(Tpath) = Av_matrix(Tpath)+1;
end
X = cumsum(Av_matrix/movea);                                                % for computing CIs and median


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Define the limits of the colorbars depending on the misfit of the accepted data
nbthicks=5;
col = round(100*(max(misfit(idefix)):(min(misfit(idefix))-max(misfit(idefix)))/nbthicks:min(misfit(idefix))))/100; % create colorbar limits for the misfit of accepted data, and round the values to 2 decimals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (1) Contrasts (n/N)mod with (n/N)nat %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure(1);

h1 = scatter(MTEMP(:),MaxnNs(:),30,INDEX(:),'filled'); hold on;             % predicted nNmod
h2 = errorbar(TypeSignal,nNnat,snNnat,'ko','MarkerEdgeColor','k','LineWidth',1.5,'MarkerSize',15); % observed values nNnat
legend([h2,h1],'nN_{nat}','nN_{mod}','Location','NorthWest');
pos1 = get(gca,'Position'); axis square; box on;

xlabel('IRSL (°C)'); ylabel('(n/N)');
xticks(TypeSignal); xticklabels(string(TypeSignal));
axis([min(TypeSignal)-20 max(TypeSignal)+20 0 1]);

%%% Define the colorbar to correspond to the misfit of the accepted data
c1 = colorbar('yticklabel',col);
set(c1,'ylim',[min(c1.Limits) max(c1.Limits)]);
set(c1,'XTick',linspace(min(c1.Limits),max(c1.Limits),nbthicks+1));
set(c1,'XTickLabel',strsplit(num2str(col)));
set(get(c1,'title'),'string','Misfit'); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (2) Modelled tT for individual signals %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure (2);

c1 = colorbar('yticklabel',col);
set(c1,'ylim',[min(c1.Limits) max(c1.Limits)]);
set(c1,'XTick',linspace(min(c1.Limits),max(c1.Limits),nbthicks+1));
set(c1,'XTickLabel',strsplit(num2str(col)));
set(get(c1,'title'),'string','Misfit'); hold on;

xlabel('Time (Ma)'); ylabel('Temperature (°C)');
axis([time_min time_max 0 Temp_max],'square');
axis ij;                                                                    % plots the Y-axis with values in reverse order
xticks = time_min:time_max/5:time_max;                                      % defines the tick-marks
xticklabels = time_max:-time_max/5:time_min;                                % defines the tick-labels
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);                       % plots the X-axis with values in reverse order
box on; hold on;
P = plot(time,Temp(idefix,:));                                              % 'spaghetti' plot
for curve = 1:movea
    set(P(curve),'Color',map(index_color(idefix(curve)),:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (3) Modelled dose response for individual signals %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f3 = figure(3);

set(gcf,'renderer','Painters')

for k = 1:nt
    subplot(2,(nt+mod(nt,2))/2,k);

    c1 = colorbar('yticklabel',col);
    set(c1,'ylim',[min(c1.Limits) max(c1.Limits)]);
    set(c1,'XTick',linspace(min(c1.Limits),max(c1.Limits),nbthicks+1));
    set(c1,'XTickLabel',strsplit(num2str(col)));
    set(get(c1,'title'),'string','Misfit'); hold on;

    xlabel('Time (Ma)'); ylabel('(n/N)');
    axis([time_min time_max 0 1],'square');
    xticks = time_min:0.2:time_max;                                         % define the tick-marks
    xticklabels = time_max:-0.2:time_min;                                   % define the tick-labels
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);                   % plot the X-axis with values in reverse order
    box on; hold on;

    P = plot(time,nNs(idefix,:,k));                                         % 'spaghetti' plot
    h = errorbar(time_max,nNnat(k),snNnat(k),'ko','MarkerFaceColor','k','LineWidth',1.5,'MarkerSize',5); % plot the natural signals
    legend(h,'nN_{nat}','Location','NorthWest');
    for curve = 1:movea
        set(P(curve),'Color',map(index_color(idefix(curve)),:));
    end
    title (sprintf('%s %d °C',cell2mat(TypeMeasurement(k)),TypeSignal(k)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (4) Cooling histories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure(4);

set(gcf,'renderer','Painters'); axis square; box on; hold on;

xlabel('Time (Ma)'); ylabel('Temperature (°C)');
axis([time_max-ToI time_max 0 tempoI],'square');
axis ij;                                                                    % plot the Y-axis with values in reverse order
xticks = time_max-ToI:time_max/10:time_max;                                 % define the tick-marks
xticklabels = time_min+ToI:-time_max/10:time_min;                           % define the tick-labels
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);                       % plot the X-axis with values in reverse order
c4 = colorbar; set(get(c4,'title'),'string','PDF'); hold on
[cs,hc] = contourf(vec_time,vec_Temp,Av_matrix./movea,100);                 % plot PDF
set(hc,'EdgeColor','none');
[cont3,h3] = contour(vec_time,vec_Temp,X,[0.025,0.975],'k','LineWidth',2);  % 95% CI = 2-sigma (CI=Confidence Interval)
[cont2,h2] = contour(vec_time,vec_Temp,X,[0.16,0.84],'g','LineWidth',2);    % 68% CI = 1-sigma
[cont1,h1] = contour(vec_time,vec_Temp,X,1,'r','LineWidth',2);              % median model
clim([0 0.1]);                                                              % limit of the colorbar, arbitrary scaled to allow good visualisation, max. value of 1

AgeOSL = Tt.AgeOSL(:,1)./1e3;
maxAgeOSL = Tt.maxAgeOSL./1e3;
plotvec=zeros(length(AgeOSL(:,1)),1)+tempoI;
h5 = plot(time_max-maxAgeOSL,plotvec,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[17 17 17]/255);
h4 = plot(time_max-AgeOSL(:,1),plotvec,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','yellow');

leg = legend([h1 h2 h3 h4 h5],'    Median','    68% CI','    95% CI','    OSL Age','    Max OSL Age','Location', 'SouthEast');


%% Save plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select output file format
print(f1,'-dpng',['./Figures/' filename '_nN_misfit_' SAR_MODEL '_' ITH_MODEL '.png']);
print(f2,'-dpng',['./Figures/' filename '_Tt_misfit_' SAR_MODEL '_' ITH_MODEL '.png']);
print(f3,'-dpng',['./Figures/' filename '_DResp_misfit_' SAR_MODEL '_' ITH_MODEL '.png']);
print(f4,'-dpng',['./Figures/' filename '_Cool_misfit_' SAR_MODEL '_' ITH_MODEL '.png']);
% print(f1,'-dsvg',['./Figures/' filename '_nN_misfit_' SAR_MODEL '_' ITH_MODEL '.svg']);
% print(f2,'-dsvg',['./Figures/' filename '_Tt_misfit_' SAR_MODEL '_' ITH_MODEL '.svg']);
% print(f3,'-dsvg',['./Figures/' filename '_DResp_misfit_' SAR_MODEL '_' ITH_MODEL '.svg']);
% print(f4,'-dsvg',['./Figures/' filename '_Cool_misfit_' SAR_MODEL '_' ITH_MODEL '.svg']);
% print(f1,'-depsc',['./Figures/' filename '_nN_misfit_' SAR_MODEL '_' ITH_MODEL '.eps']);
% print(f2,'-depsc',['./Figures/' filename '_Tt_misfit_' SAR_MODEL '_' ITH_MODEL '.eps']);
% print(f3,'-depsc',['./Figures/' filename '_DResp_misfit_' SAR_MODEL '_' ITH_MODEL '.eps']);
% print(f4,'-depsc',['./Figures/' filename '_Cool_misfit_' SAR_MODEL '_' ITH_MODEL '.eps']);


%%% Running time
tEnd = toc(tStart);
fprintf('Stage3b_PlotTt took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));