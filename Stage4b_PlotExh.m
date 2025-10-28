%%% STAGE 4b, PlotExh %%%%
%%% Plots the results of Stage 4a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rabiul Biswas, 2018       biswasrabiul@gmail.com %
% Nadja Stalder, 2020       nadjafranziska.stalder@unil.ch %
% Frédéric Herman, 2020     frederic.herman@unil.ch %
% Georgina King, 2022       georgina.king@unil.ch %
% Chloé Bouscary, 2025      chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Define figures parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'units','centimeters')
set(gca,'FontSize',6)
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

%%% Load inversionExh output (Stage4a) and original data (Stage2a) %%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_zt.mat']);
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.mat']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting parameters (TO CHANGE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define parameters for plotting
DoI   = 2;                      % depth of interest, last x [km]
ToI   = 0.3;                    % time of interest, last x [Ma]
G_cut = 100;                    % cut-off for geothermal gradient [°C/km]
time_window = ToI;              % restrict to a time-window of interest for plotting

%%% Resolution of matrix
nAvM = 1200;                    % for temporal dt = ToI/nAvM & spatial dZ = DoI/nAvM resolutions
resampling = 100;               % how many times median is resampled (min. 100)
nb_path    = 0.10;              % number of resampled path in % of total accepted path (min. 10%)
vel_interval = 10;              % define velocity smoothing interval in [kyr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Measured natural signal
nNnat = zt.nNnat;
snNnat = zt.snNnat;

%%% Inversion model output
misOUT      = zt.misfit;
Zt          = zt.temp;
time        = zt.time;
G_end       = zt.G_end;
nNsort      = zt.maxnNsort;
maxnNsort   = max(nNsort,[],8);
nt          = length(maxnNsort(1,:));
ntmin       = 1;
ntmax       = nt;
MTemp       = 1:nt;

%%% Extract (n/N)ss, fading corrected ages and maximum possible ages (2*D0) from Stage2a output
for k = 1:nt
    Kars(k,:) = [records(k).params.nNss];
    AgeOSL(k,:) = [records(k).params.AgeOSL];
    maxAgeOSL(k,:) = [records(k).params.maxAgeOSL];
end
Kars = Kars(ntmin:ntmax,:); nNss = Kars(:,1); snNss = Kars(:,2);
AgeOSL = AgeOSL(ntmin:ntmax,:)./1e3;
maxAgeOSL = maxAgeOSL(ntmin:ntmax,:)./1e3;

%%% Extract name of the measured signals
TypeMeasurement = zt.TypeMeasurement;
TypeSignal = zt.TypeSignal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity Field Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define box size
time_max=max(time);
time_min=time_max-time_window;                                              % limit matrix time to RoI using time_window
Z_max = max(Zt(1,:));
Z_min = 0;

%%% Compute PDF for the likelihood, Zt paths interpolated onto grid
Sz       = ones(length(misOUT(:,1)),nt);
[m,nt]   = size(Sz);
dt       = (time_max-time_min)/(nAvM-1);                                    % calculate time step for each part of the matrix
dZ       = (Z_max-Z_min)/(nAvM-1);                                          % -1 so that length(Zvec) = length(nAvM)
vec_time = time_min:dt:time_max;                                            % time vector of nAvM grid
Zvec     = Z_min:dZ:Z_max;                                                  % depth vector of nAvM grid

%%% Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prob = exp(-misOUT);                                                        % prob = likelihood
scale = max(prob); prob = prob/scale;                                       % normalisation by the maximum of the likelihood, max(prob)

%%% Accept all final geothermal gradient, and apply rejection algorithm
%R = rand(m,1); test = prob>R;                                              % rejection criteria = for every position of prob: 1 if true (prob>R), 0 if wrong (prob<=R)

%%% Accept only final geothermal gradient > G_cut
% **set probability very low if geothermal gradient is higher than G_cut**
% **need to take care when applying this as if the fit is poor will lead**
% **to acceptance of all of the data**
for i = 1:length(G_end); if G_end(i) > G_cut; prob(i) = 1E-16; end; end
test = prob>1E-16;                                                          % rejection criteria = for every position of prob: 1 if true (prob>1E-16), 0 if wrong (prob<=1E-16)

id=find(test);                                                              % gives the name of the position for which test is true (=1)
id=id(end:-1:1);                                                            % invert the data for the line before, to order them in decreasing order (now worse fit to best fit)
movea=length(id);                                                           % length of the selected data (nb of data that passed the rejection criteria)

max_likelihood=max(prob(id)); min_likelihood=min(prob(id));                 % bounds for the colour scheme
index_color=max(1,floor(63*(prob-min_likelihood)/(max_likelihood-min_likelihood))+1); % scaling the value between 1 and 64
MTEMP=ones(movea,1)*MTemp(1:nt); INDEX=index_color(id)*ones(1,nt);          % create index for colors
MaxnNsort=maxnNsort(id,:);

%%% Add accepted model to a matrix to compute the PDF, based on the rejection algorithm
Av_matrix = zeros(nAvM);
for k = id'
    vec_Z            = interp1(time,Zt(k,:),vec_time,'linear');
    Zpath            = (0:nAvM-1)*nAvM+round((vec_Z-Z_min)/dZ)+1;
    Av_matrix(Zpath) = Av_matrix(Zpath)+1;                                  % grid with the number of paths pathing by each cell
end
X = cumsum(Av_matrix./movea);                                               % for computing CIs and median


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (1) Final geothermal gradient by misfit %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure(1);
set(gcf,'renderer','Painters');

subplot(1,2,1); axis square; box on; hold on
scatter(misOUT, G_end, 'b', 'filled'); hold on
scatter(misOUT(id), G_end(id), 'r', 'filled');
xlabel('Misfit'); ylabel('Geothermal gradient (°C/km)');
title('All data')

subplot(1,2,2); axis square; box on; hold on
scatter(misOUT(id), G_end(id), 'r', 'filled');
xlabel('Misfit'); ylabel('Geothermal gradient (°C/km)');
ylim([0 G_cut+10]);
title('Accepted data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (2) Modelled vs. measured n/N values %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col1 = round(100*(max(misOUT(id)):(min(misOUT)-max(misOUT(id)))/6:min(misOUT)))/100;

f2 = figure(2);
set(gcf,'renderer','Painters');

%%% Right-side: Modelled vs. measured n/N values %%%
subplot(1,2,2);
xlabel('Luminescence signal'); ylabel('(n/N)'); axis square; box on; hold on;
c1 = colorbar('yticklabel',col1,'location','northoutside'); set(get(c1,'title'),'string','Misfit'); colormap('viridis');
h1 = scatter(MTEMP(:),MaxnNsort(:),30,INDEX(:),'filled');
h2 = errorbar(MTemp(1:nt),nNnat(1:nt),snNnat(1:nt),'ko','MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',1.2,'MarkerSize', 5);
h3 = errorbar(MTemp(1:nt),nNss(1:nt),snNss(1:nt),'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 5);
legend([h2,h1,h3],'nN_{nat}','nN_{mod}', 'nN_{ss}','Location','NorthWest');
xlim([0 nt+1]); ylim([0 1]);
set(gca,'xtick',1:nt,'xticklabel',string(TypeSignal));

%%% Left-side: Probability plot using nAvM %%%
subplot(1,2,1);
axis square; box on; hold on; 
c3 = colorbar('location', 'northoutside'); ylabel(c3, 'PDF'); colormap('viridis');
contourf(vec_time,Zvec,Av_matrix./movea,[0 1e-8 0.01:0.01:1],'edgecolor','none'); shading flat; % PDF, where we have data, and then 1 category every 0.5%=0.05 pdf
onesig = contour(vec_time,Zvec,X,[0.16,0.84],'g','LineWidth',1);            % 68% CI = 1-sigma
twosig = contour(vec_time,Zvec,X,[0.025,0.975],'k','LineWidth',1);          % 95% CI = 2-sigma (CI=Confidence Interval)
median = contour(vec_time,Zvec,X,1,'color','r','LineWidth',1);              % median model
clim([0 0.05]);
axis([time_max-ToI time_max 0 DoI],'square');                               % plot only region of interest
% % % axis([time_max-ToI time_max 0 DoI-3],'square');                             % plot only region of interest
ax = gca;
ax.XTick=time_max-ToI:0.05:time_max;                                        % define the tick-marks
ax.XTickLabel=(ToI*1000):-50:0;                                             % define the tick-labels
set(gca,'YDir','Reverse');
xlabel('Time (ka)'); ylabel('Depth (km)');

%%% Add text to figure
txt=['# acc. paths = ' num2str(movea)];
plott=time_max-ToI+0.30; plotd=DoI-0.2; plotd2=DoI-0.4;
text(plott,plotd,txt,'fontweight','bold','col','w');
text(plott,plotd2,filename,'fontweight','bold','col','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (3) Exhumation rate on linear scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Compute velocity of median model for all accepted paths
t_median = median(1,2:end); z_median = median(2,2:end);                     % time and depth vectors of median
[t1_median,indx_median] = unique(t_median); z1_median = z_median(indx_median); % remove double time entries and sort from lowest to highest
npt = round((time_max-time_min)*1000); t_median_intp = linspace(time_min,time_max,npt); % interpolate median on time vector that has 1 ka intervals
z_median_intp = interp1(t1_median,z1_median,t_median_intp, 'linear');       % depth vector of median in 1 ka intervals

%%% Calculate velocity at a vel_interval to smooth erosion rate
t_vel = (time_max:-(vel_interval/1000):time_max-ToI)';                      % time-median value at every vel_interval kyr
z_vel = interp1(t_median_intp,z_median_intp,t_vel);                         % depth-median value at every vel_interval kyr

Vel_median = -diff(z_vel)./diff(t_vel); Vel_time = t_vel(1:end-1);

%%% Remove NaN value from final data point
for j = length(Vel_time)-1:-1:1; if isnan(Vel_median(j)); Vel_median(j) = Vel_median(j+1); end; end

%%% Resampling median for error propagation and calculate velocity fields
% for the statistics of the median velocity. Accepted paths are resampled
% to compute medians m_1 to m_resampling and their velocities, which are
% used to estimate the confidence intervals of the median velocity:
%
% Vel_res =          m1   m2   m3   m4 ...  m98   m99   m100
%               t1 | v11  v12  v13  v14 ... v198  v199  v1100 |
%               t2 | v21  v22  v23  v24 ... v298  v299  v2100 |
%                  |                     .                    |
%                  |                     .                    |
%                  |                     .                    |
%               t6 | v61  v62  v63  v64 ... v698  v699  v6100 |
%
% sort v along rows: lowest ----------------------> highest v
% If resampled 100 times, the 90%-CI are simply the 5th and 95th columns

nb_resPath     = ceil(nb_path*movea);
z_vel_res      = zeros(length(z_vel),100);
Vel_median_res = zeros(length(Vel_median),100);

for i = 1:resampling
    Av_matrix_res = zeros(nAvM);
    R_res         = randi(movea,1,nb_resPath);                              % choose randomly 10% of path out of accepted paths
    id_res        = id(R_res,1);                                            % indices of 100 path
    for k = id_res'
        vec_Z_res = interp1(time,Zt(k,:),vec_time,'linear');                % interpolate accepted modelled Zt-path on nAvM grid
        Zpath_res = (0:nAvM-1)*nAvM+round((vec_Z_res-Z_min)/dZ)+1;
        Av_matrix_res(Zpath_res) = Av_matrix_res(Zpath_res)+1;              % add 1 to each AvMn grid cell through which the kth path goes
    end
    X_res = cumsum(Av_matrix_res/nb_resPath);                               % for computing CIs and median

    figure(34);
    title('Resample for median model');
    median_res = contour(vec_time,Zvec,X_res,1,'color',[0.8,0.0,0.2],'LineWidth',1.5); hold on;
    ax = gca;
    ax.XTick=time_max-ToI:0.05:time_max;                                    % define the tick-marks
    ax.XTickLabel=(ToI*1000):-50:0;                                         % define the tick-labels
    xlabel('Time (ka)'); ylabel('Depth (km)');
    set(gca,'YDir','Reverse');

    t_median_res = median_res(1,2:end);                                     % time vector of median
    z_median_res = median_res(2,2:end);                                     % depth vector of median
    [t1_median_res,indx_median_res] = unique(t_median_res);                 % remove double time entries and sort from lowest to highest
    z1_median_res = z_median_res(indx_median_res);
    z_median_res_intp = interp1(t1_median_res,z1_median_res,t_median_intp, 'linear'); % depth vector of median in 1 ka intervals

    z_vel_res(:,i) = interp1(t_median_intp,z_median_res_intp,t_vel);        % calculate velocity on a 5 ka or 10 ka interval to smooth erosion rate (vel_int)
    Vel_median_res(:,i) = -diff(z_vel_res(:,i))./diff(t_vel);

    %%% Manipulation of last point as last point sometime have NaN value
    for j = length(Vel_time)-1:-1:1; if isnan(Vel_median_res(j,i)); Vel_median_res(j,i) = Vel_median_res(j+1,i); end; end

    clear Av_matrix_res R_res id_res vec_Z_res Zpath_res Av_matrix_res
    clear X_res median_res cut t_median_res z_median_res t1_median_res z1_median_res Z_median_res
end

%%% Sort data
[sortedVel_res,I] = sort(Vel_median_res,2);                                 % sort velocity vector along rows
Vel_68CI_u  = sortedVel_res(:,round(0.84*resampling));                      % upper bounds of 68% CI
Vel_68CI_l  = sortedVel_res(:,round(0.16*resampling));                      % lower bounds of 68% CI
Vel_95CI_u  = sortedVel_res(:,round(0.975*resampling));                     % upper bounds of 95% CI
Vel_95CI_l  = sortedVel_res(:,round(0.025*resampling));                     % upper bounds of 95% CI
med_control = sortedVel_res(:,round(0.5*resampling));                       % median of resampled medians

nstep        = length(Vel_time)-1;
tplot(1,:)   = Vel_time(1:nstep);

Vplot_median(1,:)  = Vel_median(1:nstep);
Vplot_68CI_l(1,:)  = Vel_68CI_l(1:nstep);
Vplot_68CI_u(1,:)  = Vel_68CI_u(1:nstep);
Vplot_95CI_l(1,:)  = Vel_95CI_l(1:nstep);
Vplot_95CI_u(1,:)  = Vel_95CI_u(1:nstep);
Vplot_med_res(1,:) = med_control(1:nstep);

Vplot_median(Vplot_median<=0)  = 0.001;                                     % fix minimum values to avoid negative rates
Vplot_68CI_l(Vplot_68CI_l<=0)  = 0.001;
Vplot_68CI_u(Vplot_68CI_u<=0)  = 0.001;
Vplot_95CI_l(Vplot_95CI_l<=0)  = 0.001;
Vplot_95CI_u(Vplot_95CI_u<=0)  = 0.001;
Vplot_med_res(Vplot_68CI_u<=0) = 0.001;

%%% Plotting velocity
X_plot = [tplot,fliplr(tplot)];                                             % create continuous x value array for plotting: (10 .. 9.8 9.8 .. 10)
Y_95CI = [Vplot_95CI_l,fliplr(Vplot_95CI_u)];                               % create y values for out and then back
Y_68CI = [Vplot_68CI_l,fliplr(Vplot_68CI_u)];                               % plot filled area

col = [0.8500    0.3250    0.0980];                                         % color for median

%%% Exhumation rate on linear scale

f3 = figure(3);
set(gcf,'renderer','Painters'); axis square; box on; hold on;

h4 = fill(X_plot,Y_95CI,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha', 0.6); % 95% CI
h3 = fill(X_plot,Y_68CI,[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha', 0.8); % 68% CI
h2 = plot(tplot,Vplot_median,'Color',col,'LineWidth', 2);                   % median

plotvec=zeros(length(AgeOSL(:,1)),1);
h6 = plot(time_max-AgeOSL(:,1),plotvec,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','yellow');
h7 = plot(time_max-maxAgeOSL,plotvec,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[17 17 17]/255);
leg = legend([h2 h3 h4 h6 h7],'    Median','    68% CI','    95% CI','    OSL Age','    Max OSL Age','Location', 'NorthEast');
title(leg, ['# acc. path = ' num2str(movea)],  'FontSize', 12, 'FontWeight','Normal'); legend ('boxoff');
xlim([time_max-ToI time_max]); ylim([0 ceil(max(Y_95CI)/10)*10+10]);
ax = gca;
ax.XTick=time_max-ToI:0.05:time_max;                                        % define the tick-marks
ax.XTickLabel=(ToI*1000):-50:0;                                             % define the tick-labels
xlabel('Time (ka)'); ylabel('Exhumation rate (mm/yr)');
set(gca,'YMinorGrid','Off');
set(gca,'FontSize',14); set(gca,'Layer','top'); ax = gca; ax.LineWidth = 1.5; ax.TickLength = [0.015,0.015];

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure (4) Exhumation rate on log scale %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure(4);
set(gcf,'renderer','Painters'); axis square; box on; hold on;

h4 = fill(X_plot,Y_95CI,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha', 0.6); % 95% CI
h3 = fill(X_plot,Y_68CI,[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha', 0.8); % 68% CI
h2 = plot(tplot,Vplot_median,'Color',col,'LineWidth', 2);                   % median

plotveclog=ones(length(AgeOSL(:,1)),1).*0.1;
h6 = plot(time_max-AgeOSL(:,1),plotveclog,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','yellow');
h7 = plot(time_max-maxAgeOSL,plotveclog,'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[17 17 17]/255);
leg = legend([h2 h3 h4 h6 h7],'    Median','    68% CI','    95% CI','    OSL Age','    Max OSL Age','Location', 'NorthEast','Autoupdate','off');
title(leg, ['# acc. path = ' num2str(movea)],  'FontSize', 12, 'FontWeight','Normal'); legend ('boxoff');
xlim([time_max-ToI time_max]); ylim([0.1 10^3]);
ax = gca;
ax.XTick=time_max-ToI:0.05:time_max;                                        % define the tick-marks
ax.XTickLabel=(ToI*1000):-50:0;                                             % define the tick-labels
xlabel('Time (ka)'); ylabel('Exhumation rate (mm/yr)');
set(gca, 'YScale', 'log');
set(gca,'YMinorGrid','Off'); set(gca,'FontSize',14); set(gca,'Layer','top'); ax = gca; ax.LineWidth = 1.5; ax.TickLength = [0.015,0.015];

hold off;


%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zt.Plotres.median_erosion = [tplot; Vplot_median];
zt.Plotres.CI_erosion = [Y_95CI; Y_68CI];
zt.Plotres.median_depth = [t_median; z_median];

Geo_cut = num2str(G_cut);
save(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_zt_' Geo_cut '.mat'],'zt');

%% Save plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select output file format
print(f1,'-dpng',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Geotherm_' Geo_cut '.png']);
print(f2,'-dpng',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Depth_Ages_' Geo_cut '.png']);
print(f3,'-dpng',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_lin_' Geo_cut '.png']);
print(f4,'-dpng',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_log_' Geo_cut '.png']);
% print(f1,'-dsvg',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Geotherm_' Geo_cut '.svg']);
% print(f2,'-dsvg',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Depth_Ages_' Geo_cut '.svg']);
% print(f3,'-dsvg',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_lin_' Geo_cut '.svg']);
% print(f4,'-dsvg',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_log_' Geo_cut '.svg']);
% print(f1,'-depsc',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Geotherm_' Geo_cut '.eps']);
% print(f2,'-depsc',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Depth_Ages_' Geo_cut '.eps']);
% print(f3,'-depsc',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_lin_' Geo_cut '.eps']);
% print(f4,'-depsc',['./Figures/' filename '_' SAR_MODEL '_' ITH_MODEL '_Vel_log_' Geo_cut '.eps']);


%%% Running time
tEnd = toc(tStart);
fprintf('Stage4b_PlotExh took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));