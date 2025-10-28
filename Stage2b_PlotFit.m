%%%% STAGE 2b, PlotLxTx %%%%
%%%% Plots data fits with selected models, requires input from Stage 2a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Georgina King, 2015       georgina.king@unil.ch %
% Chloé Bouscary, 2025      chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.mat']);
nt = length(records);                                                       % number of signals


%%% Plot the figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nt
    figure(i);                                                              % one figure per signal
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Fading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Read in vectors and data for plotting
    fx = records(i).plot.xf;
    fy = records(i).plot.yf;
    fdelta = records(i).plot.syf';
    fdelta = fdelta';
    ft = records(i).plot.tf;
    fL = records(i).plot.Lf;

    fxx = [fx fliplr(fx)];
    fd = [fy-fdelta fliplr(fy+fdelta)];

    %%% Plot fading data
    subplot(1,3,1,'xScale','log');
    xlabel('t* (s)'); ylabel('IRSL intensity (a.u.)'); %ylabel('(n/N)');
    axis([1e2 1e7 0.7 1.1],'square'); box on; hold on;
    set(gca,'XTick',logspace(2,7,6),'XScale','log');

    fill(fxx,fd,[0.5 0.5 1],'EdgeColor','none','MarkerEdgeColor', 'k');
    plot(fx,fy, 'r-');
    scatter(ft,fL,[], 'r','filled', 'MarkerEdgeColor', 'k');

    %%% Add text to figure
    text(3e6,1.08,'A','fontweight','bold');
    TxF1 = sprintf('%s = %0.2f +/- %0.2f', 'log_{10}\rho', records(i).params.rhop10(1), records(i).params.rhop10(2));
    TxF2 = sprintf('%s = %0.2f +/- %0.2f %/dec.', 'g2d', records(i).params.g2d(1), abs(records(i).params.g2d(2)));
    txtFAD=strvcat(TxF1,TxF2);
    text(3e4,0.77,txtFAD,'VerticalAlignment','top');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dose Response Curve (SAR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Read in vectors and data for plotting
    x = records(i).plot.x;
    y = records(i).plot.y;
    delta = records(i).plot.sy;
    t = records(i).plot.t;
    L = records(i).plot.L;

    %%% Transform regeneratives doses from [s] to [Gy]
    for j = 1:nSAR
        labDdot_ali(j) = records(i).rawdata(j).labDdot;                     % laboratory dose rate per aliquot [Gy/s]
    end
    labDdot = mean(labDdot_ali);
    x = x * labDdot;
    t = t * labDdot;

    ok = isfinite(t);
    UFy = y./max(y);
    ffun = find(UFy==1);
    UFy(ffun(1):end)=1;

    xx = [x,fliplr(x)];
    d = [y-delta fliplr(y+delta)];
    t(t==0 | isnan(t)) = 1;

    %%% Plot dose response curve
    %{
    % On a linear scale
    subplot(1,3,2);
    xlabel('Dose (Gy)'); ylabel('IRSL intensity (a.u.)'); %ylabel('(n/N)');
    axis([1 round(max(t)+1000,-3) 0 1.1],'square'); box on; hold on;
    %}
    % %{
    % On a logarithmic scale
    subplot(1,3,2,'xScale','log');
    xlabel('Dose (Gy)'); ylabel('IRSL intensity (a.u.)'); %ylabel('(n/N)');
    axis([1 1e6 0 1.1],'square'); box on; hold on;
    set(gca,'XTick',logspace(0,6,7),'XScale','log');
    % %}

    fill(xx,d,[0.5 0.5 1],'EdgeColor','none');
    plot(x,y,'r-');
    plot(x,UFy,'k-');
    scatter(t(ok),L(ok),[],'r','filled','MarkerEdgeColor','k');
    n=scatter(t(~ok),L(~ok),[],'y','filled','MarkerEdgeColor','k');
    legend(n,{'nN_n_a_t'},'Location','SouthEast');

    %%% Add text to figure
    if strcmp(SAR_MODEL,'SSE')
        TxT1 = sprintf('%s = %0.2f +/- %0.2f', 'nN_n_a_t', records(i).params.nNnat(1), records(i).params.nNnat(2));
        TxT2 = sprintf('%s = %0.f +/- %0.f Gy', 'D_0', round(records(i).params.D0(1)), round(records(i).params.D0(2)));
        TxTl = sprintf('%s = %0.2f Gy.s^-^1', 'Dr_l_a_b', labDdot);
        TxT=strvcat(TxT1,TxT2,TxTl);
    elseif strcmp(SAR_MODEL,'GOK')
        TxT1 = sprintf('%s = %0.2f +/- %0.2f', 'nN_n_a_t', records(i).params.nNnat(1), records(i).params.nNnat(2));
        TxT2 = sprintf('%s = %0.f +/- %0.f Gy', 'D_0', round(records(i).params.D0(1)), round(records(i).params.D0(2)));
        TxT3 = sprintf('%s = %0.1f +/- %0.1f', 'a', round(records(i).params.GOK_a(1),1), round(records(i).params.GOK_a(2),1));
        TxTl = sprintf('%s = %0.3f Gy.s^-^1', 'Dr_l_a_b', labDdot);
        TxT=strvcat(TxT1,TxT2,TxT3,TxTl);
    end
    %{
    % For linear scale
    text(round(max(t)+1000,-3)-500,1.05,'B','fontweight','bold');   
    text((0+250),1.05,SAR_MODEL,'fontweight','bold');                       % dose response curve model
    text(2500,0.30,TxT,'horizontalalignment','left');
    %}
    % %{
    % For logarithmic scale
    text(3e5,1.05,'B','fontweight','bold');
    text(2,1.05,SAR_MODEL,'fontweight','bold');                             % dose response curve model
    text(2e3,0.30,TxT,'horizontalalignment','left');
    % %}


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Define colors for plot
    colmat = colormap(parula);
    index = round(linspace(1,length(colmat),nITH));

    %%% Read in vectors and data for plotting
    ithx = records(i).plot.xith;
    ithy = records(i).plot.yith;
    ithdelta = records(i).plot.syith;
    itht = records(i).plot.tith;
    ithL = records(i).plot.Lith;
    A = records(i).plot.scale;

    isoT = [records(i).rawdata(nSAR+1:nSAR+nITH).T]; ok=isfinite(isoT); isoT=isoT(ok); % isothermal temperatures [°C]

    isNAN = sum(isnan(A));
    if isNAN < length(A)

        invA = ones(length(ithy),1)*(1./A');
        d1 = (ithy-ithdelta).*invA; d1(d1<0)=0.0001;
        d2 = (ithy+ithdelta).*invA; d2(d2<0)=0.0001;
        ithy = ithy.*invA;
        invA = ones(length(itht),1)*(1./A');
        ithL = ithL.*invA;

        ithxx = [ithx; flipud(ithx)];
        ithd = [d1; flipud(d2)];

        %%% Plot isothermal decay data
        subplot(1,3,3,'xScale','log');
        xlabel('Time (s)'); ylabel('IRSL intensity (a.u.)'); %ylabel('(n/N)');
        axis([1 1e6 0 1.1],'square'); box on; hold on;
        set(gca,'XTick',logspace(0,6,7),'XScale','log');

        fill(ithxx,ithd,[0.5 0.5 1],'EdgeColor','none');
        P = plot(ithx,ithy);

        %%% Create the legend
        lgd=cell(nITH,1);
        for k = 1:length(P)
            l(k)=scatter(itht(:,k),ithL(:,k),[],colmat(index(k),:),'filled','MarkerEdgeColor','k');
            set(P(k),'Color',colmat(index(k),:));
            lgd{k}=strcat(num2str(isoT(k)),' °C');
        end
        legend(l,lgd,'Location','best');

        %%% Add text to figure
        text(3e5,1.05,'C','fontweight','bold');
        text(2,1.05,ITH_MODEL,'fontweight','bold');                         % isothermal decay model
        if strcmp(ITH_MODEL,'BTS')
            TxT1 = sprintf('%s = %0.2f +/- %0.2f eV', 'Et', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
            TxT2 = sprintf('%s = %0.3f +/- %0.3f eV', 'Eu', round(records(i).params.Eu(1),3), round(records(i).params.Eu(2),3));
            TxT3 = sprintf('%s = %0.2f +/- %0.2f s^-^1', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
            TxT=strvcat(TxT1,TxT2,TxT3);
        elseif strcmp(ITH_MODEL,'GOK')
            TxT1 = sprintf('%s = %0.2f +/- %0.2f eV', 'E', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
            TxT2 = sprintf('%s = %0.2f +/- %0.2f s^-^1', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
            TxT3 = sprintf('%s = %0.1f +/- %0.1f', 'b', round(records(i).params.GOK_b(1),1), round(records(i).params.GOK_b(2),1));
            TxT=strvcat(TxT1,TxT2,TxT3);
        elseif strcmp(ITH_MODEL,'GAUSS')
            TxT1 = sprintf('%s = %0.2f +/- %0.2f eV', 'E', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
            TxT2 = sprintf('%s = %0.2f +/- %0.2f s^-^1', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
            TxT=strvcat(TxT1,TxT2);
        end
        text(2e0,0.15,TxT,'horizontalalignment','left');

    end


    %%% Save plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,3,1);
    text(2e2,0.76,filename,'fontweight','bold');                            % sample name
    TypeMeasurement(i) = records(i).typeMeasurement;                        % type of measured signal (IRSL, OSL, etc.)
    TypeSignal(i) = records(i).typeSignal;                                  % temperature of the signal measured in [°C] (50, 100, etc.)
    text(2e2,0.73,[cell2mat(TypeMeasurement(i)) ' ' num2str(TypeSignal(i)) ' °C'],'fontweight','bold'); % signal measured

    print('-dpng',['./Figures/' filename '_' records(i).id '_' SAR_MODEL '_' ITH_MODEL '.png']);
    % print('-dsvg',['./Figures/' filename '_' records(i).id '_' SAR_MODEL '_' ITH_MODEL '.svg']);
    % print('-depsc',['./Figures/' filename '_' records(i).id '_' SAR_MODEL '_' ITH_MODEL '.eps']);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kars plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Read in data for plotting
for i = 1:nt
    nNnat(i,:) = [records(i).params.nNnat];
    snNnat(i,:) = [records(i).params.nNnat];
    nNss(i,:) = [records(i).params.nNss];
    snNss(i,:) = [records(i).params.nNss];
end

nNnat=nNnat(:,1);
snNnat=snNnat(:,2);
nNss=nNss(:,1);
snNss=snNss(:,2);

col=[0 0 1
    0 1 0
    1 1 0
    1 0 0];                                                                 % defines the colors for the Kars plot, as many raw as there are measured signals
col=col(1:nt,:);                                                            % defines the colors depending on the number of measurement signals

figure(i+1);
xlabel('(n/N)_n_a_t'); ylabel('(n/N)_s_s');
axis([0 1 0 1],'square'); box on; hold on;

%%% Plot data
fill([0 1 1],[0 0 1],[0.5 0.5 1]);                                          % fills the saturation area
m15_x=[0,1]; m15_y=[0,1.15]; plot(m15_x,m15_y,'blue');                      % at 15% of saturation limit
errorbar(nNnat,nNss,snNss,snNss,snNnat,snNnat,'ko');

%%% Plot data and generate legend
lgd=cell(nt,1);                                                             % creates the legend
for i=1:nt
    h(i)=scatter(nNnat(i),nNss(i),40,col(i,:),'filled','o','MarkerEdgeColor','k');
    lgd{i}=sprintf('%s %d °C',cell2mat(TypeMeasurement(i)),TypeSignal(i));
end
legend(h,lgd,'Location','southeast');
text(0.05,0.95,filename,'fontweight','bold');

%%% Save Kars plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select output file format
print('-dpng',['./Figures/' filename '_Kars_' SAR_MODEL '_' ITH_MODEL '.png']);
% print('-dsvg',['./Figures/' filename '_Kars_' SAR_MODEL '_' ITH_MODEL '.svg']);
% print('-depsc',['./Figures/' filename '_Kars_' SAR_MODEL '_' ITH_MODEL '.eps']);


%%% Running time
tEnd = toc(tStart);
fprintf('Stage2b_PlotFit took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));