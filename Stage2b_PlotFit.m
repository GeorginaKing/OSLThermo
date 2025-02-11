%%%% STAGE 2b, PlotLxTx %%%%
%%%% Plots data fits with selected models, requires input from Stage 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Georgina King, 2015 georgina.king@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR; 
close all; clc;

load(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.mat']); 
nt = length(records);                   %number of signals

colmat = colormap;
index = linspace(1,61,sum([ones(size(records(1).rawdata(4:3+NITL)))]));
Index = ones(12,1)*index; % This might be to change
logtick = logspace(0,9,10);
MTemp=[50 225 50 225 50 225 50 225 50 225 50 225 50 225 50 225 50 225];  %Arbitrary numbers, can be changed 

textmod = '%s: %0.2f +/- %0.2f';

for i=1:nt
	figure(i);
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8])
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % Read in vectors and data for plotting
    fx = records(i).plot.xf;
	fy = records(i).plot.yf;
	fdelta = records(i).plot.syf'; 
    fdelta = fdelta';
	ft = records(i).plot.tf;
	fL = records(i).plot.Lf;

	fxx = [fx fliplr(fx)];
	fd = [fy-fdelta fliplr(fy+fdelta)];

    % Plot fading data
	subplot(1,3,1,'xScale','log'); xlabel('t* (s)');
	ylabel('Luminescence Response (a.u.)');
	axis([1e2 1e7 0.7 1.1],'square'); box on; hold on;

    fill([fxx],[fd],[0.5 0.5 1],'EdgeColor','none','MarkerEdgeColor', 'k');
	plot(fx,fy, 'r-'); scatter(ft,fL,[], 'r','filled', 'MarkerEdgeColor', 'k');

    % Add text to figure to list sample parameters
	text(3e6,1.075,'A','fontweight','bold'); 
	name = fieldnames([records(i).params]);
	value = struct2cell([records(i).params]);
	txt = sprintf(textmod, name{4}, value{4}(1), abs(value{4}(2) - value{4}(1))); %rhop
	text(2e2,0.85,txt,'VerticalAlignment','top');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dose Response Curve %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % Read in vectors and data for plotting
    x = records(i).plot.x;
	y = records(i).plot.y;
	delta = records(i).plot.sy'; 
    delta = delta';
	t = records(i).plot.t;
	L = records(i).plot.L;

	ok = isfinite(t);
	UFy = y./max(y);
	ffun = find(UFy==1);
	UFy(ffun:end)=1;

	xx = [x, fliplr(x)];
	d = [y-delta fliplr(y+delta)];
	t(t==0 | isnan(t)) = 2; 

    % Plot dose response curve
	subplot(1,3,2,'xScale','log');
	xlabel('Time (s)'); ylabel('(n/N)');
	axis([1 1e6 0 1.2],'square'); box on; hold on;

	fill(xx,d,[0.5 0.5 1], 'EdgeColor','none');
	plot(x,y, 'r-');
	plot(x,UFy, 'k-');
	scatter(t(ok),L(ok),[], 'r','filled','MarkerEdgeColor', 'k');
	scatter(t(~ok),L(~ok),[], 'y','filled','MarkerEdgeColor', 'k');

    text(3e5,1.05,'B','fontweight','bold');
    if strcmp(SAR_MODEL,'GAUSS')
        TxT1 = sprintf(textmod, 'n/N', records(i).nNnat(1), records(i).nNnat(2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'D_0', round(records(i).params.D0(1),1), round(records(i).params.D0(2),1));
        TxT3 = sprintf('%s: %0.1f +/- %0.1f', '\sigma_{D_0}', round(records(i).params.sigmaD0(1),1), round(records(i).params.sigmaD0(2),1));
        TxT=strvcat(TxT1,TxT2,TxT3);
    elseif strcmp(SAR_MODEL,'GOK')
        TxT1 = sprintf(textmod, 'n/N', records(i).nNnat(1), records(i).nNnat(2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'D_0', round(records(i).params.D0(1),1), round(records(i).params.D0(2),1));
        TxT3 = sprintf('%s: %0.1f +/- %0.1f', 'a', round(records(i).params.GOK_a(1),1), round(records(i).params.GOK_a(2),1));
	    labddotlab = sprintf('Laboratory dose rate: %0.2f Gy s^-^1', records(i).rawdata(1).Ddot); 
        TxT=strvcat(TxT1,TxT2,TxT3);
    else
        TxT1 = sprintf(textmod, 'n/N', records(i).nNnat(1), records(i).nNnat(2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'D_0', round(records(i).params.D0(1),1), round(records(i).params.D0(2),1));
        TxT=strvcat(TxT1,TxT2);
    end
	text(3,0.90,TxT,'horizontalalignment','left');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % Read in vectors and data for plotting
    itlx = records(i).plot.xitl;
	itly = records(i).plot.yitl;
	itldelta = records(i).plot.syitl;
	itlt = records(i).plot.titl;
	itlL = records(i).plot.Litl;
	A = records(i).plot.scale;

	invA = ones(length(itly),1)*(1./A');
	d1 = (itly-itldelta).*invA; d1(d1<0)=0.0001;
	d2 = (itly+itldelta).*invA; d2(d2<0)=0.0001;
	itly = itly.*invA;
	invA = ones(length(itlt),1)*(1./A');
	itlL = itlL.*invA;

	itlxx = [itlx; flipud(itlx)]; 
	itld = [d1; flipud(d2)];

    % Plot isothermal decay data
	subplot(1,3,3,'xScale','log');
	xlabel('Time (s)');
	ylabel('Luminescence Response (a.u.)'); 
	axis([1 1e6 0 1.1],'square'); box on; hold on;

	fill(itlxx,itld,[0.5 0.5 1], 'EdgeColor','none');
	P = plot(itlx,itly);
	scatter(itlt(:),itlL(:),[], Index(:),'filled','MarkerEdgeColor','k');

    text(3e5,1.05,'C','fontweight','bold'); 
    if strcmp(ITL_MODEL,'GAUSS')
        TxT1 = sprintf(textmod, 'E', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
        TxT3 = sprintf(textmod, 'log_{10}\rho', round(records(i).params.rhop10(1),2), round(abs(records(i).params.rhop10(3)-records(i).params.rhop10(1)),2));
        TxT=strvcat(TxT1,TxT2,TxT3);
    elseif strcmp(ITL_MODEL,'BTS')
        TxT1 = sprintf(textmod, 'E', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
        TxT3 = sprintf(textmod, 'log_{10}\rho', round(records(i).params.rhop10(1),2), round(abs(records(i).params.rhop10(3)-records(i).params.rhop10(1)),2));
        TxT4 = sprintf('%s: %0.3f +/- %0.3f', 'Eu', round(records(i).params.Eu(1),3), round(records(i).params.Eu(2),3));
        TxT=strvcat(TxT1,TxT2,TxT3,TxT4);
    else
        TxT1 = sprintf(textmod, 'E', round(records(i).params.Et(1),2), round(records(i).params.Et(2),2));
        TxT2 = sprintf('%s: %0.1f +/- %0.1f', 'log_{10}s', round(records(i).params.s10(1),2), round(records(i).params.s10(2),2));
        TxT3 = sprintf(textmod, 'log_{10}\rho', round(records(i).params.rhop10(1),2), round(abs(records(i).params.rhop10(3)-records(i).params.rhop10(1)),2));
        TxT4 = sprintf('%s: %0.1f +/- %0.1f', 'b', round(records(i).params.GOK_b(1),1), round(records(i).params.GOK_b(2),1));
        TxT=strvcat(TxT1,TxT2,TxT3,TxT4);
    end
	text(2,0.2,TxT,'horizontalalignment','left');
    
   	subplot(1,3,1)
    text(2e2,1.15,[records(i).id],'fontweight','bold')
    subplot(1,3,2); text(2,1.05, [SAR_MODEL],'fontweight','bold'); 
    subplot(1,3,3); text(2,1.05, [ITL_MODEL],'fontweight','bold');

    print(['./Figures/' filename '_' records(i).id '_' SAR_MODEL '_' ITL_MODEL '.svg'],'-dsvg','-vector'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kars plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in data for plotting
nNnat = reshape([records.nNnat],2,nt)'; nNnat=nNnat(:,1);
snNnat = reshape([records.nNnat],2,nt)'; snNnat=snNnat(:,2);
KarsresAv = reshape([records.KarsresAv],2,nt)'; KarsresAv=KarsresAv(:,1);
sKarsresAv = reshape([records.KarsresAv],2,nt)'; sKarsresAv=sKarsresAv(:,2);

col=[0 1 0
	0 0 1
	1 0 0
	1 1 0
    0 1 1
    1 0 1
    0.5 0.5 0.5
    1 0.3 0.3
    0.3 1 0.3
    0.3 0.3 1
    0.8 0.2 0.2
    0.2 0.8 0.2
    0.2 0.2 0.8
    0.1 0.8 0.8
    0.8 0.1 0.8
    0.8 0.8 0.1
    0.4 0.6 0.6
    0.6 0.4 0.6
    0.6 0.6 0.4];
col=col(1:nt,:);
colorPerSample = [[255,51,51];[255,255,1];[51,255,51];[51,255,255];[51,153,255];[255,102,255];[255,178,102];[0.5,0.5,0.5];[204,0,102];...
    [255,51,51];[255,255,1];[51,255,51];[51,255,255];[51,153,255];[255,102,255];[255,178,102];[0.5,0.5,0.5];[204,0,102];]./255;

h1 = figure(i+1);
xlabel('(n/N)_n_a_t'); ylabel('(n/N)_s_s');  
axis([0 1 0 1],'square'); box on; hold on;

fill([0 1 1],[0 0 1],[0.5 0.5 1]);
p15_x=[0,1]; p15_y=[0,0.85]; plot(p15_x,p15_y,'blue');
m15_x=[0,1]; m15_y=[0,1.15]; plot(m15_x,m15_y,'blue');
for j=1:nt
errorbar(nNnat(j),KarsresAv(j),sKarsresAv(j),sKarsresAv(j),snNnat(j),snNnat(j),'ko','MarkerFaceColor',colorPerSample(j,:));
end

% Generate dummy info for plot handles "h"
h = zeros(nt/2,1);
names = {};
for j=1:nt/2
h(j) = scatter(nNnat(j),KarsresAv(j),40,'MarkerFaceColor',colorPerSample(j,:),'MarkerEdgeColor','k');
names{j} = records(j*2-1).id(:);
end

% Define legend according to handles "h"
% names = ["IRSL 50 ^oC","IRSL 225 ^oC","IRSL 50 ^oC","IRSL 225 ^oC"]
legend(h,names(:),'Location','southeast');
text(0.05,0.95,filename,'fontweight','bold');

%Select output file format
% print(['./Figures/' filename '_Kars' '_' SAR_MODEL '_' ITL_MODEL '.png'],'-dpng');
print(['./Figures/' filename '_Kars' '_' SAR_MODEL '_' ITL_MODEL '.eps'],'-depsc');
