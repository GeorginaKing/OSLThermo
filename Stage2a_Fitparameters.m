%%%% STAGE 2a, FitParameters %%%%
%%%% Fits data with selected models, requires input from Stage 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Georgina King, 2015, modified 2022 georgina.king@unil.ch %
% Arnaud Duverger, 2015 arnaud.duverger@ens.fr %
% Chloé Bouscary, 2018 chloe.bouscary@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic %records time of execution

clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR; close all; clc; 

load(['./ComputeData/' filename '.mat']);   % loads data file from Stage1
nt = length(records);                       % number of signals (e.g. IRSL50, IRSL100 etc.)
global rhop10 isoT measL Hs                 % makes parameters available to the functions used when fitting the data


na = nSAR;                                  % number of aliquots 

% Compute vectors for fitting the data with the different models
fvec = logspace(0,8,100);                   % Creates a matrix with 100 columns, generates 100 points between 10^0 and 10^8
modvec = logspace(0,6,100);                 % Creates a matrix with 100 columns, generates 100 points between 10^0 and 10^6
tvec = logspace(0,6,100);                   % Creates a matrix with 100 columns, generates 100 points between 10^0 and 10^6
natdose = linspace(0,10000,10000);          % Creates a vector of natural dose time (s) 
tmat = tvec'*ones(1,NITL);                  % NITL = Number of isothermal holding temperatures %Creates a matrix with n columns and 100 lines

%Define constants
Ma = 1e6*365*24*3600;                       % 1 Ma in seconds
ka = Ma./1000;                              % 1 ka in seconds
kb = 8.617343*10^(-5);                      % Boltzmann constant eV/K
alpha = 1; Hs = 3e15;                       % values after Huntley (2006) J. Phys. D.
rprime = linspace(0.01,2.5,100)';           % creates vector of rprime distances between 0.01 and 2.5
pr = 3.*rprime.^2.*exp(-rprime.^3.);        % calculate probability of r'=rprime

for i=1:nt  % loop through the number of traps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fading %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract fading data from mat file
	ft = [records(i).rawdata(end-na+1:end).t];   
	fL = [records(i).rawdata(end-na+1:end).L];
	ok = isfinite(ft); fx = ft(ok); fy = fL(ok);    % remove NaN values

	% fit data to determine rhop using nlinfit function
	beta01 = [1 -5.5]';                             % initial parameters, scaling parameter, rhop estimate
	[beta1,R,J,Cov,MSE] = nlinfit(fx,fy,@FfitGK,beta01);   
    
	% caclulate confidence interval
	[fYpred,fdelta] = nlpredci(@FfitGK,fvec,beta1,R,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta1 = nlparci(beta1,R,'covar',Cov,'alpha',1-.68); 
	    
	% Extract rhop
	rhop10 = beta1(2); rhop10(rhop10<=-6.5)=-7;     % added for infinite fading rates
	srhop10 = [sbeta1(2) sbeta1(4)];
    
	%Calculate g2days for reporting only (not used in further calculations)
	beta02 = [1 -0.05]';                            %Initial parameters, I, m
	gy = fy./fy(1); gx = fx./172800;
	[beta2,R,J,Cov,MSE] = nlinfit(gx,gy,@Ffit2GK,beta02);
    
	%calculate confidence interval 
	[gYpred,gdelta] = nlpredci(@Ffit2GK,fvec,beta2,R,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68); 
	sigma2 = abs(beta2-sbeta2(:,2));

	%Extract g2days
	g2d = (-100*beta2(2));
	sg2d = g2d.*(sqrt((sigma2(2)/beta2(2))^2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Dose Response Curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract data from mat file
	labDdot = records(i).rawdata(1).Ddot;           % laboratory dose rate 
    natDdot = records(i).params(1).natDdot;         % environmental dose rate
	t = [records(i).rawdata(1:na).t];               % measurement times
	L = [records(i).rawdata(1:na).L];               % luminescence data
	ok = isfinite(t); x = t(ok); y = L(ok);         % remove NaN values%%

    if SAR_model==1 % Selects SSE fit
            
        % nlinfit solution
        beta0 = [1 1000]';                          % initial parameters, scaling parameter, D0 estimation
        [beta,R,J,Cov,MSE] = nlinfit(x,y,@FsarfitGK,beta0);
    
        % confidence interval
        [Ypred,delta] = nlpredci(@FsarfitGK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));               % calculate parameter uncertainties

    elseif SAR_model==2 % Selects GOK fit
                
        % nlinfit solution
        beta0 = [1 300 2]';                         % initial parameters, scaling parameter, D0 estimation, kinetic order
        [beta,R,J,Cov,MSE] = nlinfit(x,y,@FsarfitGK_GOK,beta0);
    
        % confidence interval
        [Ypred,delta] = nlpredci(@FsarfitGK_GOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));               % calculate parameter uncertainties
        
        GOK_a=beta(3); sGOK_a=sbeta(3);
        
    end
            
	% extract D0
	D0 = beta(2).*labDdot; sD0 = sigma(2).*labDdot;  % Unfaded D0
		
    % Calculate n/N
	Ln = L(~ok); rmNA = isfinite(Ln); Ln=Ln(rmNA); AvLn = mean(Ln); sLn = std(Ln);
	nNnat = AvLn./beta(1); snNnat = nNnat*sqrt(((sLn/AvLn)^2)+(((labDdot/100)/labDdot)^2)+((delta/Ypred)^2));

    % Calculate age from computed natural dose response curve  
    LnTn = AvLn./beta(1);
    rhop = 10.^beta1(2); rhop(rhop<0) = 0; 
    natDdotGyka = natDdot*ka;
    K = Hs*exp(-rhop^-(1/3)*rprime);
    
    if SAR_model==1
    
        simLabDRC = [];
        for j=1:length(natdose)
            for k=1:length(rprime)
                labDdotka = labDdot/(24*365.25*1e3);
                simLabDRC(k,j)=pr(k)*((labDdot(1)/D0)/(labDdot(1)/D0+K(k))...
                *(1-exp(-natdose(j)*(1/D0+K(k)/labDdot(1)))));
            end
        end
        
        ComputedLxTx=sum(simLabDRC)/sum(pr);
                    
        simDRC = [];
        for j=1:length(natdose)
            for k=1:length(rprime)
                simDRC(k,j)=pr(k)*((natDdot(1)/D0)/(natDdot(1)/D0+K(k))...
                *(1-exp(-natdose(j)*(1/D0+K(k)/natDdot(1)))));
            end
        end
        
        ComputedLxTx=sum(simDRC)/sum(pr);
        
        % Calcuate De and D0 from the simulated dose response curve
        [beta2,R,J,Cov,MSE] = nlinfit(natdose,ComputedLxTx,@SarfitGK,beta0); 
    
    elseif SAR_model==2
        
        simDRC=[];
        for j=1:length(natdose)
            for k=1:length(rprime)
                simDRC(k,j)=pr(k)*(natDdot(1)/D0)/(natDdot(1)/D0+K(k))...
                *(1-(1+(1/D0+K(k)/natDdot(1)) * natdose(j) * (GOK_a-1))^(-1/(GOK_a-1)));
            end
        end
        
        ComputedLxTx=sum(simDRC)/sum(pr);
    
        % Calcuate De and D0 from the simulated dose response curve
        [beta2,R,J,Cov,MSE] = nlinfit(natdose,ComputedLxTx,@SarfitGK_GOK,beta0); 
    
    end  
    
    sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68); 
    sigma2 = abs(beta2-sbeta2(:,2)); % calculate error on D0
    
    De = abs((beta2(2))*log(1-([nNnat nNnat-snNnat nNnat+snNnat]/beta2(1))));
    
    if SAR_model==1
    De_err = De(1)*sqrt(((((De(1)-De(2)+De(3)-De(1))/2)/De(1))^2)+((sLn/AvLn)^2)+(((labDdot/100)/labDdot)^2)+((sD0/D0)^2));
    elseif SAR_model==2
    De_err = De(1)*sqrt(((((De(1)-De(2)+De(3)-De(1))/2)/De(1))^2)+((sLn/AvLn)^2)+(((labDdot/100)/labDdot)^2)+((sD0/D0)^2));
    end
    
    %%%Save parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    records(i).NatFadedDe = [De(1) De_err]; records(i).NatFadedD0= beta2(2); records(i).sNatFadedD0=sigma2(2); 
    records(i).FCorrAge = [(De(1)/natDdotGyka(1)) ((De_err)/natDdotGyka(1))]; 
    records(i).minFCorrAge = (2*beta2(2))/natDdotGyka(1);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract data from mat file
	isoT = [records(i).rawdata(na+1:na+NITL).T];                            % create a vector with the different T for isothermal decay
    measL = [records(i).rawdata(na+1).t]; measL = measL(isfinite(measL));   % remove the NaN data from measL  
	itlt = [records(i).rawdata(na+1:na+NITL).t];                            % create a matrix with all the time records of the isothermal decay
    itlL = [records(i).rawdata(na+1:na+NITL).L];                            % create a matrix with all the luminescence records of the isothermal decay
    ok = isfinite(itlt); itlx = itlt(ok); itly = itlL(ok); itlt=itlx ; itlL=itly; % remove the NaN data from the data
    itlsort = NaN(length(measL),length(isoT));                              % create matrix to facilitate data output storage
  
if ITL_model==1; % GOK model
    
    % nlinfit solution 
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 4 [Avec]]';     % initial estimates of s, Et, b  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FGOK_GK,beta02);
	
    % Confidence intervals %
	measL = tvec;                   % redefine time dimension for model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FGOK_GK,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params

elseif ITL_model==2; % GAUSS model 
   
    % nlinfit solution
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 0.1 [Avec]]';   % initial estimates of s, mu(Et), s(Et),  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FLiLiGauss,beta02);
	
    % Confidence intervals %
	measL = tvec;                   % redefine time dimension for model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FLiLiGauss,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params 

elseif ITL_model==3; % BTS model
    
	% nlinfit solution 
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 0.1 [Avec]]';   % initial estimates of s, Et, Eu,  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FLiLiGK,beta02);
	
    % Confidence intervals
	measL = tvec;                   % redefine time dimension for the model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FLiLiGK,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params  
    
end

    % extract parameters
	s10 = params(1); Et = params(2); Eu = params(3); A = params(4:end);
	ss10 = sigma3(1); sEt = sigma3(2); sEu = sigma3(3);  % Uncertainties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% Calculate (n/N)ss using Kars et al (2008) %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    natDdot = records(i).params.natDdot(1)*3600*24*365*1e3; % convert from Gy-s/Ma to Gy/ka
	rhop = 10.^rhop10; rhop(rhop<0) = 0;
	rhopU = 10.^srhop10(2); rhopU(rhopU<0) = 0;
	rhopL = 10.^srhop10(1); rhopL(rhopL<0) = 0;
	rhopbound = [rhop rhopU rhopL];
	rho = 3*alpha^3*rhopbound/(4*3.1415);
	r = rprime*(1./(4*3.1415*rho/3).^(1/3));
	tau = ((1/Hs)*2.71.^(alpha.*r))/(3600*24*365*1000);
	Ls = 1./(1+D0./(natDdot.*tau));
	Lstrap = (pr'*Ls)/sum(pr);

	%%%Save parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	records(i).params.rhop10 = [rhop10,srhop10];
	records(i).params.g2d = [g2d,sg2d];
	records(i).params.D0 = [D0,sD0];
    
    if SAR_model==2; records(i).params.GOK_a = [GOK_a, sGOK_a]; else records(i).params.GOK_a = [NaN, NaN]; end
	
    records(i).params.Et = [Et,sEt];
	
    if ITL_model==1; records(i).params.GOK_b = [Eu, sEu]; else records(i).params.GOK_b = [NaN, NaN]; end
    if ITL_model==2; records(i).params.sigmaEt = [Eu,sEu]; else records(i).params.sigmaEt = [NaN,NaN]; end
    if ITL_model==3; records(i).params.Eu = [Eu,sEu]; else records(i).params.Eu = [NaN,NaN];  end
    	
    records(i).params.s10 = [s10,ss10];

	%%%Save data for plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	records(i).plot.xf = fvec;
	records(i).plot.yf = fYpred;
	records(i).plot.syf = fdelta;
	records(i).plot.tf = ft;
	records(i).plot.Lf = fL;

	records(i).plot.x = modvec;
	records(i).plot.y = Ypred./beta(1);
	records(i).plot.sy = delta./beta(1);
	records(i).plot.t = t;
	records(i).plot.L = L./beta(1);

	records(i).plot.xitl = tmat;
	records(i).plot.yitl = reshape(IsoPred,size(tmat));
	records(i).plot.syitl = reshape(delta2,size(tmat));
	itlsort(1:length(itlt)) = itlt;
	records(i).plot.titl = itlsort;
	itlsort(1:length(itlL)) = itlL;
	records(i).plot.Litl = itlsort;
	records(i).plot.scale = A;

	records(i).nNnat = [nNnat,snNnat];
	records(i).KarsresAv = [Lstrap(1) max(abs(Lstrap(1)-Lstrap(2)),abs(Lstrap(1)-Lstrap(3)))];
end

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.mat'],'records')
toc
