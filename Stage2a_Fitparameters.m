%%%% STAGE 2a, FitParameters %%%%
%%%% Fits data with selected models, requires input from Stage 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Georgina King, 2015, modified 2022 georgina.king@unil.ch %
% Arnaud Duverger, 2015 arnaud.duverger@ens.fr %
% Chloé Bouscary, 2018 chloe.bouscary@unil.ch %
% Maxime Bernard, 2025 maxime.bernard@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary of assumptions made in this code:
%
% 1- We fit athermal fading data taking Hs = 3e15 (Huntley, 2006)  
%
% 2- We first fit isothermal decay curves assuming athermal fading
%
% 3- We fit DRC curves assuming athermal fading occurs during lab irradiation.
% athermal fading linked to highly instable traps of few seconds to hours lifetime,
% see Kars et al. (2008). Correction is IRSL_corr = IRSL_init * exp(-rhop*(ln(1.8*Hs*t)^3)) 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic %records time of execution

clearvars -except filename filenamevec NITL NITLvec SAR_model ITL_model SAR_MODEL ITL_MODEL nSAR; close all; clc; 

load(['./ComputeData/' filename '.mat']);   % loads data file from Stage1
nt = length(records);                       % number of signals (e.g. IRSL50, IRSL100 etc.)
global rhop10 isoT measL Hs                 % makes parameters available to the functions used when fitting the data
global nFa
global na
global index_aliquot
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
alpha = 1;                      % values after Huntley (2006) J. Phys. D.
rprime = linspace(0.01,2.5,100)';           % creates vector of rprime distances between 0.01 and 2.5 (Huntley, 2006)
pr = 3.*rprime.^2.*exp(-rprime.^3.);        % calculate probability of r'=rprime (Kars et al, 2008)

for i=1:nt  % loop through the number of traps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fading %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Record s_ath = s_th
    Hs = 3e15;
    % Extract fading data from mat file
	ft = [records(i).rawdata(end-nFa+1:end).t];   
	fL = [records(i).rawdata(end-nFa+1:end).L];
    index_aliquot = [];
    for k=1:nFa
        index_aliquot = [index_aliquot,k.*ones(size(records(i).rawdata(end-k-1).t))]; % Give each aliquot an id
        fL(index_aliquot==k) = fL(index_aliquot==k)./max(fL(index_aliquot==k));
    end
    ok = isfinite(ft); fx = ft(ok); fy = fL(ok);    % remove NaN values
    index_aliquot = index_aliquot(ok); 

	% fit data to determine rhop using nlinfit function
    A_aliquots = ones(1,nFa); % initial guesses on pre-exponential factor for each aliquot
	beta01 = [-5.5 A_aliquots]';     % initial parameters, scaling parameter, rhop estimate
	[beta1,R,J,Cov,MSE] = nlinfit(fx,fy,@FfitGK,beta01);   

	% calculate confidence interval
	[fYpred,fdelta] = nlpredci(@FfitGK,fvec,beta1,R,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta1 = nlparci(beta1,R,'covar',Cov,'alpha',1-.68); 

    % Estimate 2 days fading
    [fYpred_2days,fdelta_2days] = nlpredci(@FfitGK,172800,beta1,R,'covar',Cov,'alpha',.05,'predopt','curve');
	    
	% Extract rhop
	rhop10 = beta1(1); rhop10(rhop10<=-7)=-7;     % added for infinite fading rates
	srhop10 = [sbeta1(1,1) sbeta1(1,2)];
    
	%Calculate g2days for reporting only (not used in further calculations)
	beta02 = [1 -0.05]';             %Initial parameters, I, m
    a_factors = beta1(2:end)'; 
	gy = fy./a_factors(index_aliquot); gx = fx./(2*24*3600);
	[beta2,R,J,Cov,MSE] = nlinfit(gx,gy,@Ffit2GK,beta02);
    
	%calculate confidence interval 
	[gYpred,gdelta] = nlpredci(@Ffit2GK,fvec,beta2,R,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68); 
	sigma2 = abs(beta2-sbeta2(:,2));

	%Extract g2days
	g2d = (-100*beta2(2));
	sg2d = g2d.*(sqrt((sigma2(2)/beta2(2))^2));

    rhop = 10.^beta1(1); rhop(rhop<0) = 0; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Dose Response Curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract data from mat file
    global labDdot natDdot
    natDdot = records(i).params(1).natDdot;         % environmental dose rate
	t = [records(i).rawdata(1:na).t];               % measurement times
	L = [records(i).rawdata(1:na).L];               % luminescence data
    labDdot = zeros(na,1);                          % Initiate laboratory dose rates
    index_aliquot_raw = [];                         % Initiate index aliquots
    for k=1:na
	    labDdot(k) = records(i).rawdata(k).Ddot;           % laboratory dose rate used for each aliquot 
        index_aliquot_raw = [index_aliquot_raw,k.*ones(size(records(i).rawdata(k).t))]; % Give each aliquot an id
    end
	ok = isfinite(t); x = t(ok)'; y = L(ok)'; index_aliquot = index_aliquot_raw(ok);        % remove NaN values%%

    A_aliquots = 5.*ones(1,na); % initial guesses on pre-exponential factor for each aliquot

    if SAR_model==1 % Selects SSE fit
            
        % nlinfit solution
       beta0 = [300 A_aliquots]';                          % initial parameters, D0 estimation, scaling parameters 
       [beta,R,J,Cov,MSE] = nlinfit(x,y,@FsarfitGK,beta0);
    
        % confidence interval
       [Ypred,delta] = nlpredci(@FsarfitGK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
       
       sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
       sigma = abs(beta-sbeta(:,2));               % calculate parameter uncertainties

       % 	extract D0
        D0 = beta(1); sD0 = sigma(1);  % Unfaded D0
        A_factors = beta(2:end); 

        % Calculate faded Ln/Llabmax
	    Ln = L(~ok); rmNA = isfinite(Ln); Ln2=Ln(rmNA)./beta(2:end)'; 
        AvLn = mean(Ln2); sLn = std(Ln2);

        % save for residuals
        A_f = A_factors(index_aliquot_raw);
        records(i).Residuals.residuals = [x,y./exp(-rhop*log(1.8*Hs.*(0.5*t(ok)')).^3)./A_f(ok),R];
        records(i).Residuals.plot = [modvec;Ypred./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors);delta./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors)]';
        records(i).Residuals.param = [beta(2),D0,sD0,0,0];

    elseif SAR_model==2 % Selects GOK fit
                
         % nlinfit solution
        beta0 = [300 0.5 A_aliquots]';   % initial parameters, D0 estimation, kinetic order, scaling parameter for each aliquot
        [beta,R,J,Cov,MSE] = nlinfit(x,y,@FsarfitGK_GOK,beta0);
        beta=real(beta); R=real(R); Cov=real(Cov);
        % confidence interval
        [Ypred,delta] = nlpredci(@FsarfitGK_GOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));               % calculate parameter uncertainties
        GOK_a=beta(2)+1; sGOK_a=sigma(2);
        
        % 	extract D0  
        D0 = beta(1); sD0 = sigma(1); 
        A_factors = beta(3:end); 

        % save for residuals
        A_f = A_factors(index_aliquot_raw);
        records(i).Residuals.residuals = [x,y./exp(-rhop*log(1.8*Hs.*(0.5*t(ok)')).^3)./A_f(ok),R];
        records(i).Residuals.plot = [modvec;Ypred./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors);delta./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors)]';
        records(i).Residuals.param = [beta(1),D0,sD0,beta(3),0];

        % Calculate faded Ln/Llabmax
	    Ln = L(~ok); rmNA = isfinite(Ln); Ln2=Ln(rmNA)./beta(3:end)'; 
        AvLn = mean(Ln2); sLn = std(Ln2);

     elseif SAR_model==3 % Selects GAUSS fit
         % nlinfit solution
        beta0=[2 0.1 A_aliquots]'; % initial parameters, scaling parameter, D0 estimation, sigma D0 estimation
        [beta,R,J,Cov,MSE] = nlinfit(x',y',@Fsarfit_GAUSS,beta0);
        beta=real(beta); R=real(R); Cov=real(Cov);

        % confidence interval
        [Ypred,delta] = nlpredci(@Fsarfit_GAUSS,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));               % calculate parameter uncertainties

        % 	extract D0
	    D0 = 10.^(beta(1)); sD0 = (log(10)*10^(beta(1))*sigma(1));  % Unfaded D0, sD0 = CI on D0
	    sigmaD0 = (log(10)*10.^(beta(1))*beta(2)); ssigmaD0 = log(10)*10^(beta(2))*sigma(2); % sigmaD0 = rnage on D0, ssigmaD0 = CI on sigmaD0
        A_factors = beta(3:end);

        % Calculate faded Ln/Llabmax
	    Ln = L(~ok); rmNA = isfinite(Ln); Ln2=Ln(rmNA)./beta(3:end)'; 
        AvLn = mean(Ln2); sLn = std(Ln2);

        % save for residuals
        A_f = A_factors(index_aliquot_raw);
        records(i).Residuals.residuals = [x,y./exp(-rhop*log(1.8*Hs.*(0.5*t(ok)')).^3)./A_f(ok),R];
        records(i).Residuals.plot = [modvec;Ypred./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors);delta./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors)]';
        records(i).Residuals.param = [beta(1),D0,sD0,sigmaD0,ssigmaD0,beta(2),beta(3),sigma(2),sigma(3)];
    end
     
    % Unil - Calculate (n/N)nat = (L/Llabmax)
    nNnat = AvLn;
    snNnat = nNnat*sqrt(((sLn/AvLn)^2)+(((labDdot(1)/100)/labDdot(1))^2)+((delta/Ypred)^2));


    %%%%%% Calculate age from computed natural dose response curve %%%%%%%% 
    natDdotGyka = natDdot*ka;
    K = Hs*exp(-rhop^-(1/3)*rprime); % s_tun fading
    
    %%%% Simulate natural dose response curve %%%%
    %%% SSE %%%
    if SAR_model==1 
    
        % Compute (n/N)ss nat
        nss_nat = zeros(length(rprime),3);
        for ii=1:length(rprime) % loop over electron-hole distance
            m1=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-rhop.^(-1./3).*rprime(ii)).*x).^2;
            m2=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-(10^srhop10(2)).^(-1./3).*rprime(ii)).*x).^2;
            m3=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-(10^srhop10(1)).^(-1./3).*rprime(ii)).*x).^2;
            %Hs should be replaced by s from high-temperature best-fit
            nss_nat(ii,1)=fminbnd(m1,0,1).*pr(ii);
            nss_nat(ii,2)=fminbnd(m2,0,1).*pr(ii); % Upper boundary
            nss_nat(ii,3)=fminbnd(m3,0,1).*pr(ii); % lower boundary
        end
        fadfact_nat =sum(nss_nat)./sum(pr);
    
        simLabDRC = [];
        for j=1:length(natdose)
            % Loop through electron-hole distances
            for k=1:length(rprime)
                labDdotka = labDdot/(24*365.25*1e3);
                simLabDRC(k,j)=pr(k)*((labDdot(1)/D0)/(labDdot(1)/D0+K(k))...
                *(1-exp(-natdose(j)*(1/D0+K(k)/labDdot(1)))));
            end
        end
        
        ComputedLxTx=sum(simLabDRC)/sum(pr);
                    
        simDRC = [];
        for j=1:length(natdose)
            % Loop through electron-hole distances
            for k=1:length(rprime)
                simDRC(k,j)=pr(k)*((natDdot(1)/D0)/(natDdot(1)/D0+K(k))...
                *(1-exp(-natdose(j)*(1/D0+K(k)/natDdot(1)))));
            end
        end
        
        ComputedLxTx=sum(simDRC)/sum(pr);
        
        % Calcuate De and D0 from the simulated dose response curve
        [beta2,R,J,Cov,MSE] = nlinfit(natdose,ComputedLxTx,@SarfitGK,beta0(1:2)); 
        sigmaD0 = 0; ssigmaD0 = 0; a_nat = beta2(2); D0_nat = beta2(1);

    %%% GOK %%% Include fading
    elseif SAR_model==2 

        % Compute (n/N)ss nat
        nss_nat = zeros(length(rprime),3);
        for ii=1:length(rprime) % loop over electron-hole distance
            m1=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-rhop.^(-1./3).*rprime(ii)).*x).^2;
            m2=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-(10^srhop10(2)).^(-1./3).*rprime(ii)).*x).^2;
            m3=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-(10^srhop10(1)).^(-1./3).*rprime(ii)).*x).^2;
            %Hs should be replaced by s from high-temperature best-fit
            nss_nat(ii,1)=fminbnd(m1,0,1).*pr(ii);
            nss_nat(ii,2)=fminbnd(m2,0,1).*pr(ii); % Upper boundary
            nss_nat(ii,3)=fminbnd(m3,0,1).*pr(ii); % lower boundary
        end
        fadfact_nat =sum(nss_nat)./sum(pr);

        simDRC=[];
        for j=1:length(natdose)
            % Loop through electron-hole distances
            for k=1:length(rprime)
                simDRC(k,j)=pr(k)*(natDdot(1)/D0)/(natDdot(1)/D0+K(k))...
                *(1-(1+(1/D0+K(k)/natDdot(1)) * natdose(j) * (GOK_a))^(-1/(GOK_a)));
            end
        end
        
        ComputedLxTx=sum(simDRC)/sum(pr);  

        % Calcuate De and D0 from the simulated dose response curve
        beta0(3) = 1;
        [beta2,R,J,Cov,MSE] = nlinfit(natdose,ComputedLxTx,@SarfitGK_GOK,beta0(1:3)); 
        beta2 = real(beta2); R=real(R); Cov=real(Cov);
        sigmaD0 = 0; ssigmaD0 = 0;
        a_nat = beta2(3); D0_nat = beta2(1);
     
    %%% GAUSS %%%
    elseif SAR_model==3 

        D0s_distribution=linspace(beta(1) - 3.*beta(2), beta(1) + 3 .* beta(2), 30)'; % Distribution of D0s      
        pd=normpdf(D0s_distribution, beta(1), beta(2)); pd= pd ./ sum(pd);
        D0s_distribution = 10.^D0s_distribution;

         % Compute (n/N)ss nat
        nss_nat = zeros(length(rprime),3);
        nss_temp = zeros(length(D0s_distribution),3);
        for ii=1:length(rprime) % loop over electron-hole distance
            for k=1:length(D0s_distribution) %loop over D0s
                m1=@(x)(natDdot(1)./D0s_distribution(k).*(1-x)-Hs.*exp(-rhop.^(-1./3).*rprime(ii)).*x).^2;
                m2=@(x)(natDdot(1)./D0s_distribution(k).*(1-x)-Hs.*exp(-(10^srhop10(2)).^(-1./3).*rprime(ii)).*x).^2;
                m3=@(x)(natDdot(1)./D0s_distribution(k).*(1-x)-Hs.*exp(-(10^srhop10(1)).^(-1./3).*rprime(ii)).*x).^2;
                %Hs should be replaced by s from high-temperature best-fit
                nss_temp(k,1)=fminbnd(m1,0,1).*pd(k).*pr(ii);
                nss_temp(k,2)=fminbnd(m2,0,1).*pd(k).*pr(ii); % Upper boundary
                nss_temp(k,3)=fminbnd(m3,0,1).*pd(k).*pr(ii); % lower boundary
            end
            nss_nat(ii,1) = sum(nss_temp(:,1));
            nss_nat(ii,2) = sum(nss_temp(:,2));
            nss_nat(ii,3) = sum(nss_temp(:,3));
        end
        fadfact_nat =sum(nss_nat)./sum(pr);

        [RD, NatDose, D0s] = ndgrid(rprime, natdose, D0s_distribution); % 3D grid
        PD = reshape(pd, 1, 1, []); 
        
        % Vectorization 
        natDdot_D0 = natDdot(1) ./ D0s; 
        decay_term = natDdot_D0 ./ (natDdot_D0 + K) .* (1 - exp(-NatDose .* (1 ./ D0s + K ./ natDdot(1))));
        weighted_decay = PD .* decay_term; 
        
        % Sum over rprimes
        simDRC = pr .* squeeze(sum(weighted_decay, 3));
        ComputedLxTx = sum(simDRC) ./ sum(pr);
        
        % Calcuate De and D0 from the simulated dose response curve
        beta0(3) = 1;
        [beta2,R,J,Cov,MSE] = nlinfit(natdose,ComputedLxTx,@Sarfit_GAUSS,beta0(1:3)); 
        beta2(1) = 10^beta2(1); beta2(2) = (log(10)*10.^(beta2(1))*beta2(2));
        a_nat = beta2(3); D0_nat = beta2(1);
    end  
    
    sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68); 
    sigma2 = abs(beta2-sbeta2(:,2)); % calculate error on D0
    
    De = abs(D0_nat*log(1-([nNnat nNnat-snNnat nNnat+snNnat]/a_nat)));

    if SAR_model==1 % SSE
    De_err = De(1)*sqrt(((((De(1)-De(2)+De(3)-De(1))/2)/De(1))^2)+((sLn/AvLn)^2)+(((labDdot(1)/100)/labDdot(1))^2)+((sD0/D0)^2));
    elseif SAR_model==2 % GOK
    De_err = De(1)*sqrt(((((De(1)-De(2)+De(3)-De(1))/2)/De(1))^2)+((sLn/AvLn)^2)+( ( (labDdot(1)/100) / labDdot(1))^2)+((sD0/D0)^2));
    elseif SAR_model == 3 % GAUSS
    De_err = De(1)*sqrt(((((De(1)-De(2)+De(3)-De(1))/2)/De(1))^2)+((sLn/AvLn)^2)+(((labDdot(1)/100)/labDdot(1))^2)+((sD0/D0)^2)+ ((ssigmaD0/sigmaD0)^2) );
    end
    
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
  
    if ITL_model==1 % GOK model
    
    % nlinfit solution 
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 4 [Avec]]';     % initial estimates of s, Et, b  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FGOK_GK,beta02);
    params = real(params); resn = real(resn); Cov = real(Cov);
	
    % Confidence intervals %
	measL = tvec;                   % redefine time dimension for model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FGOK_GK,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params

     records(i).Residuals.residualsIso = [itlx',itly',resn];

    elseif ITL_model==2 % GAUSS model 
   
    % nlinfit solution
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 0.1 [Avec]]';   % initial estimates of s, mu(Et), s(Et),  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FLiLiGauss,beta02);
	
    % Confidence intervals %
	measL = tvec;                   % redefine time dimension for model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FLiLiGauss,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params 

    records(i).Residuals.residualsIso = [itlx',itly',resn];

    elseif ITL_model==3 % BTS model
    
	% nlinfit solution 
    Avec=ones(1,NITL);              % create a matrix of ones (one row, NITL columns)
	beta02 = [9 1.4 0.1 [Avec]]';   % initial estimates of s, Et, Eu,  
    [params,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@FLiLiGK,beta02);
	
    % Confidence intervals
	measL = tvec;                   % redefine time dimension for the model, 1 line, 100 columns
	[IsoPred,delta2] = nlpredci(@FLiLiGK,tmat(:),params,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta3 = nlparci(params,resn,'covar',Cov,'alpha',1-.68);
	sigma3 = abs(params-sbeta3(:,2)); % calculate error on params  

    records(i).Residuals.residualsIso = [itlx',itly',resn];
    end

    % extract parameters
	s10 = params(1); Et = params(2); Eu = params(3); A = params(4:end);
	ss10 = sigma3(1); sEt = sigma3(2); sEu = sigma3(3);  % Uncertainties
    
    %%% Save parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    records(i).NatFadedDe = [De(1) De_err]; 
    records(i).NatFadedD0= [beta2(1) sigma2(1)];
    if SAR_model == 3 % Gauss
        records(i).NatFadedsigmaD0= [beta2(2) sigma2(2)]; 
    else 
        records(i).NatFadedsigmaD0= [0.0 0.0];
    end
    records(i).FCorrAge = [(De(1)/natDdotGyka(1)) ((De_err)/natDdotGyka(1))]; 
    records(i).minFCorrAge = (2*beta2(1))/natDdotGyka(1);
    records(i).params.imax = mean(beta(3:end));   

    
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
	tau = ((1/Hs)*exp((alpha.*r)))/(3600*24*365*1000);
	Ls = 1./(1+D0./(natDdot.*tau));
	Lstrap = (pr'*Ls)/sum(pr);

	%%%Save parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	records(i).params.rhop10 = [rhop10,srhop10];
	records(i).params.g2d = [g2d,sg2d];
	records(i).params.D0 = [D0, sD0];
    records(i).params.sigmaD0 = [sigmaD0, ssigmaD0];
    
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
    records(i).plot.y = Ypred./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors); % Correct for athermal fading
	records(i).plot.sy = delta./exp(-rhop*log(1.8*Hs.*(0.5*modvec)).^3)./mean(A_factors);% Correct for athermal fading
	records(i).plot.t = t;
    ok = isfinite(t); L(ok) = L(ok)./exp(-rhop*log(1.8*Hs.*(0.5*t(ok))).^3);
    records(i).plot.L = L./A_factors(index_aliquot_raw)';

	records(i).plot.xitl = tmat;
	records(i).plot.yitl = reshape(IsoPred,size(tmat));
	records(i).plot.syitl = reshape(delta2,size(tmat));
	itlsort(1:length(itlt)) = itlt;
	records(i).plot.titl = itlsort;
	itlsort(1:length(itlL)) = itlL;
	records(i).plot.Litl = itlsort;
	records(i).plot.scale = A;

	records(i).nNnat = [nNnat,snNnat];
    records(i).KarsresAv = [fadfact_nat(1) max(abs(fadfact_nat(1)-fadfact_nat(3)),abs(fadfact_nat(1)-fadfact_nat(2)))];

end

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.mat'],'records')
toc

% save kinetic parameter values in excel file
records2 = records;
mat = zeros(nt,23); % There are 23 values
samples = strings(nt,1);
nNnat_record = zeros(nt,2);
headers = {'Sample Temp (°C)','Unc. Temp','Dr (Gy/kyr)','Unc. Dr','imax','rhop','Unc. rhop-','Unc. rhop+','g2d (%/decade)','Unc. g2d','D0 (Gy)','Unc. D0',...
    'sigmaD0','Unc. sD0','GOK_a','Unc. a','Et (eV)','Unc. Et','GOK_b','Unc. b','sigmaEt (eV)','Unc. sigEt','Eu (eV)','Unc. Eu','s10','Unc. s10','Age (ka)','Unc. age',...
    'minAge (kyr)','De (Gy)','Unc. De','NatFadedD0 (Gy)','Unc. NatD0 (Gy)','NatFadedSigmaD0 (Gy)','Unc. NatSigmaD0 (Gy)'};
for i=1:nt
    records2(i).params.natDdot = records(i).params.natDdot .*(1000.*3600.*24.*365);
    records2(i).params.Age = records2(i).FCorrAge;
    records2(i).params.minAge = records2(i).minFCorrAge;
    records2(i).params.De = records2(i).NatFadedDe;
    records2(i).params.NatD0 = records2(i).NatFadedD0;
    records2(i).params.NatsigmaD0 = records2(i).NatFadedsigmaD0;
    row(i,:) = table2array(struct2table(records2(i).params));
    samples(i) = string(records2(i).id);
    nNnat_record(i,:) = records2(i).nNnat;
end
tab = array2table(row,'VariableNames',headers);
tab2 = array2table(samples,'VariableNames',{'Samples'});
tab3 = array2table(nNnat_record,'VariableNames',{'n/N_nat','n/N err'});
tab4 = [tab2,tab,tab3];
writetable(tab4,['./ComputeData/' filename '_' SAR_MODEL '_' ITL_MODEL '_fitpar.xls'])
toc
