%%%% STAGE 2a, FitParameters %%%%
%%%% Fits data with selected models, requires input from Stage 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Georgina King, 2022       georgina.king@unil.ch %
% Arnaud Duverger, 2015  	arnaud.duverger@ens.fr %
% Chloé Bouscary, 2018  	chloebouscary@gmail.com %
% Maxime Bernard, 2025      maxime.bernard@unil.ch %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Summary of assumptions made in this code:
%
% 1- Athermal fading data are fitted using Hs = 3e15 (Huntley, 2006).
%
% 2- Isothermal decay curves are fitted accounting for athermal fading.
%
% 3- DRC curves are fitted assuming that athermal fading occurs during lab irradiation.
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec nFAD nFADvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '.mat']);                                   % loads data file from Stage 1
nt = length(records);                                                       % number of signals (e.g. IRSL50, IRSL100 etc.)

%%% Define global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rhop10 isoT measL Hs index_aliquotFAD index_aliquotSAR labDdot natDdot Ddot % makes parameters available to the functions used when fitting the data

%%% Compute vectors for fitting the data with the different models
fvec    =   logspace(0,8,100);                                              % creates a vector with 100 columns, generates 100 points between 10^0 and 10^8
modvec  =   logspace(0,6,100);                                              % creates a vector with 100 columns, generates 100 points between 10^0 and 10^6
tvec    =   logspace(0,6,100);                                              % creates a vector with 100 columns, generates 100 points between 10^0 and 10^6
natdose =   linspace(0,10000,10000);                                        % creates a vector of natural dose time [s]
tmat    =   tvec'*ones(1,nITH);                                             % creates a matrix with nITH columns (nb of isothermal holding temperatures) and 100 lines, with the values of tvec

%%% Define constants
ka = 1e3*365.25*24*3600;                                                    % 1 [kyr] in [s]
alpha = 1; Hs = 3e15;                                                       % alpha = constant, Hs=s = attempt-to-escape frequency [s-1] (Huntley, 2006, J. Phys.)
rprime = linspace(0.00,2.5,100)';                                           % creates vector of rprime, dimensionless distances between 0.01 and 2.5 (Huntley, 2006, J. Phys.)
pr = 3.*rprime.^2.*exp(-rprime.^3);                                         % calculates p(r'), the probability that the nearest recombination center is at a distance r'=rprime (Kars et al., 2008, R. M.)

%%% Fit the data with the different models %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nt                                                                  % loops through the number of traps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Athermal fading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Extract fading data from .mat file
    ft = [records(i).rawdata(end-nFAD+1:end).t];                            % athermal fading times [s]
    fL = [records(i).rawdata(end-nFAD+1:end).L];                            % luminescence signal intensity (a.u.), normalised

    %%% Normalise the data of each aliquot individualy
    index_aliquot = [];
    for k=1:nFAD
        index_aliquot = [index_aliquot,k.*ones(size(records(i).rawdata(end-k-1).t))]; % gives each aliquot an id
        fL(index_aliquot==k) = fL(index_aliquot==k)./max(fL(index_aliquot==k));
    end
    ok = isfinite(ft); fx = ft(ok); fy = fL(ok);                            % removes NaN values
    index_aliquotFAD = index_aliquot(ok);
    [~,~,index_aliquotFAD] = unique(index_aliquotFAD,'stable');
    index_aliquotFAD = index_aliquotFAD';

    %%% Fit data to determine rhop using nlinfit function
    A_aliquotFAD = ones(1,length(unique(index_aliquotFAD)));                % initial guesses on pre-exponential factor (scaling parameter A) for each aliquot
    beta01 = [-5.5 A_aliquotFAD]';                                          % initial guesses on parameters: rhop estimate, scaling parameters for each aliquot
    [beta1,R,~,Cov,~] = nlinfit(fx,fy,@Ffitrhop,beta01);

    %%% Calculate confidence interval
    [fYpred,fdelta] = nlpredci(@Ffitrhop,fvec,beta1,R,'covar',Cov,'alpha',.05,'predopt','curve'); % fdelta are 2-sigma (95%) error estimates on predicted data fYpred
    sbeta1 = nlparci(beta1,R,'covar',Cov,'alpha',1-.68);                    % extract lower and upper values for params, 1-sigma error estimates (68%)
    sigma1 = [abs(beta1-sbeta1(:,1)),abs(beta1-sbeta1(:,2))];               % calculates parameter lower and upper uncertainties

    %%% Extract rhop
    rhop10 = beta1(1);
    rhop10_l = sbeta1(1,1); rhop10_u = sbeta1(1,2);                         % lower and upper values of rhop10
    srhop10 = mean(sigma1(1,:));                                            % uncertainty on rhop10
    rhop10(rhop10<-7)=-7; rhop10_l(rhop10_l<-7)=-7; rhop10_u(rhop10_u<-7)=-7; srhop10(rhop10<-7)=NaN; % added for infinite fading rates
    rhop = 10.^rhop10; rhop(rhop<0) = 0;

    %%% Calculate g2days (not used in further calculations)
    beta02 = [1 -0.05]';                                                    % initial guesses on parameters: I, m
    a_factors = beta1(2:end)';                                              % scalling parameters for each aliquot
    gy = fy./a_factors(index_aliquotFAD); gx = fx./(2*24*3600);             % data corrected for the scalling parameters ; irradiation times scaled for 2 days in [s]
    [beta2,R,~,Cov,~] = nlinfit(gx,gy,@Ffitg2d,beta02);

    %%% Calculate confidence interval
    [gYpred,gdelta] = nlpredci(@Ffitg2d,fvec,beta2,R,'covar',Cov,'alpha',.05,'predopt','curve'); % gdelta are 2-sigma estimates (95%) error estimates on predicted data gYpred
    sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68);                    % extract lower and upper values for params, 1-sigma error estimates (68%)
    sigma2 = [abs(beta2-sbeta2(:,1)),abs(beta2-sbeta2(:,2))];               % calculates parameter lower and upper uncertainties

    %%% Extract g2days [%/dec.]
    g2d = -100*beta2(2);                                                    % put g2d in positive percent [%/dec.]
    sg2d = 100*abs(mean(sigma2(2,:)));                                      % uncertainty on g2d [%/dec.]


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dose Response Curve, regenerative dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Extract data from .mat file
    t = [records(i).rawdata(1:nSAR).t];                                     % doses = irradiation time [s]
    L = [records(i).rawdata(1:nSAR).L];                                     % luminescence signal intensity (a.u.), not normalised
    natDdot = records(i).params.natDdot;                                    % environmental dose rate [Gy/s]

    %%% Convert the regenerative doses [s] to [Gy] of each aliquot individualy
    index_aliquot = [];
    for k=1:nSAR
        labDdot(k) = records(i).rawdata(k).labDdot;                         % laboratory dose rate used for each aliquot [Gy/s]
        index_aliquot = [index_aliquot,k.*ones(size(records(i).rawdata(k).t))]; % gives each aliquot an id
    end
    ok = isfinite(t); x = t(ok)'; y = L(ok)';                               % removes NaN values
    index_aliquotSAR = index_aliquot(ok);
    ia = unique(index_aliquotSAR);                                          % to account for no data (NaN) for one or more aliquots
    Ddot = labDdot(ia);
    if length(index_aliquotSAR) < length(index_aliquot)
        index_aliquotSAR = index_aliquot(1:length(index_aliquotSAR));
    end
    [~,~,index_aliquotSAR] = unique(index_aliquotSAR,'stable');
    index_aliquotSAR = index_aliquotSAR';

    %%% Precompute Ln before fading
    Ln = L(~ok); Ln=Ln(~isnan(Ln)); rmNA = isfinite(Ln); Ln = Ln(rmNA);

    %%% Single saturating exponential (SSE) fit %%%
    if SAR_fittype==1 % SSE fit

        %%% nlinfit solution
        A_aliquotSAR = 5.*ones(1,length(Ln));                               % initial guesses on pre-exponential factor (scaling parameter A) for each aliquot
        beta0 = [300 A_aliquotSAR]';                                        % initial guesses on parameters: D0 estimation in [Gy], scaling parameters for each aliquot
        [beta,R,~,Cov,~] = nlinfit(x,y,@FSARfit_SSE,beta0);

        %%% Confidence intervals
        [Ypred,delta] = nlpredci(@FSARfit_SSE,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve'); % delta are 2-sigma estimates (95%) error estimates on predicted data Ypred
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);                  % extract lower and upper values for params, 1-sigma error estimates (68%)
        sigma = [abs(beta-sbeta(:,1)),abs(beta-sbeta(:,2))];                % calculates parameter lower and upper uncertainties

        %%% Extract D0 [Gy]
        D0 = beta(1);                                                       % unfaded D0 in [Gy]
        sD0 = mean(sigma(1,:));                                             % sD0 = CI on D0 [Gy]

        %%% Extract scalling parameters per aliquot
        A_factors = beta(2:end);

        %%% Calculate faded Ln/Llabmax
        Ln = Ln./A_factors';
        AvLn = mean(Ln);
        sLn = std(Ln); sLn(sLn==0) = 0.01*AvLn;

    %%% General order kinetics (GOK) fit %%%
    elseif SAR_fittype==2 % GOK fit

        %%% nlinfit solution
        A_aliquotSAR = 5.*ones(1,length(Ln));                               % initial guesses on pre-exponential factor (scaling parameter A) for each aliquot
        beta0 = [300 2 A_aliquotSAR]';                                      % initial guesses on parameters: D0 estimation in [Gy], GOK growth order, scaling parameters for each aliquot
        [beta,R,~,Cov,~] = nlinfit(x,y,@FSARfit_GOK,beta0);
        beta=real(beta); R=real(R); Cov=real(Cov);

        %%% Confidence intervals
        [Ypred,delta] = nlpredci(@FSARfit_GOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve'); % delta are 2-sigma estimates (95%) error estimates on predicted data Ypred
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);                  % extract lower and upper values for params, 1-sigma error estimates (68%)
        sigma = [abs(beta-sbeta(:,1)),abs(beta-sbeta(:,2))];                % calculates parameter lower and upper uncertainties

        %%% Extract D0 [Gy]
        D0 = beta(1);                                                       % unfaded D0 in [Gy]
        sD0 = mean(sigma(1,:));                                             % sD0 = CI on D0 [Gy]

        %%% Extract GOK growth order
        GOK_a = beta(2)+1;
        sGOK_a = mean(sigma(2,:));

        %%% Extract scalling parameters per aliquot
        A_factors = beta(3:end);

        %%% Calculate faded Ln/Llabmax
        Ln = Ln./A_factors';
        AvLn = mean(Ln);
        sLn = std(Ln); sLn(sLn==0) = 0.01*AvLn;

    end

    %%% Calculate (n/N)nat = (Ln/Llabmax)
    nNnat = AvLn;
    snNnat = nNnat*sqrt((sLn/AvLn)^2+((mean(labDdot)*2/100)/mean(labDdot))^2+(delta/Ypred)^2);

    %%% Simulate natural dose response curve %%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Athermal fading factor (s_tun fading)
    K = Hs*exp(-rhop^-(1/3)*rprime);

    %%% Single saturating exponential (SSE) fit %%%
    if SAR_fittype==1 % SSE fit

        %%% Compute (n/N)ss natfaded
        nss_nat = zeros(length(rprime),3);
        for k=1:length(rprime)                                              % loops through electron-hole distances
            m1=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-rhop.^(-1./3).*rprime(k)).*x).^2;
            m2=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-(10^rhop10_l).^(-1./3).*rprime(k)).*x).^2;
            m3=@(x)(natDdot(1)./D0.*(1-x)-Hs.*exp(-(10^rhop10_u).^(-1./3).*rprime(k)).*x).^2;
            nss_nat(k,1)=fminbnd(m1,0,1).*pr(k);
            nss_nat(k,2)=fminbnd(m2,0,1).*pr(k);                            % lower boundary
            nss_nat(k,3)=fminbnd(m3,0,1).*pr(k);                            % upper boundary
        end
        fadfact_nat = sum(nss_nat)./sum(pr);

        simDRC = [];
        for j=1:length(natdose)
            for k=1:length(rprime)                                          % loops through electron-hole distances
                simDRC(k,j) = pr(k)*((natDdot(1)/D0 )/(natDdot(1)/D0+K(k))...
                    *(1-exp(-natdose(j)*(1/D0+K(k)/natDdot(1)))));
            end
        end
        ComputedLxTx = sum(simDRC)/sum(pr);

        %%% Calculate D0_natfaded [Gy] from the simulated dose response curve
        beta0(2:end) = 1;
        Ddot = natDdot(1);
        [beta2,R,~,Cov,~] = nlinfit(natdose,ComputedLxTx,@SARfit_SSE,beta0);
        sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68);                % extract lower and upper values for params, 1-sigma error estimates (68%)
        sigma2 = [abs(beta2-sbeta2(:,1)),abs(beta2-sbeta2(:,2))];           % calculates parameter lower and upper uncertainties

        D0_natfaded = beta2(1); sD0_natfaded = mean(sigma2(1,:));
        a_natfaded = beta2(2);

    %%% General order kinetics (GOK) fit %%%
    elseif SAR_fittype==2 % GOK fit

        %%% Compute (n/N)ss natfaded
        nss_nat = zeros(length(rprime),3);
        for k=1:length(rprime) % loop over electron-hole distance
            m1=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-rhop.^(-1./3).*rprime(k)).*x).^2;
            m2=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-(10^rhop10_l).^(-1./3).*rprime(k)).*x).^2;
            m3=@(x)(natDdot(1)./D0.*(1-x).^GOK_a-Hs.*exp(-(10^rhop10_u).^(-1./3).*rprime(k)).*x).^2;
            nss_nat(k,1)=fminbnd(m1,0,1).*pr(k);
            nss_nat(k,2)=fminbnd(m2,0,1).*pr(k);                            % lower boundary
            nss_nat(k,3)=fminbnd(m3,0,1).*pr(k);                            % upper boundary
        end
        fadfact_nat =sum(nss_nat)./sum(pr);

        simDRC=[];
        for j=1:length(natdose)
            for k=1:length(rprime)
                simDRC(k,j) = pr(k)*(natDdot(1)/D0)/(natDdot(1)/D0+K(k))...
                    *(1-(1+(1/D0+K(k)/natDdot(1)) * natdose(j) * (GOK_a-1))^(-1/(GOK_a-1)));
            end
        end
        ComputedLxTx = sum(simDRC)/sum(pr);

        %%% Calculate D0_natfaded [Gy] from the simulated dose response curve
        beta0(3:end) = 1;
        [beta2,R,~,Cov,~] = nlinfit(natdose,ComputedLxTx,@SARfit_GOK,beta0);
        sbeta2 = nlparci(beta2,R,'covar',Cov,'alpha',1-.68);                % extract lower and upper values for params, 1-sigma error estimates (68%)
        sigma2 = [abs(beta2-sbeta2(:,1)),abs(beta2-sbeta2(:,2))];           % calculates parameter lower and upper uncertainties

        D0_natfaded = beta2(1); sD0_natfaded = mean(sigma2(1,:));
        a_natfaded = beta2(3);

    end

    %%% Calculate De [Gy] from the simulated dose response curve
    De = abs(D0_natfaded*log(1-([nNnat nNnat-snNnat nNnat+snNnat]/a_natfaded)));
    %De_err = sqrt((log(1-nNnat/a_natfaded)*sD0_natfaded)^2+(D0_natfaded/(a_natfaded*(1-nNnat/a_natfaded))*snNnat)^2);
    De_err = De(1)*sqrt(((mean(labDdot)*2/100)/mean(labDdot))^2+(snNnat/nNnat)^2+(sD0/D0)^2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Extract data from .mat file
    isoT = [records(i).rawdata(nSAR+1:nSAR+nITH).T];                        % creates a vector with the different temperatures T [°C] for isothermal decay
    measL = [records(i).rawdata(nSAR+1).t]; measL = measL(isfinite(measL)); % removes the NaN data from measL
    itht = [records(i).rawdata(nSAR+1:nSAR+nITH).t];                        % creates a matrix with all the time records of the isothermal decay [s]
    ithL = [records(i).rawdata(nSAR+1:nSAR+nITH).L];                        % creates a matrix with all the luminescence records of the isothermal decay
    ok = isfinite(itht); ithx = itht(ok); ithy = ithL(ok); itht=ithx ; ithL=ithy; % removes the NaN data from the data
    ithsort = NaN(length(measL),length(isoT));                              % creates matrix to facilitate data output storage

    isNAN = sum(isnan(ithL));
    if isNAN < length(ithL)

        %%% Band-tail state model, BTS %%%
        if ITH_fittype==1 % BTS model

            %%% nlinfit solution
            Avec = ones(1,nITH);                                            % creates a matrix of ones (one row, nITH columns)
            beta03 = [9 1.4 0.1 Avec]';                                     % initial estimates of s [s-1], Et [eV], Eu [eV], scalling parameter
            [beta3,R,~,Cov,~] = nlinfit(ithx,ithy,@FITH_BTSLiLi,beta03);

            %%% Confidence intervals
            measL = tvec;                                                   % redefines time dimension for the model, 1 line, 100 columns
            [IsoPred,delta2] = nlpredci(@FITH_BTSLiLi,tmat(:),beta3,R,'covar',Cov,'alpha',.05,'predopt','curve'); % delta2 are 2-sigma estimates (95%) error estimates on predicted data IsoPred
            sbeta3 = nlparci(beta3,R,'covar',Cov,'alpha',1-.68);            % extract lower and upper values for params, 1-sigma error estimates (68%)
            sigma3 = [abs(beta3-sbeta3(:,1)),abs(beta3-sbeta3(:,2))];       % calculates parameter lower and upper uncertainties

        %%% General order kinetics, GOK %%%
        elseif ITH_fittype==2 % GOK model

            %%% nlinfit solution
            Avec = ones(1,nITH);                                            % creates a matrix of ones (one row, nITH columns)
            beta03 = [9 1.4 4 Avec]';                                       % initial estimates of s [s-1], Et [eV], GOK order b, scalling parameter
            [beta3,R,~,Cov,~] = nlinfit(ithx,ithy,@FITH_GOK,beta03);

            %%% Confidence intervals
            measL = tvec;                                                   % redefines time dimension for model, 1 line, 100 columns
            [IsoPred,delta2] = nlpredci(@FITH_GOK,tmat(:),beta3,R,'covar',Cov,'alpha',.05,'predopt','curve'); % delta2 are 2-sigma estimates (95%) error estimates on predicted data IsoPred
            sbeta3 = nlparci(beta3,R,'covar',Cov,'alpha',1-.68);            % extract lower and upper values for params, 1-sigma error estimates (68%)
            sigma3 = [abs(beta3-sbeta3(:,1)),abs(beta3-sbeta3(:,2))];       % calculates parameter lower and upper uncertainties

        %%% Gauss model, GAUSS %%%
        elseif ITH_fittype==3 % GAUSS model

            %%% nlinfit solution
            Avec = ones(1,nITH);                                            % creates a matrix of ones (one row, nITH columns)
            beta03 = [9 1.4 0.1 Avec]';                                     % initial estimates of s [s-1], mu(Et) [eV], sigma(Et) [eV], scalling parameter
            [beta3,R,~,Cov,~] = nlinfit(ithx,ithy,@FITH_GAUSS,beta03);

            %%% Confidence intervals
            measL = tvec;                                                   % redefines time dimension for model, 1 line, 100 columns
            [IsoPred,delta2] = nlpredci(@FITH_GAUSS,tmat(:),beta3,R,'covar',Cov,'alpha',.05,'predopt','curve'); % delta2 are 2-sigma estimates (95%) error estimates on predicted data Isopred
            sbeta3 = nlparci(beta3,R,'covar',Cov,'alpha',1-.68);            % extract lower and upper values for params, 1-sigma error estimates (68%)
            sigma3 = [abs(beta3-sbeta3(:,1)),abs(beta3-sbeta3(:,2))];       % calculates parameter lower and upper uncertainties

        end

        %%% Extract parameters
        s10 = beta3(1); ss10 = mean(sigma3(1,:));
        Et = beta3(2); sEt = mean(sigma3(2,:));
        Eu = beta3(3); sEu = mean(sigma3(3,:));
        A = beta3(4:end);

    else
        IsoPred = NaN(size(tmat));
        delta2 = NaN(size(tmat));
        s10 = NaN; ss10 = NaN;
        Et = NaN; sEt = NaN;
        Eu = NaN; sEu = NaN;
        A = NaN;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    records(i).params.rhop10 = [rhop10,srhop10];
    records(i).params.g2d = [g2d,sg2d];                                     % in [%/decade]
    records(i).params.D0 = [D0,sD0];                                        % in [Gy]
    if SAR_fittype==2; records(i).params.GOK_a = [GOK_a,sGOK_a];            % for GOK model
    else records(i).params.GOK_a = [NaN,NaN]; end                           % for SSE and GAUSS models
    records(i).params.Et = [Et,sEt];                                        % in [eV]
    if ITH_fittype==1; records(i).params.Eu = [Eu,sEu];                     % for BTS model, in [eV]
    else records(i).params.Eu = [NaN,NaN]; end                              % for GOK and GAUSS models, in [eV]
    if ITH_fittype==2; records(i).params.GOK_b = [Eu,sEu];                  % for GOK model
    else records(i).params.GOK_b = [NaN, NaN]; end                          % for BTS and GAUSS models
    if ITH_fittype==3; records(i).params.sigmaEt = [Eu,sEu];                % for GAUSS model, in [eV]
    else records(i).params.sigmaEt = [NaN,NaN]; end	                        % for BTS and GOK models, in [eV]
    records(i).params.s10 = [s10,ss10];                                     % in [s-1]

    records(i).params.NatFadedDe = [De(1) De_err];                          % in [Gy]
    records(i).params.NatFadedD0 = D0_natfaded;                             % in [Gy]
    natDdotGyka = natDdot*ka;                                               % environmental dose rate [Gy/kyr]
    AgeOSL = records(i).params.NatFadedDe(1)./(natDdotGyka(1));
    sAgeOSL = AgeOSL.*sqrt((records(i).params.NatFadedDe(2)/records(i).params.NatFadedDe(1))^2+(natDdotGyka(2)/natDdotGyka(1))^2);
    records(i).params.AgeOSL = [AgeOSL sAgeOSL];                            % in [kyr]
    records(i).params.maxAgeOSL = (2*D0_natfaded)/natDdotGyka(1);           % in [kyr]

    records(i).SAR_model = SAR_MODEL;
    records(i).ITH_model = ITH_MODEL;

    %%% Save data for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fading (FAD)
    records(i).plot.xf      =   fvec;
    records(i).plot.yf      =   fYpred;
    records(i).plot.syf     =   fdelta;
    records(i).plot.tf      =   ft;
    records(i).plot.Lf      =   fL;

    %%% Dose response curve (SAR)
    records(i).plot.x       =   modvec;
    records(i).plot.y       =   Ypred./exp(-rhop*log(1.8*Hs.*(0.5.*modvec)).^3)./mean(A_factors);
    records(i).plot.sy      =   delta./exp(-rhop*log(1.8*Hs.*(0.5.*modvec)).^3)./mean(A_factors);
    ok = isfinite(t); L(ok) = L(ok)./exp(-rhop*log(1.8*Hs.*(0.5*t(ok))).^3);
    t = t(~isnan(L)); L = L(~isnan(L));
    if length(A_factors)>1
        records(i).plot.L   =   L./A_factors(index_aliquot(1:length(L)))';
    else
        records(i).plot.L   =   L./A_factors(index_aliquot(1:length(L)));
    end
    records(i).plot.t       =   t;

    %%% Isothermal decay (ITH)
    records(i).plot.xith    =   tmat;
    records(i).plot.yith    =   reshape(IsoPred,size(tmat));
    records(i).plot.syith   =   reshape(delta2,size(tmat));
    ithsort(1:length(itht)) =   itht;
    records(i).plot.tith    =   ithsort;
    ithsort(1:length(ithL)) =   ithL;
    records(i).plot.Lith    =   ithsort;
    records(i).plot.scale   =   A;

    %%% Saturation data (Kars plot)
    records(i).params.nNnat =   [nNnat,snNnat];
    records(i).params.nNss = [fadfact_nat(1) max(abs(fadfact_nat(1)-fadfact_nat(2)),abs(fadfact_nat(1)-fadfact_nat(3)))];

end

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.mat'],'records')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Excel file with the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

records2 = records;
samples = strings(nt,1);
nNnat_record = zeros(nt,2);
headers = {'Sample Temp (°C)','Unc. Temp','Dr (Gy/kyr)','Unc. Dr','rhop','Unc. rhop','g2d (%/decade)','Unc. g2d','D0 (Gy)','Unc. D0',...
    'GOK_a','Unc. a','Et (eV)','Unc. Et','Eu (eV)','Unc. Eu','GOK_b','Unc. b','sigmaEt (eV)','Unc. sigEt','s10','Unc. s10',...
    'De (Gy)','Unc. De','NatFadedD0 (Gy)','Age (ka)','Unc. age','maxAge (ka)',...
    'n/N_nat','Unc. n/N_nat','n/N_ss','Unc. n/N_ss'};
for i=1:nt
    records2(i).params.natDdot = records(i).params.natDdot.*ka;             % transform from [Gy/s] to [Gy/kyr]
    row(i,:) = table2array(struct2table(records2(i).params));
    samples(i) = string(records2(i).id);
end
tab = array2table(row,'VariableNames',headers);
tab2 = array2table(samples,'VariableNames',{'Samples'});
tab4 = [tab2,tab];
writetable(tab4,['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.xls'])


%%% Running time
tEnd = toc(tStart);
fprintf('Stage2a_Fitparameters took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));