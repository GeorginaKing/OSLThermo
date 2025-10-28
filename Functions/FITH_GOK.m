%% Fits ITH data to a GOK model %%
%%% Athermal detrapping during isothermal decay is assumed

% Georgina King [georgina.king@unil.ch], 2015


function out = FITH_GOK(beta,t)

global isoT measL rhop10 Hs 

%%% Create matrix of times, one column by isothermal temperature
mat=nan(length(measL),length(isoT));
mat(1:length(t))=t;

%%% Extract parameters
s=10.^beta(1);
E=beta(2);
b=beta(3);
A=beta(4:end); 
T=isoT; 
rhop = 10.^rhop10;

%%% Initialise output
out=[];

%%% Isothermal decay data
for j=1:length(mat(1,:))
    ok=isfinite(mat(:,j)); time=mat(:,j);
    
    %%% Anomalous fading
    kars = exp(-rhop*log(1.8*Hs.*(250+time(ok))).^3);

    %%% Compute the isothermal decay curves
    out = [out; A(j).*kars.*ThermaldecayGOK(time(ok),T(j),E,b,s)];
end

out = out';