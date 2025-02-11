% We assume athermal fading during isothermal decay
%%maxime bernard%%
% last update: 24/01/2025

function [out] = FLiLiGK(beta, t)

global isoT measL rhop10 Hs
mat=nan(length(measL),length(isoT));
mat(1:length(t))=t;

s=10.^beta(1);
Et=beta(2);
Eu=beta(3);
A=beta(4:end); 

T=isoT; 
rhop = 10.^rhop10;

out=[];

for j=1:length(mat(1,:));
    ok=isfinite(mat(:,j)); time=mat(:,j);
    Fading=exp(-rhop*log(1.8*Hs.*(250+time(ok))).^3);
    out=[out; A(j).*Fading.*ThermaldecayBandTailGK2(time(ok),T(j),Et,Eu,s)];
end
out = out';
