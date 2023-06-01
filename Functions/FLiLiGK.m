function [out] = FLiLiGK(beta, t)

global isoT measL rhop10
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
    kars=exp(-rhop*log(1.8*3e15.*(250+time(ok))).^3);
    out=[out; A(j).*kars.*ThermaldecayBandTailGK2(time(ok),T(j),Et,Eu,s)];
end
out = out';
