function [out] = FGOK_GK(beta, t)

global isoT measL rhop10
mat=nan(length(measL),length(isoT));
mat(1:length(t))=t;

s=10.^beta(1);
E=beta(2);
b=beta(3);
A=beta(4:end); 
T=isoT; 
rhop = 10.^rhop10;

out=[];

kB=8.617343e-5; %eV/K
K=@(x)s.*exp(-E./kB./(x+273.15));

for j=1:length(mat(1,:));
    ok=isfinite(mat(:,j)); time=mat(:,j);
    kars=exp(-rhop*log(1.8*3e15.*(250+time(ok))).^3);
    out=[out; A(j).*kars.*(1-(1-b).*K(T(j)).*time(ok)).^(1./(1-b))];
end
out;
out = out';


