%%Benny Guralnik%%
%%maxime Bernard%%
% last update: 22/01/2025

function out = Sarfit_GAUSS(beta0, x);
global rhop10 Hs natDdot
rhop=10.^rhop10;
D0=beta0(1); sigmaD0=beta0(2); a=beta0(3);

xx=linspace(D0-3.*sigmaD0,D0+3.*sigmaD0,30)';
p=normpdf(xx,D0,sigmaD0); p=p./sum(p);

% x input is natural doses (Gy)
t = x./natDdot(1).*(3600*24*365.25*1000); %time in s

Fading=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3); % include fading during lab irradiation [Huntley, 2006]

% Matrix operation for Fgrowth
[p_grid, t_grid] = ndgrid(p, x);        % Create grids for probabilities and time

% Compute the growth matrix
Fgrowth_matrix = p_grid .* a .* Fading.* (1 - exp(-t_grid ./ (10.^xx)));
Fgrowth = sum(Fgrowth_matrix, 1);

Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Fgrowth];
end


% Fgrowth=zeros(length(xx),length(x));
% for j=1:length(xx)
% Fgrowth(j,:) = p(j).*a.*(1-exp(-x./10.^xx(j)));
% end
% Fgrowth = sum(Fgrowth);
% Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out=[Fgrowth];
