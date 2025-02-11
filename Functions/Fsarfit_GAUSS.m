%%Benny Guralnik%%
%%maxime Bernard%%
% last update: 22/01/2025

function out = Fsarfit_GAUSS(beta0, x);
global rhop10
global labDdot index_aliquot Hs
rhop=10.^rhop10;
D0=beta0(1); sigmaD0=beta0(2); a=beta0(3:end); 

xx=linspace(D0-3.*sigmaD0,D0+3.*sigmaD0,30)';
p=normpdf(xx,D0,sigmaD0); p=p./sum(p);

t = x;
if length(x) <100 % assumes x is time read from L measurements else time is modvec
    idx = index_aliquot;
    Ddot = labDdot(idx);
else 
    idx = 1;
    a = mean(a(:));
    Ddot = mean(labDdot);
%     Ddot = labDdot(idx);
end

% fading
kars=exp(-rhop*log(1.8*Hs.*(0.5*t)).^3);

% Matrix operation for Fgrowth
[p_grid, t_grid] = ndgrid(p, t);        % Create grids for probabilities and time
[xx_grid, idx_grid] = ndgrid(xx, idx); % Grids for xx and index

labDdot_grid = Ddot(idx_grid);      % Align labDdot with the index
a_grid = a(idx_grid);                  % Align 'a' with the index

% Compute the growth matrix
Fgrowth_matrix = p_grid .* a_grid .* (1 - exp(-labDdot_grid ./ (10.^xx_grid) .* t_grid));
Fgrowth = kars .* sum(Fgrowth_matrix, 1);

Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Fgrowth];
end


% Fgrowth=zeros(length(x),length(xx));
% for j=1:length(xx)
% % Fgrowth = p(j).*a.*kars.*(1-exp(-x./xx(j))); 
% Fgrowth(:,j) = p(j).*a(idx).*(1-exp(-labDdot(idx)./(10.^xx(j)).*t));
% % Fgrowth(:,j) = p(j).*a.*(1-exp(-labDdot./(xx(j)).*x));
% end
% 
% Fgrowth = kars.*sum(Fgrowth,2);
% Fgrowth(isnan(Fgrowth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out=[Fgrowth];
% end
% 
% clf
% figure(1)
% semilogx(x,Fgrowth); hold on
% semilogx(x,y,'p')