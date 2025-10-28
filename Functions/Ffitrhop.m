%% Fits fading data to derive rho' after Huntley (2006) %%
%%% Solves for all aliquots at once
%%% Here no trapping and no thermal detrapping is assumed 
%%% during room temperature isothermal holding

% Georgina King [georgina.king@unil.ch], August 2015
% Maxime Bernard [maxime.bernard@unil.ch], 22/01/2025
% Chloé Bouscary [chloebouscary@gmail.com], 06/2025


function out = Ffitrhop(beta0,t)

global Hs index_aliquotFAD

rhop=10.^beta0(1);
c=beta0(2:end)'; 

if length(t) <100                                                           % assumes t is time read from L measurements else otherwise t is Fvec
    idx = index_aliquotFAD;
    if length(c)>1
        C_fac = c(idx);
    else
        idx(1:end) = c;
        C_fac = idx;
    end
else 
    C_fac = mean(c);
end

fade = C_fac.*exp(-rhop.*(log(1.8.*Hs.*t)).^3);

out = fade;