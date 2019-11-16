warning off;
rng(0);
DS = sin_DS();

maxtime = 6;
nepochs = 50;
fistypes = {'mamdani', 'sugeno'};
spreads = 2.^-((1:3)./2) % afeta o formato das MFs do FIS (não das componentes)
betas   = 1; % ?
deltas  = 10.^-((1 -2/4 :1/4: 2 -2/4)) % define o spread inicial
taus    = 10.^-((1 -1/2 :1/2: 2 +1/2)) % quando criar novas componentes
tmaxs   = length(DS).*[1 3 5] % tempo para assimilar instâncias
maxNCs  = [20 35 50] % número máximo de componentes estáveis (não espúrias)
% uniform = false; % proporções uniformes
% normalize = true;
% regValue = 0; % regularização da covariância

if maxtime <= 60
    combinations = length(betas) * length(deltas) * length(taus) ...
        * length(tmaxs) * length(maxNCs);
    fprintf('Coarse estimated time: %.2f h\n', combinations/60);
else
    dumps = testINFGMN(DS, maxtime, nepochs, fistypes, spreads, ...
        betas, deltas, taus, tmaxs, maxNCs);
    if maxtime >= 30*60
        save(['sin_' num2str(length(DS)) 'samples_' ...
            num2str(nepochs) 'nepochs_' ...
            num2str(length(spreads)) 'spreads_' ...
            [fistypes{:}] datestr(now) '.mat'], ...
            'dumps');
    end
    % fuzzyLogicDesigner(gmm.toFIS());
end
warning on;
