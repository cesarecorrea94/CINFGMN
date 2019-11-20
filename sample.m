warning off;
rng(0);
DS = sin_DS();
rangelen = [2*pi 2];

maxtime = 1.5*60*60;
nepochs = 50;
fistypes = {'mamdani', 'sugeno'};
spreads = 2.^-((1:3)./2); % afeta o formato das MFs do FIS (não das componentes)
betas   = 1; % ?
deltas = calcdelta(3)
taudists = calctau(4)
tmaxs   = length(DS).*[2   5   13] % tempo para assimilar instâncias
maxNCs  = [20 35   ] % número máximo de componentes estáveis (não espúrias)
% uniform = false; % proporções uniformes
normalize = false; % true;%default
% regValue = 0; % regularização da covariância

if maxtime <= 60
    combinations = length(betas) * length(deltas) * length(taudists) ...
        * length(tmaxs) * length(maxNCs);
    fprintf('Coarse estimated time: %.2f h\n', combinations/60 *length(spreads)/3);
else
    dumps = INFGMN_test(DS, maxtime, nepochs, fistypes, spreads, ...
        normalize, deltas, taudists, tmaxs, maxNCs);
    if maxtime >= 30*60
        save(['dumps/sin_' ...
            num2str(length(DS)) 'samples_' ...
            num2str(nepochs) 'nepochs_' ...
            num2str(length(spreads)) 'spreads_' ...
            [fistypes{:}] datestr(now) '.mat'], ...
            'dumps');
    end
    % fuzzyLogicDesigner(gmm.toFIS());
end
warning on;

function deltas = calcdelta(ndeltas) % define o spread inicial
    % deltas  = 10.^-((1 -2/3 :1/3: 2 -0/3)) % old
    deltas = zeros(1, ndeltas);
%     mflim = (0.5 :0.5/(length(deltas)-1): 1)./2; % normalized
    mflim = (1/3 :1/3/(ndeltas-1): 2/3).*(pi/2);
    for ii = 1:ndeltas
        spdmu = mf2mf([-mflim(ii) 0 mflim(ii)], 'trimf', 'gaussmf');
        deltas(ii) = sqrt(spdmu(1)) / (2*pi);
    end
end

function taudists = calctau(ndist)
    % taus    = 10.^-((1 -0/4 :1/4: 2 +3/4)) % old % quando criar novas componentes
    % taudists = (0.25 :0.7/(ndist-1): 0.95)./2 % normalized
    xs = (1/4 :1/4/(ndist-1): 2/4).*(pi/2);
    taudists = [pi/2; sin(pi/2)] - [xs; sin(xs)];
end
