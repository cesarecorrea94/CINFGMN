warning off;
rng(0);
DS = sin_DS(true);
rangelen = [1 1];%[2*pi 2];

maxtime = 1*60*60;
nepochs = 1;
fistypes = {'mamdani', 'sugeno'};
spreads = 0.5;%2.^-((1:3)./2); % afeta o formato das MFs do FIS (não das componentes)
% uniform = false; % proporções uniformes
normalize = true; % true;%default
% regValue = 0; % regularização da covariância
betas   = 1; % ?
deltas = calcdelta(12/2, normalize)
taudists = calctau(14/2, normalize)
tmaxs   = length(DS)./[2   5   13] % tempo para assimilar instâncias
maxNCs  = [20    50] % número máximo de componentes estáveis (não espúrias)

if maxtime <= 60
    combinations = length(betas) * length(deltas) * length(taudists) ...
        * length(tmaxs) * length(maxNCs);
    fprintf('Coarse estimated time: %.2f h\n', combinations/60 *length(spreads)/3);
else
    dumps = INFGMN_test(DS, maxtime, nepochs, fistypes, spreads, ...
        normalize, deltas, taudists, tmaxs, maxNCs, rangelen);
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

function deltas = calcdelta(ndeltas, normalize) % define o spread inicial
    % deltas  = 10.^-((1 -2/3 :1/3: 2 -0/3)) % old
    deltas = zeros(1, ndeltas);
    if normalize
        mflim = (0.5 :0.5/(length(deltas)-1): 1)./2; % normalized
        for ii = 1:ndeltas
            spdmu = mf2mf([-mflim(ii) 0 mflim(ii)], 'trimf', 'gaussmf');
            deltas(ii) = sqrt(spdmu(1)); % normalized
        end
    else
        mflim = (1/3 :1/3/(ndeltas-1): 2/3).*(pi/2);
        for ii = 1:ndeltas
            spdmu = mf2mf([-mflim(ii) 0 mflim(ii)], 'trimf', 'gaussmf');
            deltas(ii) = sqrt(spdmu(1)) / (2*pi);
        end
    end
end

function taudists = calctau(ndist, normalize)
    % taus    = 10.^-((1 -0/4 :1/4: 2 +3/4)) % old % quando criar novas componentes
    if normalize
        xs = (0.25 :0.7/(ndist-1): 0.95)./2; % normalized
        taudists = [xs; zeros(1, ndist)]; % normalized
    else
        xs = (1/4 :1/4/(ndist-1): 2/4).*(pi/2);
        taudists = [pi/2; sin(pi/2)] - [xs; sin(xs)];
    end
end
