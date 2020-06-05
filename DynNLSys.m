warning off;
dumpname = 'dumps/DynNLSys/27-May-2020 10:54:08.mat_SUB+0.5  +0.5 +0.25+0.125.mat_SUB+0.25  +0.25 +0.125+0.0625.mat_SUB+0.125  +0.125 +0.0625+0.03125.mat';
%%
if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname, NaN);
    else
        offset = struct( ...
            'log2delta',    1,      ...
            'log2tau',      1,      ... quando criar novas componentes
            'log2tmax',     0.5,    ... tempo para assimilar instâncias
            'log2maxNC',    0.25);  %%% número máximo de componentes estáveis (não espúrias)
        log2deltas  = log2(0.03162);    %-(offset.log2delta+  0   :offset.log2delta: 10   -offset.log2delta);
        log2taus    = log2(0.003981);   %-(offset.log2tau  +  0   :offset.log2tau:   10   -offset.log2tau);
        log2tmaxs   = log2(100);        % (offset.log2tmax +  4   :offset.log2tmax:   8.5 -offset.log2tmax);
        log2maxNCs  = log2(55);         % (offset.log2maxNC+  4   :offset.log2maxNC:  6   -offset.log2maxNC);
        normalize = true;%default

        DS = dynamic_nonlinear_system();
        warmup = 0;
        batchsizes = [500, 5, 2, 1];%[20, 10, 5, 2];
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_stats = true;
        save_fis = 200:200:3000;
        % maxFoCSize = 11??;
        
        dumpseries = INFGMN_series('DynNLSys', offset);
        dumpseries.create(DS, warmup, batchsizes, ...
            fis_types, save_stats, save_fis, normalize, ...
            log2deltas, log2taus, log2tmaxs, log2maxNCs);
    end
end
% fuzzyLogicDesigner(gmm.toFIS());

% try
%     for ii=1:100
%         fprintf('Iteração: %i\n', ii);
%         dumpseries.update(@core_function);
%         dumpseries.pso_update();
%     end
% catch ME
%     if ~strcmp(ME.message,'StopIteration')
%         rethrow(ME);
%     end
% end

%     myplotfis(['test/' dump_var ' t1400'], dump_cur.stats(end*1400/3000).fis);
%     myplotfis(['test/' dump_var ' t2400'], dump_cur.stats(end*2400/3000).fis);
%     myprint(['test/' dump_var], dump_cur);

self = dumpseries.myself();
[~,idxmin] = min([self.dumps.sugeno_RMSE].');
% print_all_fis(self.dumps(idxmin).stats, 'imgs/DynNLSys/dont');
print_all_fis(self.dumps(idxmin).stats, 'imgs/DynNLSys/merge9');
function print_all_fis(stats, folder)
    mkdir(folder);
    for jj = 1:length(stats)
        if  isnumeric(stats(jj).fis)
            continue;
        end
        myplotfis([folder '/t' num2str(stats(jj).t)], stats(jj).fis);
    end
end

function dumps_ii = core_function(self, dumps_ii)
    check_abort(self, dumps_ii);
    dumps_ii = INFGMN_series.step(self, dumps_ii);
end

function check_abort(self, dumps_ii)
    gmm = INFGMN(minmaxDS(self.DS), 'normalize', self.normalize, ...
        'delta', dumps_ii.delta, 'tau',  dumps_ii.tau, ...
        'tmax',  dumps_ii.tmax, 'spmin', dumps_ii.spmin );
    lenTrainDS = length(self.DS)-1;
    window = 2 * ceil(dumps_ii.tmax);
    warm_up = mod(lenTrainDS, window);
    n_batches = (lenTrainDS - warm_up) / window;
    assert(mod(n_batches,1)==0, 'n_batches is not integer');
    NCs = zeros(1, lenTrainDS);
    if warm_up
        [gmm, NCs(1:warm_up)] = gmm.train(self.DS(1:warm_up, :));
        pre_NCs = NCs(1:warm_up);
    end
    for ii = 1:n_batches
        tt = warm_up + ii * window;
        batch = (tt-(window-1)):tt;
        [gmm, NCs(batch)] = gmm.train( self.DS(batch, :) );
        pos_NCs = NCs(batch);
        if  any( dumps_ii.maxNC < movmean( [pre_NCs pos_NCs], window, ...
                        'Endpoints', 'discard') ) ...
                || sub_check_abort(dumps_ii, NCs(1:tt))
            delete(gmm);
            error('INFGMN:Abortar', 'maxNC < movmean || sub_check_abort.');
        end
        pre_NCs = pos_NCs;
    end
    delete(gmm);
end

function abort = sub_check_abort(dumps_ii, NCs)
    tt = length(NCs);
    window = 2 * ceil(dumps_ii.tmax);
    abort = window < mod(tt, 1000)      ... se passou do ruído
        &&  dumps_ii.maxNC < NCs(end);  ... e meu NC ainda está alto
    if ~abort && tt-window < 900 && tt >= 900
        mean_NC = mean( NCs(100:900) );
        abort = mean_NC < 18/2 ...
            ||  mean_NC > 18*2 ...
            ||  mean_NC > dumps_ii.maxNC;
    end
    if ~abort && tt-window < 1900 && tt >= 1900
        mean_NC = mean( NCs(1100:1900) );
        abort = mean_NC < 12/2 ...
            ||  mean_NC > 12*2 ...
            ||  mean_NC > dumps_ii.maxNC;
    end
    if ~abort && tt-window < 2900 && tt >= 2900
        mean_NC = mean( NCs(2100:2900) );
        abort = mean_NC < 18/2 ...
            ||  mean_NC > 18*2 ...
            ||  mean_NC > dumps_ii.maxNC;
    end
end

function dumps_ii = more_stats(dumps_ii)
    NCs = [dumps_ii.stats.NCs];
    dumps_ii.mean_NC_1  = mean( NCs(   1:1000) );
    dumps_ii.mean_NC_2  = mean( NCs(1001:2000) );
    dumps_ii.mean_NC_3  = mean( NCs(2001:3000) );
end
