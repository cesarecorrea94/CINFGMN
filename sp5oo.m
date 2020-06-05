warning off;
% dumpname = 'dumps/SP500/01-May-2020 22:07:24.mat';
%%
if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname, NaN);
    else
        offset = struct( ...
            'log2delta',    1,  ...
            'log2tau',      1,  ... quando criar novas componentes
            'log2tmax',     1,  ... tempo para assimilar instâncias
            'log2maxNC',    1); %%% número máximo de componentes estáveis (não espúrias)
        log2deltas  =  -(offset.log2delta+  0   :offset.log2delta: 10   -offset.log2delta);
        log2taus    =  -(offset.log2tau  +  0   :offset.log2tau:   12   -offset.log2tau);
        log2tmaxs   =   (offset.log2tmax +  3   :offset.log2tmax:   9   -offset.log2tmax);
        log2maxNCs  =   (offset.log2maxNC+  3   :offset.log2maxNC:  7   -offset.log2maxNC);
%         log2deltas  =  -( 8.5   -offset.log2delta   :offset.log2delta:  offset.log2delta    +8.5);
%         log2taus    =  -(10.5   -offset.log2tau     :offset.log2tau:    offset.log2tau      +10.5);
%         log2tmaxs   =   ( 3.5   -offset.log2tmax    :offset.log2tmax:   offset.log2tmax     +3.5);
%         log2maxNCs  =   ( 4.5   -offset.log2maxNC   :offset.log2maxNC:  offset.log2maxNC    +4.5);
        normalize = false;%true; %default
        
        DS = SP500();
        warmup = (length(DS)-1)/2;
        batchsizes = [10, 5, 2, 1];
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_stats = false;%true;
        save_fis = (length(DS)-1)-2800:200:(length(DS)-1);
        maxFoCSize = 7;
        
        dumpseries = INFGMN_series('SP500', offset);
        dumpseries.create(DS, warmup, batchsizes, ...
            fis_types, save_stats, save_fis, maxFoCSize, normalize, ...
            log2deltas, log2taus, log2tmaxs, log2maxNCs);
    end
end
% fuzzyLogicDesigner(gmm.toFIS());

for ii=1:100
    fprintf('Iteração: %i\n', ii);
    dumpseries.update(@step_core_function);
    dumpseries.pso_update();
end

%     myplotfis(['test/' dump_var ' t1400'], dump_cur.stats(end*1400/3000).fis);
%     myplotfis(['test/' dump_var ' t2400'], dump_cur.stats(end*2400/3000).fis);
%     myprint(['test/' dump_var], dump_cur);

% self = dumpseries.myself();
% print_all_fis(self.dumps(end).stats, 'imgs/SP500/numerge');
% function print_all_fis(stats, folder)
%     mkdir(folder);
%     for jj = 1:length(stats)
%         if  isnumeric(stats(jj).fis)
%             continue;
%         end
%         myplotfis([folder '/t' num2str(stats(jj).t)], stats(jj).fis);
%     end
% end

function dumps_ii = step_core_function(self, dumps_ii)
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
                        'Endpoints', 'discard') )
            delete(gmm);
            error('INFGMN:Abortar', 'maxNC < movmean.');
        end
        pre_NCs = pos_NCs;
    end
    delete(gmm);
    if sqrt(mean(NCs.^2)) < 7
            error('INFGMN:Abortar', 'RMS_NC < 7.');
    end
end
