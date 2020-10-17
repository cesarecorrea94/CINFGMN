warning off;
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.viaveis.batch10.batch5.batch2.batch1.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.viaveis.batch10.batch5.batch2.batch1.mat_dont.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.viaveis.batch10.batch5.batch2.batch1.vNEWdiv2.mat';
%%
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.subviaveis.batch1_merge.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.subviaveis.batch1_merge_vdiv4.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.subviaveis.batch1_merge_noweight.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.subviaveis.batch1_merge_maxMF11.mat';
% dumpname = 'dumps/DynNLSys/23-Jun-2020 00:25:26.subviaveis.batch1_merge_maxMF10.mat';
%%
dumpname = 'dumps/DynNLSys/final_merge.mat';
dumpname = 'dumps/DynNLSys/final_dont.mat';
% fodase(dumpname);
% function fodase(dumpname)

if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname);
    else
        offset = struct( ...
            'log2delta',    0.25,    ...
            'log2tau',      0.25,   ... quando criar novas componentes
            'log2tmax',     0,    ... tempo para assimilar instâncias
            'log2maxNC',    0.25);  %%% número máximo de componentes estáveis (não espúrias)
        paramstruct = struct( ...
            'log2delta',   -(offset.log2delta : offset.log2delta : 10), ... log2(0.03162)
            'log2tau',     -(offset.log2tau   : offset.log2tau   : 20), ... log2(0.003981)
            'log2tmax',     log2(100),    ...
            'log2maxNC',    (5 :offset.log2maxNC: 10)                    ... log2(55)
        );
        normalize = true; % default
        comb = ...
            length(paramstruct.log2delta) * ...
            length(paramstruct.log2tau) * ...
            length(paramstruct.log2tmax) * ...
            length(paramstruct.log2maxNC);
        fprintf('%i comb in %s\n', comb, duration(0,0, 1.25*comb));

        DS = dynamic_nonlinear_system();
        warmup = 0;
        batchsize = 20;%[500, 5, 2, 1];
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_stats = true;
        save_fis = 200:200:3000;
        maxFoCSize = 1;%7;%11;
        Smerge = 1;%0.7;
        
        dumpseries = INFGMN_series('DynNLSys');
        dumpseries.create(DS, warmup, batchsize, ...
            save_stats, fis_types, save_fis, Smerge, maxFoCSize, ...
            normalize, paramstruct, offset);
    end
end

try
%     self = dumpseries.myself();
%     if  all(isnan([self.dumps.sugeno_RMSE]))
%         dumpseries.update(@check_abort);
%         dumpseries = limpa_viaveis(dumpseries);
%     end
%     clear self;
    for ii=1:100
        fprintf('Iteração: %i\n', ii);
        dumpseries.deepen();
        dumpseries.update(@INFGMN_series.step);
    end
catch ME
    if ~strcmp(ME.message,'StopIteration')
        rethrow(ME);
    end
end
error('Stopped');

dumpname = dumpseries.dumpname;
self = load(dumpname);
% [~,index] = sortrows([self.dumps.sugeno_RMSE].'); self.dumps = self.dumps(index); clear index
[~,name,~] = fileparts(dumpname);
%print_all_fis(self, 1, ['imgs/DynNLSys/' name]);
myprint(self, 1, ['imgs/DynNLSys/' name]);

% end
%%

function dumpseries = limpa_viaveis(dumpseries)
    self = dumpseries.myself();
    viaveis = ~isnan([self.dumps.mean_NC]);
    self.dumps = self.dumps(viaveis);
    fprintf('%i abortados\n', sum(~viaveis));
    [self.dumps.cputime] = deal(NaN);
    [self.dumps.stats] = deal([]);
    [self.dumps.mean_NC] = deal(NaN);
    [self.dumps.RMS_NC] = deal(NaN);        %%
    [self.dumps.mamdani_RMSE] = deal(NaN);  %%
    [self.dumps.sugeno_RMSE] = deal(NaN);   %%
    self.doMerge = false; %% true;
    self.Smerge = 0.9;
    self.Sdeath = 0.75;
    self.maxFoCSize = 9;
    self.batchsize = 20;
    [filepath,name,ext] = fileparts(dumpseries.dumpname);
    dumpseries.dumpname = fullfile(filepath, [name '_dont.viaveis' ext]);
    dumpseries.save_myself(self, false);
end

function dumps_ii = check_abort(self, dumps_ii)
    [gmm, ~] = INFGMN_series.create_INFGMN(self, dumps_ii);
    lenTrainDS = length(self.DS)-1;
    batchsize = 500;
    warm_up = mod(lenTrainDS, batchsize);
    n_batches = (lenTrainDS - warm_up) / batchsize;
    assert(mod(n_batches,1)==0, 'n_batches is not integer');
    NCs = zeros(1, lenTrainDS);
    if warm_up
        [gmm, NCs(1:warm_up)] = gmm.train(self.DS(1:warm_up, :));
    end
    for ii = 1:n_batches
        tt = warm_up + ii * batchsize;
        batch = (tt-batchsize+1):tt;
        [gmm, NCs(batch)] = gmm.train( self.DS(batch, :) );
        if  sub_check_abort(tt, NCs, batchsize)
            delete(gmm);
            error('INFGMN:Abortar', 'abort.');
        end
    end
    delete(gmm);
    dumps_ii.mean_NC = mean(NCs, 'omitnan');
    if self.save_stats
        dumps_ii.stats = NCs;
    end
end

function abort = sub_check_abort(tt, NCs, batchsize)
    abort = false;
    ruido = 200;
    for etapa = 1:3
        if ~abort && tt-batchsize < etapa*1000 && tt >= etapa*1000
            filtro_NC = NCs( (ruido+1:1000) +(etapa-1)*1000 );
            mean_NC = mean( filtro_NC );
            std_NC  = std ( filtro_NC(2:end-1), 1 );
            abort = mean_NC < 7 ...
                ||  mean_NC > 7*4 ...
                ||  std_NC > 0;
        end
    end
end
