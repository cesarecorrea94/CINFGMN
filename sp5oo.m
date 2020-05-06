warning off;
% dumpname = 'dumps/SP500/01-May-2020 22:07:24.mat';
% %%
% if ~exist('dumpseries', 'var')
%     if exist('dumpname', 'var')
%         dumpseries = INFGMN_series(dumpname, NaN);
%     else
%         offset = struct( ...
%             'delta', 1,     ...
%             'tau', 1,       ... quando criar novas componentes
%             'tmax', 1,      ... tempo para assimilar instâncias
%             'maxNC', 1);    %%% número máximo de componentes estáveis (não espúrias)
%         log2deltas  = -(0  +offset.delta   :offset.delta:  -offset.delta+ 10);
%         log2taus    = -(0  +offset.tau     :offset.tau:    -offset.tau+   12);
%         log2tmaxs   =  (3  +offset.tmax    :offset.tmax:   -offset.tmax+   9);
%         log2maxNCs  =  (3  +offset.maxNC   :offset.maxNC:  -offset.maxNC+  7);
%         % uniform = false; % proporções uniformes
%         normalize = true; %default
%         % regValue = 0; % regularização da covariância
%         % combinations = length(log2deltas)*length(log2taus)*length(log2tmaxs)*length(log2maxNCs)
%         
%         dumpseries = INFGMN_series('SP500', offset);
%         DS = SP500();
%         batchsizes = [10, 5, 2, 1];
%         warmup = (length(DS)-1)/2;
%         dumpseries.create(DS, warmup, batchsizes, normalize, ...
%             log2deltas, log2taus, log2tmaxs, log2maxNCs);
%     end
% end
% % fuzzyLogicDesigner(gmm.toFIS());
% 
% for ii=1:100
%     fprintf('Iteração: %i\n', ii);
%     dumpseries.pso_update();
%     dumpseries.update(@check_break);
% end

print_all_fis(self.dumps(end-1).stats, 'imgs/SP500/dont');
print_all_fis(self.dumps(end  ).stats, 'imgs/SP500/merge');
function print_all_fis(stats, folder)
    mkdir(folder);
    for jj = 1:length(stats)
        if  isnumeric(stats(jj).fis)
            continue;
        end
        myplotfis([folder '/t' num2str(stats(jj).t)], stats(jj).fis);
    end
end

function bool = check_break(dumps_ii, stats)
    bool =  2 * dumps_ii.tmax < stats(end).t            ... se deu tempo pra estabilidade
        &&  (   dumps_ii.maxNC < stats(end).NCs(end)    ... e:  meu NC ainda está alto
            &&  dumps_ii.maxNC < mean([stats.RMS_NC])   ...     e a média geral foi alta
            );
end

%     if exist('dumps', 'var')==0
%         load(dumpname, 'DS', 'dumps');
%         [~,index] = sortrows([dumps.mean_NC].'); dumps = dumps(index(end:-1:1)); clear index
%         for ii = length(dumps)-1:-1:1
%             if isnan(dumps(ii+1).sugeno_RMSE)
%                 dumps(ii+1) = [];
%             elseif isnan(dumps(ii).sugeno_RMSE)
%                 dumps(ii) = [];
%             elseif dumps(ii).sugeno_RMSE > dumps(ii+1).sugeno_RMSE
%                 dumps(ii) = [];
%             elseif dumps(ii).mean_NC == dumps(ii+1).mean_NC
%                 dumps(ii+1) = [];
%             end
%         end
%     end
%     dump_var = 'dump_merge';% dont merge
%     dumpfile = ['test/' dump_var '.mat'];
%     if exist(dumpfile, 'file')~=0
%         load(dumpfile, 'dump_cur');
%     else
%         dump_cur = more_stats(INFGMN_series_step(DS, 1, dumps(1), @check_break));
%         save(dumpfile, 'dump_cur');
%     end
%     eval([dump_var '= dump_cur;']);
%     myplotfis(['test/' dump_var ' t1400'], dump_cur.stats(end*1400/3000).fis);
%     myplotfis(['test/' dump_var ' t2400'], dump_cur.stats(end*2400/3000).fis);
%     myprint(['test/' dump_var], dump_cur);

