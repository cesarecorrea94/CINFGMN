warning off;
% dumpname = 'dumps/abbaalone/21-May-2020 18:00:50.mat';
%%
if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname, NaN);
    else
        offset = struct( ...
            'log2delta', 1,     ...
            'log2tau', 2,       ... quando criar novas componentes
            'log2tmax', 1,    ... tempo para assimilar instâncias
            'log2maxNC', 1);  %%% número máximo de componentes estáveis (não espúrias)
        log2deltas  =  -(offset.log2delta+  0   :offset.log2delta:  5   -offset.log2delta);
        log2taus    =  -(offset.log2tau  +  0   :offset.log2tau:   12   -offset.log2tau);
        log2tmaxs   =   (offset.log2tmax +  3   :offset.log2tmax:   9   -offset.log2tmax);
        log2maxNCs  =   (offset.log2maxNC+  3   :offset.log2maxNC:  9   -offset.log2maxNC);
        normalize = true;%default

        DS = abalone_DS();
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_fis = false;
        doMerge = true;
        maxFoCSize = 3;
        
        dumpseries = INFGMN_series('abbaalone', offset);
        dumpseries.create_nonseries(DS, ...
            fis_types, save_fis, doMerge, maxFoCSize, normalize, ...
            log2deltas, log2taus, log2tmaxs, log2maxNCs);
    end
end
% fuzzyLogicDesigner(gmm.toFIS());

try
    for ii=1:1
        fprintf('Iteração: %i\n', ii);
        dumpseries.update(INFGMN_series.step_nonseries);
%         dumpseries.pso_update();
    end
catch ME
    if ~strcmp(ME.message,'StopIteration')
        rethrow(ME);
    end
end

%     myplotfis(['test/' dump_var ' t1400'], dump_cur.stats(end*1400/3000).fis);
%     myplotfis(['test/' dump_var ' t2400'], dump_cur.stats(end*2400/3000).fis);
%     myprint(['test/' dump_var], dump_cur);

% self = dumpseries.myself();
% print_all_fis(self.dumps(1).stats, 'imgs/DynNLSys/dont');
% print_all_fis(self.dumps(2).stats, 'imgs/DynNLSys/merge');
% function print_all_fis(stats, folder)
%     mkdir(folder);
%     for jj = 1:length(stats)
%         if  isnumeric(stats(jj).fis)
%             continue;
%         end
%         myplotfis([folder '/t' num2str(stats(jj).t)], stats(jj).fis);
%     end
% end
