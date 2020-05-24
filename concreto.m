warning off;
% dumpname = 'dumps/abbaalone/21-May-2020 14:23:49.mat';
%%
if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname, NaN);
    else
        offset = struct( ...
            'log2delta', 1,     ...
            'log2tau', 1,       ... quando criar novas componentes
            'log2tmax', 0.5,    ... tempo para assimilar instâncias
            'log2maxNC', 0.5);  %%% número máximo de componentes estáveis (não espúrias)
        log2deltas  =  -(offset.log2delta+  0   :offset.log2delta: 10   -offset.log2delta);
        log2taus    =  -(offset.log2tau  +  0   :offset.log2tau:   10   -offset.log2tau);
        log2tmaxs   =   (offset.log2tmax +  4   :offset.log2tmax:   7   -offset.log2tmax);
        log2maxNCs  =   (offset.log2maxNC+  4   :offset.log2maxNC:  7   -offset.log2maxNC);
        normalize = true;%default

        DS = concrete_DS();
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_fis = false;
        
        dumpseries = INFGMN_series('concreto', offset);
        dumpseries.create_nonseries(DS, ...
            fis_types, save_fis, normalize, ...
            log2deltas, log2taus, log2tmaxs, log2maxNCs);
    end
end
% fuzzyLogicDesigner(gmm.toFIS());

try
    for ii=1:1
        fprintf('Iteração: %i\n', ii);
        dumpseries.update_nonseries();
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
