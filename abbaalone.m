warning off;
%dumpname = 'dumps/abbaalone/29-Jun-2020 both.mat';
%%
%dumpname = 'dumps/abbaalone/07-Dec-2020.merge-test.mat';
%dumpseries = INFGMN_series(dumpname);
%self = dumpseries.myself;
%getINFGMNParam(self);
%%

if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname);
    else
        DS = abalone_DS();
        fis_types = {'sugeno'};%{'mamdani', 'sugeno'};
        save_fis = true;
        doMerge = true;
        maxFoCSize = 5;
        normalize = true;%default

        halfTrain = 0.6 * length(DS) / 2;
        offset = struct( ...
            'log2delta', 0.125,     ...
            'log2tau', 0.25,       ... quando criar novas componentes
            'log2tmax', 0, ...1,    ... tempo para assimilar instâncias
            'log2maxNC', 0); ...1);  %%% número máximo de componentes estáveis (não espúrias)
        paramstruct = struct( ...
            'log2delta',   -1.875,...-(offset.log2delta   :offset.log2delta:  5), ...
            'log2tau',     -21.5,...-(offset.log2tau     :offset.log2tau:    25), ...
            'log2tmax',    10.2912858339966,... [log2(halfTrain),  99], ...(offset.log2tmax +  3   :offset.log2tmax:   9   -offset.log2tmax), ...
            'log2maxNC',   7.79128583399659...  log2(halfTrain)-2.5    ...(offset.log2maxNC+  3   :offset.log2maxNC:  9   -offset.log2maxNC) ...
        );
        comb = ...
            length(paramstruct.log2delta) * ...
            length(paramstruct.log2tau) * ...
            length(paramstruct.log2tmax) * ...
            length(paramstruct.log2maxNC);
        fprintf('%i comb in %s\n', comb, duration(0,0, 2*comb));
        
        dumpseries = INFGMN_series('abbaalone');
        dumpseries.create_nonseries(DS, ...
            fis_types, save_fis, doMerge, maxFoCSize, ...
            normalize, paramstruct, offset);
    end
end

try
    dumpseries.update(@INFGMN_series.step_nonseries);
catch ME
    if ~strcmp(ME.message,'StopIteration')
        rethrow(ME);
    end
end
error('Stopped');

%     myprint(['test/' dump_var], dump_cur);

xola = {'dont','merge'};
for aa=1:2
    print_all_fis(load([dumpname '_' xola{aa} '.mat']), ...
        ['imgs/abbaalone/' xola{aa}]);
end
function print_all_fis(self, folder)
    ranges = minmaxDS(self.DS);
    for ii = 1:1%0
        subfolder = [folder '/' num2str(ii)];
        mkdir(subfolder);
        stats = self.dumps(ii).stats;
        for jj = 1:length(stats)
            myplotfis(stats(jj).fis, 1, ranges, ...
                [subfolder '/' num2str(jj) '_NC' num2str(stats(jj).NC)]);
        end
    end
end
