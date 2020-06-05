% maxNCs = [20 35 50];
% ndeltas = 5;
% tab = table('Size', [length(maxNCs)*ndeltas 3], ...
%     'VariableTypes', ["double", "double", "double"], ...
%     'VariableNames', {'maxNC',  'delta',  'tau'});
% dim = 2;
% index = 0;
% for deltai = 0:ndeltas
% for maxNC = maxNCs
%     index = index+1;
%     deltamin = 10^(-2-deltai/ndeltas);
%     % dado uma cov mínima (nenhuma componente terá cov menor)
% 	% com NC igualmente espaçadas
%     % no melhor dos casos:
% 	tau = taufordist(...
% 		ones(1,dim)./(2*maxNC^(1/dim)), ... % no ponto menos denso
% 		deltamin, 1)
%     tab(index,:) = {maxNC, deltamin, tau};
% end
% end
% tab


% gmm = INFGMN(minmaxDS(DS), 'delta', 0.0919111184034876, ...
%     'tau', 0.0361844409534152, ...
%     'tmax', 2340, 'spmin', 46.8, 'normalize', false );
% for i=1:10
%     gmm.train(DS);
%     gmm.setFisType('mamdani');
%     y = gmm.recall(DS(:,1));
%     err = dataset2mat(DS(:, end)) - dataset2mat(y(:, end));
%     RMSEml = sqrt(mean(err.^2));
%     gmm.setFisType('sugeno');
%     y = gmm.recall(DS(:,1));
%     err = dataset2mat(DS(:, end)) - dataset2mat(y(:, end));
%     RMSEts = sqrt(mean(err.^2));
%     disp([RMSEml RMSEts]);
% end
% fuzzyLogicDesigner(gmm.toFIS());

%%

% dumpname = 'dumps/SP500/01-May-2020 22:07:24.copy.mat';
% dumpname = 'dumps/SP500/01-May-2020 22:07:24.mat';
% self = load(dumpname);
% 
% offset = struct( ...
%     'log2delta', 1,     ...
%     'log2tau', 1,       ... quando criar novas componentes
%     'log2tmax', 1,      ... tempo para assimilar instâncias
%     'log2maxNC', 1);    %%% número máximo de componentes estáveis (não espúrias)
% self.offset = struct2table(offset);
% 
% save(dumpname, '-struct', 'self');
% 
% dumpname = 'dumps/SP500/01-May-2020 22:07:24.mat';
% silfa = load(dumpname);
% fildes = fieldnames(silfa.dumps);
% for ii=1:length(fildes)
%     fild=fildes{ii};
%     fprintf('%s %i\n', fild, all([silfa.dumps.(fild)]==[self.dumps.(fild)]));
% end

%%

fixFolderDumps('dumps/SP500');

function fixFolderDumps(dumpFolder)
    allmatfiles = dir([dumpFolder '/*.mat']);
    for ii = 1:length(allmatfiles)
        thisFile = fullfile(dumpFolder, allmatfiles(ii).name);
        try
            self = load(thisFile);
            if isfield(self, 'self')
                self = self.self;
            end
            self = fixFieldsSeries(self);
            save(thisFile, '-struct', 'self');
        catch ME
            if strcmp(ME.message, 'aborta')
                disp(['Jumping ' thisFile]);
            elseif strcmp(ME.message, 'Integrity is OK')
                disp([thisFile ': ' ME.message]);
            else
                disp([thisFile ': ' ME.message]);
%                 rethrow(ME);
            end
        end
        clear self;
    end
end

function self = fixFieldsSeries(self)
    try
        need_fix = false;
        check_integrity(self);
    catch
        need_fix = true;
    end
    if ~need_fix
        error('Integrity is OK');
    end
    if ~isfield(self, 'fis_types')
        self.fis_types = {'sugeno'};
    end
    if ~isfield(self, 'save_fis')
        lenTrainDS = length(self.DS)-1;
        self.save_fis = (lenTrainDS-2800):200:lenTrainDS;
    end
    %%
    paramfilds = {'delta', 'tau', 'tmax', 'maxNC'};
    assert(all(ismember(paramfilds, fieldnames(self.dumps))), 'dumps dont have params');
    log2filds = strcat('log2', paramfilds);
    if ~any(isfield(self.dumps, log2filds))
        if  isfield(self, 'pso')
            self = rmfield(self, 'pso');
            assert(  isfield(self, 'offset') && isstruct(self.offset), 'offset from pso is not a struct');
        else,assert(~isfield(self, 'offset'), 'offset is already a field without log2filds on dumps');
        end
        for ii = 1:length(log2filds)
            log2value = log2([self.dumps.(paramfilds{ii})]);
            log2cell = num2cell(log2value);
            [self.dumps.(log2filds{ii})] = deal(log2cell{:});
            log2unique = unique(log2value);
            log2diff = unique(log2unique(1:end-1) - log2unique(2:end));
            if length(log2diff) ~= 1
                disp('log2fields are unevenly spaced');
                error('aborta');
            end
            if isfield(self, 'offset')
                assert(self.offset.(log2filds{ii}) == abs(log2diff), ...
                    'offset from pso differs from the computed');
            else,self.offset.(log2filds{ii}) = abs(log2diff);
            end
        end
        endy = length(fieldnames(self.dumps));
        self.dumps = orderfields(self.dumps, [endy-3:endy, 1:endy-4]);
        self.offset = struct2table(self.offset);
    end
    assert(all(isfield(self.dumps, log2filds)), 'log2fields are partially included on dumps');
    assert(isfield(self, 'offset'), 'offset is not a field of self');
    assert(istable(self.offset), 'offset is still a struct');
    assert(~isfield(self, 'pso'), 'pso is still a field from self');
    %%
    fild = 'batchsize';
    if isfield(self.dumps, fild)
        value = unique([self.dumps.(fild)]);
        if length(value) ~= 1
            disp(['More than 1 ' fild ' value on dumps']);
            error('aborta');
        end
        if isfield(self, [fild 's'])
            value = [value self.([fild 's'])];
        end
        if length(value) ~= length(unique(value))
            disp([fild ' is not unique']);
            error('aborta');
        end
        if ~issorted(value, 'descend')
            disp([fild ' is not sorted']);
            error('aborta');
        end
        self.([fild 's']) = value;
        self.dumps = rmfield(self.dumps, fild);
    end
    assert(isfield(self, 'batchsizes'), 'batchsizes is not a field of self');
    if isfield(self, 'initial_filter')
        if self.initial_filter
            self.batchsizes = self.batchsizes(2:end);
        end
        self = rmfield(self, 'initial_filter');
    end
    %%
    dumpdiff = {'normalize', 'save_stats', 'warmup'};
    for ii = 1:length(dumpdiff)
        fild = dumpdiff{ii};
        if isfield(self.dumps, fild)
            value = unique([self.dumps.(fild)]);
            if length(value) ~= 1
                disp(['More than 1 ' fild ' value']);
                error('aborta');
            end
            if ~isfield(self, fild),           self.(fild) = value;
            elseif length(self.(fild)) ~= 1 || self.(fild) ~= value
                disp([fild ' is already a field of self, and differ from dumps']);
                error('aborta');
            end
            self.dumps = rmfield(self.dumps, fild);
        end
    end
    assert(isfield(self, 'normalize'), 'normalize is not a field of self');
    if ~isfield(self, 'save_stats')
        self.save_stats = true;
    end
    if ~isfield(self, 'warmup')
        self.warmup = 0;
    end
    %%
    if isfield(self.dumps, 'mergeFS')
        assert(~isfield(self.dumps, 'Smerge'), 'Smerge is already a field on dumps');
        value = num2cell(1 - [self.dumps.mergeFS]./2);
        [self.dumps.Smerge] = deal(value{:});
        self.dumps = rmfield(self.dumps, 'mergeFS');
    end
    if ~isfield(self.dumps, 'Smerge')
        [self.dumps.Smerge] = deal(1);
    end
    %%
    if isfield(self.dumps, 'update_needed')
        idxs = [self.dumps.update_needed];
        [self.dumps(idxs).progress] = deal(NaN);
        [self.dumps(idxs).cputime] = deal(NaN);
        self.dumps = rmfield(self.dumps, 'update_needed');
    end
    %%
    check_integrity(self);
end

function check_integrity(self)
    selffields = {'DS','dumps','offset','normalize','warmup','batchsizes',...
        'fis_types','save_stats','save_fis'};
    selfdiff = setdiff(selffields, fieldnames(self));
    assert(isempty(selfdiff), ['Fields ' selfdiff{:} ' lacked on self']);
    selfdiff = setdiff(fieldnames(self), selffields);
    assert(isempty(selfdiff), ['Fields ' selfdiff{:} ' remained on self']);
    assert(istable(self.offset), 'offset is still a struct');
    %%
    paramfilds = {'delta', 'tau', 'tmax', 'maxNC'};
    log2filds = strcat('log2', paramfilds);
    dumpfields = [log2filds paramfilds {'spmin','progress','cputime','stats',...
        'mean_NC','RMS_NC','mamdani_RMSE','sugeno_RMSE','Smerge'}];
    dumpdiff = setdiff(dumpfields, fieldnames(self.dumps));
    assert(isempty(dumpdiff), ['Fields ' dumpdiff{:} ' lacked on dumps']);
    dumpdiff = setdiff(fieldnames(self.dumps), dumpfields);
    assert(isempty(dumpdiff), ['Fields ' dumpdiff{:} ' remained on dumps']);
end
