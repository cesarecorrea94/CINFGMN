
classdef INFGMN_series < handle
    
    properties
%         DS;
%         dumps;
        dumpname;
    end
    
    methods(Static)
        function dumps = hola(dumps)
            dumps.cputime       =   NaN(height(dumps), 1);
            dumps.stats         =  cell(height(dumps), 1);
            dumps.mean_NC       =   NaN(height(dumps), 1);
            dumps.RMS_NC        =   NaN(height(dumps), 1);
            dumps.mamdani_RMSE  =   NaN(height(dumps), 1);
            dumps.sugeno_RMSE   =   NaN(height(dumps), 1);
            dumps = table2struct(dumps);
        end
        function dumps = hola_nonseries(dumps)
            dumps.cputime       =   NaN(height(dumps), 1);
            dumps.stats         =  cell(height(dumps), 1);
            dumps.mean_NC       =   NaN(height(dumps), 1);
            dumps.mamdani_RMSE  =   NaN(height(dumps), 1);
            dumps.sugeno_RMSE   =   NaN(height(dumps), 1);
            dumps = table2struct(dumps);
        end
    end
    
    methods
        
        function obj = INFGMN_series(subfolder)
            if isfile(subfolder)
                obj.dumpname = subfolder;
            else
                folder = ['dumps/' subfolder];
                mkdir(folder);
                obj.dumpname = [folder '/' datestr(now) '.mat'];
            end
        end
        
        function self = myself(obj)
            self = load(obj.dumpname);
        end
        
        function self = save_myself(obj, self, overwrite)
            if (isfile(obj.dumpname) && ~overwrite)
                disp([obj.dumpname ' already exists.' ...
                    ' Lets proceed with that file.']);
            else
                ff = waitbar(0, 'Please wait.');
                pause(0.5);
                waitbar(1/3, ff, 'Please wait...');
                pause(0.5);
                waitbar(2/3, ff, 'Saving myself.');
                [filepath,name,ext] = fileparts(obj.dumpname);
                tmpname = fullfile(filepath, ['.TMP_' name ext]);
                save(tmpname, '-struct', 'self');
                [successful,msg] = movefile(tmpname, obj.dumpname);
                if ~successful, error(msg); end
                waitbar(1, ff, 'Done.');
                close(ff);
            end
        end
        
        function obj = create(obj, DS, warmup, batchsize, ...
                save_stats, fis_types, save_fis, ...
                doMerge, Smerge, Sdeath, maxFoCSize, ...
                normalize, paramstruct, offset)
            self.DS = DS;
            self.warmup = warmup;
            self.batchsize = batchsize;
            self.fis_types = fis_types;
            self.save_stats = save_stats;
            self.save_fis = save_fis;
            self.doMerge = doMerge;
            self.Smerge = Smerge;
            self.Sdeath = Sdeath;
            self.maxFoCSize = maxFoCSize;
            self.normalize = normalize;
            self.offset = struct2table(offset);
            particles = INFGMN_series.get_particles(paramstruct);
            self.dumps = INFGMN_series.hola(particles);
            fprintf('%i new particles with batchsize = %i\n', length(self.dumps), self.batchsize);
            obj.save_myself(self, true);
        end

        function obj = create_nonseries(obj, DS, ...
                fis_types, save_fis, doMerge, maxFoCSize, ...
                normalize, paramstruct, offset)
            self.recallTest = true;
            self.Smerge = 0.9;
            self.DS = DS;
            self.fis_types = fis_types;
            self.save_fis = save_fis;
            self.doMerge = doMerge;
            self.maxFoCSize = maxFoCSize;
            self.normalize = normalize;
            self.offset = struct2table(offset);
            particles = INFGMN_series.get_particles(paramstruct);
            self.dumps = INFGMN_series.hola_nonseries(particles);
            fprintf('%i new particles\n', length(self.dumps));
            obj.save_myself(self, true);
        end

        %% update
        
        function update(obj, core_step_function)
            self = obj.myself();
            uptodate = sum(~isnan([self.dumps.cputime]));
            len_dumps = length(self.dumps);
            fprintf('%i combinations to be updated\n', len_dumps - uptodate);
            if len_dumps == uptodate, return; end
            fprintf(' step/comb.:   step  time   | estimated remain time\n');
            sumsteptime = sum([self.dumps.cputime], 'omitnan');
            next_save = cputime + 2*60;
            for ii = 1:len_dumps
                if ~isnan(self.dumps(ii).cputime)
                    continue
                end
                if cputime > next_save
                    saveref = cputime;
                    obj.save_myself(self, true);
                    savetime = cputime - saveref;
                    next_save = cputime + 5*60 + savetime * 20;
                end
                timeref = cputime;
                err_msg = '';
                try
                    self.dumps(ii) = core_step_function(self, self.dumps(ii));
                catch ME
                    if strcmp(ME.message, 'Ill Conditioned Covariance.') ...
                            || strcmp(ME.identifier, 'INFGMN:Abortar')
                        err_msg = [' - ' ME.message];
                    else
                        disp(ME);
                        rethrow(ME);
                    end
                end
                self.dumps(ii).cputime = cputime - timeref;
                sumsteptime = sumsteptime + self.dumps(ii).cputime;
                uptodate = uptodate + 1;
                
                fprintf('% 5i/% 5i: % 10s sec | % 15s sec %s\n', ...
                    uptodate, len_dumps, ...
                    duration(0,0, self.dumps(ii).cputime, 'Format', 'mm:ss.SS' ), ...
                    duration(0,0, (len_dumps - uptodate) * (sumsteptime / uptodate), ...
                        'Format', 'dd:hh:mm:ss' ), ...
                    err_msg);
            end
            fprintf('Elapsed time: %s sec\n', duration(0,0, sumsteptime, 'Format', 'dd:hh:mm:ss') );
            obj.save_myself(self, true);
        end
        
    end
    
    methods(Static)
        %% step
        function dumps_ii = step_nonseries(self, dumps_ii)
            comb_test = combnk(1:5, 2);
            n_stats = size(comb_test,1);
            stats = struct( ...
                'cputime', NaN, ...
                'NC', cell(1, n_stats), ...
                'mamdani_RMSE', NaN, ...
                'sugeno_RMSE', NaN, ...
                'fis', NaN );
            nsamples = size(self.DS,1);
            cellDS = cell(1, 5);
            for i_split = 1:5
                cellDS{i_split} = self.DS(round(nsamples*(i_split-1)/5)+1:round(nsamples*i_split/5), :);
            end
            for i_comb = 1:size(comb_test,1)
                comb_train = setdiff(1:5, comb_test(i_comb,:));
                train = vertcat(cellDS{comb_train});
                rng(0);
                train = train(randperm(size(train,1)), :);
                if self.recallTest
                    test  = vertcat(cellDS{comb_test(i_comb,:)});
                else,test = train;
                end
                expected_output = dataset2mat(test(:, end));
                % %
                gmm = INFGMN_series.create_INFGMN(self, dumps_ii);
                timeref = cputime;
                gmm.train( train );
                % gmm.train( train );
                gmm.setMergeStats(self.doMerge, ...
                    self.Smerge, 0.8, self.maxFoCSize);
                for type = self.fis_types
                    gmm.setFisType(type{1});
                    output = gmm.recall( test(:, 1:end-1) );
                    output = dataset2mat(output);
                    stats(i_comb).([type{1} '_RMSE']) = sqrt(mean((output - expected_output).^2));
                end
                stats(i_comb).NC = gmm.modelSize();
                if self.save_fis
                    stats(i_comb).fis = gmm.toFIS();
                end
                stats(i_comb).cputime = cputime - timeref;
                delete(gmm);
            end % train,test
            dumps_ii.mean_NC = mean([stats.NC], 'omitnan');
            dumps_ii.mamdani_RMSE = mean([stats.mamdani_RMSE],  'omitnan');
            dumps_ii.sugeno_RMSE  = mean([stats.sugeno_RMSE],   'omitnan');
            dumps_ii.stats = stats;
        end
        
        function dumps_ii = step(self, dumps_ii)
%             if ~isfield(self, 'Smerge')
%                 self.Smerge = 0.7;
%             end
            gmm = INFGMN_series.create_INFGMN(self, dumps_ii);
            lenTrainDS = length(self.DS)-1;
            warm_up = self.warmup + mod((lenTrainDS - self.warmup), self.batchsize);
            n_batches = (lenTrainDS - warm_up) / self.batchsize;
            assert(mod(n_batches,1)==0, 'n_batches is not integer');
            stats = struct( ...
                'cputime', NaN, ...
                't', NaN, ...
                'NCs', cell(1, n_batches), ...
                'expected_output', NaN, ...
                'mamdani_err', NaN, ...
                'sugeno_err', NaN, ...
                'fis', NaN );
            if warm_up
                [gmm, warmup_NCs] = gmm.train(self.DS(1:warm_up, :));
            end
            gmm.setMergeStats(self.doMerge, ...
                self.Smerge, self.Sdeath, self.maxFoCSize);
            for ii = 1:n_batches
                timeref = cputime;
                tt = warm_up + ii * self.batchsize;
                stats(ii).t = tt;
                batch = (tt-(self.batchsize-1)):tt;
                [gmm, stats(ii).NCs] = gmm.train( self.DS(batch, :) );
                stats(ii).expected_output = dataset2mat(self.DS(tt+1, end));
                for type = self.fis_types
                    gmm.setFisType(type{1});
                    output = gmm.recall(self.DS(tt+1, 1:end-1));
                    output = dataset2mat(output);
                    stats(ii).([type{1} '_err']) = output - stats(ii).expected_output;
                end
                if  self.save_stats && ismember(tt, self.save_fis) ...
                        && self.warmup == 0 && self.batchsize == 1
                    stats(ii).fis = gmm.toFIS();
                end
                stats(ii).cputime = cputime - timeref;
            end
            if warm_up
                stats(1).NCs = [warmup_NCs stats(1).NCs];
            end
            dumps_ii.mean_NC    =        mean(  [stats.NCs],    'omitnan');
            dumps_ii.RMS_NC     =   sqrt(mean(  [stats.NCs].^2, 'omitnan'));
            dumps_ii.mamdani_RMSE = sqrt(mean(  [stats.mamdani_err].^2, 'omitnan'));
            dumps_ii.sugeno_RMSE  = sqrt(mean(  [stats.sugeno_err] .^2, 'omitnan'));
            if self.save_stats
                dumps_ii.stats = stats;
            end
            delete(gmm);
        end
        
    end
    
    methods
        function explore(obj)
            self = obj.myself();
            if  any(isnan([self.dumps.cputime]))
                disp('there are dumps to update yet');
                return;
            end
            self.offset{:,:} = self.offset{:,:}/2;
            params = self.offset.Properties.VariableNames;
            self.dumps = INFGMN_series.get_best(self.dumps);
            best = struct2table(self.dumps, 'AsArray', true);
            best = best(:, params);
            particles = best;
            for irow = 1:height(best)
                particles = [ particles; ...
                    array2table( best{irow, params} ...
                        + self.offset{:, params} .* permn([-1, 0, 1], width(self.offset)), ...
                        'VariableNames', params) ...
                    ];
            end
            particles = unique(particles, 'rows');
%             particles = setdiff(particles, best, 'rows');
%             self.dumps = [ self.dumps; INFGMN_series.hola( particles ) ];
            self.dumps = INFGMN_series.hola( particles );
            fprintf('%i new particles\n', sum(isnan([self.dumps.cputime])) );
            [filepath,name,ext] = fileparts(obj.dumpname);
            obj.dumpname = fullfile(filepath, [name '.explore' ext]);
            obj.save_myself(self, false);
        end
        
        function deepen(obj)
            self = obj.myself();
            if  any(isnan([self.dumps.cputime]))
                disp('there are dumps to update yet');
                return;
            end
            for divisor = [2, 3, 2.5, 3.5, NaN]
                if isnan(divisor)
                    error('StopIteration');
                elseif mod(self.batchsize, divisor) == 0
                    self.batchsize = self.batchsize / divisor;
                    break;
                end
            end
            best = INFGMN_series.filter_best_params( self );
            self.dumps = INFGMN_series.hola( best(1:height(best)/divisor, :) );
            [filepath,name,ext] = fileparts(obj.dumpname);
            obj.dumpname = fullfile(filepath, [name '.batch' num2str(self.batchsize) ext]);
            obj.save_myself(self, false);
        end
        
        function self = etapa_final(~, self)
%             self.dumps = self.dumps([self.dumps.progress]==1);
            self.dumps = INFGMN_series.get_best(self.dumps);
            assert(~isempty(self.save_fis), 'save_fis is empty');
            if self.warmup ~= 0 || self.batchsize ~= 1 || ~self.save_stats ...
                    || ~isempty(setdiff({'mamdani', 'sugeno'}, self.fis_types))
                self.save_stats = true;
                self.fis_types = {'mamdani', 'sugeno'};
                self.warmup = 0;
                self.batchsize = 1;
                [self.dumps.cputime ] = deal(NaN);
            end
            clone = self.dumps;
            [clone.cputime] = deal(NaN);
            self.dumps = [ self.dumps; clone ];
            fprintf('%i particles on final step\n', length(self.dumps));
        end
        
    end

    methods(Static)
        
        function particles = get_particles(paramstruct)
            %% https://www.mathworks.com/matlabcentral/answers/433424-how-to-get-the-combinations-of-elements-of-two-arrays
            C = struct2cell(paramstruct)';
%             C = cellfun(@transpose,C,'UniformOutput',false);
            D = C;
            [D{:}] = ndgrid(C{:});
            Z = cell2mat(cellfun(@(m)m(:),D,'uni',0));
            particles = array2table(Z, 'VariableNames', fieldnames(paramstruct));
            particles = particles(randperm(height(particles)), :);
        end
        
        function [igmm, maxNC] = create_INFGMN(self, dumps_ii)
            log2spmin = dumps_ii.log2tmax - dumps_ii.log2maxNC;
            igmm = INFGMN(minmaxDS(self.DS), ...
                'normalize', self.normalize, ...
                'delta',        2^ dumps_ii.log2delta, ...
                'tau',          2^ dumps_ii.log2tau, ...
                'tmax', round(  2^ dumps_ii.log2tmax ), ...
                'spmin',        2^          log2spmin ...
                );
            maxNC = 2^ dumps_ii.log2maxNC;
        end
        
        function foremost = get_best(dumps)
            J = [dumps.sugeno_RMSE];
            viaveis = ~(isinf(J) | isnan(J));
            if any(viaveis)
                dumps = dumps(viaveis);
                J = J(viaveis);
%                 [~,index] = min(J);
                [~, indexes] = sort(J);
                foremost = dumps(indexes);
            else
                error('There is no good dump');
            end
        end
        
        function best = filter_best_params(self)
            params = self.offset.Properties.VariableNames;
            best = INFGMN_series.get_best(self.dumps);
            best = struct2table(best, 'AsArray', true);
            best = best(:, params);
        end
    end
end
