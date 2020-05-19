
classdef INFGMN_series < handle
    
    properties
%         DS;
%         dumps;
        dumpname;
    end
    
    methods(Static)
        function dumps = hola(dumps, mergeFS)
            dumps.delta         = 2 .^ dumps.log2delta;
            dumps.tau           = 2 .^ dumps.log2tau;
            dumps.tmax          = 2 .^ dumps.log2tmax;
            dumps.maxNC         = 2 .^ dumps.log2maxNC;
            dumps.spmin         = 2 .^ (dumps.log2tmax - dumps.log2maxNC);
            dumps.cputime       = zeros(height(dumps), 1);
            dumps.progress      =  ones(height(dumps), 1);
            dumps.update_needed =  true(height(dumps), 1);
            dumps.stats         =  cell(height(dumps), 1);
            dumps.mean_NC       =   NaN(height(dumps), 1);
            dumps.RMS_NC        =   NaN(height(dumps), 1);
            dumps.mamdani_RMSE  =   NaN(height(dumps), 1);
            dumps.sugeno_RMSE   =   NaN(height(dumps), 1);
            dumps.mergeFS       = false(height(dumps), 1) | mergeFS;
            dumps = table2struct(dumps);
        end
    end
    
    methods
        
        function obj = INFGMN_series(subfolder, offset)
            if isstruct(offset)
                self.offset = struct2table(offset);
                folder = ['dumps/' subfolder];
                mkdir(folder);
                obj.dumpname = [folder '/' datestr(now) '.mat'];
                obj.save_myself(self, false);
            elseif isfile(subfolder)
                obj.dumpname = subfolder;
            else
                throw(MException(['MATLAB:' 'INFGMN_series'], ...
                    [subfolder ' file does not exist.']));
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
                pause(1);
                waitbar(1/3, ff, 'Please wait...');
                pause(1);
                waitbar(2/3, ff, 'Saving myself.');
                save(obj.dumpname, '-struct', 'self');
                waitbar(1, ff, 'Done.');
                close(ff);
            end
        end
        
        %% INFGMN_series
        
        function obj = create(obj, DS, warmup, batchsizes, initial_filter, ...
                fis_types, save_stats, save_fis, normalize, ...
                delta_list, tau_list, tmax_list, maxNC_list)
            self = obj.myself();
            self.DS = DS;
            self.warmup = warmup;
            self.batchsizes = batchsizes;
            self.fis_types = fis_types;
            self.save_stats = save_stats;
            self.save_fis = save_fis;
            self.normalize = normalize;
            self.initial_filter = initial_filter;
            %% https://www.mathworks.com/matlabcentral/answers/433424-how-to-get-the-combinations-of-elements-of-two-arrays
            C = {delta_list', tau_list', tmax_list', maxNC_list'};
            D = C;
            [D{:}] = ndgrid(C{:});
            Z = cell2mat(cellfun(@(m)m(:),D,'uni',0));
            particles = array2table(Z, ...
                'VariableNames', {'log2delta', 'log2tau', 'log2tmax', 'log2maxNC'});
            self.dumps = INFGMN_series.hola(particles, false);
            fprintf('%i new particles with batchsize = %i\n', length(self.dumps), self.batchsizes(1));
            obj.save_myself(self, true);
        end

        function self = update(obj, check_break)
            self = obj.myself();
            count_updates = sum([self.dumps.update_needed]);
            fprintf('%i combinations to be updated\n', count_updates);
            fprintf(' step/comb.:  step time  | estimated remain time\n');
            count_i = 0;
            sumsteptime = 0;
            next_save = cputime + 60;
            for ii = 1:length(self.dumps)
                if ~self.dumps(ii).update_needed
                    continue
                end
                if cputime > next_save
                    saveref = cputime;
                    obj.save_myself(self, true);
                    savetime = cputime - saveref;
                    next_save = cputime + 60 + savetime * 20;
                end
                try
                    self.dumps(ii) = obj.step(self, self.dumps(ii), self.DS, check_break);
                catch ME
                    if (strcmp(ME.message,'Ill Conditioned Covariance.'))
                        self.dumps(ii).progress = 0;
                        self.dumps(ii).update_needed = false;
                        count_updates = count_updates - 1;
                        disp([ME.message ' Proceeding to the next iteration.']);
                        continue
                    else
                        disp(ME);
                        rethrow(ME);
                    end
                end
                count_i = count_i + 1;
                sumsteptime = sumsteptime + self.dumps(ii).cputime;
                fprintf('% 5i/% 5i: % 7.2f sec | % 8.2f min\n', ...
                    count_i, count_updates, self.dumps(ii).cputime, ...
                    sumsteptime / count_i * (count_updates - count_i) / 60);
            end
            fprintf('Elapsed time: %.2f min\n', sumsteptime/60);
            obj.save_myself(self, true);
        end
        
        function dumps_ii = step(~, self, dumps_ii, DS, check_break)
            gmm = INFGMN(minmaxDS(DS), 'normalize', self.normalize, ...
                'delta', dumps_ii.delta, 'tau',  dumps_ii.tau, ...
                'tmax',  dumps_ii.tmax, 'spmin', dumps_ii.spmin );
            gmm.setMergeFS(dumps_ii.mergeFS); %%
            warm_up = self.warmup + mod(((length(DS)-1)-self.warmup), self.batchsizes(1));
            n_batches = ((length(DS)-1)-warm_up) / self.batchsizes(1);
            assert(mod(n_batches,1)==0, 'n_batches is not integer');
            stats = struct( ...
                't', cell(1, n_batches), ...
                'NCs', [], ...
                'expected_output', NaN, ...
                'mamdani_output', NaN, 'mamdani_err', NaN, ...
                'sugeno_output', NaN, 'sugeno_err', NaN, ...
                'cputime', NaN, ...
                'fis', NaN );
            window = 2 * ceil(dumps_ii.tmax);
            batches_in_a_window = ceil(window / self.batchsizes(1));
            if warm_up
                [gmm, warmup_NCs] = gmm.train(DS(1:warm_up, :));
            end
            for ii = 1:n_batches
                timeref = cputime;
                tt = warm_up + ii * self.batchsizes(1);
                stats(ii).t = tt;
                [gmm, stats(ii).NCs] = gmm.train( ...
                    DS((tt-(self.batchsizes(1)-1)):tt, :) );
                stats(ii).expected_output = dataset2mat(DS(tt+1, end));
                for type = self.fis_types
                    gmm.setFisType(type{1});
                    output = gmm.recall(DS(tt+1, 1:end-1));
                    output = dataset2mat(output);
                    stats(ii).([type{1} '_output']) = output;
                    stats(ii).([type{1} '_err']) = output - stats(ii).expected_output;
                end
                if  self.save_stats && ismember(tt, self.save_fis) ...
                        && self.warmup == 0 && self.batchsizes(1) == 1
                    stats(ii).fis = gmm.toFIS();
                end
                stats(ii).cputime = cputime - timeref;
                if  ii >= 1 + batches_in_a_window ...
                        && any( dumps_ii.maxNC < movmean( ...
                                [stats((ii-batches_in_a_window):ii).NCs], window, ...
                                'Endpoints', 'discard') ) ...
                        || check_break(dumps_ii, stats(1:ii))
                    stats = stats(1:ii);
                    dumps_ii.progress = tt / (length(DS)-1);
                    break
                end
            end
            if warm_up
                stats(1).NCs = [warmup_NCs stats(1).NCs];
            end
            dumps_ii.update_needed = false;
            dumps_ii.cputime    = sum( [stats.cputime] );
            dumps_ii.mean_NC    =        mean(  [stats.NCs],    'omitnan');
            dumps_ii.RMS_NC     =   sqrt(mean(  [stats.NCs].^2, 'omitnan'));
            dumps_ii.mamdani_RMSE = sqrt(mean(  [stats.mamdani_err].^2, 'omitnan'));
            dumps_ii.sugeno_RMSE  = sqrt(mean(  [stats.sugeno_err] .^2, 'omitnan'));
            if self.save_stats
                dumps_ii.stats = stats;
            end
        end
        
        %% PSO
        
        function pso_update(obj)
            self = obj.myself();
            if  any([self.dumps.update_needed])
                disp('there are dumps to update yet');
                return;
            elseif length(self.batchsizes) == 1
                if any([self.dumps.mergeFS])
                    error('StopIteration');
                end
                obj.etapa_final();
                return;
            elseif self.initial_filter
                self.batchsizes(1) = [];
                self.initial_filter = false;
                [self.dumps([self.dumps.progress]==1).update_needed] = deal(true);
                obj.save_myself(self, true);
                return;
            end
            self.batchsizes(1) = [];
            self.offset{:,:} = self.offset{:,:}/2;
            params = self.offset.Properties.VariableNames;
            best = struct2table(INFGMN_series.get_best(self.dumps), ...
                'AsArray', true);
            particles = array2table( best{:,params} ...
                + self.offset{:,:} .* permn([-1, 0, 1], width(self.offset)), ...
                'VariableNames', params);
            self.dumps = INFGMN_series.hola( particles, false );
            fprintf('%i new particles with batchsize = %i\n', length(self.dumps), self.batchsizes(1));
            obj.dumpname = [obj.dumpname '_PSO' ...
                num2str(table2array(self.offset), '%+g') '.mat'];
            obj.save_myself(self, false);
        end
        
        function etapa_final(obj)
            self = obj.myself();
%             self.dumps = self.dumps([self.dumps.progress]==1);
            self.dumps = INFGMN_series.get_best(self.dumps);
            assert(~isempty(self.save_fis), 'save_fis is empty');
            if self.warmup ~= 0 || self.batchsizes ~= 1 || ~self.save_stats
                self.save_stats = true;
                self.warmup = 0;
                self.batchsizes = 1;
                [self.dumps.update_needed] = deal(true);
            end
            clone = self.dumps;
            [clone.mergeFS] = deal(true);
            [clone.update_needed] = deal(true);
            self.dumps = [ self.dumps; clone ];
            obj.dumpname = [obj.dumpname '_final.mat'];
            obj.save_myself(self, false);
        end
        
    end

    methods(Static)
        
        function JCV = get_JCV(last_dump)
            J = [last_dump.sugeno_RMSE]';
%             J = J .* (2 .^ (1 - min(1, [last_dump.mean_NC]'/7)));
            J(isnan(J)) = Inf;
            CV = [last_dump.RMS_NC]' ./ [last_dump.maxNC]';
            CV(isnan(CV)) = Inf;
            CV([last_dump.progress] < 1) = Inf;
            JCV = table(J, CV, 'VariableNames', {'J', 'CV'});
        end
        
        function best = get_best(dumps)
            JCV = INFGMN_series.get_JCV(dumps);
            viaveis = (JCV.CV < 1);
            if any(viaveis)
                dumps_viaveis = dumps(viaveis);
                JCV = JCV(viaveis, :);
                [~,index] = min(JCV.J);
                best = dumps_viaveis(index);
            else
                [~,index] = min(JCV.CV);
                best = dumps(index);
            end
        end
        
    end
end
