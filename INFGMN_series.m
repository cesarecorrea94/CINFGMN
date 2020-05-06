
classdef INFGMN_series < handle
    
    properties
%         DS;
%         dumps;
        dumpname;
    end
    
    methods(Static)
        function dumps_ii = hola(dumps_ii, normalize, ...
                warmup, batchsize, savestats, mergeFS)
            dumps_ii{:,:} = 2.^(dumps_ii{:,:});
            dumps_ii.delta          = max(0, min(1, dumps_ii.delta));
            dumps_ii.tau            = max(0, min(1, dumps_ii.tau));
            dumps_ii.tmax           = max(1, dumps_ii.tmax);
            dumps_ii.maxNC          = max(1, dumps_ii.maxNC);
            dumps_ii.spmin          = dumps_ii.tmax ./ dumps_ii.maxNC;
            dumps_ii.normalize      =  true(height(dumps_ii), 1) & normalize;
            dumps_ii.warmup         =  ones(height(dumps_ii), 1) * warmup;
            dumps_ii.batchsize      =  ones(height(dumps_ii), 1) * batchsize;
            dumps_ii.cputime        = zeros(height(dumps_ii), 1);
            dumps_ii.progress       =  ones(height(dumps_ii), 1);
            dumps_ii.update_needed  =  true(height(dumps_ii), 1);
            dumps_ii.save_stats     =  true(height(dumps_ii), 1) & savestats;
            dumps_ii.stats          =  cell(height(dumps_ii), 1);
            dumps_ii.mean_NC        =   NaN(height(dumps_ii), 1);
            dumps_ii.RMS_NC         =   NaN(height(dumps_ii), 1);
            dumps_ii.mamdani_RMSE   =   NaN(height(dumps_ii), 1);
            dumps_ii.sugeno_RMSE    =   NaN(height(dumps_ii), 1);
            dumps_ii.mergeFS        = false(height(dumps_ii), 1) | mergeFS;
            dumps_ii = table2struct(dumps_ii);
        end
    end
    
    methods
        
        function obj = INFGMN_series(subfolder, offset)
            if isstruct(offset)
                self.offset = offset;
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
        
        function obj = create(obj, DS, warmup, batchsizes, normalize, ...
                delta_list, tau_list, tmax_list, maxNC_list)
            self = obj.myself();
            self.DS = DS;
            self.warmup = warmup;
            self.batchsizes = batchsizes;
            self.normalize = normalize;
            %% https://www.mathworks.com/matlabcentral/answers/433424-how-to-get-the-combinations-of-elements-of-two-arrays
            C = {delta_list', tau_list', tmax_list', maxNC_list'};
            D = C;
            [D{:}] = ndgrid(C{:});
            Z = cell2mat(cellfun(@(m)m(:),D,'uni',0));
            particles = array2table(Z, ...
                'VariableNames', {'delta', 'tau', 'tmax', 'maxNC'});
            self.dumps = INFGMN_series.hola(particles, self.normalize, ...
                self.warmup, self.batchsizes(1), false, false );
            fprintf('new particles with batchsize = %i\n', self.batchsizes(1));
            self.batchsizes(1) = [];
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
                    self.dumps(ii) = obj.step(self.dumps(ii), self.DS, check_break);
                catch ME
                    if (strcmp(ME.message,'Ill Conditioned Covariance.'))
                        self.dumps(ii).progress = 0;
                        self.dumps(ii).update_needed = false;
                        count_updates = count_updates - 1;
                        disp([ME.message ' Proceeding to the next iteration.']);
                        continue
                    else
                        ME
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
        
        function dumps_ii = step(~, dumps_ii, DS, check_break)
            gmm = INFGMN(minmaxDS(DS), 'normalize', dumps_ii.normalize, ...
                'delta', dumps_ii.delta, 'tau', dumps_ii.tau, ...
                'tmax', dumps_ii.tmax, 'spmin', dumps_ii.spmin );
            gmm.setMergeFS(dumps_ii.mergeFS); %%
            n_batches = floor((length(DS)-1)/dumps_ii.batchsize);
            stats = struct( ...
                't', cell(1, n_batches), ...
                'NCs', [], 'mean_NC', NaN, 'RMS_NC', NaN, ...
                'expected_output', NaN, ...
                'mamdani_output', NaN, 'mamdani_err', NaN, ...
                'sugeno_output', NaN, 'sugeno_err', NaN, ...
                'cputime', NaN, ...
                'fis', NaN );
            window = 2 * ceil(dumps_ii.tmax);
            batches_in_a_window = ceil(window/dumps_ii.batchsize);
            resto = mod((length(DS)-1), dumps_ii.batchsize);
            if resto, gmm.train(DS(1:resto, :));
            end
            for ii = 1:n_batches
                timeref = cputime;
                tt = resto + ii * dumps_ii.batchsize;
                stats(ii).t = tt;
                [gmm, stats(ii).NCs] = gmm.train( ...
                    DS((tt-(dumps_ii.batchsize-1)):tt, :) );
                stats(ii).mean_NC =     mean( [stats(ii).NCs] );
                stats(ii).RMS_NC = sqrt(mean( [stats(ii).NCs].^2 ));
                if tt > dumps_ii.warmup
                    stats(ii).expected_output = dataset2mat(DS(tt+1, end));
                    %% alteração no 'for' a ser removida posteriormente
                    for type = {'sugeno'}%{'mamdani', 'sugeno'} %% dumps_ii.FSTypes
                        gmm.setFisType(type{1});
                        output = gmm.recall(DS(tt+1, 1:end-1));
                        output = dataset2mat(output);
                        stats(ii).([type{1} '_output']) = output;
                        stats(ii).([type{1} '_err']) = output - stats(ii).expected_output;
                    end
                    if  dumps_ii.warmup == 0 && dumps_ii.batchsize == 1 ...
                            && tt > 12000 && mod(tt-resto, 100) == 0
                        stats(ii).fis = gmm.toFIS();
                    end
                end
                stats(ii).cputime = cputime - timeref;
                if  ii >= 1 + batches_in_a_window ...
                        && any(movmean( [stats((ii-batches_in_a_window):ii).NCs], window, ...
                                'Endpoints', 'discard') > dumps_ii.maxNC ) ...
                        || check_break(dumps_ii, stats(1:ii))
                    stats = stats(1:ii);
                    dumps_ii.progress = tt / (length(DS)-1);
                    break
                end
            end
            dumps_ii.update_needed = false;
            dumps_ii.cputime    = sum( [stats.cputime] );
            dumps_ii.mean_NC    =       mean( [stats.mean_NC]  , 'omitnan');
            dumps_ii.RMS_NC     = sqrt( mean( [stats.RMS_NC].^2, 'omitnan'));
            dumps_ii.mamdani_RMSE = sqrt(mean( [stats.mamdani_err].^2, 'omitnan'));
            dumps_ii.sugeno_RMSE  = sqrt(mean( [stats.sugeno_err] .^2, 'omitnan'));
            if dumps_ii.save_stats
                dumps_ii.stats = stats;
            end
        end
        
        %% PSO
        
        function pso_update(obj)
            self = obj.myself();
            if  any([self.dumps.update_needed])
                disp('there are dumps to update yet');
            else
                if  isfield(self, 'pso')
                    [self.pso, self.dumps] = self.pso.update(self.dumps);
                    obj.save_myself(self, true);
                else
                    best = obj.get_best(self.dumps);
                    [self.pso, self.dumps] = MyPSO(best, self.offset, ...
                        self.normalize, self.warmup, self.batchsizes);
                    obj.dumpname = [obj.dumpname '_PSO' ...
                        num2str(table2array(self.pso.center), '%+g') '.mat'];
                    obj.save_myself(self, false);
                end
            end
        end
        
        function best = get_best(~, dumps)
            JCV = MyPSO.get_JCV(dumps);
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
