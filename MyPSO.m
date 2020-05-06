
classdef MyPSO < handle
    properties
        center;
        c_stats;
        offset;
        particle;
        
        normalize;
        warmup;
        batchsizes;
    end
    
    methods
        function [pso, dumps] = MyPSO(best_dump, offset, normalize, warmup, batchsizes)
            pso.normalize = normalize;
            pso.warmup = warmup;
            pso.batchsizes = batchsizes;
            best = struct2table(best_dump, 'AsArray', true);
            pso.c_stats = MyPSO.get_JCV(best);
            params = fieldnames(offset);
            pso.offset = struct2table(offset);
            pso.offset{:,:} = pso.offset{:,:}/2;
            pso.center = array2table( ...
                log2(best{:, params}), ...
                'VariableNames', params);
            pso.particle = array2table( pso.center{:,:} ...
                + pso.offset{:,:} .* permn([-1, 1], width(pso.center)), ...
                'VariableNames', params);
            dumps = pso.concat_particles(best_dump);
        end
        
        function [pso, dumps] = update(pso, dumps)
            if length(pso.batchsizes) == 1 && isnan(pso.batchsizes)
                throw(MException(['MATLAB:' 'MyPSO' ':StopIteration'], ...
                    'pso.batchsizes is empty'));
            end
            comb = height(pso.particle);
            c_update = dumps(end-comb);             p_update = dumps(end-(comb-1):end);
            MyPSO.check_same(c_update, pso.center); MyPSO.check_same(p_update, pso.particle);
            pso.c_stats = MyPSO.get_JCV(c_update);  p_stats = MyPSO.get_JCV(p_update);
            clear c_update p_update;
            for i=1:comb
                if  pso.have_improved(p_stats(i,:))
                    pso.center  = pso.particle(i,:);
                    pso.c_stats = p_stats(i,:);
                end
            end
            if ~isempty(pso.batchsizes)
                pso.offset{:,:} = pso.offset{:,:}/2;
                pso.particle = array2table(pso.center{:,:} ...
                    + pso.offset{:,:} .* permn([-1, 1], width(pso.center)), ...
                    'VariableNames', pso.center.Properties.VariableNames);
                dumps = pso.concat_particles(dumps);
            else
                dumps = [dumps; ...
                    INFGMN_series.hola( pso.center, pso.normalize, 0, 1, true, false ); ...
                    INFGMN_series.hola( pso.center, pso.normalize, 0, 1, true, true ) ...
                    ];
                pso.batchsizes = NaN;
            end
        end
    end

    methods(Access = private)
        function dumps = concat_particles(pso, dumps)
            dumps = [dumps; ...
                INFGMN_series.hola( pso.center, pso.normalize, ...
                pso.warmup, pso.batchsizes(1), true, false ); ...
                INFGMN_series.hola( pso.particle, pso.normalize, ...
                pso.warmup, pso.batchsizes(1), true, false ) ...
                ];
            fprintf('new particles with batchsize = %i\n', pso.batchsizes(1));
            pso.batchsizes(1) = [];
        end
        
        function bool = have_improved(pso, stats)
            if pso.c_stats.CV >= 1 || stats.CV >= 1
                bool = stats.CV < pso.c_stats.CV;
            else
                bool = stats.J  < pso.c_stats.J;
            end
        end
    end

    methods(Static)
        function JCV = get_JCV(last_dump)
            J = [last_dump.sugeno_RMSE]';
            J = J .* (2 .^ (1 - min(1, [last_dump.mean_NC]'/7)));
            J(isnan(J)) = Inf;
            CV = [last_dump.RMS_NC]' ./ [last_dump.maxNC]';
            CV(isnan(CV)) = Inf;
            CV([last_dump.progress] < 1) = Inf;
            JCV = table(J, CV, 'VariableNames', {'J', 'CV'});
        end
        
        function check_same(p_update, my_particle)
            same_particle = struct2table(p_update);
            same_particle = same_particle(:, my_particle.Properties.VariableNames);
            same_particle{:,:} = log2(same_particle{:,:});
            if ~all(my_particle{:,:} == same_particle{:,:})
                throw(MException(['MATLAB:' 'MyPSO' ':expectedEqualParams'], ...
                    'it should update the stats of my_particle'));
            end
        end
    end
end
