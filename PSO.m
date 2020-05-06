
classdef PSO < handle
    properties
        w;
        psi1;
        psi2;
        
        best_stats;
        best_p;
        best_i;
        h;
    end
    
    methods
        function pso = PSO(w, psi1, psi2, np, params)
            pso.w = w;
            pso.psi1 = psi1;
            pso.psi2 = psi2;
            
            pso.best_p = array2table( ...
                zeros(np, length(params)), ...
                'VariableNames', params);
            pso.h = pso.best_p;
            pso.best_stats = array2table( ...
                Inf(height(pso.best_p), 2), ...
                'VariableNames', {'J','CV'});
            pso.best_i = 1;
        end
        
        function [pso, dumps] = update(pso, dumps)
            nparam = width(pso.best_p);
            comb = height(pso.best_p);
            last_dump = dumps(end-(comb-1):end);
            p = struct2table(last_dump);
            p = p(:, pso.best_p.Properties.VariableNames);
            p{:,:} = log(p{:,:});
            stats = table( ...
                [last_dump.sugeno_RMSE]', ...
                [last_dump.RMS_NC]' ./ [last_dump.maxNC]', ...
                'VariableNames', {'J', 'CV'} ...
            );
            clear last_dump;
            for i=1:comb
                if  have_improved(pso.best_stats(i,:), stats(i,:))
                    pso.best_stats{i,:} = stats{i,:};
                    pso.best_p{i,:}     = p{i,:};
                    if have_improved(pso.best_stats(pso.best_i,:), stats(i,:))
                        pso.best_i = i;
                    end
                end
            end
            for i=1:comb
                pso.h{i, :} = pso.w * pso.h{i, :} ...
                    + pso.psi1 * rand(1, nparam) ...
                        .* (pso.best_p{ i , : } - p{i,:}) ...
                    + pso.psi2 * rand(1, nparam) ...
                        .* (pso.best_p{pso.best_i,:} - p{i,:});
                p{i, :} = p{i, :} + pso.h{i, :};
            end
            p = hola(p, fieldnames(dumps));
            unchanged = all(pso.h{:, :}==0, 2);
            if any(unchanged)
                p(unchanged).update_needed = false;
            end
            dumps = [dumps, p];
        end

    end
end

function bool = have_improved(best_pi, pi)
    if best_pi.CV >= 1 || pi.CV >= 1
        bool = pi.CV < best_pi.CV;
    else
        bool = pi.J < best_pi.J;
    end
end

function pi = hola(pi, fields)
    pi{:,:} = exp(pi{:,:});
    pi.normalize    =  true(height(pi), 1);
    pi.delta    = max(0, min(1, pi.delta));
    pi.tau      = max(0, min(1, pi.tau));
    pi.tmax     = max(1, pi.tmax);
    pi.maxNC    = max(1, pi.maxNC);
    pi.spmin    = pi.tmax ./ pi.maxNC;
    pi.cputime      = zeros(height(pi), 1);
    pi.progress     =  ones(height(pi), 1);
    pi.update_needed = true(height(pi), 1);
    pi.stats        =  cell(height(pi), 1);
    pi.mean_NC      =   NaN(height(pi), 1);
    pi.RMS_NC       =   Inf(height(pi), 1);
    pi.mamdani_RMSE =   Inf(height(pi), 1);
    pi.sugeno_RMSE  =   Inf(height(pi), 1);
    newfieldnames = setdiff(fields, pi.Properties.VariableNames);
    for field_i=1:length(newfieldnames)
        pi.(newfieldnames{field_i}) = cell(height(pi), 1);
    end
    pi = table2struct(pi(:, fields))';
end
