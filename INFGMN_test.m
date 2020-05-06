
function dumps = INFGMN_test(DS, maxtime, nepochs, fistypes, spreads, ...
        normalize, deltas, taudists, tmaxs, maxStableNC, rangelen)
    %% constants
    combinations = length(deltas) * length(taudists) ...
        * length(tmaxs) * length(maxStableNC);
    equalsteptime = maxtime / combinations;
    %% alocate dumps struct
    dumps = struct( ...
            'normalize', NaN, ...
            'delta', NaN, ...
            'tau', NaN, ...
            'tmax', NaN, ...
            'spmin', NaN, ...
            'maxStableNC', NaN, ...
            'cputime', cell(1, combinations), ...
            'stats', NaN, ...
            'mamdani_RMSE', Inf, ...
            'mamdani_NC', Inf, ...
            'sugeno_RMSE', Inf, ...
            'sugeno_NC', Inf );
    %% begin test
    fprintf(' step/comb.:  step time  | estimated remain time\n');
    dumps_i = 0;
    sumsteptime = 0;
    estimatedmeansteptime = equalsteptime;
    for delta_i = deltas
    for tau_d = taudists
    tau_i = taufordist(tau_d', delta_i, rangelen);
    for tmax_i = tmaxs
    for maxNC = maxStableNC
        dumps_i = dumps_i + 1;
        params = dumps(dumps_i);
        params.normalize = normalize;
        params.delta = delta_i;
        params.tau = tau_i;
        params.tmax = tmax_i;
        params.spmin = tmax_i / maxNC;
        params.maxStableNC = maxNC;
        maxsteptime = 2 * equalsteptime - estimatedmeansteptime;
        if maxsteptime < 1
            throw(MException(['MATLAB:' self.name ':infeasibleCombinations'], ...
                'Infeasible combinations number for the "maxtime"'));
        end
%         maxsteptime = sqrt(equalsteptime * maxsteptime);
        params = INFGMN_step(DS, params, nepochs, fistypes, spreads, maxsteptime, maxNC);
        dumps(dumps_i) = params;
        sumsteptime = sumsteptime + params.cputime;
        estimatedmeansteptime = (sumsteptime ...
            + equalsteptime * (combinations - dumps_i)) / combinations;
        fprintf('% 5i/% 5i: % 7.2f sec | % 8.2f min\n', ...
            dumps_i, combinations, params.cputime, ...
            estimatedmeansteptime * (combinations - dumps_i) / 60);
    end
    end
    end
    end
    fprintf('Elapsed time: %.2f min\n', sumsteptime/60);
end

function params = INFGMN_step(DS, params, nepochs, fistypes, spreads, maxsteptime, maxNC)
    %% constants
    timeref = cputime;
    %% begin step
    gmm = INFGMN(minmaxDS(DS), 'delta', params.delta, 'tau', params.tau, ...
        'tmax', params.tmax, 'spmin', params.spmin, 'normalize', params.normalize );

    params.stats = struct(...
        'epochs', cell(1, nepochs), ...
        'mamdani_spread', NaN, ...
        'mamdani_RMSE', Inf, ...
        'sugeno_spread', NaN, ...
        'sugeno_RMSE', Inf, ...
        'NC', NaN );
    aging = 0;
    for epoch = 1:nepochs
        gmm.train(DS);
        stats = params.stats(epoch);
        stats.epochs = epoch;
        stats.NC = gmm.modelSize();
        prop = stats.NC / params.maxStableNC;
        aging = aging + prop * log2(prop+1);
        for type = fistypes
            gmm.setFisType(type{1});
            for spread_i = spreads
                gmm.setSpread(spread_i);
                y = gmm.recall(DS(:,1));
                err = dataset2mat(DS(:, end)) - dataset2mat(y(:, end));
                RMSE = sqrt(mean(err.^2));
                if RMSE < stats.([type{1} '_RMSE'])
                    stats.([type{1} '_RMSE']) = RMSE;
                    stats.([type{1} '_spread']) = spread_i;
                end
            end

            if stats.([type{1} '_RMSE']) < params.([type{1} '_RMSE']) ...
                    && stats.NC <= maxNC
                params.([type{1} '_RMSE']) = stats.([type{1} '_RMSE']);
                params.([type{1} '_NC'])   = stats.NC;
                aging = 0;
            end
        end
        params.stats(epoch) = stats;
        params.cputime = cputime - timeref;
        consumedtime = params.cputime / maxsteptime;
        regValue = 2^(consumedtime^2 -1);
        if regValue * epoch > 2 * sqrt(nepochs) ...
                && regValue * aging > sqrt(nepochs)
            break;
        end
    end
    params.stats = params.stats(1:epoch);
end
