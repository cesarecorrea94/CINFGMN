
function params = INFGMN_series_step(DS, batchsize, max_NC, params)
    gmm = INFGMN(minmaxDS(DS), 'normalize', params.normalize, ...
        'delta', params.delta, 'tau', params.tau, ...
        'tmax', params.tmax, 'spmin', params.spmin );
    stats = struct(...
        't', cell(1, floor((length(DS)-1)/batchsize)), ...
        'NC', NaN, 'mean_NC', NaN, ...
        'expected_output', NaN, ...
        'mamdani_output', NaN, ...
        'mamdani_err', NaN, ...
        'sugeno_output', NaN, ...
        'sugeno_err', NaN, ...
        'cputime', NaN );
    for tt = batchsize:batchsize:length(DS)-1
        timeref = cputime;
        sumNC = 0;
        for ss = 1+tt-batchsize:tt
            gmm.train(DS(ss, :));
            sumNC = sumNC + gmm.modelSize();
        end
        ii = tt/batchsize;
        stats(ii).t = tt;
        stats(ii).NC = gmm.modelSize();
        stats(ii).mean_NC = sumNC / batchsize;
        stats(ii).expected_output = dataset2mat(DS(tt+1, end));
        for type = {'mamdani', 'sugeno'}
            gmm.setFisType(type{1});
            output = gmm.recall(DS(tt+1, 1:end-1));
            output = dataset2mat(output);
            stats(ii).([type{1} '_output']) = output;
            stats(ii).([type{1} '_err']) = output - stats(ii).expected_output;
        end
        stats(ii).cputime = cputime - timeref;
        if mod(tt, 1000) > params.tmax * sqrt(2) && ... se passou do ruído
                stats(ii).NC > max_NC && ... e meu NC ainda está alto
                stats(ii).mean_NC > max_NC && ... e a média do último batch foi alta
                mean([stats(1:ii).mean_NC]) > max_NC % e a média geral foi alta
            stats = stats(1:ii);
            break
        end
    end
    params.cputime = sum( [stats.cputime] );
    params.progress = tt / (length(DS)-1);
    params.stats = stats;
    %% mean NC
    params.mean_NC  = mean( [stats.mean_NC] );
    t1000 = [stats.t] <= 1000;
    t3000 = [stats.t] > 2000;
    t2000 = ~(t1000 | t3000);
    params.mean_NC_1  = mean( [stats(t1000).mean_NC] );
    params.mean_NC_2  = mean( [stats(t2000).mean_NC] );
    params.mean_NC_3  = mean( [stats(t3000).mean_NC] );
    %% RMSE
    params.mamdani_RMSE = sqrt(mean( [stats.mamdani_err] .^2 ));
    params.sugeno_RMSE  = sqrt(mean( [stats.sugeno_err] .^2 ));
end
