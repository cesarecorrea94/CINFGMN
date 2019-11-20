
function dumps = INFGMN_series_update(dumpname, batchsize)
    load(dumpname); % DS, dumps
    combinations = length(dumps);
    nbatchs = floor((length(DS)-1)/batchsize);
    count_updates = 0;
    for ii = 1:combinations
        if needs_update(dumps, ii, nbatchs)
            count_updates = count_updates + 1;
        end
    end
    fprintf('%i combinations to be updated\n', count_updates);
    fprintf(' step/comb.:  step time  | estimated remain time\n');
    lastsave = cputime;
    count_i = 0;
    sumsteptime = 0;
    for ii = 1:combinations
        if ~(needs_update(dumps, ii, nbatchs))
            continue
        end
        if cputime - lastsave > 60
            mysave(dumpname, DS, dumps);
            lastsave = cputime;
        end
        dumps(ii) = INFGMN_series_step(DS, batchsize, ...
            dumps(ii).maxStableNC, dumps(ii));
        count_i = count_i + 1;
        sumsteptime = sumsteptime + dumps(ii).cputime;
        fprintf('% 5i/% 5i: % 7.2f sec | % 8.2f min\n', ...
            count_i, count_updates, dumps(ii).cputime, ...
            sumsteptime / count_i * (count_updates - count_i) / 60);
    end
    fprintf('Elapsed time: %.2f min\n', sumsteptime/60);
    mysave(dumpname, DS, dumps);
end

function mysave(dumpname, DS, dumps)
    % Para evitar corrompimento em caso de interrupção da execução,
    % somente após o processo de save estar completo,
    save([dumpname '~'], 'DS', 'dumps');
    % que então, o arquivo temporário recebe o nome final.
    movefile([dumpname '~'], dumpname);
    % Enquanto isso, saves anteriores permanecem seguros.
end

function boolean = needs_update(dumps, ii, nbatchs)
    boolean = length(dumps(ii).stats) < nbatchs ...
            && dumps(ii).progress > 0.99 ...
            && dumps(ii).mean_NC   < dumps(ii).maxStableNC ...
            && dumps(ii).mean_NC_1 < dumps(ii).maxStableNC ...
            && dumps(ii).mean_NC_2 < dumps(ii).maxStableNC ...
            && dumps(ii).mean_NC_3 < dumps(ii).maxStableNC ...
            && dumps(ii).mean_NC   > 5 ...
...            && dumps(ii).mean_NC   < 30 ...25 ...
            && dumps(ii).mean_NC_1 > 5 ...10 ...
...            && dumps(ii).mean_NC_1 < 30 ...25 ...
            && dumps(ii).mean_NC_2 > 5 ...
...            && dumps(ii).mean_NC_2 < 25 ...20 ...
            && dumps(ii).mean_NC_3 > 10 ...
...            && dumps(ii).mean_NC_3 < 40 ...25 ...
            ;
end
