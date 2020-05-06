warning off;
dumpname = '/home/cesar/INFGMNthesis/dumps/DynNLSys_1x2288.mat';
%%
if ~exist('dumpseries', 'var')
    if exist('dumpname', 'var')
        dumpseries = INFGMN_series(dumpname);
    else
        deltas  = 2.^-(0.5:1.5:9.5);
        % deltas  = 10.^-((1 -0/2 :1/2: 2 +0/2));
        taus    = 2.^-(1.5:1.5:10.5); % quando criar novas componentes
        % taus    = 10.^-((1 -0/2 :1/2: 2 +0/2)); % quando criar novas componentes
        tmaxs   = 2.^(6:8); % tempo para assimilar instâncias
%         tmaxs   = [50 100 200 250]; % tempo para assimilar instâncias
        maxNCs  = 2.^(4:6); % número máximo de componentes estáveis (não espúrias)
%         maxNCs  = [25 30 40 55]; % número máximo de componentes estáveis (não espúrias)
        % uniform = false; % proporções uniformes
        normalize = true;%default
        % regValue = 0; % regularização da covariância

        dumpseries = INFGMN_series(NaN);
        dumpseries.create('DynNLSys', dynamic_nonlinear_system(), ...
            normalize, deltas, taus, tmaxs, maxNCs);
    end
end
% fuzzyLogicDesigner(gmm.toFIS());

for ii=1:100
    fprintf('Iteração: %i\n', ii);
    dumpseries.update(@check_break);
    dumpseries.pso_update();
end

%     if exist('dumps', 'var')==0
%         load(dumpname, 'DS', 'dumps');
%         [~,index] = sortrows([dumps.mean_NC].'); dumps = dumps(index(end:-1:1)); clear index
%         for ii = length(dumps)-1:-1:1
%             if isnan(dumps(ii+1).sugeno_RMSE)
%                 dumps(ii+1) = [];
%             elseif isnan(dumps(ii).sugeno_RMSE)
%                 dumps(ii) = [];
%             elseif dumps(ii).sugeno_RMSE > dumps(ii+1).sugeno_RMSE
%                 dumps(ii) = [];
%             elseif dumps(ii).mean_NC == dumps(ii+1).mean_NC
%                 dumps(ii+1) = [];
%             end
%         end
%     end
%     dump_var = 'dump_merge';% dont merge
%     dumpfile = ['test/' dump_var '.mat'];
%     if exist(dumpfile, 'file')~=0
%         load(dumpfile, 'dump_cur');
%     else
%         dump_cur = more_stats(INFGMN_series_step(DS, 1, dumps(1), @check_break));
%         save(dumpfile, 'dump_cur');
%     end
%     eval([dump_var '= dump_cur;']);
%     myplotfis(['test/' dump_var ' t1400'], dump_cur.stats(end*1400/3000).fis);
%     myplotfis(['test/' dump_var ' t2400'], dump_cur.stats(end*2400/3000).fis);
%     myprint(['test/' dump_var], dump_cur);

function bool = check_break(dumps_ii, stats)
    bool =  2 * dumps_ii.tmax < mod(stats(end).t, 1000) ... se passou do ruído
        &&  (   dumps_ii.maxNC < stats(end).NC          ... e:  meu NC ainda está alto
            ||  dumps_ii.maxNC < mean([stats.RMS_NC])   ...     ou a média geral foi alta
            );
end

function dumps_ii = more_stats(dumps_ii)
    t1000 = [dumps_ii.stats.t] <= 1000;
    t3000 = [dumps_ii.stats.t] > 2000;
    t2000 = ~(t1000 | t3000);
    dumps_ii.mean_NC_1  = mean( [dumps_ii.stats(t1000).mean_NC] );
    dumps_ii.mean_NC_2  = mean( [dumps_ii.stats(t2000).mean_NC] );
    dumps_ii.mean_NC_3  = mean( [dumps_ii.stats(t3000).mean_NC] );
end

function bool = needs_update(dumps_ii, nbatchs)
    bool =  length(dumps_ii.stats) < nbatchs ...
        &&  dumps_ii.progress > 0.99 ...
        &&  (   isnan(dumps_ii.mean_NC) ...
            ||  dumps_ii.mean_NC   < dumps_ii.maxNC ...
            &&  dumps_ii.mean_NC_1 < dumps_ii.maxNC ...
            &&  dumps_ii.mean_NC_2 < dumps_ii.maxNC ...
            &&  dumps_ii.mean_NC_3 < dumps_ii.maxNC ...
            &&  dumps_ii.mean_NC   > 5 ...
...            &&  dumps_ii.mean_NC   < 30 ...25 ...
            &&  dumps_ii.mean_NC_1 > 5 ...10 ...
...            &&  dumps_ii.mean_NC_1 < 30 ...25 ...
            &&  dumps_ii.mean_NC_2 > 5 ...
...            &&  dumps_ii.mean_NC_2 < 25 ...20 ...
            &&  dumps_ii.mean_NC_3 > 10 ...
...            &&  dumps_ii.mean_NC_3 < 40 ...25 ...
            );
end

function check_update_needed(dumpname, batchsize)
    load(dumpname, 'DS', 'dumps');
    nbatchs = floor((length(DS)-1)/batchsize);
    for ii = 1:length(dumps)
        dumps(ii) = more_stats(dumps(ii));
        dumps(ii).update_needed = needs_update(dumps(ii), nbatchs);
    end
    save(dumpname, 'DS', 'dumps');
end
