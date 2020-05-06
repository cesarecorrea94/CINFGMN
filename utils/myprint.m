
function myprint(prefixname, dumps)
    fig = figure;
    for ii=1:length(dumps)
        myprintstep(fig, [prefixname '(' num2str(ii) ')'], dumps(ii));
    end
    close(fig);
end

function myprintstep(fig, prefixname, dump)
        %% Number of rules
    %     tiledlayout(2,1)
    %     ax1 = nexttile;
        ax1 = newaxis(fig);
    %     hold(ax1, 'on');
        plot([dump.stats.NC]);
        plotlabels(ax1, 'Time Instance t', 'Number of rules');
        saveas(fig, [prefixname '_NC'], 'png');
        %% RMSE
    %     ax1 = nexttile;
        ax1 = newaxis(fig);
    %     hold(ax1, 'on');
        err = zeros(length(dump.stats)/30,30);
        fistype = {'mamdani', 'sugeno'};
        fiscolor = {'b:', 'r--'};
        for fis_i = 1:2
            err(:) = [dump.stats.([fistype{fis_i} '_err'])];
            err(:) = [0 err(1:end-1)];
            plot((0:30).*100, [NaN movavg(err)], fiscolor{fis_i}, 'DisplayName', [fistype{fis_i} ' RMSE'], 'LineWidth', 1.2);
%             plot([0 100:30*100], [NaN zmovavg(err)], fiscolor{fis_i}, 'DisplayName', [fistype{fis_i} ' RMSEz'], 'LineWidth', 1.2);
        end
        legend;
        plotlabels(ax1, 'Time Instance t', 'RMSE');
        saveas(fig, [prefixname  '_RMSE'], 'png');
    %     saveas(fig, [prefixname '_NC_RMSE'], 'png');
        %% output
        range = length(dump.stats)/15;
    %     tiledlayout(2,1)
        for center = [length(dump.stats)/3, length(dump.stats)*2/3]
    %         ax1 = nexttile;
            ax1 = newaxis(fig);
    %         hold(ax1, 'on');
            x = center-range/2:center+range/2;
            stats = dump.stats(x);
            plot(x, [stats.expected_output],  'k-', 'DisplayName', 'actual');
            plot(x, [stats. mamdani_output], 'b-.', 'DisplayName', 'ML output', 'LineWidth', 1.2);
            plot(x, [stats.  sugeno_output], 'r--', 'DisplayName', 'TS output', 'LineWidth', 1.2);
            legend;
            plotlabels(ax1, 'Time Instance t', 'y(t+1)');
            saveas(fig, [prefixname '_output' num2str(center)], 'png');
        end
    %     saveas(fig, [prefixname '_outputs'], 'png');
end

function ret = movavg(err)
    ret = sqrt(mean(err.^2, 1));
end

function ret = zmovavg(err)
    window = 100;
    wgtsqrerr = (err(:).^2) ./window;
    ret = zeros(1, length(wgtsqrerr)-(window-1));
    for ii=1:length(wgtsqrerr)-(window-1)
        ret(ii) = sqrt(sum(wgtsqrerr(ii:ii+(window-1))));
    end
end
