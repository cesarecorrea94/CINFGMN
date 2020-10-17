
function myprint(self, topN, folder)
    for ii=1:topN
        subfolder = [folder '/' num2str(ii)];
        mkdir(subfolder);
        myprintstep(subfolder, self.dumps(ii));
    end
end

function myprintstep(subfolder, dump)
        fig = figure('Position', [0,0, 4*50*4, 50*4]);
        position = [0.1300, 0.1100+0.060, 0.7750, 0.8150-0.030];
        %% output
        range = length(dump.stats)*2*200/3000;
        for center = length(dump.stats)*[1000, 2000]/3000
            ax1 = newaxis(fig);
            ax1.Position = position;
            myticks = -1.5:0.5:2.5;
            yticks(myticks);
            ylim(ax1, minmax(myticks));
            x = center-range/2:center+range/2;
            expected_output = [dump.stats(x).expected_output];
%             mamdani= expected_output + [dump.stats(x).mamdani_err];
            sugeno = expected_output + [dump.stats(x). sugeno_err];
            plot(x, expected_output,  '-', 'Color', '#777', 'DisplayName', 'saída real', 'LineWidth', 1);
%             plot(x, mamdani, 'b-.', 'DisplayName', 'saída ML', 'LineWidth', 1.5);
%             plot(x,  sugeno, 'r--', 'DisplayName', 'saída TS', 'LineWidth', 1.5);
            plot(x,  sugeno, 'k-',  'DisplayName', 'saída TS', 'LineWidth', 1.5);
            hid = plot([center center], [-1.5 2.5],  'r-', 'LineWidth', 1);
            set(get(get(hid,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            if center == 1000
                legend('Location', 'southeast');
            else,legend('Location', 'southwest');
            end
            plotlabels(ax1, 'Instância de Tempo t', 'y(t+1)');
            saveas(fig, [subfolder '/output_t' num2str(center)], 'png');
        end
        %% Number of rules
        ax1 = newaxis(fig);
        ax1.Position = position;
        xticks(0:200:3000);
        NCs = [dump.stats.NCs];
        plot(NCs, 'k-', 'LineWidth', 1);
        plot([1000 1000 nan 2000 2000], [0 ceil(max(NCs)/5)*5 nan 0 ceil(max(NCs)/5)*5],  'r-', 'LineWidth', 1);
        plotlabels(ax1, 'Instância de Tempo t', 'Número de regras');
        saveas(fig, [subfolder '/NCs'], 'png');
        %% RMSE
        ax1 = newaxis(fig);
        ax1.Position = position;
        xticks(0:5:30);
        myticks = 0:0.03:0.21;
        yticks(myticks);
        ylim(ax1, minmax(myticks));
        err = zeros(length(dump.stats)/30,30);
%         fistype = {'mamdani', 'sugeno'};
%         fiscolor = {'b:', 'r--'};
%         for fis_i = 1:2
            err(:) = [dump.stats.sugeno_err];
%             err(:) = [dump.stats.([fistype{fis_i} '_err'])];
            err(:) = [0 err(1:end-1)];
            hola = movavg(err);
            plot((0:30), [NaN hola], 'k-', 'DisplayName', 'Erro de aprendizado online', 'LineWidth', 1);
%             plot((0:30), [NaN hola], fiscolor{fis_i}, 'DisplayName', [fistype{fis_i} ' RMSE'], 'LineWidth', 1.5);
%         end
        hid = plot([10 10 nan 20 20], [minmax(myticks) nan minmax(myticks)],  'r-', 'LineWidth', 1);
        set(get(get(hid,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend('Location', 'northwest');
        plotlabels(ax1, 'Instância de Tempo t x 100', 'RMSE');
        saveas(fig, [subfolder '/RMSE'], 'png');
        close(fig);
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
