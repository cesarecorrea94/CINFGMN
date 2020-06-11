
function myplotfis(fis, ncols, rangesDS, filename)
%     fis.Inputs = [fis.Inputs fis.Outputs];
    nrows = ceil(length(fis.Inputs)/ncols);
    numPoints = 200;
    fig = figure('visible','off', ...
        'Position', 2*numPoints*[0,0,ncols,nrows/4]);
    tiledlayout(nrows, ncols);%('flow');
    for inp = 1:length(fis.Inputs)
        ranges = double(rangesDS(:,fis.Inputs(inp).Name))';
        [~, proportion] = mapminmax(ranges, -1, 1);
        ax = nexttile;
        hold(ax, 'on');
        xlim(ax, round(ranges, 2));
        xticks(round(mapminmax('reverse', -1:0.5:1, proportion), 2));
        [xOut,yOut] = plotmf(fis, 'input', inp, numPoints);
        for mf = 1:size(xOut,2)
            plot(mapminmax('reverse', xOut(:,mf), proportion), yOut(:,mf));
        end
        title(ax, replace(fis.Inputs(inp).Name, '_', ' '));
    end
    saveas(fig, filename, 'png');
    close(fig);
end

function oldmyplotfis(prefixname, fis)
    fig = figure;
    fisvar = {fis.Inputs, fis.Outputs};
    vartype = {'input', 'output'};
    for type=1:1
        for inp = 1:length(fisvar{type})
            ax = newaxis(fig);
            [xOut,yOut] = plotmf(fis, vartype{type}, inp);
            for mf=1:size(xOut,2)
                plot(xOut(:,mf), yOut(:,mf));%,  'k-', 'DisplayName', 'mf1');
            end
            saveas(fig, [prefixname ' (' convertStringsToChars(fisvar{type}(inp).Name) ')'], 'png');
        end
    end
    close(fig);
end
