
function myplotfis(prefixname, fis)
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
