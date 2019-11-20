
function range = minmaxDS(DS)
    DSmat = dataset2mat(DS);
    range = array2table(minmax(DSmat')', ...
        'VariableNames', DS.Properties.VarNames);
    range = table2dataset(range);
end
