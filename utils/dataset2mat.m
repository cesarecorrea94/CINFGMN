
function DSmat = dataset2mat(DS)
    DScell = dataset2cell(DS);
    DSmat = cell2mat(DScell(2:end,:));
end
