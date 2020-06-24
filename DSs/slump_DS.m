function DS = slump_DS()
    rng(0);
    DS = dataset('File','slump_test.data.csv', ...
        'ReadVarNames',true,'ReadObsNames',false,'Delimiter',',');
    DS = DS(randperm(size(DS,1)), 2:end-2);
end
