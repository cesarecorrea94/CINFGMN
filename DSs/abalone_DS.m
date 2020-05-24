function DS = abalone_DS()
    rng(0);
    DS = dataset('File','abalone.data.csv','ReadVarNames',false,'ReadObsNames',false,'Delimiter',',');
    DS.Properties.VarNames = {'Sex', ...
        'Length','Diameter','Height','Whole_weight', ...
        'Shucked_weight','Viscera_weight','Shell_weight','Rings'};
    DS = DS(randperm(size(DS,1)), 2:end);
end
