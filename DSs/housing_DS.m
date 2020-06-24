function DS = housing_DS()
    rng(0);
    DS = dataset('File','housing.data.csv','ReadVarNames',false,'ReadObsNames',false,'Delimiter',',');
    DS.Properties.VarNames = {'CRIM','ZN','INDUS','CHAS','NOX', ...
        'RM','AGE','DIS','RAD','TAX','PTRATIO','B','LSTAT','MEDV'};
    DS = DS(randperm(size(DS,1)), [1:3 5:end]);
end
