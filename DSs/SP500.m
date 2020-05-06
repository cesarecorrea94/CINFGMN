function DS = SP500()
    X = dataset('File','^GSPC.csv','ReadVarNames',true,'ReadObsNames',false,'Delimiter',',');
    DS = zeros(length(X)-5, 7);
    DS(:,2) = X(1:end-5, 'Close');
    DS(:,3) = X(2:end-4, 'Close');
    DS(:,4) = X(3:end-3, 'Close');
    DS(:,5) = X(4:end-2, 'Close');
    DS(:,6) = X(5:end-1, 'Close');
    t       = X(5:end-1, 'Date');
    DS(:,7) = X(6:end-0, 'Close');
    DS = mat2dataset(DS, 'VarNames', [{'t'} ...
                {'y(t-4)'} {'y(t-3)'} {'y(t-2)'} {'y(t-1)'} {'y(t)'} {'y(t+1)'}]);
    DS(:,'t') = t;
    DS(:,'t') = [];
end