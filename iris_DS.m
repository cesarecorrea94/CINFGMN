
function DS = iris_DS()
    iris = load('iris_dataset.mat');
    DS = mat2dataset([iris.irisInputs' iris.irisTargets'], 'VarNames', ...
        {'sepal length','sepal width','petal length','petal width', ...
        'setosa','virginica','versicolor'});
end
