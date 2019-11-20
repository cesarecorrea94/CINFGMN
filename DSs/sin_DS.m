
function DS = sin_DS()
    x = 2*pi*randperm(180)/180;
    DS = mat2dataset([x' sin(x)'], 'VarNames', {'x', 'sin(x)'});
end
