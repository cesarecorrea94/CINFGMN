
function DS = sin_DS()
    x = 2*pi*(0:100)/100;
    DS = mat2dataset([x' sin(x)'], 'VarNames', {'x', 'sin(x)'});
end
