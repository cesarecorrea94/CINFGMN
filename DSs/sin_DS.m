
function DS = sin_DS(permute)
    n = 360;
    if permute
        x = 2*pi*randperm(n)/n;
    else, x = 2*pi.*(1:n)./n;
    end
    DS = mat2dataset([x' sin(x)'], 'VarNames', {'x', 'sin(x)'});
end
