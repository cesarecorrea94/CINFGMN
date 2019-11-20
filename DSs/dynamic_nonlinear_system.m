function DS = dynamic_nonlinear_system()
    n = 3001;
	u = zeros(1,n);
	y = zeros(1,n+1);
	for t = 1:n
		u(t) = ufunc(t);
		y(t+1) = y(t)/(1+y(t)^2) + u(t)^3 + f(t);
	end
	DS = mat2dataset([u' y(1:end-1)' y(2:end)'], 'VarNames', {'u(t)', 'y(t)', 'y(t+1)'});
end

function res = ufunc(t)
	res = sin(2*pi*t/100);
end

function res = f(t)
    if (1001 <= t) && (t <= 2000)
        res = 1;
    else, res = 0;
    end
end