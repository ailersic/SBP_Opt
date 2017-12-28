n = 4;
m = 40;

di = linspace(0.2, 0.3, m+2);
di = di(2:end-1);
erf = zeros(m, 1);

err_func = @(d) h_norm_err(n, d);

for i=1:m
    erf(i) = err_func(di(i));
    disp([num2str(100*i/m), '%'])
end

plot(di, erf);
axis equal