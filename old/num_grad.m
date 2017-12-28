function gradf = num_grad(f, x0)
    gradf = zeros(length(x0), 1);
    h = 1e-5;
    
    for i=1:length(x0)
        xp1 = x0;
        xp1(i) = xp1(i) + h;
        xm1 = x0;
        xm1(i) = xm1(i) - h;
        
        gradf(i) = (f(xp1) - f(xm1))/(2*h);
    end
end
