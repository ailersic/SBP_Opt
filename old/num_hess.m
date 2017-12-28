function hessf = num_hess(f, x0)
    hessf = zeros(length(x0), length(x0));
    h = 1e-5;
    
    for i=1:length(x0)
        xp1 = x0;
        xp1(i) = xp1(i) + h;
        xm1 = x0;
        xm1(i) = xm1(i) - h;
        hessf(:,i) = (num_grad(f, xp1) - num_grad(f, xm1))/(2*h);
    end
end
