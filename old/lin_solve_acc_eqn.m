function err = lin_solve_acc_eqn(n, d)
    p = 3;

    h = sym('h', [ceil(n/2), 1]);
    q = sym('q', [ceil(n/2)*floor(n/2), 1]);

    if size(d) ~= [floor((n-2)/2), 1]
        disp('Wrong shape of d vector!')
        D = -1;
        return
    end
    
    for i=1:ceil(n/2)
        assumeAlso(h(i) > 0);
    end
    
    %% Generate matrices and vectors

    H = sym(eye(n));

    for i=1:ceil(n/2)
        H(i, i) = h(i);
        H(n-i+1, n-i+1) = h(i);
    end

    Q = sym(zeros(n,n));
    Q(1, 1) = -1/2;
    Q(n, n) = 1/2;

    k = 1;
    for i=1:ceil(n/2)
        for j=i+1:n-i+1
            Q(i, j) = q(k);
            Q(n-j+1, n-i+1) = q(k);
            Q(j, i) = -q(k);
            Q(n-i+1, n-j+1) = -q(k);
            k = k + 1;
        end
    end
    
    % x is on domain [-1,1]
    if mod(n, 2) == 1
        x = [-1; d-ones(length(d),1); 0; ones(length(d),1)-flipud(d); 1];
    else
        x = [-1; d-ones(length(d),1);    ones(length(d),1)-flipud(d); 1];
    end

    At = sym(zeros(n*p, length(h) + length(q)));
    bt = sym(zeros(n*p, 1));

    for j=1:p
        [A, b] = equationsToMatrix(j.*H*(x.^(j-1)) == Q*(x.^j), [h; q]);
        At(n*(j-1)+1:n*j, :) = A;
        bt(n*(j-1)+1:n*j) = b;
    end
    
    [A, b] = equationsToMatrix(Q*ones(n, 1) == zeros(n, 1), [h; q]);
    At = [At; A];
    bt = [bt; b];
    
    hqsol = At\bt;

    hsol = hqsol(1:floor(n/2));
    qsol = hqsol(floor(n/2)+1:end);

    Hsol = double(subs(H, h, hsol));
    Qsol = double(subs(Q, q, qsol));

    [~,spd] = chol(Hsol);
    spd = (spd == 0);
    
    D1 = inv(Hsol)*Qsol;

    e = D1*x.^(p+1) - (p+1)*x.^p;
    err = double(e'*Hsol*e)
end