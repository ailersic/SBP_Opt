function D = acc_eqn_free_d(n, d)
    p = n - 1 - floor((n-2)/2);

    h = sym('h', [ceil(n/2), 1]);
    q = sym('q', [ceil(n/2)*floor(n/2), 1]);

    if size(d) ~= [floor((n-2)/2), 1]
        disp('Wrong shape of d vector!')
        D = -1;
        return
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

    %% Solve accuracy equations

    eqns = sym(zeros(p*length(x), 1));

    for j=1:p
        eqns((j-1)*length(x)+1:j*length(x)) = Q*x.^j - j*H*x.^(j-1);
    end

    sol = solve(eqns);
    D = subs(inv(H)*Q, sol);
end