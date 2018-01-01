function [Hsol, Qsol, D1] = gen_d1(n)
%% Generate an SBP first derivative operator
    % n = number of nodes
    % m = number of each starting offset in global search
    % obj = objective function, can be 'specrad', 'acceqn', or 'weight'
    % srwgt = coefficient for specrad function if obj is 'weight'
    % aewgt = coefficient for acceqn function if obj is 'weight'

    disp(['Finding SBP operator for n = ', num2str(n)])

    h = sym('h', [ceil(n/2), 1]);
    q = sym('q', [ceil(n/2)*floor(n/2), 1]);
    d = sym('d', [floor((n-2)/2), 1]);

    %% Generate matrices and vectors

    disp('Generating symbolic H, Q matrices')

    % Define diagonal matrix H
    H = sym(eye(n));

    for i=1:ceil(n/2)
        H(i, i) = h(i);
        H(n-i+1, n-i+1) = h(i);
    end

    % Define bisymmetric matrix Q
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
    
    p = ceil(n/2); % Max degree of accuracy equations to be solved

    disp(['Solving accuracy equations to degree ', num2str(p)])

    At = sym(zeros(n*p, length(h) + length(q)));
    bt = sym(zeros(n*p, 1));

    % Define accuracy equations and store them in matrix At and vector bt
    for j=1:p
        [A, b] = equationsToMatrix(j.*H*(x.^(j-1)) == Q*(x.^j), [h; q]);
        At(n*(j-1)+1:n*j, :) = A;
        bt(n*(j-1)+1:n*j) = b;
    end
    [A, b] = equationsToMatrix(Q*ones(n, 1) == zeros(n, 1), [h; q]);
    At = [At; A];
    bt = [bt; b];

    % Solve linear system for solution vector [h; q]
    hqsol = At\bt;

    % Take h and q values and plug into respective matrices
    hsol = hqsol(1:ceil(n/2));
    qsol = hqsol(ceil(n/2)+1:end);

    Hsol = subs(H, h, hsol);
    Qsol = subs(Q, q, qsol);

    % Solve for derivative operator
    D1 = inv(Hsol)*Qsol;
end