function [dsol, srsol, aesol, D1sol] = sbp_d1_opt(n, m, obj, srwgt, aewgt)
%% Generate an optimized SBP first derivative operator
    % n = number of nodes
    % m = number of each starting offset in global search
    % obj = objective function, can be 'specrad', 'acceqn', or 'weight'
    % srwgt = coefficient for specrad function if obj is 'weight'
    % aewgt = coefficient for acceqn function if obj is 'weight'

    disp(['Finding SBP operator for n = ', num2str(n), ' and obj = ', obj])
    if strcmp(obj, 'weight')
        disp(['With specrad and acceqn weights ', num2str(srwgt), ', ', ...
              num2str(aewgt)]);
    end

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

    %% Minimize objective function

    A = inv(Hsol)*(Qsol + diag([1; zeros(n-1, 1)]));
    specrad = @(d_) max(abs(eig(double(subs(A, d, d_)))));

    e = D1*x.^(p+1) - (p+1)*x.^p;
    acceqn = @(d_) double(subs(e'*abs(Hsol)*e, d, d_));
    % I know taking the absolute value of Hsol should be unnecessary, but
    % floating point error is a bastard. I enforce that each element of
    % Hsol is > 0 in nonlin_con anyway.

    if strcmp(obj, 'specrad')
        err = @(d_) specrad(d_);
    elseif strcmp(obj, 'acceqn')
        err = @(d_) acceqn(d_);
    elseif strcmp(obj, 'weight')
        err = @(d_) srwgt*specrad(d_) + aewgt*acceqn(d_);
    else
        disp('ERROR: No objective function')
        return
    end

    % Define parameters for optimization
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', ... 
     'Display', 'final', 'OptimalityTolerance', 1e-8, ...
     'SpecifyObjectiveGradient', false, 'CheckGradients', false, ... 
     'stepTolerance', 1e-14', 'MaxFunctionEvaluations', 100000, ... 
     'MaxIterations', 500, 'ScaleProblem', 'obj-and-constr', ...
     'ConstraintTolerance', 0, 'Display', 'off');

    % Try a spread of starting offsets, store the converged offsets and
    % their minimized objective function values in vectors
    dsolvec = [];
    errvvec = [];
    d0 = gen_d0(m, 1, length(d));

    disp(['Starting global search with ', num2str(length(d0)), ...
          ' initial guesses']);

    % Try each set of offsets
    for i=1:size(d0, 2)
        [dsol, errv] = fmincon(err, d0(:,i), [], [], [], [], [], [], ...
                               @(d_) nonlin_con(Hsol, d, d_), options);
        dsolvec = [dsolvec, dsol];
        errvvec = [errvvec, errv];
        disp([num2str(100*i/size(d0, 2)), '%'])
    end

    % Choose the best solution and present it
    [errv, i] = min(errvvec);
    dsol = dsolvec(:,i);
    srsol = specrad(dsol);
    aesol = acceqn(dsol);
    D1sol = double(subs(D1, d, dsol));

    disp(['Final offset(s): [', num2str(dsol'), ...
          '] with spectral radius: ', num2str(srsol), ...
          ', truncation error: ', num2str(aesol), ...
          ', and weighted error: ', num2str(errv)])
end