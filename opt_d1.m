function [dsol, srsol, aesol, D1sol] = opt_d1(H, Q, D1, m, obj, srwgt, aewgt, desvar, id)
%% Minimize objective function for given D1 as a function of deltas
    n = length(D1);
    p = ceil(n/2); % Max degree of accuracy equations to be solved
    
    d = sym('d', [floor((n-2)/2), 1]);
    
    % These warnings get really annoying
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');

    % x is on domain [-1,1]
    if mod(n, 2) == 1
        x = [-1; d-ones(length(d),1); 0; ones(length(d),1)-flipud(d); 1];
    else
        x = [-1; d-ones(length(d),1);    ones(length(d),1)-flipud(d); 1];
    end

    A = inv(H)*(Q + diag([1; zeros(n-1, 1)]));
    specrad = @(d_) max(abs(eig(double(subs(A, d, d_)))));

    e1 = D1*x.^(p+1) - (p+1)*x.^p;
    acceqn = @(d_) double(subs(e1'*abs(H)*e1, d, d_));
    
    e2 = D1*x.^(p+2) - (p+2)*x.^(p+1);
    e12 = e1 + e2;
    acceqn2 = @(d_) double(subs(e12'*abs(H)*e12, d, d_));
    % I know taking the absolute value of Hsol should be unnecessary, but
    % floating point error is a bastard. I enforce that each element of
    % Hsol is > 0 in nonlin_con anyway.
    
    if strcmp(obj, 'specrad')
        err = @(d_) specrad(d_);
        ceqf = [];
    elseif strcmp(obj, 'acceqn')
        err = @(d_) acceqn(d_);
        ceqf = [];
    elseif strcmp(obj, 'acceqn2')
        err = @(d_) acceqn2(d_);
        ceqf = [];
    elseif strcmp(obj, 'desspecrad')
        err = @(d_) acceqn(d_);
        ceqf = @(d_) specrad(d_);
    elseif strcmp(obj, 'desspecrad2')
        err = @(d_) acceqn2(d_);
        ceqf = @(d_) specrad(d_);
    elseif strcmp(obj, 'desacceqn')
        err = @(d_) specrad(d_);
        ceqf = @(d_) acceqn(d_);
    elseif strcmp(obj, 'desacceqn2')
        err = @(d_) specrad(d_);
        ceqf = @(d_) acceqn2(d_);
    elseif strcmp(obj, 'weight')
        err = @(d_) srwgt*specrad(d_) + aewgt*acceqn(d_);
        ceqf = [];
    elseif strcmp(obj, 'weight2')
        err = @(d_) srwgt*specrad(d_) + aewgt*acceqn2(d_);
        ceqf = [];
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
    d0 = init_deltas(m, 1, length(d));

    disp(['Starting global search with ', num2str(length(d0)), ...
          ' initial guesses']);

    % Try each set of offsets
    for i=1:size(d0, 2)
        [dsol, errv] = fmincon(err, d0(:,i), [], [], [], [], [], [], ...
                       @(d_) nonlin_con(H, d, d_, ceqf, desvar), options);
        dsolvec = [dsolvec, dsol];
        errvvec = [errvvec, errv];
        disp(['Global search in process ', num2str(id), ': ', ...
              num2str(100*i/size(d0, 2)), '%'])
    end

    % Choose the best solution and present it
    [errv, i] = min(errvvec);
    dsol = dsolvec(:,i);
    srsol = specrad(dsol);
    if strcmp(obj(end), '2')
        aesol = acceqn2(dsol);
    else
        aesol = acceqn(dsol);
    end
    D1sol = double(subs(D1, d, dsol));

    disp(['Final offset(s): [', num2str(dsol'), ...
          '] with spectral radius: ', num2str(srsol), ...
          ', truncation error: ', num2str(aesol), ...
          ', and weighted error: ', num2str(errv)])