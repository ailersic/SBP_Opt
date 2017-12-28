clear all
clc

n = 5;                  % Number of nodes
p = ceil(n/2);          % Degree of accuracy equations
m = 4;                  % Number of each starting delta in global search
objfunc = 'specrad';    % Objective function, can be 'specrad', 'acceqn'

h = sym('h', [ceil(n/2), 1]);
q = sym('q', [ceil(n/2)*floor(n/2), 1]);
d = sym('d', [floor((n-2)/2), 1]);

%% Generate matrices and vectors

disp('Generating symbolic H, Q matrices')

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

disp(['Solving accuracy equations to degree ', num2str(p)])

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

hsol = hqsol(1:ceil(n/2));
qsol = hqsol(ceil(n/2)+1:end);

Hsol = subs(H, h, hsol);
Qsol = subs(Q, q, qsol);

D1 = inv(Hsol)*Qsol;

%% Minimize objective function

if strcmp(objfunc, 'specrad')
    A = inv(Hsol)*(Qsol + diag([1; zeros(n-1, 1)]));
    err = @(d_) max(abs(eig(double(subs(A, d, d_)))));
elseif strcmp(objfunc, 'acceqn')
    e = D1*x.^(p+1) - (p+1)*x.^p;
    err = @(d_) double(subs(e'*Hsol*e, d, d_));
else
    disp('ERROR: No objective function')
    return
end

options = optimoptions(@fmincon, 'Algorithm', 'interior-point', ... 
 'Display', 'final', 'OptimalityTolerance', 1e-8, ...
 'SpecifyObjectiveGradient', false, 'CheckGradients', false, ... 
 'stepTolerance', 1e-14', 'MaxFunctionEvaluations', 100000, ... 
 'MaxIterations', 500, 'ScaleProblem', 'obj-and-constr', ...
 'ConstraintTolerance', 0, 'Display', 'off');

dsolvec = [];
errvvec = [];

d0 = gen_d0(m, 1, length(d));

disp(['Starting global search with ', num2str(length(d0)), ...
      ' initial guesses']);

for i=1:size(d0, 2)
    [dsol, errv] = fmincon(err, d0(:,i), [], [], [], [], [], [], ...
                           @(d_) nonlin_con(Hsol, d, d_), options);
    dsolvec = [dsolvec, dsol];
    errvvec = [errvvec, errv];
end

[errv, i] = min(errvvec);
dsol = dsolvec(:,i);

disp(['Final offset(s): [', num2str(dsol'), '] with error: ', ...
      num2str(errv)])
