%% Define symbolic variables
clear all
clc

n = 6;
p = 3;%n - ceil((n-2)/2);

h = sym('h', [floor(n/2), 1]);
q = sym('q', [ceil(n/2)*floor(n/2), 1]);

%% Generate matrices and vectors

H = sym(eye(n));

for i=1:floor(n/2)
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

d = sym('d', [floor((n-2)/2), 1]);

% x is on domain [-1,1]
if mod(n, 2) == 1
    x = [-1; d-ones(length(d),1); 0; ones(length(d),1)-flipud(d); 1];
else
    x = [-1; d-ones(length(d),1);    ones(length(d),1)-flipud(d); 1];
end

sec = ceil(n/2);

At = sym(zeros(sec*p, length(h) + length(q)));
bt = sym(zeros(sec*p, 1));
eqnsTot = [];

for j=1:p
    eqns = j.*H*(x.^(j-1)) == Q*(x.^j);
    eqnsTot = [eqnsTot; eqns];
    [A, b] = equationsToMatrix(eqns(1:sec), [h; q]);
    At(sec*(j-1)+1:sec*j, :) = A;
    bt(sec*(j-1)+1:sec*j) = b;
end

hqsol = At\bt;

hsol = hqsol(1:floor(n/2));
qsol = hqsol(floor(n/2)+1:end);

Hsol = subs(H, h, hsol);
Qsol = subs(Q, q, qsol);

D1 = inv(Hsol)*Qsol;

%%

e = D1*x.^(p+1) - (p+1)*x.^p;
err = @(di) double(subs(e'*Hsol*e, d, di));

%%

err = @(di) max(abs(eig(double(subs(D1, d, di)))));

%%

d0 = [0.3453;0.9980];

%%

d0 = [0.6;0.8];

%%

%options = optimoptions(@fmincon,'StepTolerance', 1e-12);
for i=1:1
    df = fmincon(err, d0, [1 -1; -1 0; 0 1], [0; -0.1; 0.9])%, [],[],[],[],[],options)
    err(df)
    d0 = df;
end

%%

np = 1000;

d1 = linspace(0.001, 0.999, np);
d2 = linspace(0.001, 0.999, np);
errm = zeros(np, np);
k = 0;

for i=1:np-1
    for j=i+1:np
        errm(i,j) = err([d1(i); d2(j)]);
        k = k + 1;
    end
    disp(2*k/(np^2 - ceil(np/2)))
end

surf(d1, d2, errm)
axis([0 1 0 1 0 3])