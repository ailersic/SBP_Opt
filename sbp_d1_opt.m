%% Define symbolic variables
clear all
clc

n = 6;
p = n - 1;

d = sym('d', [floor((n-2)/2), 1]);
h = sym('h', [ceil(n/2), 1]);
q = sym('q', [ceil(n/2)*floor(n/2), 1]);

for i=2:floor((n-2)/2)
    assume(d(i-1) <= d(i));
end
assume(d(1) >= 0)
assume(d(end) <= 0.5);

disp([num2str(p*n), ' equations and ', num2str(floor((n-2)/2) + ceil(n/2) + ceil(n/2)*floor(n/2)), ' unknowns'])

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

% x is on domain [0,1]
if mod(n, 2) == 1
    x = [0; d; 0.5; 1.-flipud(d); 1];
else
    x = [0; d; 1.-flipud(d); 1];
end

%% Solve accuracy equations

eqns = sym(zeros(p*length(x), 1));

for j=1:p
    eqns((j-1)*length(x)+1:j*length(x)) = (Q*x.^j == j*H*x.^(j-1));
end

sol = vpasolve(eqns);

H = subs(H, sol);
Q = subs(Q, sol);
x = subs(x, sol);
