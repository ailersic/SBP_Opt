clear all
clc

n = 8;              % Number of nodes
m = 6;              % Number of each starting delta in global search
obj = 'weight';     % Objective function: 'specrad', 'acceqn', or 'weight'
wgtincr = 0.05;     % Increment in spread of weights

srwgt = 0:wgtincr:0.4;
aewgt = 1 - srwgt;
dsols = zeros(floor((n-2)/2), length(srwgt));
srsols = zeros(1, length(srwgt));
aesols = zeros(1, length(srwgt));
D1sols = zeros(n, n, length(srwgt));

[H, Q, D1] = gen_d1(n);

disp('Starting optimization loop')

% srwgt and aewgt are weights for specrad error and acceqn error
parfor i = 1:length(srwgt)
    disp(['Using specrad and acceqn weights ', num2str(srwgt(i)), ...
          ', ', num2str(aewgt(i))])
    [dsols(:, i), srsols(i), aesols(i), D1sols(:, :, i)] = ...
        opt_d1(H, Q, D1, m, obj, srwgt(i), aewgt(i));
    disp(['Optimization loop: ', num2str(100*i/length(srwgt)), '%'])
end

scatter(srsols, aesols, 'o')
xlabel('Spectral Radius')
ylabel('Truncation Error')
title(['Truncation Error over Spectral Radius for n = ', num2str(n)])