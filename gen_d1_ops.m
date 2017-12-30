clear all
clc

n = 6;              % Number of nodes
m = 6;              % Number of each starting delta in global search
obj = 'weight';     % Objective function: 'specrad', 'acceqn', or 'weight'
wgtincr = 0.01;     % Increment in spread of weights

srwgt = 0:wgtincr:1;
aewgt = 1 - srwgt;
dsols = zeros(floor((n-2)/2), length(srwgt));
srsols = zeros(1, length(srwgt));
aesols = zeros(1, length(srwgt));
D1sols = zeros(n, n, length(srwgt));

% srwgt and aewgt are weights for specrad error and acceqn error
for i = 1:length(srwgt)
    [dsols(:, i), srsols(i), aesols(i), D1sols(:, :, i)] = ...
        sbp_d1_opt(n, m, obj, srwgt(i), aewgt(i));
    disp(['Weight loop: ', num2str(100*i/length(srwgt)), '%'])
end

scatter(srsols, aesols, 'o')
xlabel('Spectral Radius')
ylabel('Truncation Error')
title(['Truncation Error over Spectral Radius for n = ', num2str(n)])