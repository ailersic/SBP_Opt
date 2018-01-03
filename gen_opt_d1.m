clear all
clc

n = 5;              % Number of nodes
m = 6;              % Number of each starting delta in global search
obj = 'desacceqn';  % Objective function:
                    %  'acceqn':     minimize H norm truncation error for
                    %                next order accuracy equations
                    %  'specrad':    minimize spectral radius for matrix
                    %                A = H^-1(Q + diag([1, zeros(n-1, 1)]))
                    %  'desacceqn':  minimize spectral radius for a desired
                    %                value of truncation error
                    %  'desspecrad': minimize truncation error for a
                    %                desired value of spectral radius
                    %  'weight':     minimize a weighted sum of spectral
                    %                radius and truncation error

wgtlow = 0.0;       % Lower limit for spread of weights
wgthigh = 1.0;      % Upper limit for spread of weights
wgtincr = 0.02;     % Increment in spread of weights

desvarlow = 0.0;    % Lower limit for spread of desired variable value
desvarhigh = 1.1;   % Upper limit for spread of desired variable value
desvarincr = 0.1;   % Increment in spread of desired variable value

srwgt = wgtlow:wgtincr:wgthigh;
aewgt = 1 - srwgt;
desvar = desvarlow:desvarincr:desvarhigh;

dsols = zeros(floor((n-2)/2), length(srwgt));
srsols = zeros(1, length(srwgt));
aesols = zeros(1, length(srwgt));
D1sols = zeros(n, n, length(srwgt));

warning('off', 'symbolic:mldivide:RankDeficientSystem');

% Generate H, Q, and D1 matrices as explicit functions of d offsets
[H, Q, D1] = gen_d1(n);

disp('Starting optimization loop')

% Choose what to iterate over based on objective function
if strcmp(obj, 'desacceqn') || strcmp(obj, 'desspecrad')
    k = length(desvar);
    srwgt = zeros(1, k);
    aewgt = zeros(1, k);
elseif strcmp(obj, 'weight')
    k = length(srwgt);
    desvar = zeros(1, k);
else
    k = 1;
end

% Iterate over weight or desired variable, store results in vectors
parfor i = 1:k
    [dsols(:, i), srsols(i), aesols(i), D1sols(:, :, i)] = ...
        opt_d1(H, Q, D1, m, obj, srwgt(i), aewgt(i), desvar(i));
    disp(['Optimization loop: ', num2str(100*i/k), '%'])
end

% Pareto plot of truncation error over spectral radius
scatter(srsols, aesols, 'o')
xlabel('Spectral Radius')
ylabel('Truncation Error')
title(['Truncation Error over Spectral Radius for n = ', num2str(n)])