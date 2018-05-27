clear all
clc

n = 5;              % Number of nodes
m = 6;              % Number of each starting delta in global search
obj = 'desacceqn2';  % Objective function:
                    %  'acceqn':     minimize H norm truncation error for
                    %                next order accuracy equations
                    %  'acceqn2':    minimize H norm truncation error for
                    %                next 2 orders accuracy equations
                    %  'specrad':    minimize spectral radius for matrix
                    %                A = H^-1(Q + diag([1, zeros(n-1, 1)]))
                    %  'desacceqn':  minimize spectral radius for a desired
                    %                value of truncation error
                    %  'desacceqn2': minimize spectral radius for a desired
                    %                value of 2 term truncation error
                    %  'desspecrad': minimize truncation error for a
                    %                desired value of spectral radius
                    %  'desspecrad2':minimize 2 term truncation error for a
                    %                desired value of spectral radius
                    %  'weight':     minimize a weighted sum of spectral
                    %                radius and truncation error
                    %  'weight2':    minimize a weighted sum of spectral
                    %                radius and 2 term truncation error

wgtlow = 0.0;       % Lower limit for spread of weights
wgthigh = 1.0;      % Upper limit for spread of weights
wgtincr = 0.5;     % Increment in spread of weights

desvarlow = 0.6503;    % Lower limit for spread of desired variable value
desvarhigh = 1.8654*20;   % Upper limit for spread of desired variable value
desvarincr = (1.8654*20-0.6503)/100;  % Increment in spread of desired variable value

srwgt = wgtlow:wgtincr:wgthigh;
aewgt = 1 - srwgt;
desvar = desvarlow:desvarincr:desvarhigh;

warning('off', 'symbolic:mldivide:RankDeficientSystem');

% Generate H, Q, and D1 matrices as explicit functions of d offsets
[H, Q, D1] = gen_d1(n);

disp('Starting optimization loop')

% Choose what to iterate over based on objective function
if strcmp(obj, 'desacceqn') || strcmp(obj, 'desspecrad') || ...
        strcmp(obj, 'desacceqn2') || strcmp(obj, 'desspecrad2')
    k = length(desvar);
    srwgt = zeros(1, k);
    aewgt = zeros(1, k);
elseif strcmp(obj, 'weight') ||  strcmp(obj, 'weight2')
    k = length(srwgt);
    desvar = zeros(1, k);
else
    k = 1;
end

dsols = zeros(floor((n-2)/2), k);
srsols = zeros(1, k);
aesols = zeros(1, k);
D1sols = zeros(n, n, k);

% Iterate over weight or desired variable, store results in vectors
% I would use the parallel toolbox, but it crashes for whatever reason
%parfor i = 1:k
%    thread = getCurrentTask();
%    id = thread.ID;
for i = 1:k
    thread = 0;
    id = 0;
    
    disp(['Process ', num2str(id), ' started iteration i/k = ', ...
          num2str(i), '/', num2str(k)])
    [dsols(:, i), srsols(i), aesols(i), D1sols(:, :, i)] = ...
        opt_d1(H, Q, D1, m, obj, srwgt(i), aewgt(i), desvar(i), id);
    disp(['Process ', num2str(id), ' finished iteration i/k = ', ...
          num2str(i), '/', num2str(k)])
end

% Pareto plot of truncation error over spectral radius
figure
scatter(srsols, aesols, 'o')
xlabel('Spectral Radius')
ylabel('Truncation Error')
title(['Truncation Error over Spectral Radius for n = ', num2str(n)])