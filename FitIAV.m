clear;clc;

% load experimental data
data_m_nk = xlsread('data_to_fit.xlsx', 2, 'B13:E15');
data_m = data_m_nk([1 2], :);
data_nk = data_m_nk([1 3], :);

data_v = xlsread('data_to_fit.xlsx', 2, 'B17:I18');

data_t_te = xlsread('data_to_fit.xlsx', 2, 'B20:C22');
data_t = data_t_te([1 2], :);
data_te = data_t_te([1 3], :);

% calculate cell numbers with H = 10 ^ 7
data_m(2, :) = data_m(2, :) .* 10 ^ 5;
data_nk(2, :) = data_nk(2, :) .* 10 ^ 5;
data_t(2, :) = data_t(2, :) .* 10 ^ 5;
data_te(2, :) = data_te(2, :) .* 10 ^ 5;

data_v(2, :) = 10 .^ data_v(2, :);

FitFcn = @(parFit)ODE_fit(data_v, data_m, data_nk, data_t, data_te, parFit);

% initial parameters
par_bound = xlsread('par_IAV.xlsx');
lb = par_bound(:, 1);
ub = par_bound(:, 2);
par0 = par_bound(:, 3);
% translate initial par of odes to parFit0
lbFit = lb;
ubFit = ub;
parFit0 = par0;
for i = [(1:15), (22:37)]
    lbFit(i) = Safe_log10(lb(i));
    ubFit(i) = Safe_log10(ub(i));
    parFit0(i) = Safe_log10(par0(i));
end

% optimize using fminsearch
% options = optimset('PlotFcns',@optimplotfval);
% [p, fminres] = fminsearch(FitFcn, parFit0, options);

% optimize using simulannealbnd
ini_tem = ones(1, length(par0)) * 100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf}, ...
          'InitialTemperature', ini_tem, 'TemperatureFcn', {@temperatureexp}, ...
          'HybridFcn', {@fmincon}, 'AnnealingFcn', {@newPar}, 'MaxIterations', 1000);
par_est = simulannealbnd(FitFcn, parFit0, lbFit, ubFit, options);

options_fmin = optimset('PlotFcns',@optimplotfval);
[p, fminres] = fminsearch(FitFcn, par_est, options_fmin);

save('best_par.txt', 'p', '-ascii');
