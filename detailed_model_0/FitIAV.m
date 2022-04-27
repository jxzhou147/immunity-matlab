clear;clc;

% load experimental data
data_m_nk = xlsread('data_to_fit.xlsx', 2, 'B13:E15');
data_m = data_m_nk([1 2], :);
data_nk = data_m_nk([1 3], :);

data_v = xlsread('data_to_fit.xlsx', 2, 'B17:I18');

data_t_te = xlsread('data_to_fit.xlsx', 2, 'B20:G22');
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
for i = [(1:12), (17:27)]
    lbFit(i) = Safe_log10(lb(i));
    ubFit(i) = Safe_log10(ub(i));
    parFit0(i) = Safe_log10(par0(i));
end

% optimize using simulannealbnd
hybridopts = optimoptions('fminunc', 'Display','iter');

ini_tem = ones(1, length(par0)) * 100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf}, ...
          'InitialTemperature', ini_tem, 'TemperatureFcn', {@temperatureexp}, ...
          'HybridFcn', {@fmincon, hybridopts}, 'AnnealingFcn', {@newPar}, 'MaxIterations', 1000);
par_est = simulannealbnd(FitFcn, parFit0, lbFit, ubFit, options);

% optimize using fminsearch
options_fmin = optimset('PlotFcns',@optimplotfval);
[p, fminres] = fminsearch(FitFcn, par_est, options_fmin);

par_name = char('beta', 'gamma', 'eta', 'alpha_MI ', 'alpha_KI ', 'alpha_TI ', ...
        'alpha_MV ', 'c_IM', 'c_M', 'c_IK', 'c_MK', 'c_MT', 'n_IM', ...
        'n_IK', 'n_MK', 'n_MT', 'K_IM', 'K_M', 'K_IK', 'K_MK', 'K_MT', 'K_T', ...
        'K_TI', 'd_I', 'd_M', 'd_K', 'd_T');
file = fopen('best_par.txt', 'w');
for i = 1:length(par_name)
fprintf(file, '%s, %f\n', par_name(i, :), p(i));
end
fclose(file);
