clear;clc;

% load experimental data
data_h = xlsread('data_to_fit.xlsx', 3, 'B1:E3');
data_m = xlsread('data_to_fit.xlsx', 3, 'B5:E7');
data_neu = xlsread('data_to_fit.xlsx', 3, 'B9:E11');
data_nk = xlsread('data_to_fit.xlsx', 3, 'B13:E15');
data_v = xlsread('data_to_fit.xlsx', 3, 'B17:I19');
data_t = xlsread('data_to_fit.xlsx', 3, 'B21:G23');
data_te = xlsread('data_to_fit.xlsx', 3, 'B25:G27');
data_il6 = xlsread('data_to_fit.xlsx', 3, 'B29:E31');
data_ccl2 = xlsread('data_to_fit.xlsx', 3, 'B33:E35');

data_v(2:3, :) = 10 .^ data_v(2:3, :);


% initial parameters
par_bound = importdata('par_IAV_clock_bound.txt');
lb = par_bound.data(64:77, 1);
ub = par_bound.data(64:77, 2);

par = importdata('fit_IAV_ZT11.txt');
par_IAV_fixed = par.data(1:63);
parFit0 = par.data(64:77);

% par_bound = importdata('par_IAV_clock_bound.txt');
% lb = par_bound.data([68; 71], 1);
% ub = par_bound.data([68; 71], 2);
% 
% par = importdata('fit_IAV_ZT11.txt');
% par_IAV_fixed = par.data([1:67 69:70 72:77]);
% parFit0 = par.data([68 71]);

FitFcn = @(parFit)Fit_err_clock_IAV( ...
    data_h, data_m, data_neu, data_nk, data_v, data_t, data_te, data_il6, data_ccl2, par_IAV_fixed, parFit);

% optimize using simulannealbnd
hybridopts = optimoptions('fminunc', 'Display','iter','MaxIterations',10);

ini_tem = ones(1, length(parFit0)) * 100;
options = optimoptions('simulannealbnd','PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf}, ...
          'InitialTemperature', ini_tem, 'TemperatureFcn', {@temperatureexp}, ...
          'HybridFcn', {@fmincon, hybridopts}, 'AnnealingFcn', {@newPar}, 'MaxIterations', 1000);
par_est = simulannealbnd(FitFcn, parFit0, lb, ub, options);


% optimize using fminsearch
options_fmin = optimset('PlotFcns',@optimplotfval);
[p, fminres] = fminsearch(FitFcn, par_est, options_fmin);

par_name = char('beta', 'gamma', 'eta', 'a_NH', 'a_MI', 'a_NI', 'a_KI', ...
        'a_TI', 'a_NV', 'c_IL6_M', 'c_IM', 'c_DM', 'c_CCL2_M', 'c_IN', 'c_CXCL5_N', ...
        'c_M_IL6', 'c_K_IL6', 'c_N_IL6', 'c_I_IL6', 'c_M_IL10', 'c_M_CCL2', ...
        'c_I_CCL2', 'c_M_CXCL5', 'c_I_CXCL5', 'c_IK', 'c_CCL2_K', ...
        'c_MT', 'K_IL6_M', 'K_IM', 'K_DM', 'K_CCL2_M', 'K_IN', 'K_CXCL5_N', 'K_IL10_IL6', ...
        'K_IL10_CCL2', 'K_IL10_CXCL5 ', 'K_IK', 'K_CCL2_K', 'K_T', 'K_MT', 'n_IM', ...
        'n_DM', 'n_CCL2_M', 'n_IN', 'n_CXCL5_N', 'n_IK', 'n_CCL2_K', 'n_MT', 'd_H', 'd_I', ...
        'd_V', 'd_M', 'd_N', 'd_IL6', 'd_IL10', 'd_CCL2', 'd_CXCL5', 'd_K', 'd_T', ...
        'b_H', 'b_M', 'b_N', 'b_K', ...
        'K_B_M', 'K_BMAL1_IL6', 'c_P_K', 'K_P_K', ...
        'K_REV_V', 'K_REV_IL6', 'K_REV_IL10', 'K_REV_CCL2', 'K_REV_CXCL5', ...
        'c_ROR_IL6', 'K_ROR_IL6', 'c_IL6_REV', 'K_IL6_REV', 'n_IL6_REV');

par_fitted = [par_IAV_fixed; p];

file = fopen('fitted_par_IAV_clock.txt', 'w');

for i = 1:length(par_fitted)
    fprintf(file, '%s, %f\n', par_name(i, :), par_fitted(i));
end

fclose(file);
