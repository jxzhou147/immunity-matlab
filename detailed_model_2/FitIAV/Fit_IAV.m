clear;clc;

% load experimental data
data_h = xlsread('data_to_fit.xlsx', 3, 'B1:E3');
data_m = xlsread('data_to_fit.xlsx', 3, 'B5:E7');
data_mono = xlsread('data_to_fit.xlsx', 3, 'B9:E11');
data_neu = xlsread('data_to_fit.xlsx', 3, 'B13:E15');
data_nk = xlsread('data_to_fit.xlsx', 3, 'B17:E19');
data_v = xlsread('data_to_fit.xlsx', 3, 'B21:I23');
data_t = xlsread('data_to_fit.xlsx', 3, 'B25:G27');
data_te = xlsread('data_to_fit.xlsx', 3, 'B29:G31');
data_il1b = xlsread('data_to_fit.xlsx', 3, 'B33:E35');
data_ccl2 = xlsread('data_to_fit.xlsx', 3, 'B37:E39');

data_v(2:3, :) = 10 .^ data_v(2:3, :);

parBase = [350; 80; 15; 40];

FitFcn = @(parFit)Fit_err_IAV( ...
    data_h, data_m, data_mono, data_neu, data_nk, data_v, data_t, data_te, data_il1b, data_ccl2, parFit, parBase);

% initial parameters
par_bound = importdata('par_IAV_bound.txt');
lb = par_bound.data(1:60, 1);
ub = par_bound.data(1:60, 2);

parFit0 = importdata('fitted_par_IAV.txt');
parFit0 = parFit0.data(1:60);


% optimize using simulannealbnd
% hybridopts = optimoptions('fminunc', 'Display','iter','MaxIterations',10);
% 
% ini_tem = ones(1, length(parFit0)) * 500;
% options = optimoptions('simulannealbnd','PlotFcns',...
%           {@saplotbestx,@saplotbestf,@saplotx,@saplotf}, ...
%           'InitialTemperature', ini_tem, 'TemperatureFcn', {@temperatureexp}, ...
%           'AnnealingFcn', {@newPar}, 'MaxIterations', 3000);
% %           'HybridFcn', {@fmincon, hybridopts}, 'AnnealingFcn', {@newPar}, 'MaxIterations', 500);
% par_est = simulannealbnd(FitFcn, parFit0, lb, ub, options);


% optimize using fminsearch
options_fmin = optimset('PlotFcns',@optimplotfval);
[p, fminres] = fminsearch(FitFcn, parFit0, options_fmin);

par_name = char('beta', 'gamma', 'eta', 'a_NH', 'a_MH', 'a_MI', 'a_NI', 'a_KI', ...
        'a_TI', 'a_MD', 'a_NV', 'c_IL1b_M', 'c_IL1b_Mono', 'c_IL1b_N', 'c_IL1b_K', 'c_IM', 'c_CCL2_Mono', 'c_N', 'c_CXCL5_N', ...
        'c_M_IL1b', 'c_I_IL1b','c_D_IL1b', 'c_M_IL10', 'c_M_CCL2', ...
        'c_I_CCL2', 'c_N_CXCL5', 'c_I_CXCL5', 'c_K', 'c_MK', 'c_IK', 'c_MT', ...
        'K_IL1b_M', 'K_IL1b_Mono', 'K_IL1b_N', 'K_IL1b_K', 'K_MI', 'K_I_M', 'K_D_IL1b', 'K_IL10_IL1b', 'K_IL10_CCL2', 'K_IL10_CXCL5 ', 'K_T', 'K_MT', 'n_I_M', 'n_D_IL1b', 'n_MT', ...
        'd_H', 'd_I', 'd_V', 'd_M0', 'd_M', 'd_Mono', 'd_N', 'd_IL1b', 'd_IL10', 'd_CCL2', 'd_CXCL5', 'd_K', 'd_T', ...
        'Mono_tot', 'b_H', 'b_M0', 'b_N', 'b_K');
file = fopen('fitted_par_IAV.txt', 'w');
for i = 1:length(p)
    fprintf(file, '%s, %f\n', par_name(i, :), p(i));
end
for i = 1:4
    fprintf(file, '%s, %f\n', par_name(i + length(p), :), parBase(i));
end

fclose(file);
