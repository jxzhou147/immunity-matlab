clear;clc;

% load experimental data
% data_h = xlsread('data_to_fit.xlsx', 3, 'B1:E3');
% data_m = xlsread('data_to_fit.xlsx', 3, 'B5:E7');
% data_mono = xlsread('data_to_fit.xlsx', 3, 'B9:E11');
% data_neu = xlsread('data_to_fit.xlsx', 3, 'B13:E15');
% data_nk = xlsread('data_to_fit.xlsx', 3, 'B17:E19');
% data_v = xlsread('data_to_fit.xlsx', 3, 'B21:I23');
% data_t = xlsread('data_to_fit.xlsx', 3, 'B25:G27');
% data_te = xlsread('data_to_fit.xlsx', 3, 'B29:G31');
% data_il1b = xlsread('data_to_fit.xlsx', 3, 'B33:E35');
% data_ccl2 = xlsread('data_to_fit.xlsx', 3, 'B37:E39');
%     
% data_v(2:3, :) = 10 .^ data_v(2:3, :);

% load immunity parameters
par_IAV = importdata('par_base_IAV_values.txt');
par_IAV = par_IAV.data;

% load inflammation parameters
par_infla = importdata('par_infla.txt');
par_infla = par_infla.data;

y0 = importdata('init_base.txt');
y0 = y0.data;

tmax = 1000;
tspan = 0:1:tmax;

[t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV);

infla = zeros(tmax+1, 1);
for i = 1:(tmax+1)
    infla(i) = inflammation(y(i, 5), y(i, 6), y(i, 7), par_infla);
end

%%

var_name = ["H", "If", "V", "M_0", "M", "Mono", "N", "IL1b", "IL10", "CCL2", "CXCL5", "K", "T", "T_E"];
figure;
for i = 1:14
    subplot(4, 4, i);
    plot(tspan, y(:, i)); hold on;
    ylabel(var_name(i)); hold off;
end

subplot(4, 4, 15);
plot(tspan, infla); hold on;
xlabel('Time (h)'); ylabel('Inflammation'); hold off;
% legend('Only c_{M-IL10} oscillates', 'Both oscillate', 'Only c_{M-CCL2} oscillates', 'Location',[0.8 0.2 0 0]);
% sgtitle('Circadian control of c_{M-CCL2} and c_{M-IL10}');