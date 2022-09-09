% fit relation between weight change wc and survival curve
clear;clc

% load experimental data
data1 = xlsread('data_to_fit.xlsx', 4, 'B6:C15');
data2 = xlsread('data_to_fit.xlsx', 4, 'F6:G15');

data3 = xlsread('data_to_fit.xlsx', 4, 'B24:C31');
data4 = xlsread('data_to_fit.xlsx', 4, 'F24:G31');
data5 = xlsread('data_to_fit.xlsx', 4, 'J24:K30');
data6 = xlsread('data_to_fit.xlsx', 4, 'N24:O30');

data = cat(1, data1, data2, data3, data4, data5, data6);
data(:, 1) = abs(data(:, 1));

FitFcn = @(par)Fit_err_survival_wc(data, par);

par0 = [10; 3];

% optimize using fminsearch
options_fmin = optimset('PlotFcns',@optimplotfval);
[p, fminres] = fminsearch(FitFcn, par0, options_fmin);

% plot fitted curve vs experimental data
figure;
wc = (1:0.1:50);
s = (1 - wc .^ p(2) ./ (p(1) ^ p(2) + wc .^ p(2)));
plot(wc, s, 'LIneWidth', 2); hold on;
scatter(data(:, 1), data(:, 2));
xlabel('weight loss (%)'); ylabel('survival rate'); box on;
set(gca, 'FontSize', 26);
hold off;