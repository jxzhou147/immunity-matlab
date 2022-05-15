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

% load fitted parameter set
parFit = importdata('fitted_par_IAV_clock.txt');
parFit = parFit.data;

% translate parFit to par in the odes
log_par_ind = [1:40 49:59 62 69 71:72];
par_IAV = parFit;
for i = log_par_ind
    par_IAV(i) = 10 .^ parFit(i);
end

par_clock = load('par_clock.csv');

% solve odes
tmax = 400;
tspan = 0:1:tmax;
y0 = zeros(25, 1);
y0(16) = 40;
y0(13) = par_IAV(74);
y0(17) = par_IAV(75);
y0(18) = par_IAV(76);
y0(23) = par_IAV(77);
y0(24) = 10;


% For ZT23
t_IAV_23 = 115;
[t, y_23] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23);

% For ZT11
t_IAV_11 = 127;
[t, y_11] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_11);

y_23 = real(y_23);
y_11 = real(y_11);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_h(1, :) + t_IAV_23, data_h(2, :), 'b');   hold on;
plot(tspan, y_23(:, 13), 'b', 'LineWidth', 2); hold on;
scatter(data_h(1, :) + t_IAV_11 - 12, data_h(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 13), 'r', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('H (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, y_23(:, 14), 'b', 'LineWidth', 2); hold on;
plot(tspan - 12, y_11(:, 14), 'r', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('I (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :) + t_IAV_23, Safe_log10(data_v(2, :)), 'b');   hold on;
plot(tspan, Safe_log10(y_23(:, 16)), 'b', 'LineWidth', 2);  hold on;
scatter(data_v(1, :) + t_IAV_11 - 12, Safe_log10(data_v(3, :)), 'r');   hold on;
plot(tspan - 12, Safe_log10(y_11(:, 16)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :) + t_IAV_23, data_m(2, :), 'b');   hold on;
plot(tspan, y_23(:, 17), 'b', 'LineWidth', 2);  hold on;
scatter(data_m(1, :) + t_IAV_11 - 12, data_m(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 17), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('M (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_neu(1, :) + t_IAV_23, data_neu(2, :), 'b');   hold on;
plot(tspan, y_23(:, 18), 'b', 'LineWidth', 2);  hold on;
scatter(data_neu(1, :) + t_IAV_11 - 12, data_neu(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 18), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Neu (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :) + t_IAV_23, data_nk(2, :), 'b');   hold on;
plot(tspan, y_23(:, 23), 'b', 'LineWidth', 2);  hold on;
scatter(data_nk(1, :) + t_IAV_11 - 12, data_nk(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 23), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('NK (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_t(1, :) + t_IAV_23, data_t(2, :), 'b');   hold on;
plot(tspan, y_23(:, 24), 'b', 'LineWidth', 2);  hold on;
scatter(data_t(1, :) + t_IAV_11 - 12, data_t(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 24), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('T (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :) + t_IAV_23, data_te(2, :), 'b');   hold on;
plot(tspan, y_23(:, 25), 'b', 'LineWidth', 2);  hold on;
scatter(data_te(1, :) + t_IAV_11 - 12, data_te(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 25), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('TE (cells)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_il6(1, :) + t_IAV_23, data_il6(2, :), 'b');   hold on;
plot(tspan, y_23(:, 19), 'b', 'LineWidth', 2);  hold on;
scatter(data_il6(1, :) + t_IAV_11 - 12, data_il6(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 19), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('IL6 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_ccl2(1, :) + t_IAV_23, data_ccl2(2, :), 'b');   hold on;
plot(tspan, y_23(:, 21), 'b', 'LineWidth', 2);  hold on;
scatter(data_ccl2(1, :) + t_IAV_11 - 12, data_ccl2(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 21), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('CCL2 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, y_23(:, 20), 'b', 'LineWidth', 2);  hold on;
plot(tspan - 12, y_11(:, 20), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, y_23(:, 22), 'b', 'LineWidth', 2);  hold on;
plot(tspan - 12, y_11(:, 22), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;