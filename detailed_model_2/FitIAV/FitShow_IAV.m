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
data_il6 = xlsread('data_to_fit.xlsx', 3, 'B33:E35');
data_ccl2 = xlsread('data_to_fit.xlsx', 3, 'B37:E39');

data_v(2:3, :) = 10 .^ data_v(2:3, :);

% load fitted parameter set ZT11
parFit_11 = importdata('fitted_par_IAV.txt');
parFit_11 = parFit_11.data;

% translate parFit to par in the odes
log_par_ind = [1:34 36:48];
par_IAV_11 = parFit_11;
for i = log_par_ind
    par_IAV_11(i) = 10 .^ parFit_11(i);
end

% % load fitted parameter set ZT23
% parFit_23 = importdata('fit_ZT23.txt');
% parFit_23 = parFit_23.data;
% 
% % translate parFit to par in the odes
% log_par_ind = [1:40 49:59];
% par_IAV_23 = parFit_23;
% for i = log_par_ind
%     par_IAV_23(i) = 10 .^ parFit_23(i);
% end


% solve odes
tmax = 250;
tspan = 0:1:tmax;
y0 = zeros(14, 1);
y0(4) = 40;
y0(1) = par_IAV_11(49);
y0(5) = par_IAV_11(50);
y0(7) = par_IAV_11(51);
y0(12) = par_IAV_11(52);
y0(13) = 10;


% For ZT11
[t, y_11] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV_11);
y_11 = real(y_11);

% For ZT23
% [t, y_23] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV_23);
% y_23 = real(y_23);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_h(1, :), data_h(2, :), 'b');   hold on;
scatter(data_h(1, :), data_h(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 1), 'b', 'LineWidth', 2); hold on;
plot(tspan, y_11(:, 1), 'r', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('H (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 2), 'b', 'LineWidth', 2); hold on;
plot(tspan, y_11(:, 2), 'r', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('I (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :), Safe_log10(data_v(2, :)), 'b');   hold on;
scatter(data_v(1, :), Safe_log10(data_v(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_23(:, 4)), 'b', 'LineWidth', 2);  hold on;
plot(tspan, Safe_log10(y_11(:, 4)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 3), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 3), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('D (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), data_m(2, :), 'b');   hold on;
scatter(data_m(1, :), data_m(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 5), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 5), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('M (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_mono(1, :), data_mono(2, :), 'b');   hold on;
scatter(data_mono(1, :), data_mono(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 6), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 6), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Mono (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_neu(1, :), data_neu(2, :), 'b');   hold on;
scatter(data_neu(1, :), data_neu(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 7), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 7), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Neu (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), data_nk(2, :), 'b');   hold on;
scatter(data_nk(1, :), data_nk(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 12), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 12), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('NK (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% scatter(data_t(1, :), data_t(2, :), 'b');   hold on;
% plot(tspan, y_23(:, 13), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 13), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('T (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), data_te(2, :), 'b');   hold on;
scatter(data_te(1, :), data_te(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 14), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 14), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('TE (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_il6(1, :), data_il6(2, :), 'b');   hold on;
scatter(data_il6(1, :), data_il6(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 8), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 8), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('IL6 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_ccl2(1, :), data_ccl2(2, :), 'b');   hold on;
scatter(data_ccl2(1, :), data_ccl2(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 10), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 10), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('CCL2 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 9), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 9), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 11), 'b', 'LineWidth', 2);  hold on;
plot(tspan, y_11(:, 11), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;