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

% load fitted parameter set ZT23
parFit = importdata('fitted_par_IAV.txt');
parFit = parFit.data;

% translate parFit to par in the odes
log_par_ind = [1:43 47:60];
par_IAV = parFit;
for i = log_par_ind
    par_IAV(i) = 10 .^ parFit(i);
end

% load inflammation parameters
par_infla = importdata('par_infla.txt');
par_infla = par_infla.data;

% solve odes
numV0 = 20;


tmax = 10000;
tspan = 0:1:tmax;
y_all = zeros(tmax+1, 18, numV0);   % 1:15, immune vars; 16, weight change; 17, inflammation; 18, survival
y0_1 = par_IAV(61);
y0_5 = par_IAV(62);
y0_8 = par_IAV(63);
y0_13 = par_IAV(64);

survival_min = zeros(numV0, 1);



p = parpool(20);

parfor i = 1:numV0
    
    infla = zeros(tmax+1, 1);
    survival = zeros(tmax+1, 1);
    
    par_var = zeros(numV0, 4);
%     par_var = 10 ^ -4.6;
    par_var(i, 1) = 0.2;
    par_var(i, 2) = 0.0005 * 10;
    par_var(i, 3) = 0.0005 * 10;
    par_var(i, 4) = 0.0005 * 10;
    
    y0 = zeros(16, 1);
    y0(4) = 40;
    y0(1) = y0_1;
    y0(5) = y0_5;
    y0(8) = y0_8;
    y0(13) = y0_13;
    y0(14) = 10;
    y0(4) = 100 * (i - 1) + 10;

    [t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV, par_var(i, :), [49, 24, 25, 27], par_infla);
    y = real(y);
    
    for j = 1:(tmax+1)
        infla(j) = inflammation(y(j, 6), y(j, 7), y(j, 8), par_infla);
        survival(j) = (1 - y(j, 16) .^ 5 ./ (24 ^ 5 + y(j, 16) .^ 5));
    end
    y_all(:, :, i) = [y infla survival];
    
    survival_min(i) = min(survival);
end
delete(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
dq = jet(numV0);

figure;
xSize = 30; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,3,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_h(1, :), data_h(2, :), 'b');   hold on;
scatter(data_h(1, :), data_h(3, :), 'r');   hold on;

for i = 1:numV0
    plot(tspan, y_all(:, 1, i), 'color', dq(i, :), 'LineWidth', 2); hold on;
end

xlabel('Time (h)'); ylabel('H (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,3,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 2), 'b', 'LineWidth', 2); hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 2, i), 'color', dq(i, :), 'LineWidth', 2); hold on;
end
xlabel('Time (h)'); ylabel('I (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,3,3); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :), Safe_log10(data_v(2, :)), 'b');   hold on;
scatter(data_v(1, :), Safe_log10(data_v(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_23(:, 4)), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
%     plot(tspan, Safe_log10(y_all(:, 4, i)), 'color', dq(i, :), 'LineWidth', 2);  hold on;
    plot(tspan, y_all(:, 4, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :), Safe_log10(data_v(2, :)), 'b');   hold on;
scatter(data_v(1, :), Safe_log10(data_v(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_23(:, 4)), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, Safe_log10(y_all(:, 4, i)), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 3), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 3, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('D (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), data_m(2, :), 'b');   hold on;
scatter(data_m(1, :), data_m(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 5), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 5, i) + y_all(:, 6, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('M0 (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 5), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 6, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('M (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_mono(1, :), data_mono(2, :), 'b');   hold on;
scatter(data_mono(1, :), data_mono(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 7), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 7, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
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
% plot(tspan, y_23(:, 8), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 8, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('Neu (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), data_nk(2, :), 'b');   hold on;
scatter(data_nk(1, :), data_nk(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 13), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 13, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
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
% plot(tspan, y_23(:, 14), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 14, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('T (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), data_te(2, :), 'b');   hold on;
scatter(data_te(1, :), data_te(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 15), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 15, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('TE (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% scatter(data_il1b(1, :), data_il1b(2, :), 'b');   hold on;
% scatter(data_il1b(1, :), data_il1b(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 9), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 9, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('IL1b (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_ccl2(1, :), data_ccl2(2, :), 'b');   hold on;
scatter(data_ccl2(1, :), data_ccl2(3, :), 'r');   hold on;
% plot(tspan, y_23(:, 11), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 11, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('CCL2 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 10), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 10, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% plot(tspan, y_23(:, 12), 'b', 'LineWidth', 2);  hold on;
for i = 1:numV0
    plot(tspan, y_all(:, 12, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 30; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,3,1); hold on; set(gca,'Fontsize',26); box on;
for i = 1:numV0
    plot(tspan, y_all(:, 17, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('inflammation'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,3,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:numV0
    plot(tspan, y_all(:, 16, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('weight loss'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,3,3); hold on; set(gca,'Fontsize',26); box on;
for i = 1:numV0
    plot(tspan, y_all(:, 18, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('survival rate'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% xSize = 30; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
set(gca,'Fontsize',26); box on;
% plot(par_var(:, 1), survival_min, 'color', 'r', 'LineWidth', 2);  hold on;
plot(survival_min, 'color', 'r', 'LineWidth', 2);  hold on;
xlabel('parameter'); ylabel('survival rate'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;