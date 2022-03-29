clear;clc;

% load experimental data
data_m = xlsread('data_to_fit.xlsx', 3, 'B1:E3');
data_nk = xlsread('data_to_fit.xlsx', 3, 'B5:E7');
data_v = xlsread('data_to_fit.xlsx', 3, 'B9:I11');
data_t = xlsread('data_to_fit.xlsx', 3, 'B13:G15');
data_te = xlsread('data_to_fit.xlsx', 3, 'B17:G19');

% calculate cell numbers with H = 10 ^ 7
data_m(2, :) = data_m(2, :) .* 10 ^ 5;
data_nk(2, :) = data_nk(2, :) .* 10 ^ 5;
data_t(2, :) = data_t(2, :) .* 10 ^ 5;
data_te(2, :) = data_te(2, :) .* 10 ^ 5;

data_v(2, :) = 10 .^ data_v(2, :);

data_m(3, :) = data_m(3, :) .* 10 ^ 5;
data_nk(3, :) = data_nk(3, :) .* 10 ^ 5;
data_t(3, :) = data_t(3, :) .* 10 ^ 5;
data_te(3, :) = data_te(3, :) .* 10 ^ 5;

data_v(3, :) = 10 .^ data_v(3, :);

% load fitted parameter set
parFit = importdata('test_of_best_par_IAV_clock.txt');
parFit = parFit.data;

% translate parFit to par in the odes
par_IAV = parFit;
for i = [(1:12), (17:27)]
    par_IAV(i) = 10 .^ parFit(i);
end

par_clock = load('par_clock.csv');

% solve odes
tmax = 400;
tspan = 0:1:tmax;
y0 = zeros(19, 1);
y0(13) = 10 ^ 7;
y0(14) = 40;
y0(16) = 10 ^ 6;
y0(18) = 1000;

% For ZT23
t_IAV_23 = 115;
[t, y_23] = ode45(@ODE_Clock_IAV_try, tspan, y0, [], par_clock, par_IAV, t_IAV_23);

% For ZT11
t_IAV_11 = 127;
[t, y_11] = ode45(@ODE_Clock_IAV_try, tspan, y0, [], par_clock, par_IAV, t_IAV_11);

y_23 = real(y_23);
y_11 = real(y_11);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, Safe_log10(y_23(:, 13)), 'b', 'LineWidth', 2); hold on;
plot(tspan - 12, Safe_log10(y_11(:, 13)), 'r', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('H (log_{10} cells)');
set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, Safe_log10(y_23(:, 14)), 'b', 'LineWidth', 2);  hold on;
plot(tspan - 12, Safe_log10(y_11(:, 14)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('I (log_{10} cells)');
set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :) + t_IAV_23, Safe_log10(data_v(2, :)), 'b');   hold on;
plot(tspan, Safe_log10(y_23(:, 15)), 'b', 'LineWidth', 2);  hold on;
scatter(data_v(1, :) + t_IAV_11 - 12, Safe_log10(data_v(3, :)), 'r');   hold on;
plot(tspan - 12, Safe_log10(y_11(:, 15)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% scatter(data_te(1, :) + t_IAV_23, Safe_log10(data_te(2, :)), 'b');   hold on;
% plot(tspan, Safe_log10(y_23(:, 19)), 'b', 'LineWidth', 2);  hold on;
% scatter(data_te(1, :) + t_IAV_11, Safe_log10(data_te(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_11(:, 19)), 'r', 'LineWidth', 2);  hold on;
scatter(data_te(1, :) + t_IAV_23, Safe_log10(data_te(2, :)), 'b');   hold on;
plot(tspan, Safe_log10(y_23(:, 19)), 'b', 'LineWidth', 2);  hold on;
scatter(data_te(1, :) + t_IAV_11 - 12, Safe_log10(data_te(3, :)), 'r');   hold on;
plot(tspan - 12, Safe_log10(y_11(:, 19)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('CD8T_E (log_{10} cells)');
set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :) + t_IAV_23, data_m(2, :), 'b');   hold on;
plot(tspan, y_23(:, 16), 'b', 'LineWidth', 2);  hold on;
scatter(data_m(1, :) + t_IAV_11 - 12, data_m(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 16), 'r', 'LineWidth', 2);  hold on;
% scatter(data_m(1, :) + t_IAV_23, Safe_log10(data_m(2, :)), 'b');   hold on;
% plot(tspan, Safe_log10(y_23(:, 16)), 'b', 'LineWidth', 2);  hold on;
% scatter(data_m(1, :) + t_IAV_11, Safe_log10(data_m(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_11(:, 16)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('M (cells)');
set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'Fontsize', 26, 'linewidth', 2);
%set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :) + t_IAV_23, data_nk(2, :), 'b');   hold on;
plot(tspan, y_23(:, 17), 'b', 'LineWidth', 2);  hold on;
scatter(data_nk(1, :) + t_IAV_11, data_nk(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 17), 'r', 'LineWidth', 2);  hold on;
% scatter(data_nk(1, :) + t_IAV_23, Safe_log10(data_nk(2, :)), 'b');   hold on;
% plot(tspan, Safe_log10(y_23(:, 17)), 'b', 'LineWidth', 2);  hold on;
% scatter(data_nk(1, :) + t_IAV_11, Safe_log10(data_nk(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_11(:, 17)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('NK (cells)');
set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'Fontsize', 26, 'linewidth', 2);
%set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :) + t_IAV_23, data_m(2, :), 'b');   hold on;
plot(tspan, y_23(:, 16), 'b', 'LineWidth', 2);  hold on;
scatter(data_m(1, :) + t_IAV_11 - 12, data_m(3, :), 'r');   hold on;
plot(tspan - 12, y_11(:, 16), 'r', 'LineWidth', 2);  hold on;
% scatter(data_m(1, :) + t_IAV_23, Safe_log10(data_m(2, :)), 'b');   hold on;
% plot(tspan, Safe_log10(y_23(:, 16)), 'b', 'LineWidth', 2);  hold on;
% scatter(data_m(1, :) + t_IAV_11, Safe_log10(data_m(3, :)), 'r');   hold on;
% plot(tspan, Safe_log10(y_11(:, 16)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('M (cells)');
set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'Fontsize', 26, 'linewidth', 2);
%set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :) + t_IAV_23, Safe_log10(data_v(2, :)), 'b');   hold on;
plot(tspan, Safe_log10(y_23(:, 15)), 'b', 'LineWidth', 2);  hold on;
scatter(data_v(1, :) + t_IAV_11 - 12, Safe_log10(data_v(3, :)), 'r');   hold on;
plot(tspan - 12, Safe_log10(y_11(:, 15)), 'r', 'LineWidth', 2);  hold on;
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
% hold on;
% subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
% % scatter(data_t(1, :) + t_IAV_23, Safe_log10(data_t(2, :)), 'b');   hold on;
% % plot(tspan, Safe_log10(y_23(:, 18)), 'b', 'LineWidth', 2);  hold on;
% % scatter(data_t(1, :) + t_IAV_11, Safe_log10(data_t(3, :)), 'r');   hold on;
% % plot(tspan, Safe_log10(y_11(:, 18)), 'r', 'LineWidth', 2);  hold on;
% scatter(data_t(1, :) + t_IAV_23, data_t(2, :), 'b');   hold on;
% plot(tspan, y_23(:, 18), 'b', 'LineWidth', 2);  hold on;
% scatter(data_t(1, :) + t_IAV_11, data_t(3, :), 'r');   hold on;
% plot(tspan - 12, y_11(:, 18), 'r', 'LineWidth', 2);  hold on;
% xlabel('Time (h)'); ylabel('CD8T (log_{10} cells)');
% %set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
% hold on;
% 
% subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
% % scatter(data_te(1, :) + t_IAV_23, Safe_log10(data_te(2, :)), 'b');   hold on;
% % plot(tspan, Safe_log10(y_23(:, 19)), 'b', 'LineWidth', 2);  hold on;
% % scatter(data_te(1, :) + t_IAV_11, Safe_log10(data_te(3, :)), 'r');   hold on;
% % plot(tspan, Safe_log10(y_11(:, 19)), 'r', 'LineWidth', 2);  hold on;
% scatter(data_te(1, :) + t_IAV_23, data_te(2, :), 'b');   hold on;
% plot(tspan, y_23(:, 19), 'b', 'LineWidth', 2);  hold on;
% scatter(data_te(1, :) + t_IAV_11, data_te(3, :), 'r');   hold on;
% plot(tspan - 12, y_11(:, 19), 'r', 'LineWidth', 2);  hold on;
% xlabel('Time (h)'); ylabel('Effector CD8T (log_{10} cells)');
% % set(gca, 'XTick', [100:250:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
% hold off;