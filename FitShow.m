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

% load fitted parameter set
parFit = importdata('best_par.txt');
parFit = parFit.data;

% translate parFit to par in the odes
par = parFit;
for i = [(1:12), (17:27)]
    par(i) = 10 .^ parFit(i);
end

% solve odes
tmax = 250;
tspan = 0:1:tmax;
y0 = [10^7, 40, 0, 0, 0, 1000, 0];
[t, y] = ode45(@ODE_IAV, tspan, y0, [], par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, Safe_log10(y(:, 1)), 'k', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('H (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, Safe_log10(y(:, 2)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('I (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :), Safe_log10(data_v(2, :)));   hold on;
plot(tspan, Safe_log10(y(:, 3)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), Safe_log10(data_m(2, :)));   hold on;
plot(tspan, Safe_log10(y(:, 4)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('M (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), Safe_log10(data_nk(2, :)));   hold on;
plot(tspan, Safe_log10(y(:, 5)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('NK (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_t(1, :), Safe_log10(data_t(2, :)));   hold on;
plot(tspan, Safe_log10(y(:, 6)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('CD8T (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), Safe_log10(data_te(2, :)));   hold on;
plot(tspan, Safe_log10(y(:, 7)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('Effector CD8T (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;