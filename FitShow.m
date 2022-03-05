% load experimental data
data_m_nk = xlsread('data_to_fit.xlsx', 2, 'B13:E15');
data_m = data_m_nk([1 2], :);
data_nk = data_m_nk([1 3], :);

data_v = xlsread('data_to_fit.xlsx', 2, 'B17:I18');

data_t_te = xlsread('data_to_fit.xlsx', 2, 'B20:C22');
data_t = data_t_te([1 2], :);
data_te = data_t_te([1 3], :);

% calculate cell numbers with H = 10 ^ 7
data_m(2, :) = data_m(2, :) .* 10 ^ 5;
data_nk(2, :) = data_nk(2, :) .* 10 ^ 5;
data_t(2, :) = data_t(2, :) .* 10 ^ 5;
data_te(2, :) = data_te(2, :) .* 10 ^ 5;

data_v(2, :) = 10 .^ data_v(2, :);

% load fitted parameter set
parFit = load('best_par.txt');

% translate parFit to par in the odes
par = parFit;
for i = [(1:15), (22:37)]
    par(i) = 10 .^ parFit(i);
end

% solve odes
tmax = 250;
tspan = 0:1:tmax;
y0 = [10^7, 40, 0, 0, 0, 0, 0, 0];
[t, y] = ode45(@ODE_IAV, tspan, y0, [], par);

% Calcute the actual number of immune cells
% by adding the baseline cell numbers
par(34) = 10000;
par(35) = 10000;
par(36) = 10000;
par(37) = 0;
y(:, 5) = y(:, 5) + par(34);
y(:, 6) = y(:, 6) + par(35);
y(:, 7) = y(:, 7) + par(36);
y(:, 8) = y(:, 8) + par(37);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, log10(y(:, 1)), 'k', 'LineWidth', 2); hold on;
xlabel('Time (h)'); ylabel('H (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, log10(y(:, 2)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('I (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_v(1, :), log10(data_v(2, :)));   hold on;
plot(tspan, log10(y(:, 3)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, log10(y(:, 4)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('D (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), log10(data_m(2, :)));   hold on;
plot(tspan, log10(y(:, 5)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('M (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), log10(data_nk(2, :)));   hold on;
plot(tspan, log10(y(:, 6)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('NK (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_t(1, :), log10(data_t(2, :)));   hold on;
plot(tspan, log10(y(:, 7)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('CD8T (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), log10(data_te(2, :)));   hold on;
plot(tspan, log10(y(:, 8)), 'k', 'LineWidth', 2);
xlabel('Time (h)'); ylabel('Effector CD8T (log_{10} cells)');
set(gca, 'XTick', [0:100:200], 'XLim', [-5 260], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;