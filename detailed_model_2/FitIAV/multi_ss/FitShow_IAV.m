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

% load fitted parameter set ZT11
par_base = importdata('multi_ss_par_256.txt');
par_base = par_base(1, :);
par_consider_idx = (1:58);

% translate parFit to par in the odes
% log_par_ind = [1:39 43:54];
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end

% load initial values
y0s = importdata('lhs_init.txt');
y0s = y0s(1:1000, :);
[m, n] = size(y0s);

tspan = 0:1:10000;

y = zeros(length(tspan), n, m);

p = parpool(20);
tic;
parfor i = 1:m
    y_i = solve_odes(par_base, par_consider_idx, par_base, y0s(i, :), tspan);
    if size(y_i) == [length(tspan) n]
        y(:, :, i) = y_i;
    else
        y(:, :, i) = 404;
    end
end
toc;
disp(['run time:', num2str(toc)]);
delete(p);

% solve odes
% tmax = 50000;
% tspan = 0:1:tmax;
% y0 = zeros(14, 1);
% y0(3) = 0;
% y0(1) = par_base(55);
% y0(4) = par_base(56);
% y0(7) = par_base(57);
% y0(12) = par_base(58);
% y0(13) = 10;
% 
% [t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_base);
% y = real(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(data_h(1, :), data_h(2, :), 'b');   hold on;
scatter(data_h(1, :), data_h(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 1, i), 'r', 'LineWidth', 2); hold on;
end
xlabel('Time (h)'); ylabel('H (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, y(:, 2, i), 'r', 'LineWidth', 2); hold on;
end
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
for i = 1:m
    plot(tspan, Safe_log10(y(:, 3, i)), 'r', 'LineWidth', 2);  hold on;
end
    xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), data_m(2, :), 'b');   hold on;
scatter(data_m(1, :), data_m(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 4) + y(:, 5, i), 'r', 'LineWidth', 2);  hold on;
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
for i = 1:m
    plot(tspan, y(:, 5, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('M (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_mono(1, :), data_mono(2, :), 'b');   hold on;
scatter(data_mono(1, :), data_mono(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 6, i), 'r', 'LineWidth', 2);  hold on;
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
for i = 1:m
    plot(tspan, y(:, 7, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('Neu (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), data_nk(2, :), 'b');   hold on;
scatter(data_nk(1, :), data_nk(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 12, i), 'r', 'LineWidth', 2);  hold on;
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
for i = 1:m
    plot(tspan, y(:, 13, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('T (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), data_te(2, :), 'b');   hold on;
scatter(data_te(1, :), data_te(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 14, i), 'r', 'LineWidth', 2);  hold on;
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
for i = 1:m
    plot(tspan, y(:, 8, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('IL1b (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_ccl2(1, :), data_ccl2(2, :), 'b');   hold on;
scatter(data_ccl2(1, :), data_ccl2(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 10, i), 'r', 'LineWidth', 2);  hold on;
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
for i = 1:m
    plot(tspan, y(:, 9, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, y(:, 11, i), 'r', 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;