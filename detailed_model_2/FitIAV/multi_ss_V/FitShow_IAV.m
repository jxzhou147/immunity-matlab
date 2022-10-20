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
par_base = importdata('xppaut\par_cmccl2_bifurcation.txt');
par_base = par_base.data;

% translate parFit to par in the odes
% log_par_ind = [1:38 42:53];
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end
% 
% par_consider = importdata('tight_parset\lhs_par.txt');
% par_consider = par_consider(1, :);
% par_consider_idx = [14, 20, 21, 33, 40, 46, 50];

% par_base = importdata('wide_parset_data\multi_ss_par_488_in_bound.txt');
% par_base = par_base(1, :);
par_consider = par_base;
par_consider_idx = (1:57);


% load initial values
y0s = importdata('tight_parset\lhs_init.txt');
% y0s = importdata('wide_parset_data\lhs_init.txt');
y0s = y0s(1:100, :);
[m, n] = size(y0s);

% y0s_low = importdata('tight_parset\ss_init_low.txt');
% [m_low, ~] = size(y0s_low);

% m = 7;
% n = 14;
% 
% y0_base = importdata('init_base.txt');
% y0_base = y0_base.data;
% y0s = zeros(m, n);
% for i = 1:m
%     y0s(i, :) = y0_base;
%     y0s(i, 3) = 10 * (i - 1);
% end

tspan = 0:1:10000;

y = zeros(length(tspan), n, m);
TF_arr = false(m, 1);

% y_low = zeros(length(tspan), n, m_low);
% TF_arr_low = false(m_low, 1);

% ss_init = zeros(m, n);  % steady state values(n var) of each initial for a given parameter set

p = parpool(20);
tic;
parfor i = 1:m
    y_i = solve_odes(par_base, par_consider_idx, par_consider, y0s(i, :), tspan);
    if size(y_i) == [length(tspan) n]
        y(:, :, i) = y_i;
%         ss_init(i, :) = y_i(end, :);
    else
        TF_arr(i) = true;
    end
end

% write ss_init to a file
% file_id = fopen('tight_parset\ss_init.txt', 'w');
% for i = 1:m
%     fprintf(file_id, '%d ', i);
%     fprintf(file_id, '%f ', ss_init(i, :));
%     fprintf(file_id, '\n');
% end

y(:, :, TF_arr) = [];
% ss_init(TF_arr, :) = [];
[~, ~, m] = size(y);

% parfor i = 1:m_low
%     y_i = solve_odes(par_base, par_consider_idx, par_base, y0s_low(i, :), tspan);
%     if size(y_i) == [length(tspan) n]
%         y_low(:, :, i) = y_i;
% %         ss_init(i, :) = y_i(end, :);
%     else
%         TF_arr_low(i) = true;
%     end
% end
% y_low(:, :, TF_arr_low) = [];
% % ss_init(TF_arr, :) = [];
% [~, ~, m_low] = size(y_low);

toc;
disp(['run time:', num2str(toc)]);
delete(p);

%% 
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 1, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('H (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, y(:, 2, i), 'r', 'LineWidth', 2); hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 2, i), 'b', 'LineWidth', 2); hold on;
% end
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
%     plot(tspan, Safe_log10(y(:, 3, i)), 'r', 'LineWidth', 2);  hold on;
    plot(tspan,y(:, 3, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 3, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('Virus pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_m(1, :), data_m(2, :), 'b');   hold on;
scatter(data_m(1, :), data_m(3, :), 'r');   hold on;
for i = 1:m
%     plot(tspan, y(:, 4) + y(:, 5, i), 'r', 'LineWidth', 2);  hold on;
    plot(tspan, y(:, 4, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 4, i), 'b', 'LineWidth', 2); hold on;
% end
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 5, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('M (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_mono(1, :), data_mono(2, :), 'b');   hold on;
scatter(data_mono(1, :), data_mono(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 6, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 6, i), 'b', 'LineWidth', 2); hold on;
% end
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 7, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('Neu (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_nk(1, :), data_nk(2, :), 'b');   hold on;
scatter(data_nk(1, :), data_nk(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 12, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 12, i), 'b', 'LineWidth', 2); hold on;
% end
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 13, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('T (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_te(1, :), data_te(2, :), 'b');   hold on;
scatter(data_te(1, :), data_te(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 14, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 14, i), 'b', 'LineWidth', 2); hold on;
% end
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 8, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('IL1b (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(data_ccl2(1, :), data_ccl2(2, :), 'b');   hold on;
scatter(data_ccl2(1, :), data_ccl2(3, :), 'r');   hold on;
for i = 1:m
    plot(tspan, y(:, 10, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 10, i), 'b', 'LineWidth', 2); hold on;
% end
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
% for i = 1:m_low
%     plot(tspan, y_low(:, 9, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, y(:, 11, i), 'r', 'LineWidth', 2);  hold on;
end
% for i = 1:m_low
%     plot(tspan, y_low(:, 11, i), 'b', 'LineWidth', 2); hold on;
% end
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;