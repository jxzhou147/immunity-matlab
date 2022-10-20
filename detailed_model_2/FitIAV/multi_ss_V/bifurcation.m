% plot bifurcation diagram
clear;clc
% fix 1 parameter set and let V be the consider par
par_base = importdata('tight_parset\par_base.txt');
% par_base = par_base(2, :);
par_base = par_base.data;

% translate parFit to par in the odes
log_par_ind = [1:38 42:53];
for i = log_par_ind
    par_base(i) = 10 .^ par_base(i);
end

par_set_num = 20;
par_consider_idx = 20;
par_consider = zeros(par_set_num, 1);
multi_ss_bool = false(par_set_num, 1);
multi_ss = zeros(2, 14, par_set_num);

p = parpool(20);
tic;
for i = 1:par_set_num
    par_consider(i) = (i - 1) * 5e-2;
    [multi_ss_bool(i), multi_ss(:, :, i)] = if_multi_ss(par_base, par_consider_idx, par_consider(i));
end
delete(p);
toc;

disp(['run time:', num2str(toc)]);
disp(['number of multi_ss: ', num2str(sum(multi_ss_bool))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 1, :)), 'r', par_consider, squeeze(multi_ss(2, 1, :)), 'b');
xlabel('Time (h)'); ylabel('H (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 2, :)), 'r', par_consider, squeeze(multi_ss(2, 2, :)), 'b');
xlabel('Time (h)'); ylabel('I (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 3, :)), 'r', par_consider, squeeze(multi_ss(2, 3, :)), 'b');
xlabel('Time (h)'); ylabel('Virus (log_{10} pfu)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 4, :)) + squeeze(multi_ss(1, 5, :)), 'r', par_consider, squeeze(multi_ss(2, 4, :)) + squeeze(multi_ss(2, 5, :)), 'b');
xlabel('Time (h)'); ylabel('M0 (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 5, :)), 'r', par_consider, squeeze(multi_ss(2, 5, :)), 'b');
xlabel('Time (h)'); ylabel('M (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 6, :)), 'r', par_consider, squeeze(multi_ss(2, 6, :)), 'b');
xlabel('Time (h)'); ylabel('Mono (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 7, :)), 'r', par_consider, squeeze(multi_ss(2, 7, :)), 'b');
xlabel('Time (h)'); ylabel('Neu (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 12, :)), 'r', par_consider, squeeze(multi_ss(2, 12, :)), 'b');
xlabel('Time (h)'); ylabel('NK (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 13, :)), 'r', par_consider, squeeze(multi_ss(2, 13, :)), 'b');
xlabel('Time (h)'); ylabel('T (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 14, :)), 'r', par_consider, squeeze(multi_ss(2, 14, :)), 'b');
xlabel('Time (h)'); ylabel('TE (10^4/ml)');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 8, :)), 'r', par_consider, squeeze(multi_ss(2, 8, :)), 'b');
xlabel('Time (h)'); ylabel('IL1b (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 10, :)), 'r', par_consider, squeeze(multi_ss(2, 10, :)), 'b');
xlabel('Time (h)'); ylabel('CCL2 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 9, :)), 'r', par_consider, squeeze(multi_ss(2, 9, :)), 'b');
xlabel('Time (h)'); ylabel('IL10 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, squeeze(multi_ss(1, 11, :)), 'r', par_consider, squeeze(multi_ss(2, 11, :)), 'b');
xlabel('Time (h)'); ylabel('CXCL5 (pg/ml)'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;