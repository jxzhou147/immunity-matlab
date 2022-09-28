% plot distribution of variables at different parameters
multi_ss = importdata('merged_multi_ss_488.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 3), multi_ss(:, 4)); hold on;
xlabel('H high'); ylabel('H low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 5), multi_ss(:, 6)); hold on;
xlabel('I high'); ylabel('I low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 7), multi_ss(:, 8)); hold on;
xlabel('V high'); ylabel('V low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 9) + multi_ss(:, 11), multi_ss(:, 10) +multi_ss(:, 12)); hold on;
xlabel('M0 high'); ylabel('M low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 11), multi_ss(:, 12)); hold on;
xlabel('M high'); ylabel('M low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 13), multi_ss(:, 14)); hold on;
xlabel('Mono high'); ylabel('Mono low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 15), multi_ss(:, 16)); hold on;
xlabel('Neu high'); ylabel('Neu low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 25), multi_ss(:, 26)); hold on;
xlabel('NK high'); ylabel('NK low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 27), multi_ss(:, 28)); hold on;
xlabel('T high'); ylabel('T low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 29), multi_ss(:, 30)); hold on;
xlabel('TE high'); ylabel('TE low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 17), multi_ss(:, 18)); hold on;
xlabel('IL1b high'); ylabel('IL1b low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 21), multi_ss(:, 22)); hold on;
xlabel('CCL2 high'); ylabel('CCL2 low');
hold off;

% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 19), multi_ss(:, 20)); hold on;
xlabel('IL10 high'); ylabel('IL10 low');
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
scatter(multi_ss(:, 23), multi_ss(:, 24)); hold on;
xlabel('CXCL5 high'); ylabel('CXCL5 low');
hold off;