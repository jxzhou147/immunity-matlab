clear;clc

% load clock model parameters
par_clock = load('par_clock.csv');

% solve odes
tspan = 0:1:500;
y0 = zeros(12, 1);

[t, y] = ode45(@ODE_Clock, tspan, y0, [], par_clock);

% plot figures clock mRNAs
figure;
xSize = 20; X=xSize; ySize = 60;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;

subplot(3, 2, 1); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 1), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('Per mRNA');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 2); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 2), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('Cry mRNA');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 3); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 3), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('Rev mRNA');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445], 'YLim', [0 15]);
hold on;

subplot(3, 2, 4); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 4), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('Ror mRNA');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 5); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 5), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('Bmal1 mRNA');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445], 'YLim', [0 2.5]);
hold off;

% plot figures clock proteins
figure;
xSize = 20; X=xSize; ySize = 60;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;

subplot(3, 2, 1); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 6), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('PER');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 2); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 7), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('CRY');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 3); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 8), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('REV');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold on;

subplot(3, 2, 4); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 9), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('ROR');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445], 'YLim', [0.4 0.55]);
hold on;

subplot(3, 2, 5); hold on; set(gca, 'FontSize', 26); box on;
plot(tspan, y(:, 10), 'k', 'LineWidth', 2); hold on;
xlabel('Circadian Time (h)'); ylabel('BMAL1');
xticks([403 415 427 439]); xticklabels({'0', '12', '0', '12'});
set(gca, 'XLim', [397 445]);
hold off;