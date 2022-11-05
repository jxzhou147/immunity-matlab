clear;clc;

% load experimental data
% data_h = xlsread('data_to_fit.xlsx', 3, 'B1:E3');
% data_m = xlsread('data_to_fit.xlsx', 3, 'B5:E7');
% data_mono = xlsread('data_to_fit.xlsx', 3, 'B9:E11');
% data_neu = xlsread('data_to_fit.xlsx', 3, 'B13:E15');
% data_nk = xlsread('data_to_fit.xlsx', 3, 'B17:E19');
% data_v = xlsread('data_to_fit.xlsx', 3, 'B21:I23');
% data_t = xlsread('data_to_fit.xlsx', 3, 'B25:G27');
% data_te = xlsread('data_to_fit.xlsx', 3, 'B29:G31');
% data_il1b = xlsread('data_to_fit.xlsx', 3, 'B33:E35');
% data_ccl2 = xlsread('data_to_fit.xlsx', 3, 'B37:E39');
%     
% data_v(2:3, :) = 10 .^ data_v(2:3, :);

% load immunity and clock-immunity parameters
par_IAV = importdata('par_base_IAV_clock_values.txt');
par_IAV = par_IAV.data;

% translate parFit to par in the odes
% log_par_ind = [1:38 42:53];
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end

% load clock parameters
par_clock = load('par_clock.csv');

% load inflammation parameters
par_infla = importdata('par_infla.txt');
par_infla = par_infla.data;

y0 = zeros(27, 1);
y0_IAV = importdata('init_base.txt');
y0_IAV = y0_IAV.data;
for i = 13:26
    y0(i) = y0_IAV(i - 12);
end

tmax = 1000;
tspan = 0:1:tmax;

% par_IAV(64) = 0.02;
% par_IAV(65) = 10000;
% t_IAV_11 = 127;
% [~, y1] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_11, par_infla);
% infla1 = zeros(tmax+1, 1);
% for i = 1:(tmax+1)
%     infla1(i) = inflammation(y1(i, 17), y1(i, 18), y1(i, 19), par_infla);
% end

par_IAV(68) = 0.02;
par_IAV(69) = 10000;
t_IAV_23 = 115;
[~, y2] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23, par_infla);
infla2 = zeros(tmax+1, 1);
for i = 1:(tmax+1)
    infla2(i) = inflammation(y2(i, 17), y2(i, 18), y2(i, 19), par_infla);
end

% par_IAV(64) = 0.02;
par_IAV(68) = 10000;
par_IAV(69) = 0.02;
% par_clock(8) = 0.6;
t_IAV_23 = 115;
[~, y3] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23, par_infla);
infla3 = zeros(tmax+1, 1);
for i = 1:(tmax+1)
    infla3(i) = inflammation(y3(i, 17), y3(i, 18), y3(i, 19), par_infla);
end
 
% par_clock(8) = 0.4;
par_IAV(68) = 0.02;
par_IAV(69) = 0.02;
t_IAV_23 = 115;
[~, y4] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23, par_infla);
infla4 = zeros(tmax+1, 1);
for i = 1:(tmax+1)
    infla4(i) = inflammation(y4(i, 17), y4(i, 18), y4(i, 19), par_infla);
end

%% 
% clock_name = ["Per", "Cry", "Rev", "Ror", "Bmal1", "PER", "CRY", "REV", "ROR", "BMAL1", "PER-CRY", "CLOCK-BMAL1"];
% figure;
% for i = 1:12
%     subplot(3, 4, i);
%     plot(tspan, y1(:, i), 'r', 'LineWidth', 1); hold on;
%     plot(tspan, y2(:, i), 'b', 'LineWidth', 1); hold on;
% %     plot(tspan, y3(:, i), 'r', 'LineWidth', 1); hold on;
% %     plot(tspan, y4(:, i), 'b', 'LineWidth', 1); hold on;
%     xlim([400 500]);
%     xlabel('Time (h)');
%     ylabel(clock_name(i));
% end
% sgtitle('Time evolution of clock components');
% hold off;

% figure;
% xSize = 30; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
% hold on;
% subplot(1, 3, 1);
% plot(par_IAV(69) ./ (par_IAV(69) + y2(:, 8)) * par_IAV(20), 'k', 'LineWidth', 2); hold on;
% plot([100 200], [0.8 0.8], 'g', 'LineStyle', '--', 'LineWidth', 2); hold off;
% xlim([100 200]); ylim([0.6 1.3]); xlabel('Time (h)'); ylabel('c_{M-CCL2}'); title('Normal Clock');
% subplot(1, 3, 2);
% plot(par_IAV(69) ./ (par_IAV(69) + y3(:, 8)) * par_IAV(20), 'r', 'LineWidth', 2); hold on;
% plot([100 200], [0.8 0.8], 'g', 'LineStyle', '--', 'LineWidth', 2); hold off;
% xlim([100 200]); ylim([0.6 1.3]); xlabel('Time (h)'); ylabel('c_{M-CCL2}'); title('Clock Knock Out');
% subplot(1, 3, 3);
% plot(par_IAV(69) ./ (par_IAV(69) + y4(:, 8)) * par_IAV(20), 'b', 'LineWidth', 2); hold on;
% plot([100 200], [0.8 0.8], 'g', 'LineStyle', '--', 'LineWidth', 2); hold off;
% xlim([100 200]); ylim([0.6 1.3]); xlabel('Time (h)'); ylabel('c_{M-CCL2}'); title('Smaller Amplitude');

var_name = ["H", "If", "V", "M_0", "M", "Mono", "N", "IL1b", "IL10", "CCL2", "CXCL5", "K", "T", "T_E"];
figure;
for i = 13:26
    subplot(4, 4, i - 12);
%     plot(tspan, y(:, i), 'k', 'LineWidth', 1); hold on;
%     plot(tspan, y1(:, i), 'r', 'LineWidth', 1);
%     plot(tspan - t_IAV_11, y1(:, i), 'r', 'LineWidth', 1); hold on;
    plot(tspan - t_IAV_23, y2(:, i), 'r', 'LineWidth', 1); hold on;
    plot(tspan - t_IAV_23, y3(:, i), 'b', 'LineWidth', 1);
    plot(tspan - t_IAV_23, y4(:, i), 'k', 'LineWidth', 1);
%     xlim([0 500]);
    xlabel('Time (h)');
    ylabel(var_name(i - 12)); hold off;
end
subplot(4, 4, 15);
% plot(tspan - t_IAV_11, infla1, 'r', 'LineWidth', 1); hold on;
plot(tspan - t_IAV_23, infla2, 'r', 'LineWidth', 1); hold on;
plot(tspan - t_IAV_23, infla3, 'b', 'LineWidth', 1);
plot(tspan - t_IAV_23, infla4, 'k', 'LineWidth', 1);
xlabel('Time (h)'); ylabel('Inflammation'); hold off;
legend('Only c_{M-IL10} oscillates', 'Only c_{M-CCL2} oscillates', 'Both oscillate', 'Location',[0.8 0.2 0 0]);
% legend('ZT11', 'ZT23', 'Location', [0.8 0.2 0 0]);
% legend('Normal Clock', 'Clock Knock Out', 'Smaller Amplitude', 'Location', [0.8 0.2 0 0]);
% sgtitle('Immune response of different infection time with intact clock');
% sgtitle('The circadian clock and immunity interact in both directions');
% sgtitle('c_{M-CCL2} is controlled by different clocks');
sgtitle('Circadian control of c_{M-CCL2} and c_{M-IL10}');