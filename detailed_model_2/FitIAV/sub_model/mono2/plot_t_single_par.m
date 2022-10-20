
% solve and plot results of OED_Mt

% load base parameters
par_base = importdata('par_M_init.txt');
par_base = par_base.data;
% translate to par in the odes
% log_par_ind = (1:7);
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end

tmax = 100;
tspan = 0:1:tmax;
y0 = [0; 0];

par_base(6) = 1;
[t, y_min] = ode15s(@ODE_Mt, tspan, y0, [], par_base, 0);
y_min = real(y_min);

t_lag = 12;
par_base(6) = 3;
[t, y_max] = ode15s(@ODE_Mt, tspan, y0, [], par_base, t_lag);
y_max = real(y_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, y_min(:, 1), 'b', 'LineWidth', 2); hold on;
plot(tspan(t_lag+1:end) - t_lag, y_max(t_lag+1:end, 1), 'r', 'LineWidth', 2);
xlabel('t (h)'); ylabel('Mono');
% legend('c_{I-CCL2} = 1', 'c_{I-CCL2} = 3');
legend('t_0 = 0', 't_0 = 12')
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(tspan, y_min(:, 2), 'b', 'LineWidth', 2); hold on;
plot(tspan(t_lag+1:end) - t_lag, y_max(t_lag+1:end, 2), 'r', 'LineWidth', 2);
xlabel('t (h)'); ylabel('CCL2');
% legend('c_{I-CCL2} = 1', 'c_{I-CCL2} = 3');
legend('t_0 = 0', 't_0 = 12')
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;