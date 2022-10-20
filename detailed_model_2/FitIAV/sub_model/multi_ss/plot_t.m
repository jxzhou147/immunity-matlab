clear;clc

% plot the time evolution of the system

% load base parameters
par_base = importdata('par_M_init.txt');
par_base = par_base.data;
% translate to par in the odes
log_par_ind = (1:7);
for i = log_par_ind
    par_base(i) = 10 .^ par_base(i);
end

% load hls parameters
par_consider_idx = (1:8);
par_consider = importdata('multi_ss_par.txt');
par_consider = par_consider(10, :);
% par_consider = par_base;

% replace parameters considered with par_consider
tmp = 1;
for i = par_consider_idx
    par_base(i) = par_consider(tmp);
    tmp = tmp + 1;
end

% load initial values
y0s = importdata('lhs_init.txt');
y0s = y0s(1:100, :);
[m, n] = size(y0s);

% solve odes
tmax = 10000;
tspan = 0:1:tmax;
ys = zeros(tmax+1, n, m);

p = parpool(20);
tic;
parfor i = 1:m
    [t, y] = ode15s(@ODE_M, tspan, y0s(i, :), [], par_base);
    y = real(y);
    ys(:, :, i) = y;
end
toc;
disp(['run time:', num2str(toc)]);
delete(p);

% plot figuires
dq = jet(m);

figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, ys(:, 1, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('Mono'); hold on;
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
for i = 1:m
    plot(tspan, ys(:, 2, i), 'color', dq(i, :), 'LineWidth', 2);  hold on;
end
xlabel('Time (h)'); ylabel('CCL2');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;

figure;
for i = 1:m
    plot(ys(:, 1, i), ys(:, 2, i), 'color', 'k'); hold on;
end
hold off;