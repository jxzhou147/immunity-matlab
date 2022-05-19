clear;clc

y_V = load('y_V.txt');
y_M = load('y_M.txt');
y_K = load('y_K.txt');

tspan = 0:250;

prediction_interval = 95;
p_V = prctile(y_V, [50 - prediction_interval / 2, 50 + prediction_interval / 2], 1);
p_M = prctile(y_M, [50 - prediction_interval / 2, 50 + prediction_interval / 2], 1);
p_K = prctile(y_K, [50 - prediction_interval / 2, 50 + prediction_interval / 2], 1);

t_conf = [tspan tspan(end:-1:1)];

V_conf = [p_V(1, :), p_V(2, end:-1:1)];
M_conf = [p_M(1, :), p_M(2, end:-1:1)];
K_conf = [p_K(1, :), p_K(2, end:-1:1)];

figure;
set(gca, 'FontSize', 20);
p = fill(t_conf, log10(V_conf), 'r');
p.FaceColor = [191, 220, 238] / 255;
p.EdgeColor = 'none';
hold on;
plot(tspan, log10(mean(y_V)), 'Color', [217, 83, 25] / 255)
legend({'95% region', 'mean value'});
xlabel('Time (hour)'); ylabel('log(V)')
hold off;

figure;
set(gca, 'FontSize', 20);
p = fill(t_conf, M_conf, 'r');
p.FaceColor = [191, 220, 238] / 255;
p.EdgeColor = 'none';
hold on;
plot(tspan, mean(y_M), 'Color', [217, 83, 25] / 255)
legend({'95% region', 'mean value'});
xlabel('Time (hour)'); ylabel('M');
hold off;

figure;
set(gca, 'FontSize', 20);
p = fill(t_conf, K_conf, 'r');
p.FaceColor = [191, 220, 238] / 255;
p.EdgeColor = 'none';
hold on;
plot(tspan, mean(y_K), 'Color', [217, 83, 25] / 255)
legend({'95% region', 'mean value'});
xlabel('Time (hour)'); ylabel('K');
hold off;